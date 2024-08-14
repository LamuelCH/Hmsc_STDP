# Install and load the dplyr tidyr package
#install.packages("dplyr")
#install.packages("tidyr")
#Load necessary library
  library(tidyr)
  library(dplyr)
  library(terra)
  library(ggplot2)
  library(ggmap)
  library(maps)
  library(mapdata)
  
  library(stringr)
  library(sp)
  library(rgdal)
  library(pracma)
  library(geosphere)
  
  #library(sf)
  #library(raster)
  #library(geodata)
  
##########################################################################
# STEP 1: SETTING UP WORKING ENVIRONMENT ####
# 1.1. set a seed number for replicability
set.seed(1)

# 1.2. Also setup the working directory of the project
setwd(".")
localDir = "."
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir, "models")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir) #create the folder if not existed yet

# 1.3. optional, we can assign the number of cores for computation,
# recommend to leave at least one core for bg task, so we set to use all all core -1
options(mc.cores = parallel::detectCores()-1)

###########################################################################
# STEP 2: SETTING UP STUDY AREA AND SPATIAL SCALE ####
# 2.1. Restrict study area to TAS only 
aoi <- ext(225231.8, 629809.5, 5151752.7, 5660841.6) #make sure it is GDA 94 ZONE 55

# 2.2. Also set up the target coordinate reference system and spatial scale
target_crs <- "EPSG:32755" #UTM 55S
target_resolution = 1000 #5000x5000m i.e. 25km sq, change it if necessary 

###########################################################################
# STEP 3: EMPTY RASTER 1KM ####
# 3.1. Create an empty raster layer for projection raster projection
raster <- rast (ext = aoi,
                crs = target_crs,
                resolution = target_resolution)

# 3.2. Input the Tasmania mask raster and project to relevant crs and resolution
env_tasMask <- rast("data/env_tasMask_1km.tif") %>% terra::project(., raster)

###########################################################################
# STEP 4: OCCURRENCE DATA WRANGLING ####
# 4.1. Read raw data
raw = read.csv("data/raw/STDP_2012-2023.csv")

# 4.2. Defining focal species and fields we want to contain in the dataset
focal_species = c(unique(raw$SPECIES))
focal_fields = c("EASTING", "NORTHING", "LONGITUDE", "LATITUDE", "SITE", "LOCATION_FOREIGN_ID", "ACTIVITY_DATE", "SPECIES")

# 4.3. convert the raw data set into standard input of HMSC
# Filter the data frame to include only the specified species and drop NA
occ = raw %>%
  filter(SPECIES %in% focal_species) %>%                            # include only focal species
  mutate(SPECIES = ifelse(SPECIES == "", "Absence", SPECIES)) %>%   #replace empty string with Absence
  filter(SITE != "") %>%                                            # drop out empty site
  dplyr::select(all_of(focal_fields)) %>% 
  drop_na()                                                         # drop all row with NA value

# rename the columns for better readability
occ = data.frame(x = occ$EASTING,
                 y = occ$NORTHING,
                 lon = occ$LONGITUDE,
                 lat = occ$LATITUDE,
                 site = occ$SITE,
                 trapID = occ$LOCATION_FOREIGN_ID,
                 date = occ$ACTIVITY_DATE,
                 species = occ$SPECIES)

# visualize the data so to have a glimpse on how species distributed
# recall online map registration
register_stadiamaps("1794b86f-c7ad-4aab-82b4-23ae3e2f039d", write = TRUE)

# create plots for site and species respectively
plot.SITE = qmplot(lon, lat, data = occ, zoom = 9, maptype = "stamen_terrain", darken = c(0.5, "black"), color = site)
plot.SPECIES = qmplot(lon, lat, data = occ, zoom = 9, maptype = "stamen_terrain", darken = c(0.5, "black"), color = species)

ggsave("site.png", plot = plot.SITE, device = "png", path = resultDir)
ggsave("species.png", plot = plot.SPECIES, device = "png", path = resultDir)


# Summarize species presence/absence at each site, as we are working on Easting and Northing, drop out the lon lat and only group by E N
summary_df <- occ %>%
  group_by(x, y, site) %>%
  summarise(
    TD = ifelse(any(species == "Tasmanian Devil"), 1, 0),
    EQ = ifelse(any(species == "Eastern Quoll"), 1, 0),
    STQ = ifelse(any(species == "Spotted-tailed Quoll"), 1, 0)
  )

# Check for any duplicated x and y coordinates
duplicate_indices <- duplicated(summary_df[c("x", "y")]) | duplicated(summary_df[c("x", "y")], fromLast = TRUE)
unique(duplicate_indices) # should only contain "FALSE"

###########################################################################
# STEP 5: CREATE DF MATCH WITH XDATA RESOLUTION
# 5.1. Convert summary_df to a SpatVector
summary_vect <- vect(summary_df, geom = c("x", "y"), crs = crs(env_tasMask))

# 5.2. Extract coordinate
coords <- geom(summary_vect)[, c("x", "y")]

# 5.3.Extract cell numbers from the Tasmania mask raster
cell_numbers <- cellFromXY(env_tasMask, coords)

# 5.4. Add these cell numbers as a new column in the summary_df
summary_df$cell <- cell_numbers

# 5.5. Aggregate by the cell numbers
coarse_summary_df <- summary_df %>%
  group_by(cell) %>%
  summarise(
    x = mean(x),
    y = mean(y),
    site = paste(unique(site), collapse = ", "),
    TD = max(TD),
    EQ = max(EQ),
    STQ = max(STQ)
  ) %>%
  ungroup() %>%          # Ungroup the data to avoid issues with grouping
  select(-cell)          # Drop the 'cell' column

##########################################################################
# STEP 6: CREATE FIXED EFFECT DATAFRAME ####
# 6.1. Recall covariates data from directory 
env_wc2 = rast("data/env_wc2_1km.tif")
env_forestCover = rast("data/env_forestCover_1km.tif")
env_foliageCover = rast("data/env_foliage_1km.tif")
env_der = rast("data/env_der_1km.tif")
env_road = rast("data/env_roadDensity_10km_1km.tif")

env.stack = c(env_wc2,
               env_forestCover,
               env_foliageCover,
               env_der,env_road) %>% 
  terra::project(., env_tasMask)

# 6.2. extract value from the survey point
env.values <- terra::extract(env.stack, data.frame(coarse_summary_df$x, coarse_summary_df$y))

##########################################################################
# STEP 7 : CREATE FINAL DATAFRAME ####
#Combine the occurrence data and env covariate, drop NA value, then save seperately 
df <- cbind(coarse_summary_df, env.values) %>% 
  drop_na()

df <- df[, colSums(df != 0) > 0]

# Remove existing ID columns if present
df <- df %>% dplyr::select(-starts_with("ID"))

# Add a new ID column with row numbers
df <- df %>%
  mutate(ID = seq_len(nrow(.))) %>%
  dplyr::select(ID, everything())

# Arrange the coordinates according to the leading eigenvalue
X = cbind(df$x, df$y)
pc = prcomp(X)                                                                  #PCA finds the directions (principal components) in which the data varies the most.
proj = X %*% pc$rotation[,1]
optOrder = rank(proj)

final_df = df[optOrder,]

# Check if the final dataframe contain dulplicate coordinates
{duplicates <- final_df %>%
    group_by(x, y) %>%
    filter(n() > 1) %>%
    ungroup()
  
  if (nrow(duplicates) > 0) {
    print("There are duplicate coordinates in the dataframe.")
  } else {
    print("There are no duplicate coordinates in the dataframe.")
  }
}

# Renumber the ID column to be continuous
final_df$ID <- seq_len(nrow(final_df))

# Reorder columns
final_df <- final_df %>%
  ungroup() %>% 
  dplyr::select(ID, everything())

summary(final_df)

# Create the plot
ggplot(final_df, aes(x = x, y = y, color = ID)) +
  geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Spatial Points Ordered by PCA",
       x = "Easting",
       y = "Northing",
       color = "Row Number") +
  theme_minimal()

################################################
#save env_values as environment.csv for modelling

write.csv(final_df, file = "data/data.csv", row.names = F)

#environment_csv <- final_df[, c(1,55:77)]
#write.csv(environment_csv, file = "data/atlas_environment_5km.csv", row.names = F)

#species_csv <- final_df[, c(1:54)]
#write.csv(x = species_csv, file = "data/atlas_species_5km.csv", row.names = F)

################################################
"DATA READY FOR HMSC"
