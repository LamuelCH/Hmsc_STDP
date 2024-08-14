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
aoi.pre = ext(140, 150, -45, -38) #use it to pre crop CMIP data before projection to reduce computation burden
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

##########################################################################
# STEP 6: CREATE FIXED EFFECT DATAFRAME ####
# List of new names for the bioclimatic variables
bio_names <- paste0("bio", 1:19)

# SSP126 2021-2040

# SSP585 2021-2040, there are 12 GCMs available in total
{
ssp585_2021.access_1km = rast("data/raw/CMIP6/2021-2040/ssp585/wc2.1_30s_bioc_ACCESS-CM2_ssp585_2021-2040.tif") %>% crop(aoi.pre) %>% terra::project(., env_tasMask)*env_tasMask
ssp585_2021.CMCC_1km = rast("data/raw/CMIP6/2021-2040/ssp585/wc2.1_30s_bioc_CMCC-ESM2_ssp585_2021-2040.tif") %>% crop(aoi.pre) %>% terra::project(., env_tasMask) *env_tasMask
ssp585_2021.EC_1km = rast("data/raw/CMIP6/2021-2040/ssp585/wc2.1_30s_bioc_EC-Earth3-Veg_ssp585_2021-2040.tif") %>% crop(aoi.pre) %>% terra::project(., env_tasMask) *env_tasMask
ssp585_2021.FIO_1km = rast("data/raw/CMIP6/2021-2040/ssp585/wc2.1_30s_bioc_FIO-ESM-2-0_ssp585_2021-2040.tif") %>% crop(aoi.pre) %>% terra::project(., env_tasMask) *env_tasMask
ssp585_2021.GISS_1km = rast("data/raw/CMIP6/2021-2040/ssp585/wc2.1_30s_bioc_GISS-E2-1-G_ssp585_2021-2040.tif") %>% crop(aoi.pre) %>% terra::project(., env_tasMask) *env_tasMask
ssp585_2021.HADGEM3_1km = rast("data/raw/CMIP6/2021-2040/ssp585/wc2.1_30s_bioc_HadGEM3-GC31-LL_ssp585_2021-2040.tif") %>% crop(aoi.pre) %>% terra::project(., env_tasMask) *env_tasMask
ssp585_2021.INM_1km = rast("data/raw/CMIP6/2021-2040/ssp585/wc2.1_30s_bioc_INM-CM5-0_ssp585_2021-2040.tif") %>% crop(aoi.pre) %>% terra::project(., env_tasMask) *env_tasMask
ssp585_2021.IPSL_1km = rast("data/raw/CMIP6/2021-2040/ssp585/wc2.1_30s_bioc_IPSL-CM6A-LR_ssp585_2021-2040.tif") %>% crop(aoi.pre) %>% terra::project(., env_tasMask) *env_tasMask
ssp585_2021.MICROC6_1km = rast("data/raw/CMIP6/2021-2040/ssp585/wc2.1_30s_bioc_MIROC6_ssp585_2021-2040.tif") %>% crop(aoi.pre) %>% terra::project(., env_tasMask)*env_tasMask
ssp585_2021.MPI_1km = rast("data/raw/CMIP6/2021-2040/ssp585/wc2.1_30s_bioc_MPI-ESM1-2-HR_ssp585_2021-2040.tif") %>% crop(aoi.pre) %>% terra::project(., env_tasMask) *env_tasMask
ssp585_2021.MRI_1km = rast("data/raw/CMIP6/2021-2040/ssp585/wc2.1_30s_bioc_MRI-ESM2-0_ssp585_2021-2040.tif") %>% crop(aoi.pre) %>% terra::project(., env_tasMask) *env_tasMask
ssp585_2021.UKESM1_1km = rast("data/raw/CMIP6/2021-2040/ssp585/wc2.1_30s_bioc_UKESM1-0-LL_ssp585_2021-2040.tif") %>% crop(aoi.pre) %>% terra::project(., env_tasMask) *env_tasMask
}

# Renaming the layers for each raster
{names(ssp585_2021.access_1km) <- bio_names
names(ssp585_2021.CMCC_1km) <- bio_names
names(ssp585_2021.EC_1km) <- bio_names
names(ssp585_2021.FIO_1km) <- bio_names
names(ssp585_2021.GISS_1km) <- bio_names
names(ssp585_2021.HADGEM3_1km) <- bio_names
names(ssp585_2021.INM_1km) <- bio_names
names(ssp585_2021.IPSL_1km) <- bio_names
names(ssp585_2021.MICROC6_1km) <- bio_names
names(ssp585_2021.MPI_1km) <- bio_names
names(ssp585_2021.MRI_1km) <- bio_names
names(ssp585_2021.UKESM1_1km) <- bio_names
}

{writeRaster(ssp585_2021.access_1km, "data/CMIP6/2021-2040/ssp585/ssp585_2021_ACCESS_1km.tif", overwrite = TRUE)
writeRaster(ssp585_2021.CMCC_1km, "data/CMIP6/2021-2040/ssp585/ssp585_2021_CMCC_1km.tif", overwrite = TRUE)
writeRaster(ssp585_2021.EC_1km, "data/CMIP6/2021-2040/ssp585/ssp585_2021_EC_1km.tif", overwrite = TRUE)
writeRaster(ssp585_2021.FIO_1km, "data/CMIP6/2021-2040/ssp585/ssp585_2021_FIO_1km.tif", overwrite = TRUE)
writeRaster(ssp585_2021.GISS_1km, "data/CMIP6/2021-2040/ssp585/ssp585_2021_GISS_1km.tif", overwrite = TRUE)
writeRaster(ssp585_2021.HADGEM3_1km, "data/CMIP6/2021-2040/ssp585/ssp585_2021_HADGEM3_1km.tif", overwrite = TRUE)
writeRaster(ssp585_2021.INM_1km, "data/CMIP6/2021-2040/ssp585/ssp585_2021_INM_1km.tif", overwrite = TRUE)
writeRaster(ssp585_2021.IPSL_1km, "data/CMIP6/2021-2040/ssp585/ssp585_2021_IPSL_1km.tif", overwrite = TRUE)
writeRaster(ssp585_2021.MICROC6_1km, "data/CMIP6/2021-2040/ssp585/ssp585_2021_MICROC6_1km.tif", overwrite = TRUE)
writeRaster(ssp585_2021.MPI_1km, "data/CMIP6/2021-2040/ssp585/ssp585_2021_MPI_1km.tif", overwrite = TRUE)
writeRaster(ssp585_2021.MRI_1km, "data/CMIP6/2021-2040/ssp585/ssp585_2021_MRI_1km.tif", overwrite = TRUE)
writeRaster(ssp585_2021.UKESM1_1km, "data/CMIP6/2021-2040/ssp585/ssp585_2021_UKESM1_1km.tif", overwrite = TRUE)}

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