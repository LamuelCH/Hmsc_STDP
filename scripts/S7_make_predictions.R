# Step 1: MAKING SPATIAL PREDICTION ####
library(terra)
library(magrittr)
library(dplyr)
library(Hmsc)

##########################################################################
# STEP 2: SETTING UP WORKING ENVIRONMENT ####
# 2.1. set a seed number for replicability
set.seed(1)

# 2.2. Also setup the working directory of the project
setwd(".")
localDir = "."
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir, "models")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir) #create the folder if not existed yet

# 2.3. optional, we can assign the number of cores for computation,
# recommend to leave at least one core for bg task, so we set to use all all core -1
options(mc.cores = parallel::detectCores()-1)
nParallel=parallel::detectCores()-1

# 2.4. load data and models into the environment
load("models/models_thin_10_samples_250_chains_4CMIP.Rdata") 
data = read.csv(file.path(dataDir, "data.csv"),stringsAsFactors=TRUE)

##################################################################################################
# STEP 3: CONSTRUCT DATA for the ENTIRE STUDY AREA ####
# 3.1. Rasterize data into RStudio and stack them
env_wc2 <- rast("data/env_wc2_1km.tif") #replace it with CMIP6 data for climate projection 
env_der <- rast("data/env_der_1km.tif") # Depth of Regolith

env_stack <- c(env_wc2,env_der) #stack all raster

# 3.2. Convert raster stack into dataframe and extract coordinates
grid_1km <- terra::as.data.frame(env_stack, xy = TRUE) %>%
  na.omit(.) %>% # whole map of Tasmania, dropping out NA 
  arrange(x) # [optional?] sort by x axis

xy.grid <- as.matrix(cbind(grid_1km$x, grid_1km$y)) # extract coordinate information

# 3.3. Select the covariates for prediction, this should be matched with the trained model
XData.grid = grid_1km %>% dplyr::select(bio1, bio3, bio5, bio7, bio8, bio12, bio15, Regolith)



##################################################################################################
# STEP 4: SPATIAL PREDICTION ####
# 4.1.1 Making conditional spatial prediction using the new information ####
# 4.1.2 Create Y matrix for conditional prediction
Yc = as.matrix(data[,c(5:7)]) 
Yc = apply(Yc, MARGIN = 2,FUN = as.numeric) # make sure the matrix is numeric instead of integer

# 3.5. Define new study design and spatial random level for the study area
studyDesign.new = data.frame(sample = as.factor(1:nrow(XData.grid))) # sample number will be the total grid number of the study area 
rL.new = HmscRandomLevel(sData = xy.grid, sMethod = "NNGP", longlat = FALSE) 
rL.new = setPriors(rL.new, nfMin=2, nfMax=3)

system.time(
Pred.Yc = predict(models$mFULL, 
               XData= XData.grid,
               studyDesign=studyDesign.new, 
               ranLevels = list(sample=rL_new), 
               Yc= Yc, 
               predictEtaMeanField = TRUE, # Incorporate the uncertainty and possible variability of the random effects into the predictions.
               expected=TRUE) 
)

save(Pred.Yc, file = "results/predYc_mFULL_thin_1_Samples_250_CMIP")

# 4.1.2 [Alternative] Making spatial prediction without considering species association ####
# Prepare Gradient
Gradient =  prepareGradient(models$mFULL,
                            XDataNew = XData.grid, 
                            sDataNew = list(sample = xy.grid))

# Run your prediction
{print(date())
Pred.Y <- predict(models$mFULL, 
                  Gradient = Gradient, 
                  predicEtaMeanField = TRUE, 
                  expected = TRUE, 
                  nParallel = nParallel)
print(date())}

save(Pred.Y, file = "results/predY_mFULL_thin_100_Samples_250_CMIP")

##################################################################################################

