# This script was made seperately to show prediction of the HMSC model
# Step 1: SETTING UP WORKING ENVIRONMENT ####
library(terra)
library(magrittr)
library(dplyr)
library(Hmsc)

set.seed(1) 

# 1.1. Also setup the working directory of the project
setwd(".")
localDir = "."
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir, "models")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir) #create the folder if not existed yet

# 1.2. optional, we can assign the number of cores for computation,
# recommend to leave at least one core for bg task, so we set to use all all core -1
options(mc.cores = parallel::detectCores()-1)

# 1.3. load data and models into the environment
load("results/predY_mFULL_thin_10_Samples_250_CMIP") 
data = read.csv(file.path(dataDir, "data.csv"),stringsAsFactors=TRUE)

# Step 2: EXPLORE PREDICTION OBJECT ####
class(Pred.Y) #list
length(Pred.Y) #length = 1000, with 250 samples for four chains
dim(Pred.Y[[1]]) #Each prediction is a matrix with dimensions 66057 x 3
head(Pred.Y[[1]]) #sampling unit filled with occurrence probability

EpredY=Reduce("+",Pred.Y)/length(Pred.Y)
dim(EpredY)

# STEP 5: ILLUSTRATE PREDICTION ####
# 5.1. Construct the entire Study area
env_wc2 <- rast("data/env_wc2_1km.tif") #replace it with CMIP6 data for climate projection 
env_der <- rast("data/env_der_1km.tif") # Depth of Regolith

env_stack <- c(env_wc2,env_der) #stack all raster

# 5.3. Convert raster stack into dataframe and extract coordinates
grid_1km <- terra::as.data.frame(env_stack, xy = TRUE) %>%
  na.omit(.) %>% # whole map of Tasmania, dropping out NA 
  arrange(x) # [optional?] sort by x axis

xy.grid <- as.matrix(cbind(grid_1km$x, grid_1km$y)) # extract coordinate information

# 5.4. Select the covariates for prediction, this should be matched with the trained model
XData.grid = grid_1km %>% dplyr::select(bio1, bio3, bio5, bio7, bio8, bio12, bio15, Regolith)

# 5.5. Create dataframe for a normal prediction
TD = EpredY[,1] #Tasmania Devil
EQ = EpredY[,2] #Eastern Quoll
STQ = EpredY[,3] #Spotted-tail quoll
S=rowSums(EpredY) #Species richness
xy = grid_1km[,1:2]

mapData = data.frame(xy, STQ, TD, EQ, S, XData.grid, stringsAsFactors = TRUE)

# We will use the ggplot function from the ggplot2 package, so let's load the data
library(ggplot2)

ggplot(data = mapData, aes(x=x, y=y, color=STQ))+geom_point(size=0.01) + ggtitle(expression(italic("Dasyurus maculatas")))+ scale_color_gradient(low="gray", high="red") + coord_equal()
ggplot(data = mapData, aes(x=x, y=y, color=TD))+geom_point(size=0.01) + ggtitle(expression(italic("Sarcophilus harrisii")))+ scale_color_gradient(low="gray", high="red") + coord_equal()
ggplot(data = mapData, aes(x=x, y=y, color=EQ))+geom_point(size=0.01) + ggtitle(expression(italic("Dasyurus viverrinus")))+ scale_color_gradient(low="gray", high="red") + coord_equal()
ggplot(data = mapData, aes(x=x, y=y, color=FC))+geom_point(size=0.01) + ggtitle(expression(italic("Felis catus")))+ scale_color_gradient(low="gray", high="red") + coord_equal()

# This prediction is reassuringly very similar to that based on the single-species model of Chapter 5.7
# We next plot predicted species richness, which is highest in Southern Finland

ggplot(data = mapData, aes(x=x, y=y, color=S))+geom_point(size=0.01) + ggtitle("Species richness")+ scale_color_gradient(low="grey", high="red") + coord_equal()

############################################################################################