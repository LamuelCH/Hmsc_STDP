# SPATIAL PREDICTION ###########################################################################################
library(terra)
library(magrittr)
library(dplyr)
library(Hmsc)

load("models/models_thin_10_samples_250_chains_4CMIP.Rdata")
data = read.csv("data/data.csv")

# Prepare known condition Yc for condiional prediction

# Prepare Xdata for the whole Study area
# Recall covariate raster
env_wc2 <- rast("data/env_wc2_1km.tif")
env_der <- rast("data/env_der_1km.tif")

env_stack <- c(env_wc2,env_der) #stack all covariates

grid_1km <- terra::as.data.frame(env_stack, xy = TRUE) %>% na.omit(.) %>% arrange(x) #whole map of Tasmania, spatial sorting by x axis

xy.grid <- as.matrix(cbind(grid_1km$x, grid_1km$y)) #extract spatial information

XData.grid = grid_1km %>% dplyr::select(bio1, bio3, bio5, bio7, bio8, bio12, bio15, Regolith)


# Prepare Gradient
Gradient.mFULL <- prepareGradient(models$mFULL, XDataNew = Xdata.grid, sDataNew = list(ID = xy.grid))
Gradient.mSINGLE <- prepareGradient(models$mSINGLE, XDataNew = Xdata.grid, sDataNew = list(ID = xy.grid))

# We are now ready to compute the posterior predictive distribution
nParallel=parallel::detectCores()-1

{# mFULL spatial prediction
  print(date())
  print(paste0("model = ",names(models)[1]))
  predY.mFULL = predict(models$mFULL, Gradient=Gradient.mFULL, predictEtaMean = TRUE, expected = FALSE, nParallel=nParallel)
  save(predY.mFULL, file = "results/predY_mFULL_thin_1_Samples_1000_PA")
  print(date())
  
  
  # mSINGLE spatial prediction
  print(date())
  print(paste0("model = ",names(models)[2]))
  predY.mSINGLE = predict(models$mSINGLE, Gradient=Gradient.mSINGLE, predictEtaMean = TRUE, expected = FALSE, nParallel=nParallel)
  save(predY.mSINGLE, file = "results/predY_mSINGLE_thin_1_Samples_1000_PA")
  print(date())
}

# EXPLORE PREDICTION OBJECT###########################################################################################
# Let's explore the prediction object for the FULL model.
class(predY.full)
class(predY.single)

# It is a list... 
length(predY.full)
# ...of length 1000, four chains with 250 samples from each

dim(predY.full[[1]])
# Each prediction is a matrix with dimensions 66057 x 3, as there are 66057 prediction locations (pixels) and 3 species.

head(predY.full[[1]])
# Each matrix is filled in with occurrence probabilities
# We may simply by ignoring parameter uncertainty and just looking at 
# the posterior mean prediction. 
EpredY=Reduce("+",predY.single)/length(predY.single)
dim(EpredY)
# EpredY is a 66057 x 4 matrix of posterior mean occurrence probabilities

# ILLUSTRATE PREDICTION ###########################################################################################
# The next step is to post-process the predictions to those community features
# that we wish to illustrate over the prediction space. With the script below,
# we derive from the predictions the occurrence probability of STQ (species number 1),
# the species richness
# We also include data on habitat type and climatic conditions to the dataframe mapData that
# includes all the information we need to visualize the predictions as maps

STQ = EpredY[,1]
TD = EpredY[,2]
EQ = EpredY[,3]
FC = EpredY[,4]
S=rowSums(EpredY) #Species richness
xy = grid_1km[,1:2]

bio1 = Xdata.grid$bio1
bio12 = Xdata.grid$bio12
bio4 = Xdata.grid$bio4
bio15 = Xdata.grid$bio15
roadLength = Xdata.grid$roadLength
Regolith = Xdata.grid$Regolith
forestCover = Xdata.grid$forestCover
foliageCover = Xdata.grid$foliageCover

mapData=data.frame(xy,STQ,#TD,EQ,FC,
                   S,
                   bio1, bio12, bio4, bio15, roadLength, Regolith, forestCover, foliageCover, 
                   stringsAsFactors=TRUE)

# We will use the ggplot function from the ggplot2 package, so let's load the data
library(ggplot2)

# We first plot variation in the habitat and climatic conditions on which the predictions are based on.

ggplot(data = mapData, aes(x=x, y=y, color=quartz))+geom_point(size=2) + ggtitle("quartz") + scale_color_gradient(low="blue", high="red") + coord_equal()
ggplot(data = mapData, aes(x=x, y=y, color=bio1))+geom_point(size=2) + ggtitle("bio1") + scale_color_gradient(low="blue", high="red") + coord_equal()
ggplot(data = mapData, aes(x=x, y=y, color=bio12))+geom_point(size=2) + ggtitle("bio12") + scale_color_gradient(low="blue", high="red") + coord_equal()
ggplot(data = mapData, aes(x=x, y=y, color=Regolith))+geom_point(size=2) + ggtitle("regolith") + scale_color_gradient(low="blue", high="red") + coord_equal()

# We then exemplify prediction for the focal species

ggplot(data = mapData, aes(x=x, y=y, color=STQ))+geom_point(size=0.65) + ggtitle(expression(italic("Dasyurus maculatas")))+ scale_color_gradient(low="gray", high="red") + coord_equal()
ggplot(data = mapData, aes(x=x, y=y, color=TD))+geom_point(size=0.01) + ggtitle(expression(italic("Sarcophilus harrisii")))+ scale_color_gradient(low="gray", high="red") + coord_equal()
ggplot(data = mapData, aes(x=x, y=y, color=EQ))+geom_point(size=0.01) + ggtitle(expression(italic("Dasyurus viverrinus")))+ scale_color_gradient(low="gray", high="red") + coord_equal()
ggplot(data = mapData, aes(x=x, y=y, color=FC))+geom_point(size=0.01) + ggtitle(expression(italic("Felis catus")))+ scale_color_gradient(low="gray", high="red") + coord_equal()

# This prediction is reassuringly very similar to that based on the single-species model of Chapter 5.7
# We next plot predicted species richness, which is highest in Southern Finland

ggplot(data = mapData, aes(x=x, y=y, color=S))+geom_point(size=2) + ggtitle("Species richness")+ scale_color_gradient(low="grey", high="red") + coord_equal()

############################################################################################