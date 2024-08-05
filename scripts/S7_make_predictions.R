# Set the base directory using your favorite method
# setwd("...")

##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (BEGINNING)
##################################################################################################
#	INPUT. the Fitted models.

#	OUTPUT. Predictions over environmental gradients (for highest RUN of S2) in the file
# "results/predictions.pdf".
##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (END)
##################################################################################################


##################################################################################################
# MAKE THE SCRIPT REPRODUCIBLE (BEGINNING)
##################################################################################################
set.seed(1)
##################################################################################################
## MAKE THE SCRIPT REPRODUCIBLE (END)
##################################################################################################


##################################################################################################
# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (BEGINNING)
##################################################################################################
species.list = NULL #one example species shown for each model,
#selected as prevalence closest to 0.5 (probit models) or most abundant species (other models)
trait.list = NULL #community weighted mean shown for all traits
env.list = NULL #predictions constructed over all environmental gradients
##################################################################################################
# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (END)
##################################################################################################

##################################################################################################
# CHANGE DEFAULT OPTIONS BY REMOVING COMMENT AND SETTING VALUE (BEGINNING)
# NOTE THAT THIS IS THE ONLY SECTION OF THE SCRIPT THAT YOU TYPICALLY NEED TO MODIFY
##################################################################################################
#use species.list to select which species are used as examples for which predictions are shown
#species.list should be a list of length the number of models. 
#for each element provide either 0 (use default) or a vector of species indices
species.list = list()
species.list[[1]] = 0
species.list[[2]] = 0
species.list[[3]] = c(1,2)

#use trait.list to select for which traits predictions for community weighted mean traits are shown
#trait.list should be a list of length the number of models. 
#for each element provide either 0 (use default) or a vector of trait indices
#see models[[j]]$trNames to see which trait each index corresponds to
#trait.list = list()
#trait.list[[1]] = c(2,10)
#trait.list[[2]] = 0

#use env.list to select over which environmental gradients predictions are generated
#env.list should be a list of length the number of models. 
#for each element provide either 0 (use default) or a vector of environmental variables
#env.list = list()
#env.list[[1]] = 0
#env.list[[2]] = c("tree","decay")
##################################################################################################
# CHANGE DEFAULT OPTIONS BY REMOVING COMMENT AND SETTING VALUE (END)
# NOTE THAT THIS IS THE ONLY SECTION OF THE SCRIPT THAT YOU TYPICALLY NEED TO MODIFY
##################################################################################################

##################################################################################################
# SET DIRECTORIES (BEGINNING)
##################################################################################################
localDir = "."
modelDir = file.path(localDir, "models")
resultDir = file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)
##################################################################################################
# SET DIRECTORIES (END)
##################################################################################################


library(Hmsc)
library(ggplot2)

samples_list = c(1000)
thin_list = c(1)
nst = length(thin_list)
nChains = 4

for (Lst in nst:1) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  filename = file.path(modelDir,paste("models_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),
                                      ".Rdata",sep = ""))
  if(file.exists(filename)){break}
}




if(file.exists(filename)){
  load(filename)
  nm = length(models)
  modelnames = names(models)
  
  if (is.null(species.list) || length(species.list) < nm) {
    species.list = vector("list", nm)
    for (j in 1:nm) species.list[[j]] = list() # Use an empty list or an appropriate default value
  }
  
  if (is.null(trait.list) || length(trait.list) < nm) {
    trait.list = vector("list", nm)
    for (j in 1:nm) trait.list[[j]] = list() # Use an empty list or an appropriate default value
  }
  
  if (is.null(env.list) || length(env.list) < nm) {
    env.list = vector("list", nm)
    for (j in 1:nm) env.list[[j]] = list() # Use an empty list or an appropriate default value
  }
  
#  if(is.null(species.list)){
#    species.list = list()
#    for(j in 1:nm) species.list[[j]] = 0
#  }
#  if(is.null(trait.list)){
#    trait.list = list()
#    for(j in 1:nm) trait.list[[j]] = 0
#  }
#  if(is.null(env.list)){
#    env.list = list()
#    for(j in 1:nm) env.list[[j]] = 0
#  }
  
  pdf(file= file.path(resultDir,"predictions.pdf"))
  for(j in 1:nm){
    m = models[[j]]
    if(all(env.list[[j]]==0)){
      if(m$XFormula=="~."){
        covariates = colnames(m$XData)
      } else {
        covariates = all.vars(m$XFormula)
      }
    } else {
      covariates = env.list[[j]]
    }
    ex.sp = which.max(colMeans(m$Y,na.rm = TRUE)) #most common species as example species
    if(m$distr[1,1]==2){
      ex.sp = which.min(abs(colMeans(m$Y,na.rm = TRUE)-0.5))
    }
    if(!all(species.list[[j]])==0){
      ex.sp = species.list[[j]]
    }
    if(length(covariates)>0){
      for(k in 1:(length(covariates))){
        covariate = covariates[[k]]
        Gradient = constructGradient(m,focalVariable = covariate)
        Gradient2 = constructGradient(m,focalVariable = covariate,non.focalVariables = 1)
        predY = predict(m, Gradient=Gradient, expected = TRUE)  
        predY2 = predict(m, Gradient=Gradient2, expected = TRUE)  
        par(mfrow=c(2,1))
        pl = plotGradient(m, Gradient, pred=predY, yshow = 0, measure="S", showData = TRUE, 
                          main = paste0(modelnames[j],": summed response (total effect)"))
        if(inherits(pl, "ggplot")){
          print(pl + labs(title=paste0(modelnames[j],": summed response (total effect)")))
        }
        pl = plotGradient(m, Gradient2, pred=predY2, yshow = 0, measure="S", showData = TRUE, 
                          main = paste0(modelnames[j],": summed response (marginal effect)"))
        if(inherits(pl, "ggplot")){
          print(pl + labs(title=paste0(modelnames[j],": summed response (marginal effect)")))
        }
        for(l in 1:length(ex.sp)){
          par(mfrow=c(2,1))
          pl = plotGradient(m, Gradient, pred=predY, yshow = if(m$distr[1,1]==2){c(-0.1,1.1)}else{0}, measure="Y",index=ex.sp[l], showData = TRUE, 
                            main = paste0(modelnames[j],": example species (total effect)"))
          if(inherits(pl, "ggplot")){
            print(pl + labs(title=paste0(modelnames[j],": example species (total effect)")))
          }
          pl = plotGradient(m, Gradient2, pred=predY2, yshow = if(m$distr[1,1]==2){c(-0.1,1.1)}else{0}, measure="Y",index=ex.sp[l], showData = TRUE, 
                            main = paste0(modelnames[j],": example species (marginal effect)"))
          if(inherits(pl, "ggplot")){
            print(pl + labs(title=paste0(modelnames[j],": example species (marginal effect)")))
          }
        }
        if(m$nt>1){
          traitSelection = 2:m$nt
          if(!all(trait.list[[j]]==0)) traitSelection = trait.list[[j]]
          for(l in traitSelection){
            par(mfrow=c(2,1))
            pl = plotGradient(m, Gradient, pred=predY, measure="T",index=l, showData = TRUE,yshow = 0,
                              main = paste0(modelnames[j],": community weighted mean trait (total effect)"))
            if(inherits(pl, "ggplot")){
              print(pl + labs(title=paste0(modelnames[j],": community weighted mean trait (total effect)")))
            }
            pl = plotGradient(m, Gradient2, pred=predY2, measure="T",index=l, showData = TRUE, yshow = 0,
                              main = paste0(modelnames[j],": community weighted mean trait (marginal effect)"))
            if(inherits(pl, "ggplot")){
              print(pl + labs(title=paste0(modelnames[j],": community weighted mean trait (marginal effect)")))
            }
          }
        }
      }
    }
  }
  dev.off()
}


# SPATIAL PREDICTION ###########################################################################################
library(terra)
library(magrittr)
library(dplyr)
library(Hmsc)
load("models/models_thin_1_samples_1000_chains_4.Rdata")

# MAKING PREDICTION FOR THE ENTIRE TASMANIA ####
# Construct whole map of Tasmania
env_wc2 <- rast("data/env_wc2_1km.tif")
env_forestCover <- rast("data/env_forestCover_1km.tif")
env_foliageCover <- rast("data/env_foliage_1km.tif")
env_der <- rast("data/env_der_1km.tif")
env_road <- rast("data/env_roadDensity_10km_1km.tif")

env_stack <- c(env_wc2,env_forestCover,env_foliageCover,env_der,env_road) #stack all covariates

grid_1km <- terra::as.data.frame(env_stack, xy = TRUE) %>% na.omit(.) %>% arrange(x) #whole map of Tasmania, spatial sorting by x axis
xy.grid <- as.matrix(cbind(grid_1km$x, grid_1km$y))

Xdata.grid <- data.frame(bio1 = grid_1km$bio1,
                         bio12 = grid_1km$bio12,
                         bio4 = grid_1km$bio4,
                         bio15 = grid_1km$bio15,
                         roadLength = grid_1km$roadLength,
                         Regolith = grid_1km$Regolith,
                         forestCover = grid_1km$forestCover,
                         foliageCover = grid_1km$foliageCover)

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