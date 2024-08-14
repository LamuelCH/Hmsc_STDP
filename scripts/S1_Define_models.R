# Set the base directory
# setwd("...")

##################################################################################################
# INPUT AND OUTPUT OF THIS SCRIPT (BEGINNING)
##################################################################################################
#	INPUT. Original datafiles of the case study, placed in the data folder.

#	OUTPUT. Unfitted models, i.e., the list of Hmsc model(s) that have been defined but not fitted yet,
# stored in the file "models/unfitted_models.RData".
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
# LOAD PACKAGES (BEGINNING)
##################################################################################################
library(Hmsc)
library(GGally)
library(tidyverse)
library(ggplot2)
library(ape)
##################################################################################################
# LOAD PACKAGES (END)
##################################################################################################


##################################################################################################
# SET DIRECTORIES (BEGINNING)
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
#options(mc.cores = parallel::detectCores()-1)

###########################################################################
# SET DIRECTORIES (END)
##################################################################################################


##################################################################################################
# READ AND SELECT SPECIES DATA (BEGINNING)
##################################################################################################
#data1 = read.csv(file.path(dataDir, "atlas_data.csv"),stringsAsFactors=TRUE)
#data = read.csv(file.path(dataDir, "atlas_data_5km.csv"),stringsAsFactors=TRUE)
data = read.csv(file.path(dataDir, "data.csv"),stringsAsFactors=TRUE)
# Create the plot to confirm if the order is good 
ggplot(data, aes(x = x, y = y, color = ID)) +
  geom_point(size = 2) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = "Longitude",
       y = "Latitude",
       color = "Row Number") +
  theme_minimal()

##################################################################################################
Y = as.matrix(data[,c(5:7)])

Y.single = as.matrix(data$STQ)
colnames(Y.single) = "Dasyurus.maculatus"

# It is important to ensure the presence/absence are presented in numeric but not integer,
# which will potentially cause error when predicting gradients.
Y = apply(Y,MARGIN = 2,FUN = as.numeric)

Y.single = apply(Y.single, MARGIN = 2,FUN = as.numeric)


hist(Y)
##################################################################################################
# 06 DEFINE XY MATRIX ###############################################################################################
# The data contains also the x- and y-coordinates of each camera trap. 
# We will store these as the xy-matrix to be able to fit a spatial model
xy = as.matrix(cbind(data$x,data$y)) 
rownames(xy)=data$ID #set row name to ID name
colnames(xy)=c("x-coordinate","y-coordinate") # rename the column

##################################################################################################


# 07 INSPECT HMSC DATA #################################################################################################
# It is always a good idea to eyeball the data. Let us first have a look at the community data.

dim(Y)
dim(Y.single)

head(Y)
head(Y.single)##################################################################################################


##################################################################################################
data$ID = as.factor(data$ID) #convert to FACTOR
#data$effort = as.factor(data$effort)

XData = data %>% dplyr::select(bio1, bio3, bio5, bio7, bio8, bio12, bio15, Regolith)

ggcorr(XData, label =TRUE)

##################################################################################################

###############################################################################################

##################################################################################################
# SELECT COMMON SPECIES (END)
##################################################################################################


##################################################################################################
# SET UP THE MODEL (BEGINNING)
##################################################################################################
# Spatial random effect
studyDesign = data.frame(sample = as.factor(1:nrow(data))) #SITE = XData$site) #effort = XData$effort
rL.spatial = HmscRandomLevel(sData = xy, sMethod = "NNGP", longlat = FALSE)
rL.spatial = setPriors(rL.spatial, nfMin=2, nfMax=3) #limit the model to one latent variables for visualization purpose

# RANDOM EFFECT STRUCTURE, HERE Site (hierarchical study design)
# and optionally id, if we are interested in species associations at that level
#rL.site = HmscRandomLevel(unit = levels(studyDesign$site))
#rL.id = HmscRandomLevel(sData=xy, sMethod = "NNGP", nNeighbours = 100)
#rL.id = setPriors(rL.id, nfMin=1, nfMax=3) #limit the model to one latent variables for visualization purpose

XFormula = as.formula(paste('~', paste("poly(",colnames(XData),",", "degree = 2", ",", "raw = TRUE)", collapse = '+')))

# CONSTRUCT THE MODELS.

# PRESENCE-ABSENCE MODEL FOR INDIVIDUAL SPECIES (COMMON ONLY)
mFULL = Hmsc(Y=Y, XData = XData, XFormula=XFormula,
             distr="probit", studyDesign=studyDesign,
             ranLevels=list(sample = rL.spatial))

mSINGLE = Hmsc(Y=Y.single, XData = XData, XFormula=XFormula,
               distr="probit", studyDesign=studyDesign, ranLevels=list(sample=rL.spatial))

##################################################################################################
# SET UP THE MODEL (END)
##################################################################################################


##################################################################################################
# COMBINING AND SAVING MODELS (START)
##################################################################################################
models = list(mFULL) #mSINGLE) #mSPACE)
names(models) = c("mFULL")#, "mSINGLE")#"mSPACE")
save(models, file = file.path(modelDir, "unfitted_models_CMIP.RData"))
##################################################################################################
# COMBINING AND SAVING MODELS (END)
##################################################################################################


##################################################################################################
# TESTING THAT MODELS FIT WITHOUT ERRORS (START)
##################################################################################################
system.time(
for(i in 1:length(models)){
  print(i)
  sampleMcmc(models[[i]],samples=2)
}
)
##################################################################################################
# TESTING THAT MODELS FIT WITHOUT ERRORS (END)
##################################################################################################
