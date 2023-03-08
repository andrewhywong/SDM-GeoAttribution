####
# # To run JAVA-related model need to download x64bit java and install package 'rJava':
# install.packages('rJava')
library(rJava)
####
# devtools::install_github("meeliskull/prg/R_package/prg")

library(biomod2)
library(prg)
library(precrec) # AUC_ROC; PR curve
library(ecospat) # Boyce index
library(dismo) # Maxent
library(randomForest) # RF model
library(scales)
library(ggplot2)
library(doParallel) # parallel processing
library(stringr) # dataframe processing
# library(disdat)
# library(PresenceAbsence)
# library(rgbif)

## raster dependency:
library(raster)
library(sf)
library(rgeos)
library(sp)

###################################################################################
###################################################################################
###################################################################################
## set your working directory to the parent folder
setwd('C:/Users/adwhy/Box/NYMPHS_2023/')
## Select models to use: 
modelList = c('biomod2','rfdown','rf','brt','maxent')
###################################################################################
## study area extent
studyArea="POLYGON((-90 43, -125 43,-125 25, -90 25 , -90 43))"

## read-in de-corred bioclimate data, note diff version
# clim_alt <- stack('./env_variables/clim_alt_TXCA_decorr.tif')
clim_alt <- stack('./env_variables/clim_alt_2.1_TXCA_decorr.tif')

## change name to meaningful; remember update for new CATX study area
names(clim_alt) <- c("clim_alt.2" , "clim_alt.3",  "clim_alt.7" , "clim_alt.8" , "clim_alt.9"  ,
                     "clim_alt.10" ,"clim_alt.14", "clim_alt.15" ,"clim_alt.16" ,"clim_alt.18" ,
                     "clim_alt.19", "SRTM")
# plot(clim_alt$SRTM)

# ## base mask and a func to remove duplicated entries, for spaital thinning:
# base_mask <- raster('./env_variables/all_unique_mask.tif')
# RemoveDups <- function(df, column) {
#   dups = duplicated(df[, column], fromLast = TRUE)
#   df   = df[!dups, ]
# }

## occurrences data input
occdataLoc = './data_beeCSV/occdata_rmFull_GEE_CLEANED/' 
occimport = list.files(occdataLoc,"*.csv",full.names = T)
occSymbol = list.files(occdataLoc,"*.csv",full.names = F)
symbolList = sub('\\.csv$', '', occSymbol) 

## make sure the imported are the same with symbol name:
# str_sub(sapply(str_split(occimport, "/"),`[[`,4)
identical(str_sub(occSymbol,1,-5),symbolList)


## get target bees 
beeSelected = read.csv('./data_beeCSV/bee_selected/selected_Bees_stack_vall.csv',
                       na.strings=c("","NA"))
beeSelected[beeSelected$Dataset=='TX2',]$Selected.Bees <- str_sub(
  beeSelected[beeSelected$Dataset=='TX2',]$Selected.Bees,2,-1)
beeSelected$Selected.Bees <- paste0(beeSelected$Dataset,'_',beeSelected$Selected.Bees)

# if need to reduce bees to those we focused: in current study
reducedBees <- read.csv("./data_beeCSV/bee_selected/selectedBees_LAAN81.csv")
reducedBees <- as.vector(reducedBees$x)
beeSelected <- beeSelected[beeSelected$Selected.Bees %in% reducedBees,]

## double check bee duplicates
length(unique(beeSelected$Selected.Bees)) == length(beeSelected$Selected.Bees)
beeSelected$Selected.Bees[which(duplicated(beeSelected$Selected.Bees))]

## get all unique species on bees
beeVec <- as.vector(matrix(t(beeSelected[,4:ncol(beeSelected)])))
beeVec <- unique(na.omit(beeVec))
any(!beeVec %in% symbolList)
setdiff(beeVec,symbolList)

## species that are already run need to be excluded; otherwise:
# truebeeVec = 'DAPU3'
truebeeVec = beeVec
truebeeVec

###################################################################################
## call in functions
funcDir = './GA_code/functions_tune'
source(paste0(funcDir,'/Fun_TuneMaxent.R'))
source(paste0(funcDir,'/Fun_trainTest_gcloud.R'))
source(paste0(funcDir,'/Fun_allModel.R'))
###################################################################################

###################################################################################
## 1) specify objective, either 'evaluate': 80%/20% train/test partitioning
## or 'full': use all data to generate suitability map;
## 2) specify absenseMode to either 'background' for large sampling
## or 'pseudoAbsence' using k-means clustering pseudo absence
## check working directory before running
objective = 'full' # full or evaluate
absenseMode = 'background' #'background' 'pseudoAbsence'
addon = 'wc2'
## import presampled 50000 when absense mode is background
bgdf <- read.csv('./env_variables/bg_CATX_50000.csv')


## read in species scientific names to get unique seed for each running:
fullSci = read.csv('./data_beeCSV/All_uni_symbols_SciName.csv')

## species that has less than 30 occ triger model errors
lessthan30occsp <- c() #train/test failed: c('LAAN81', 'OXARR','BULBI'); full failed: 'LAAN81'
errorSpecies <- c() # store error species

#### config parallel running model-wise:
n.cores <- detectCores()
cl2 <- makeCluster(8)
registerDoParallel(cl2)
stTime = Sys.time()

## set working directory to the target model folder for biomod2 export
setwd(paste0('./SPmodels_',objective,'_',absenseMode,'_',addon))

## check again
objective
absenseMode
getwd()

# remove temp raster at each loop
tmp_dir <- paste0(tempdir(),'/raster')
plot(clim_alt$clim_alt.2)
dev.off() # clear png device before running

for (s in truebeeVec[53:length(truebeeVec)]){  # run 37ROOF  
  # match('ROOF',truebeeVec)
  # s = truebeeVec[53]
  n <- match(s,fullSci$Symbol) # seeds
  print(paste0('----------Modeling species ',n,s,' to ',objective,'----------'))
  
  # import formatted occ data
  occdata = read.csv(paste0('.',occimport)[which(symbolList==s)]) 
  occdata$filter_uniq <- NULL
  
  # # if too few data, go to next loop;
  if (objective=='evaluate' & nrow(occdata[occdata$occ==1,]) < 30){
    print(paste0(s," have less than 30 occ"))
    lessthan30occsp <- append(lessthan30occsp,s)
    next
  }

  # call func to create a list that stores needed files
  paraList <- traintestPrep(objective,occdata,bgdf,clim_alt)
  
  # parallel can be applied to species-wise to reach maximum speed, but not stable
  all.runs <- foreach(i = 1:length(modelList), 
                      .packages = c("biomod2","prg","disdat",
                                    "precrec","ggplot2","ecospat",
                                    "dismo","randomForest")) %dopar% {
                                      stacked <- train_mod(
                                        paraList$traindf,paraList$Xnorm_traindf,
                                        paraList$testdf,paraList$Xnorm_testdf,
                                        paraList$traindf_xy,paraList$normed_stack,
                                        modelList[[i]],s,objective,
                                        errorSpecies
                                      )
                                    }
  
  ## output model results:
  for (m in 1:length(all.runs)){ # length(all.runs)
    if (objective=='evaluate'){
      # export predictions info
      write.csv(all.runs[[m]][[2]],
                paste0('./model_outputs/',
                       all.runs[[m]][[1]],'/', s, '.csv'),
                row.names = FALSE)
    } else if (objective=='full'){
      # export predicted probability maps
      writeRaster(all.runs[[m]][[2]],
                  paste0('./model_outputs/',
                         all.runs[[m]][[1]],'/', s, '_pred_map.tif'),
                  overwrite=T)
    }
  }
  
  # clear up temp files to release space
  files <- list.files(tmp_dir, full.names = T)
  file.remove(files)
  
}
stopCluster(cl2)
errorSpecies
edTime = Sys.time()
edTime - stTime

############  Code #1 ends  ############  
