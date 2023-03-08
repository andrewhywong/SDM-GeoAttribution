#### quantile of pred map
# install.packages('rgdal')
# allpacks <- installed.packages()[,1]
# lapply(allpacks, require, character.only = TRUE)
library(rgdal)
library(raster)
library(stringr)
library(stringi)
library(data.table)
library(doParallel)
# rm(list = ls())

###################################################################################
###################################################################################
###################################################################################
## set your working directory to the root folder first
setwd('C:/Users/adwhy/Box/NYMPHS_2023/')

## get target bees 
beeSelected = read.csv('./data_beeCSV/bee_selected/selected_Bees_stack_vall.csv',na.strings=c("","NA"))
beeSelected[beeSelected$Dataset=='TX2',]$Selected.Bees <- str_sub(
  beeSelected[beeSelected$Dataset=='TX2',]$Selected.Bees,2,-1)
beeSelected$Selected.Bees <- paste0(beeSelected$Dataset,'_',beeSelected$Selected.Bees)
any(duplicated(beeSelected$Selected.Bees))

########################################################
# if need to reduce bees to those we focused: in current study
reducedBees <- read.csv("./data_beeCSV/bee_selected/selectedBees_LAAN81.csv")
reducedBees <- as.vector(reducedBees$x)
beeSelected <- beeSelected[beeSelected$Selected.Bees %in% reducedBees,]
########################################################

## specify current mode
objective = 'full'
absenseMode = 'background' #'background' 'pseudoAbsence'

## set working directory to the target model folder for biomod2 export
setwd(paste0('./SPmodels_',objective,'_',absenseMode,'/model_outputs/'))


## available models
modelslist <- list.dirs('.', full.names = FALSE, recursive = FALSE)
modelslist


# ## get selected bees + their species list, and xy loc 
# beeSelected <- read.csv('C:/Users/adwhy/OneDrive - The University of Texas at Austin/nymphs_TGIS/data_beeCSV/bee_selected/selected_Bees_stack_v2.csv',
#                        na.strings=c("","NA"))
# # update beeID concatenate
# # remove the X before bee name of TX2
# beeSelected[beeSelected$Dataset=='TX2',]$Selected.Bees <- str_sub(beeSelected[beeSelected$Dataset=='TX2',]$Selected.Bees,2,-1)
# beeSelected$Selected.Bees <- paste0(beeSelected$Dataset,'_',beeSelected$Selected.Bees)
# # MAKE SURE BEE NAME NOT DUPLICATED
# any(duplicated(beeSelected$Selected.Bees))

# occdata source from a occ pool:
occdataLoc = 'C:/Users/adwhy/OneDrive - The University of Texas at Austin/nymphs_TGIS/data_beeCSV/occdata_rmFull_GEE_CLEANED/' # refined occ data, use a random n = 100
occimport = list.files(occdataLoc,"*.csv",full.names = T)
occSymbol = list.files(occdataLoc,"*.csv",full.names = F)
symbolList = sub('\\.csv$', '', occSymbol)
symbolList

### diff between eval and maps?
setdiff(
        str_sub(list.files('C:/Users/adwhy/OneDrive - The University of Texas at Austin/nymphs_TGIS/nymphs_eval/stack_outputs/biomod2',"*.csv",full.names = F),1,-5),
        str_sub(list.files('C:/Users/adwhy/OneDrive - The University of Texas at Austin/nymphs_TGIS/nymphs_maps/biomod2',"*.tif",full.names = F),1,-14) )

## rescale rasters
rescale0to1 <- function(rasterForCalculation){
  if (class(rasterForCalculation) != "RasterLayer"){
    warning("Supplied argument is not a raster./n", sep = "")
    return(NULL)
  }
  rescaledRaster <- (rasterForCalculation - rasterForCalculation@data@min)/(rasterForCalculation@data@max - rasterForCalculation@data@min)
  return(rescaledRaster)
}

# # setdiff(tifs,list.files('C:/Users/adwhy/OneDrive - The University of Texas at Austin/nymphs_TGIS/nymphs_maps/biomod2',"*.tif",full.names = F))
# 
### beeXspecies list
beeXspList <- list(read.csv('C:/Users/adwhy/OneDrive - The University of Texas at Austin/nymphs_TGIS/data_beeCSV/CA_NATURAL.csv'),
               read.csv('C:/Users/adwhy/OneDrive - The University of Texas at Austin/nymphs_TGIS/data_beeCSV/CA_URBAN.csv'),
               read.csv('C:/Users/adwhy/OneDrive - The University of Texas at Austin/nymphs_TGIS/data_beeCSV/TX_NATURAL.csv'),
               read.csv('C:/Users/adwhy/OneDrive - The University of Texas at Austin/nymphs_TGIS/data_beeCSV/TX_NATURAL_2.csv'))
# eliminate duplicated bees (input raw file error)
# note that duplicate only for numeric
beeXspList_Xdupli <- list()
n <- 0
for (each in beeXspList){
  n <- n+1
  if (is.character(each[,2])){
    each <- each[!stri_duplicated(as.character(each$BeeID)),]
  }else{
    each <- each[!duplicated(each$BeeID),]
  }
  beeXspList_Xdupli[[n]] <- each
}
beeXspList_Xdupli[[1]]$BeeID <- paste0('CAN_',beeXspList_Xdupli[[1]]$BeeID)
beeXspList_Xdupli[[2]]$BeeID <- paste0('CAU_',beeXspList_Xdupli[[2]]$BeeID)
beeXspList_Xdupli[[3]]$BeeID <- paste0('TX1_',beeXspList_Xdupli[[3]]$BeeID)
beeXspList_Xdupli[[4]]$BeeID <- paste0('TX2_',beeXspList_Xdupli[[4]]$BeeID)



#### 1. 1 bee situation
targetVersion <- 'bee_wise_SS' # version of the joint suitability method

## UPDATE FINDXY before start loop
findXY <- read.csv('../../data_beeCSV/beeXY.csv')
findXY$BeeID <- paste0(findXY$State,'_',findXY$BeeID)
findXY <- findXY[!duplicated(findXY[,1]),]

for(i in modelslist[2:5]){ # RUN NOTES: 
  # i='biomod2'
  mapDir <- paste0('./',i,'/',targetVersion)
  if(!file.exists(mapDir)){
    dir.create(file.path(mapDir))
    print("The directory is created")
  }
  print(i)
  
  findXY$SSprob <- NA
  findXY$UpperBee <- NA
  findXY$TotalCells <- NA
  findXY$nSp <- NA
  
  # for each bee:
  for (j in 1:nrow(beeSelected)){ # issue bee: 
    # j=89
    beeName <- beeSelected$Selected.Bees[j]
    print(j)
    thisBee <- as.character(as.vector(beeSelected[j,4:ncol(beeSelected)]))
    thisBee <- thisBee[!is.na(thisBee)]
    
    # # Weighting of environmental space: 
    # if ()

    # if don't do weight:
    tifs <- paste0(thisBee,'_pred_map.tif')
    mapList <- stack(paste0('./',i,'/',tifs))
    
    # scale sum
    scaledSumTIF <- calc(mapList, sum)/length(thisBee)
    
    # rescale 0-1
    scaledSumTIF <- rescale0to1(scaledSumTIF)
    
    # BEE LOC PROB
    XY <- cbind(findXY[findXY$BeeID==beeName,]$Lon,findXY[findXY$BeeID==beeName,]$Lat)
    bee.mean = extract(scaledSumTIF, XY)
    
    # write out this ss raster and prob
    writeRaster(scaledSumTIF,paste0('./',i,'/',targetVersion,'/',beeName,'.tif'),overwrite = T)
    
    # record the prob at location, and the cell No. > or equal to bee.mean
    
    findXY[findXY$BeeID==beeName,]$SSprob = bee.mean
    findXY[findXY$BeeID==beeName,]$UpperBee = length(scaledSumTIF[scaledSumTIF>=bee.mean])
    findXY[findXY$BeeID==beeName,]$nSp = nlayers(mapList)
    
    }
  ## Total non-na pixel
  findXY$TotalCells = length(scaledSumTIF[!is.na(scaledSumTIF)])
  findXY <- findXY[!rowSums(is.na(findXY)),]
  write.csv(findXY,paste0('./',i,'/',targetVersion,'/bees_',i,'.csv'))
}



