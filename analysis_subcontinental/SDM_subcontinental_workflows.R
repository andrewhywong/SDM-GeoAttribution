###############################################################################################################
################################################  Code Part 1  ################################################
###############################################################################################################

# # To run JAVA-related model need to download x64bit java and install package 'rJava':
# install.packages('rJava')
library(rJava)
####

#### import required libraries
library(biomod2)
library(prg)
library(randomForest)
library(scales)
library(ggplot2)
library(doParallel) # parallel processing
library(stringr) # dataframe processing

#### raster dependency
library(raster)
library(sf)
library(rgeos)
library(sp)

##############################################################################################################
## set your working directory to the parent folder
setwd('/Box/')
## Select models to test:
modelList = c('biomod2','rfdown','rf','brt','maxent')
##############################################################################################################
## a subcontinental study area extent
studyArea="POLYGON((-90 43, -125 43,-125 25, -90 25 , -90 43))"

## read-in bioclimate data, note different versions
clim_alt <- stack('./env_variables/clim_alt_TXCA_decorr.tif')

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
beeSelected = read.csv('./data_beeCSV/bee_selected/selected_Bees_stack_vall.csv',na.strings=c("","NA"))
beeSelected[beeSelected$Dataset=='TX2',]$Selected.Bees <- str_sub(
  beeSelected[beeSelected$Dataset=='TX2',]$Selected.Bees,2,-1)
beeSelected$Selected.Bees <- paste0(beeSelected$Dataset,'_',beeSelected$Selected.Bees)

## if need to reduce bees to those we focused: in current study
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
funcDir = './functions_tune'
source(paste0(funcDir,'/Fun_TuneMaxent.R'))
source(paste0(funcDir,'/Fun_trainTest_gcloud.R'))
source(paste0(funcDir,'/Fun_allModel.R'))
###################################################################################

###################################################################################
## 1) specify objective, either 'evaluate': 80%/20% train/test partitioning, 
## or 'full': use all data to generate suitability map;
## 2) specify absenseMode to either 'background' for large sampling
## or 'pseudoAbsence' using k-means clustering pseudo-absence

## check working directory before running
objective = 'full' # full or evaluate
absenseMode = 'background' #'background' or 'pseudoAbsence'

## import presampled 50000 when absense mode is background
bgdf <- read.csv('./bg_CATX_50000.csv')


## read in species scientific names to get unique seed for each running:
fullSci = read.csv('./All_uni_symbols_SciName.csv')

## species that has less than 30 occ will likely cause model errors
lessthan30occsp <- c() 
errorSpecies <- c() # store error species

#### config parallel running model-wise:
n.cores <- detectCores()
cl <- makeCluster(n.cores)
registerDoParallel(cl)
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
stopCluster(cl)
errorSpecies
edTime = Sys.time()
edTime - stTime

###############################################################################################################
################################################  Code Part 2  ################################################
###############################################################################################################

targetVersion <- 'bee_wise_SS' # version of the joint suitability method

## UPDATE FINDXY before start loop
findXY <- read.csv('../../data_beeCSV/beeXY.csv')
findXY$BeeID <- paste0(findXY$State,'_',findXY$BeeID)
findXY <- findXY[!duplicated(findXY[,1]),]

for(i in modelslist[2:5]){ 
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
    
    # suitability at the bee sampling location
    XY <- cbind(findXY[findXY$BeeID==beeName,]$Lon,findXY[findXY$BeeID==beeName,]$Lat)
    bee.mean = extract(scaledSumTIF, XY)
    
    # write out this ss raster and prob
    writeRaster(scaledSumTIF,paste0('./',i,'/',targetVersion,'/',beeName,'.tif'),overwrite = T)
    
    # record the suitability at location, and the cell No. > or equal to bee.mean
    findXY[findXY$BeeID==beeName,]$SSprob = bee.mean
    findXY[findXY$BeeID==beeName,]$UpperBee = length(scaledSumTIF[scaledSumTIF>=bee.mean])
    findXY[findXY$BeeID==beeName,]$nSp = nlayers(mapList)
    
  }
  ## Total non-na pixel
  findXY$TotalCells = length(scaledSumTIF[!is.na(scaledSumTIF)])
  findXY <- findXY[!rowSums(is.na(findXY)),]
  write.csv(findXY,paste0('./',i,'/',targetVersion,'/bees_',i,'.csv'))
}


###############################################################################################################
################################################  Code Part 3  ################################################
###############################################################################################################

##### 3.1 Pixel-wise rate of identfying bees #####
thres_list <- c(0.9700,0.9900,0.9950,0.9990,0.9995,0.9999) # change to have different %tiles
source('../../GA_code/functions_post/pixel_HitRate.R')

# used time <30min
cl <- registerDoParallel(makeCluster(detectCores()))
thresStack <- foreach(i = modelslist, .packages = c('raster')
                      ) %dopar% {stacked <- thresfun(i,targetVersion,thres_list)}
stopCluster(cl)

##########################################################################################################
## 2. Polygon version: normal buffer
buf_size <- 0.05 # 0.05 degrees ~= 5.5km 
thres_list <- c(0.9700,0.9900,0.9950,0.9990,0.9995,0.9999) # change to have different %tiles

chullsFolder <- '../nymphs_beeplots_nbuffer/' # note parent folder is at upper loc
targetVersion <- 'bee_wise_SS'
shp_outpath <- 'shp_n_5.5km' # to make diff output file names
csvTail <- '_n_5km'

modelslist

stTime = Sys.time()
for(i in modelslist[2:5]){ # run notes: run 3:4 during worldcup
  # i = 'biomod2'
  print(i)
  maps <- list.files(paste0('./',i,'/',targetVersion,'/'),"*.tif$",full.names = T)
  percs <- read.csv(paste0('./',i,'/',targetVersion,'/','bees_',i,'.csv'))
  ori_ncol <- ncol(percs)
  
  ## IF CODE INTERUPTED: 
  # percs <- read.csv('./',i,'/',targetVersion,'/','bees_',i,'_gunary_Polygon_a01','.csv'))
  # percs <- percs[,2:ncol(percs)] # row.name =F
  
  ## prepare new cols
  # poly hitrate for cols
  percs[,paste0("hit",thres_list*1e4)] <- NA
  # area of total polygon AOIs
  percs[,paste0("hit",thres_list*1e4,'_tarea')] <- NA
  # distance to the Closest polygon, if point X fall in
  percs[,paste0("hit",thres_list*1e4,'_dtoc')] <- NA
  # distance to the maximum suitability pixel
  percs$dtoMax <- NA 
  
  for (tmap in maps){
    # visBack <- match('CAN_SJGP_9',str_sub(maps,49,-5)) # to visit back and plot
    # tmap = maps[visBack]
    # tmap = maps[1]
    
    print(match(tmap,maps))
    tmap = raster(tmap)
    
    # bee loc
    XY <- cbind(percs[percs$BeeID==names(tmap),]$Lon,percs[percs$BeeID==names(tmap),]$Lat)
    XY_sf <- st_as_sf(SpatialPointsDataFrame(coordinates(XY),as.data.frame(XY),
                                             proj4string = CRS("+init=epsg:4326")))
    
    for (thres in thres_list){
      # thres = thres_list[1]
      print(thres)
      # get upper %tile cells only 
      tmapTest = tmap
      tmapTest[tmapTest[] < quantile(tmapTest,probs=thres,na.rm=T)[[1]]] = NA # get top xx %tile points
      tmapTest[!is.na(tmapTest[])] = 1
      
      tmapPoints = rasterToPoints(tmapTest)[,1:2]
      tmapPoints = tmapPoints[!(duplicated(tmapPoints[,1:2]) | duplicated(tmapPoints[,1:2], fromLast = TRUE)), ]
      
      p <- st_buffer(st_point(cbind(0, 0)), buf_size)
      l <- replicate(nrow(tmapPoints), p, simplify = FALSE)
      for (tp in seq_len(nrow(tmapPoints))) {
        l[[tp]] <- l[[tp]] + tmapPoints[tp, ]
      }
      chulls <- st_union(st_sfc(l))
      st_crs(chulls) <- CRS("+init=epsg:4326")
      
      if (!st_is_valid(chulls)){
        print('fixed chulls geometry')
        chulls <- st_make_valid(chulls)
      }
      # plot(chulls)
      
      ##############################################################################################################
      st_write(chulls,
               paste0(chullsFolder,shp_outpath,'/',i,'.',
                      names(tmap),'.',thres,'.','Xholes','.shp'),
               quiet=T,delete_dsn=T)
      ##############################################################################################################
      
      # check xy coord within 1km distance to AOI
      IsWithin <- st_within(XY_sf, chulls,sparse=F)
      shortestDist <- st_distance(XY_sf, chulls) # note this func returns 0 when inside
      
      # hit or not
      percs[percs$BeeID==names(tmap),][,ori_ncol+match(thres,thres_list)] <- ifelse(IsWithin[,1],1,0)
      # area of total AOIs
      percs[percs$BeeID==names(tmap),][,ori_ncol+length(thres_list)+match(thres,thres_list)] <- 
        as.numeric(st_area(chulls)) #m2; km2 in 4_bee_6 code
      # dist to nearest polygon/point; note on change of condition of ifelse()
      percs[percs$BeeID==names(tmap),][,ori_ncol+length(thres_list)*2+match(thres,thres_list)] <- 
        as.numeric(shortestDist)
      
    }
    
    # dist to maximum prob pixel; avoid redundant calc
    percs[percs$BeeID==names(tmap),][,ncol(percs)] <- 
      min(as.numeric(
        st_distance(XY_sf,
                    st_as_sf(as.data.frame(xyFromCell(tmap,which.max(tmap))),
                             coords =c('x','y'),crs=CRS("+init=epsg:4326"))) 
      ))
    
  }
  ##############################################################################################################
  # write added cols version out to same directory, with a new name
  write.csv(percs,paste0('./',i,'/',targetVersion,'/','bees_',i,'_roundBuffer_',thres,csvTail,'.csv'))
  # write.csv(broken_list,paste0(getwd(),'/',i,'/',targetVersion,'/','bees_',i,'_roundBuffer_',thres,csvTail,'.csv'))
  ##############################################################################################################
}

edTime = Sys.time()
edTime - stTime

###############################################################################################################
################################################  Code Part 4  ################################################
###############################################################################################################

## Raw files used in Part 4 were generated in previous code parts (CSVs)
for(i in modelslist){
  # i = 'biomod2'
  if (i=='biomod2'){
    modelBees <- read.csv(paste0('./',i,'/',targetVersion,'/bees_',i,'.csv')) 
    modelBees$model <- i
  }else{
    t <- read.csv(paste0('./',i,'/',targetVersion,'/bees_',i,'.csv'))
    t$model <- i
    modelBees <- rbind(modelBees,t)
  }
}

modelBees <- modelBees[modelBees$BeeID %in% reducedBees,]

## add a col for % upper cells/total cells
modelBees$seachPercent <- 100-(modelBees$UpperBee/modelBees$TotalCells)*100

## specify plotting 
level_order  <- factor(modelBees$model, level = c(modelslist))
modelBees$model_level_order <- level_order
nsp_order <- factor(modelBees$nSp, levels = seq(1,20)[-11]) 
orderbynSp = unique(modelBees[order(modelBees$nSp),]$BeeID)
beeID_ordered  <- factor(modelBees$BeeID,level = orderbynSp) # order by NSP

t <- ggplot(modelBees, aes(x=nsp_order, y=seachPercent))+
  # geom_point(alpha = 0.5, size = 2) +
  geom_boxplot(outlier.colour="red", outlier.shape=4,outlier.size=1) +
  facet_grid(~model_level_order, scales="free_x")  +
  geom_hline(yintercept=99.9, linetype="dashed", color = "red") +
  # geom_text(aes(2,99.9,label = 99.9, vjust = -1))+
  geom_hline(yintercept=99.5, linetype="dashed", color = "red") +
  # geom_errorbar(aes(x=level_order,ymin = mean-sd, ymax = mean+sd), width = 0.2)+
  xlab("No. of species identified on an object") +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab(paste0("SSR Score")) +
  ylim(90,100)

# when we specify x-axis variable inside the aesthetics function aes(). 
# reorder() function sorts the carriers by mean values by default
modelBees$GeographicArea #,fill=GeographicArea
t <- ggplot(modelBees, aes(x=beeID_ordered, y=seachPercent)) +  # bee genus level ploting is intesesting too
  # geom_point(alpha = 0.5, size = 2) +
  geom_boxplot(outlier.colour="red", outlier.shape=4,outlier.size=1) +
  facet_grid(.~State, scales="free_x") + # "free_x"
  geom_hline(yintercept=99.9, linetype="dashed", color = "red") +
  # geom_text(aes(2,99.9,label = 99.9, vjust = -1))+
  geom_hline(yintercept=99.5, linetype="dashed", color = "red") +
  # geom_errorbar(aes(x=level_order,ymin = mean-sd, ymax = mean+sd), width = 0.2)+
  xlab("Bees") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab(paste0("Search Score")) 

################################################################################################################################################
#### Freq table of %tiles captured, POLYGON-based
targetVersion <- 'bee_wise_SS' 
##### change csv_tail!!
for(i in modelslist){
  # i = 'biomod2'
  if (i=='biomod2'){
    modelBees <- read.csv(
      paste0('./',i,'/',targetVersion,'/bees_',i,'_roundBuffer_0.9999_n_5km.csv')
    )
    modelBees$model <- i
    # mTiles <- modelBees[,c(10,19:24)] %>% group_by(State) %>% summarize_all(sum)
    # mTiles <- as.data.frame(t(colSums(modelBees[19:24])))
    # mTiles$model <- i
  }else{
    t <- read.csv(
      paste0('./',i,'/',targetVersion,'/bees_',i,'_roundBuffer_0.9999_n_5km.csv')
    )
    t$model <- i
    modelBees <- rbind(modelBees,t)
  }
}

# melt df to plot hit
mTiles <- modelBees[,c(10,15:ncol(modelBees))]
# check rows at hit9900 to hit9999, they should be binary
nrow(mTiles[!rowSums(is.na(mTiles[,c(2:7)])), ])
nrow(mTiles[mTiles$State=='CAN',])
colnames(mTiles)

# change unit to km2
# change area unit to km2, and distance unit to km
mTiles[,c(8:13)] <- mTiles[,c(8:13)] / 10^6
mTiles[,c(14:(ncol(mTiles)-1))] <- mTiles[,c(14:(ncol(mTiles)-1))] / 10^3

# if use ..count.. in ggplot:
mTiles_a0 <- mTiles[,c(1:7,ncol(mTiles))] # area plotting X
colnames(mTiles_a0)[2:7] <- c('97.00%', '99.00%' ,'99.50%' ,'99.90%', '99.95%', '99.99%')

# store to polygon-hitrate
polyHitrate <- mTiles_a0 
polyHitrate$base <- '% captured by buffer AOI'

################################################################################################
# # 1. area of that polygon if point fall in
# mTiles_a1 <- mTiles[,c(1,7:11,ncol(mTiles)-1)]
# colnames(mTiles_a1)[2:6] <- c('99.00%' ,'99.50%' ,'99.90%', '99.95%', '99.99%')

##### 2.area of total polygon AOIs
mTiles_a2 <- mTiles[,c(1,8:13,ncol(mTiles))]
# mTiles_a2_avg <- mTiles_a2 %>% group_by(State,model) %>% summarize_all(mean)
colnames(mTiles_a2)[2:7] <- c('97.00%' ,'99.00%' ,'99.50%' ,'99.90%', '99.95%', '99.99%')

##### 3.Distance to the Closest polygon, if point X fall in
##### update dtoc cols condition on hits cols, since I forgot to add in alpha algorithm
distmTiles <- mTiles
for (colNo in 14:19){ # change all dist cols
  distmTiles[,colNo] <- ifelse(distmTiles[,colNo-12]==0,distmTiles[,colNo],NA)
}
mTiles_d1 <- distmTiles[,c(1,14:19,ncol(distmTiles))]
colnames(mTiles_d1)[2:7] <- c('97.00%','99.00%' ,'99.50%' ,'99.90%', '99.95%', '99.99%')

# # 4. Average Distance to all AOI polygons (even fall in)
# mTiles_d2 <- mTiles[,c(1,22:ncol(mTiles)-1)]
# colnames(mTiles_d2)[2:6] <- c('99.00%' ,'99.50%' ,'99.90%', '99.95%', '99.99%')

# # 5. distance to highest pixel
# mTiles_dmax <- mTiles[,c(1,ncol(mTiles)-1,ncol(mTiles)-2)]
# mTiles_a2_avg <- mTiles_dmax %>% group_by(State,model) %>% summarize_all(mean)
# colnames(mTiles_dmax)[3] <- c('Distance to Pixel')
# colnames(mTiles_a2_avg)[3] <- c('Mean Distance to Pixel')

dfToMelt <- mTiles_a0
mTiles_melted <- melt(dfToMelt, id.vars=c('State','model'))

################################################################################################
# to cont'd
# # penalized prediction error
# ggplot(mTiles_melted[mTiles_melted$model=='biomod2',], aes(x = value)) + 
#   facet_grid(State~variable)+  # , scales="free_x" State~variable for all except freq
#   geom_histogram(aes(y = stat(density)*5))+  #, bins = 8 / sum(count)
#   scale_y_continuous(labels = scales::percent)
m1 <- mTiles_melted[mTiles_melted$model=='rfdown' &
                      mTiles_melted$variable=='99.99%',] #mTiles_melted$model=='rfdown'& 

grp_plots <- by(m1, m1$State, function(sub){
  print(sub$State[1])
  totalwithNA <- length(sub$State)
  print(totalwithNA)
  sub <- sub[!rowSums(is.na(sub)),] # na=captured; rm NA to plot ONLY THE UNCAPTURED 
  print(length(sub$State))
  print(max(sub$value))
  p <- ggplot(sub, aes(x=value)) + 
    geom_histogram(aes(y =  ..count..),breaks=seq(0,600,50)) + # #..count../totalwithNA facet_wrap wont give sum count
    # geom_text(stat="bin", vjust = 1, size = 1.5, colour = "white",breaks = seq(0,800,30))  + 
    # scale_y_continuous(labels = scales::percent,limits=c(0,0.5)) +
    scale_y_continuous(limits=c(0,40)) +
    # facet_grid(State~variable) + 
    # ggtitle(sub$State[[1]]) +
    # theme(plot.title = element_text(hjust = 0.5))+
    # xlab('Prediction Error (km)') + 
    # ylab(paste0('Missing Count/Total Count')) +
    ylab("")+
    xlab("")+
    theme(aspect.ratio=1) + # make fixed 1:1
    # theme(axis.text.x = element_blank(),axis.text.y = element_blank())+ #
    theme(plot.margin=unit(c(0.01,0.01,0.01,0.01),"cm"))
  # {if(sub$State[1]=='TX2')theme(axis.text.y = element_blank())
  #   else theme(axis.text.x = element_blank(),axis.text.y = element_blank())}
  
  # if (!sub$State[1]=='TX2'){
  #   p + theme(axis.text.x = element_blank(),axis.text.y = element_blank())
  # }else{
  #   p + theme(axis.text.y = element_blank())
  # }
})
grid.arrange(grobs = grp_plots, nrow=4)

png(filename = paste0(concatwd,'XhitfreqDistBand_9999_rfdown',".png"),# HittedArea_9900_9950_AOI
    width = 5000, height = 5000,res = 600) #10000 4000 600 # 3000 2000 300
grid.arrange(grobs = grp_plots, nrow=4)
dev.off()

# cont'd
mTiles_melted <- mTiles_melted[!rowSums(is.na(mTiles_melted)),]
unique(mTiles_melted$variable)
submtiles <- mTiles_melted[mTiles_melted$variable %in% 
                             c('97.00%','99.00%','99.50%' ,'99.90%', '99.95%', '99.99%'),] #'97.00%','99.00%','99.50%' ,'99.90%', '99.95%', '99.99%'
level_order  <- factor(submtiles$model,
                       level = c("brt","maxent","rf","rfdown","Ensemble","biomod2"))

max(submtiles$value)

t <- ggplot(submtiles, 
            aes(x=level_order, y=value)) + 
  # geom_point(alpha = 0.5, size = 2) +
  # geom_boxplot(outlier.colour="red", outlier.shape=4,outlier.size=1) +
  geom_bar(stat="identity") +
  facet_grid(State~variable) + # , scales="free_x" State~variable for all except freq
  xlab("Models") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1))+
  ylab('No. of times one of AOIs captured the target') +
  # 'Averaged Distance (km) from an object to all AOIs'
  # bquote("Area ("~km^2~") of the AOI where the object is identified")
  # 'Distance (km) from an object (if not captured) to its nearest AOI'
  # No. of times one of AOIs captured the target; 
  # coord_flip() +
  # geom_text(aes(label = round(value,1)), vjust = 1, size = 1.5, colour = "white") +
  # ylim(0,600) +
  theme(aspect.ratio=1) # make fixed 1:1


##############################################################################################################################
#### Freq table of %tiles captured, pixel-based
for(i in modelslist){
  # i = 'biomod2'
  if (i=='biomod2'){
    modelBees <- read.csv(paste0('./',i,'/',targetVersion,'/bees_',i,'_recursPercTiles_Pixel.csv'))
    modelBees$model <- i
    # mTiles <- modelBees[,c(10,19:24)] %>% group_by(State) %>% summarize_all(sum)
    # mTiles <- as.data.frame(t(colSums(modelBees[19:24])))
    # mTiles$model <- i
  }else{
    t <- read.csv(paste0('./',i,'/',targetVersion,'/bees_',i,'_recursPercTiles_Pixel.csv'))
    t$model <- i
    modelBees <- rbind(modelBees,t)
  }
}

# same process as the AOI-based method
mTiles <- modelBees[,c(10,15:ncol(modelBees))]
# check rows at hit9900 to hit9999, they should be binary
nrow(mTiles[!rowSums(is.na(mTiles[,c(2:7)])), ])
nrow(mTiles[mTiles$State=='CAN',])
colnames(mTiles)

# if use ..count.. in ggplot:
mTiles_a0 <- mTiles[,c(1:7,ncol(mTiles))] # area plotting X
colnames(mTiles_a0)[2:7] <- c('97.00%', '99.00%' ,'99.50%' ,'99.90%', '99.95%', '99.99%')

pixelHitrate <- mTiles_a0
pixelHitrate$base <- '% captured by pixel'
###############################################--------------------------############################

mTiles <- melt(mTiles, id.vars=c('State','model'))
mTiles <- mTiles[mTiles$value==1,]
level_order  <- factor(mTiles$model, 
                       level = c("brt","maxent","rf","biomod2")) #,"rfdown","Ensemble"
t <- ggplot(mTiles, aes(x=level_order)) +  #, y=value
  # geom_point(alpha = 0.5, size = 2) +
  geom_bar() + #stat="identity"
  # The dot-dot notation (`..count..`) was deprecated in ggplot2 3.4.0.
  # Please use `after_stat(count)` instead.
  geom_text(aes(label = after_stat(count)),stat = "count", vjust = 1, size = 1.5, colour = "white") +
  facet_grid(State~variable) + # , scales="free_x"
  xlab("Models") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab(paste0("Percentage of captureing objects with an AOI, pixel-based")) +
  theme(aspect.ratio=1)  # make fixed 1:1 

#### stack polyhitrate and pixelhitrate and compare
hitrate <- rbind(polyHitrate,pixelHitrate)
hitrate <- melt(hitrate,id.vars=c('State','model','base'))
hitrate<- hitrate[hitrate$value==1,]

level_order  <- factor(hitrate$model,
                       level = c(modelslist))
phit <- ggplot(hitrate, aes(fill=base,x=level_order)) +  
  # geom_point(alpha = 0.5, size = 2) +
  geom_bar(position = position_dodge(preserve = "single")) + 
  facet_grid(State~variable) + # , scales="free_x"
  xlab("Models") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab('No. of successfully identified objects') +
  ylim(0,30)+
  geom_text(aes(label = ..count..),
            stat = "count", vjust = 1, size = 1.5, colour = "black",
            position = position_dodge(width = 0.5)) +
  theme(aspect.ratio=1)+ 
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1))

