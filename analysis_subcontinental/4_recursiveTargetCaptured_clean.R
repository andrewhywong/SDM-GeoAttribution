#### quantile of pred map
# install.packages('alphahull')
library(doParallel)
library(rgdal)
library(raster)
library(ggplot2)
library(stringr)
library(geosphere)
library(concaveman)
library(alphahull)
library(rgeos)
library(sf)
library(maptools)
library(exactextractr)
# rm(list = ls())

##########################################################################################################
##########################################################################################################
##########################################################################################################

## set your working directory to the root folder first
setwd('C:/Users/adwhy/Box/NYMPHS_2023/')
targetVersion <- 'bee_wise_SS' # version of the joint suitability method


## specify current mode
objective = 'full'
absenseMode = 'background' #'background' 'pseudoAbsence'

## set working directory to the target model folder for biomod2 export
setwd(paste0('./SPmodels_',objective,'_',absenseMode,'/model_outputs/'))

## available models
modelslist <- list.dirs('.', full.names = FALSE, recursive = FALSE)
modelslist

# ####
# A loop to assess how many times we fall into 99%, then onto 99.5%, 99.95%, 99.99%... - out of
# total number of times we at least hit 95% (depending on how well models perform generally);
# The output will be a a frequency table with rows of bees and cols of %tile
# Also, an average distance from target point to %tile AOI can be calculated

##### 1. Pixel-wise version. #####
thres_list <- c(0.9700,0.9900,0.9950,0.9990,0.9995,0.9999) # change to have different %tiles, according to area size of the site
source('../../GA_code/functions_post/pixel_HitRate.R')

# used time <30min
n.cores <- detectCores() 
cl2 <- makeCluster(n.cores)
registerDoParallel(cl2)
thresStack <- foreach(i = modelslist,
                    .packages = c('raster')) %dopar% {stacked <- thresfun(i,targetVersion,thres_list)}
stopCluster(cl2)

##########################################################################################################
##########################################################################################################
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
  # percs <- read.csv(paste0(getwd(),'/',i,'/',targetVersion,'/','bees_',i,'_gunary_Polygon_a01','.csv'))
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


