# pixel-based hit rate function
thresfun <- function(mod,targetVersion,thres_list){
  # mod='biomod2'
  print(mod)
  maps <- list.files(paste0('./',mod,'/',targetVersion,'/'),"*.tif$",full.names = T)
  percs <- read.csv(paste0('./',mod,'/',targetVersion,'/','bees_',mod,'.csv'))
  ori_ncol <- ncol(percs)
  
  # prepare thresholds to columns
  percs[,paste0("hit",thres_list*1e4)] <- NA
  
  for (tmap in maps){
    # tmap = maps[1]
    print(match(tmap,maps))
    
    tmap = raster(tmap)
    
    # try raster to point + convex/cave hull directly since the contour line is not accurate enough
    for (thres in thres_list){
      # condition of %tile
      percs[percs$BeeID==names(tmap),][,ori_ncol+match(thres,thres_list)] <- ifelse(
        percs[percs$BeeID==names(tmap),]$SSprob >= quantile(tmap,probs=thres,na.rm=T)[[1]],1,0
      )
    }

  }
  
  # write added cols version out to same directory, with a new name
  write.csv(percs,paste0('./',mod,'/',targetVersion,'/','bees_',mod,
                         '_recursPercTiles_Pixel.csv'))
  return(percs)
}