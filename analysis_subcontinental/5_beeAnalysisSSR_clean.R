#### quantile of pred map
# install.packages('gridExtra')
library(rgdal)
library(raster)
library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)
library(sf)
library(stringr)
library(tidyr)

#####
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      median = median(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  # data_sum <- rename(data_sum, replace = c("mean" = varname))
  return(data_sum)
}

##########################################################################################################
##########################################################################################################
##########################################################################################################

concatwd = './nymphs_beeplots/'

## set your working directory to the root folder first
setwd('C:/Users/adwhy/Box/NYMPHS_2023/')
targetVersion <- 'bee_wise_SS' # version of the joint suitability method

########################################################
# if need to reduce bees to those we focused: in current study
reducedBees <- read.csv("./data_beeCSV/bee_selected/selectedBees_LAAN81.csv")
reducedBees <- as.vector(reducedBees$x)
########################################################

## specify current mode
objective = 'full'
absenseMode = 'background' #'background' 'pseudoAbsence'

## set working directory to the target model folder for biomod2 export
setwd(paste0('./SPmodels_',objective,'_',absenseMode,'/model_outputs/'))

## available models
modelslist <- list.dirs('.', full.names = FALSE, recursive = FALSE)
modelslist


## read in bee scores csv 
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


########################################################
#### 1 bee standard 
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
# png(filename = paste0(concatwd,'nSp_models',".png"),width = 5000, height = 2000,res = 300)
plot(t)
# dev.off()


# when we specify x-axis variable inside the aesthetics function aes(). 
# reorder() function sorts the carriers by mean values by default.s
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
  # ylim(90,100)

# png(filename = paste0(concatwd,'species_state',".png"),
# width = 10000, height = 4000,res = 600) 
plot(t)  
# dev.off()

#### after examined new results, import old results and stack boxplots #### 
for(i in modelslist){
  # i = 'biomod2'
  if (i=='biomod2'){
    modelBees_old <- read.csv(paste0('./nymphs_maps','/',i,'/',targetVersion,'/bees_',i,'.csv')) #_2beeCombo;
    modelBees_old$model <- i
  }else{
    t <- read.csv(paste0('./nymphs_maps','/',i,'/',targetVersion,'/bees_',i,'.csv'))
    t$model <- i
    modelBees_old <- rbind(modelBees_old,t)
  }
}

modelBees_old <- modelBees_old[modelBees_old$BeeID %in% modelBees$BeeID,]
modelBees_old$seachPercent <- 100-(modelBees_old$UpperBee/modelBees_old$TotalCells)*100

modelBees_old$addHuman <- 'Env'
modelBees$addHuman <- 'Env+Human'

stacked_modelBees <- rbind(modelBees_old,modelBees)

## levels 
level_order  <- factor(stacked_modelBees$model, 
                       level = c("brt","maxent","rf","rfdown","Ensemble","biomod2"))
# stacked_modelBees$model_level_order <- level_order

# leg_order <- factor(stacked_modelBees$nSp_gr, levels = c("5", "10", "20"))
nsp_order <- factor(stacked_modelBees$nSp, levels = seq(1,20)[-11]) # change accordingly
orderbynSp = unique(stacked_modelBees[order(stacked_modelBees$nSp),]$BeeID)
beeID_ordered  <- factor(stacked_modelBees$BeeID, 
                         level = orderbynSp) # order by NSP

t <- ggplot(stacked_modelBees, aes(x=beeID_ordered, y=seachPercent,
                                   colour=addHuman)) +  # bee genus level ploting is intesesting too
  # geom_point(alpha = 0.5, size = 2) +
  geom_boxplot(lwd=0.5, #0.3 for overall
               outlier.colour="grey", outlier.shape=4,outlier.size=0.1) +
  facet_grid(.~State, scales="free_x") + # "free_x"
  geom_hline(yintercept=99.9, linetype="dashed", color = "grey") +
  # geom_text(aes(2,99.9,label = 99.9, vjust = -1))+
  geom_hline(yintercept=99.5, linetype="dashed", color = "grey") +
  # geom_errorbar(aes(x=level_order,ymin = mean-sd, ymax = mean+sd), width = 0.2)+
  xlab("Bees") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab(paste0("Search Score")) +
  ylim(95,100)
  # theme(panel.grid.major = element_blank(), 
  #       panel.grid.minor = element_blank(),
  #       axis.line = element_line(colour = "grey"),
  #       axis.ticks=element_line(colour = "grey"))

# png(filename = paste0(concatwd,'model_weighted',".png"),
#     width = 3000, height = 2000,res = 300) # 10000 4000 600 for bigger plots
png(filename = paste0(concatwd,'bees_SSR_4x40_zoomin',".png"),
    width = 10000, height = 6000,res = 600) #10000 4000 600 # 3000 2000 300
plot(t)  
dev.off()

################################################################################################################################################
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
polyHitrate <- mTiles_a0 # for plot dodge 
polyHitrate$base <- '% captured by buffer AOI'

# # 1. area of that polygon if point fall in
# mTiles_a1 <- mTiles[,c(1,7:11,ncol(mTiles)-1)]
# colnames(mTiles_a1)[2:6] <- c('99.00%' ,'99.50%' ,'99.90%', '99.95%', '99.99%')

##### 2.area of total polygon AOIs; take average for this? 
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
# # add a penalize to area; for each quantile draw the histogram?
# mTiles_meanArea <- mTiles[,c(1,7:11,ncol(mTiles)-1)] %>% group_by(State,model) %>% summarize_all(mean)


# # 4. Average Distance to all AOI polygons (even fall in)
# mTiles_d2 <- mTiles[,c(1,22:ncol(mTiles)-1)]
# colnames(mTiles_d2)[2:6] <- c('99.00%' ,'99.50%' ,'99.90%', '99.95%', '99.99%')

# # 5. distance to highest pixel
# mTiles_dmax <- mTiles[,c(1,ncol(mTiles)-1,ncol(mTiles)-2)]
# mTiles_a2_avg <- mTiles_dmax %>% group_by(State,model) %>% summarize_all(mean)
# colnames(mTiles_dmax)[3] <- c('Distance to Pixel')
# colnames(mTiles_a2_avg)[3] <- c('Mean Distance to Pixel')


# melt the df to plot ################################################
dfToMelt <- mTiles_a0
mTiles_melted <- melt(dfToMelt, id.vars=c('State','model'))

# to cont'd
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #
# # penalized prediction error
# ggplot(mTiles_melted[mTiles_melted$model=='biomod2',], aes(x = value)) + 
#   facet_grid(State~variable)+  # , scales="free_x" State~variable for all except freq
#   geom_histogram(aes(y = stat(density)*5))+  #, bins = 8 / sum(count)
#   scale_y_continuous(labels = scales::percent)
m1 <- mTiles_melted[mTiles_melted$model=='rfdown' &
                      mTiles_melted$variable=='99.99%',] #mTiles_melted$model=='rfdown'& 
# m1$value <- round(m1$value,0)
# max(m1$value)
# length(m1$State[m1$State=='CAN'])

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

# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## # 
# FOR COMPARISON BETWEEN AOI-AND PIXEL-CAPTURED, JUST NEED TO MAKE CROSS-EQUAL
# ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## #

# cont'd
mTiles_melted <- mTiles_melted[!rowSums(is.na(mTiles_melted)),]
unique(mTiles_melted$variable)
# c('99.00%' ,'99.50%' ,'99.90%', '99.95%', '99.99%')
# 'hit9900_farea' ,'hit9950_farea', 'hit9990_farea', 'hit9995_farea', 'hit9999_farea'
# 'hit9900_tarea' ,'hit9950_tarea' 'hit9990_tarea' ,'hit9995_tarea', 'hit9999_tarea'
# 'hit9900_dtoc' ,'hit9950_dtoc' ,'hit9990_dtoc', 'hit9995_dtoc' ,'hit9999_dtoc'
# 'hit9900_avd' ,'hit9950_avd', 'hit9990_avd', 'hit9995_avd' ,'hit9999_avd'
submtiles <- mTiles_melted[mTiles_melted$variable %in% 
                             c('97.00%','99.00%','99.50%' ,'99.90%', '99.95%', '99.99%'),] #'97.00%','99.00%','99.50%' ,'99.90%', '99.95%', '99.99%'
level_order  <- factor(submtiles$model,
                       level = c("brt","maxent","rf","rfdown","Ensemble","biomod2"))
# level_order  <- factor(mTiles_melted$model,
#                        level = c("brt","maxent","rf","rfdown","Ensemble","biomod2"))

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

t
# png(filename = paste0(concatwd,'model_weighted',".png"),
#     width = 3000, height = 2000,res = 300) # 10000 4000 600 for bigger plots
png(filename = paste0(concatwd,'DtoClosestAOI_9995_9999_AOI',".png"),
    width = 5000, height = 5000,res = 600) #10000 4000 600 # 3000 2000 300
plot(t)  
dev.off()


##############################################################################################################################
##############################################################################################################################
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
  # ylim(90,100)
  theme(aspect.ratio=1)  # make fixed 1:1 

# png(filename = paste0(concatwd,'model_weighted',".png"),
#     width = 3000, height = 2000,res = 300) # 10000 4000 600 for bigger plots
concatwd = "E:\\nymphs_backup\\nymphs_beeplots\\"
png(filename = paste0(concatwd,'ptilesXRegionsXModels_WSS_chulls_10000k',".png"), #100k; 500k; 2000k is good
    width = 5000, height = 5000,res = 600) #10000 4000 600 # 3000 2000 300
plot(t)  
dev.off()

# # Horizontal bar plot
# p + coord_flip()

#### after examined new results, import old results and stack boxplots#### #### #### #### 
for(i in modelslist){
  # i = 'biomod2'
  if (i=='biomod2'){
    modelBees_old <- read.csv(paste0('./SPmodels_3_10_full/model_outputs','/',i,'/',targetVersion,'/bees_',i,'.csv')) #_2beeCombo;
    modelBees_old$model <- i
  }else{
    t <- read.csv(paste0('./SPmodels_3_10_full/model_outputs','/',i,'/',targetVersion,'/bees_',i,'.csv'))
    t$model <- i
    modelBees_old <- rbind(modelBees_old,t)
  }
}

modelBees_old <- modelBees_old[modelBees_old$BeeID %in% modelBees$BeeID,]
modelBees_old$seachPercent <- 100-(modelBees_old$UpperBee/modelBees_old$TotalCells)*100

modelBees_old$addHuman <- 'Env'
modelBees$addHuman <- 'Env+Human'

stacked_modelBees <- rbind(modelBees_old,modelBees)

################################################################################
#### stack polyhitrate and pixelhitrate; 08022022+across env/human comparisons
################################################################################
# ## export hit rate file to graphs loc
# write.csv(pixelHitrate,'./nymphs_beeplots_nbuffer/pixelHitrate_Env.csv',
#           row.names = F)

########################################
#### stack polyhitrate and pixelhitrate and compare
hitrate <- rbind(polyHitrate,pixelHitrate)
hitrate <- melt(hitrate,id.vars=c('State','model','base'))
# hitrate[hitrate == "polygon"] <- '% captured by AOI'
# hitrate[hitrate == "Pixel-based"] <- '% captured by pixel'
hitrate<- hitrate[hitrate$value==1,]
  
# level_order  <- factor(hitrate$model,
#                        level = c("brt","rf")) # GEE
level_order  <- factor(hitrate$model,
                       level = c(modelslist)) # SUB
phit <- ggplot(hitrate, aes(fill=base,x=level_order)) +  #, y=value
  # geom_point(alpha = 0.5, size = 2) +
  geom_bar(position = position_dodge(preserve = "single")) + #stat="identity" #position = position_dodge(preserve = "single")
  facet_grid(State~variable) + # , scales="free_x"
  xlab("Models") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab('No. of successfully identified objects') +
  ylim(0,30)+
  # geom_text(aes(label = value,y=Pos),size = 3, position = position_dodge(width = 0.001)) +
  geom_text(aes(label = ..count..),
            stat = "count", vjust = 1, size = 1.5, colour = "black",
            position = position_dodge(width = 0.5)) +
  theme(aspect.ratio=1)+ # make fixed 1:1
  theme(legend.position="none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1))

phit

png(filename = paste0('../results/HitRate_dodge_buffer5km.png'),# HittedArea_9900_9950_AOI
    width = 6000, height = 5000,res = 600) #10000 4000 600 # 3000 2000 300
plot(phit)  
dev.off()



