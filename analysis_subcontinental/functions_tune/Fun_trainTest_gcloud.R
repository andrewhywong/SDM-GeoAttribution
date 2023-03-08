#### func to make train/test data set 
traintestPrep <- function(objective,occdata,bgdata,clim_alt){ 
  
  occdf <- occdata[occdata$occ==1,c('x','y',"occ")]
  bgdf_test <- occdata[occdata$occ==0,c('x','y',"occ")]
  
  if (absenseMode == 'background'){
    ## use pseudoabsence-method to evaluate, 
    ## but use large background points to model
    bgdf <- bgdata
      
  }else if(absenseMode == 'pseudoAbsence'){
    bgdf <- occdata[occdata$occ==0,c('x','y',"occ")]
    # if (!nrow(occdf)==nrow(bgdf)){print(paste0('diff in occ and bg: ',nrow(bgdf)-nrow(occdf)))}
    
  }


  ## split for train and test in 80/20. 
  # 12122022 update bg/occ already been sampled in GEE
  # add a conditional to make full 100% occ for prediction (see do.full.model=T in biomod2)
  if (objective=='evaluate'){
    
    set.seed(78731) 
    bg_train_id <- sample(seq_len(nrow(bgdf_test)), size = floor(0.8 * nrow(bgdf_test)))
    # identical(bg_train_id_test,bg_train_id)
    bg_train <- bgdf_test[bg_train_id, ]
    bg_test <- bgdf_test[-bg_train_id, ]
    
    set.seed(78731) 
    occ_train_id <- sample(seq_len(nrow(occdf)), size = floor(0.8 * nrow(occdf)))
    occ_train <- occdf[occ_train_id, ]
    occ_test <- occdf[-occ_train_id, ]
    
    # bind occ & bg
    occbg_train <- rbind(occ_train,bg_train)
    occbg_test <- rbind(occ_test,bg_test)
    # any(is.na(occbg_train)) # NA produced by raster ext
    
    # extract env data and append
    coordinates(occbg_train)= ~ x + y
    coordinates(occbg_test)= ~ x + y
    
    # extract env var 
    ext_train = raster::extract(clim_alt, occbg_train)
    ext_test = raster::extract(clim_alt, occbg_test)
    
    occbgEnv_train <- cbind(occbg_train,ext_train)
    occbgEnv_test <- cbind(occbg_test,ext_test)
    
    # onlyenv test
    occbg_test_env <- occbgEnv_test[,2:ncol(occbgEnv_test)]
    
    # after extraction, coerce back to df
    occbgEnv_traindf <- occbgEnv_train
    occbgEnv_traindf <- as.data.frame(occbgEnv_traindf)
    occbgEnv_testdf <- as.data.frame(occbg_test_env)

    # any(is.na(occbgEnv_testdf)) # NARM here
    occbgEnv_traindf <- occbgEnv_traindf[!rowSums(is.na(occbgEnv_traindf)),]
    occbgEnv_testdf <- occbgEnv_testdf[!rowSums(is.na(occbgEnv_testdf)),]
    
    ## drop xy 
    occbgEnv_traindfXxy <- subset(occbgEnv_traindf, select=-c(x,y))
    occbgEnv_testdfXxy <- subset(occbgEnv_testdf, select=-c(x,y))
    
    # normalize the covariates
    # *notice: not all the models are fitted on normalized data in
    # the main analysis! Please check the main text.
    final_train <- occbgEnv_traindfXxy
    env_test <- occbgEnv_testdfXxy
    
    # X normalized version of train
    final_train_Xnorm <- final_train
    env_test_Xnorm <- env_test
    
    clim_alt_norm <- stack()
    
    for(v in colnames(final_train[2:ncol(final_train)])){
      # print(v)
      meanv <- mean(final_train[,v])
      sdv <- sd(final_train[,v])
      final_train[,v] <- (final_train[,v] - meanv) / sdv
      env_test[,v] <- (env_test[,v] - meanv) / sdv
      clim_alt_norm <- stack(clim_alt_norm, (clim_alt[[v]] - meanv) / sdv) #normalize stack too for pred mapping
    }
    
    # get a trueLabel of test data for metric calculation
    trueLabel <- as.data.frame(occbgEnv_test)[!rowSums(is.na(as.data.frame(occbgEnv_test))),]$occ
    ttList <- list("traindf" = final_train,
                   "Xnorm_traindf" = final_train_Xnorm,
                   "traindf_xy"=occbgEnv_traindf, 
                   "testdf" = env_test,
                   "Xnorm_testdf" = env_test_Xnorm,
                   "trueLabel" = trueLabel,
                   "normed_stack" = clim_alt_norm)
    
  } else if (objective=='full'){
    
    # predict the map use all data available
    occ_train <- occdf
    bg_train <- bgdf
    
    # bind occ & bg
    occbg_train <- rbind(occ_train,bg_train)

    # any(is.na(occbg_train)) # NA produced by raster ext
    
    # extract env data and append
    coordinates(occbg_train)= ~ x + y
    
    # extract env var 
    ext_train = raster::extract(clim_alt, occbg_train)
    occbgEnv_train <- cbind(occbg_train,ext_train)
    
    # after extraction, coerce back to df
    occbgEnv_traindf <- occbgEnv_train
    occbgEnv_traindf <- as.data.frame(occbgEnv_traindf)
    occbgEnv_traindf <- occbgEnv_traindf[!rowSums(is.na(occbgEnv_traindf)),]
    
    ## drop xy 
    occbgEnv_traindfXxy <- subset(occbgEnv_traindf, select=-c(x,y))
    
    # normalize the covariates (exept vegsys which is categorical)
    # *notice: not all the models are fitted on normalized data in
    # the main analysis! Please check the main text.
    final_train <- occbgEnv_traindfXxy

    # X normalized version of train
    final_train_Xnorm <- final_train
    
    clim_alt_norm <- stack()
    
    for(v in colnames(final_train[2:ncol(final_train)])){
      # print(v)
      meanv <- mean(final_train[,v])
      sdv <- sd(final_train[,v])
      final_train[,v] <- (final_train[,v] - meanv) / sdv
      clim_alt_norm <- stack(clim_alt_norm, (clim_alt[[v]] - meanv) / sdv) #normalize stack too for pred mapping
    }
    
    # get a trueLabel of test data for metric calculation
    trueLabel <- as.data.frame(occbgEnv_train)[!rowSums(is.na(as.data.frame(occbgEnv_train))),]$occ
    ttList <- list("traindf" = final_train,
                   "Xnorm_traindf" = final_train_Xnorm,
                   "traindf_xy"=occbgEnv_traindf, 
                   # "trueLabel" = trueLabel,
                   "normed_stack" = clim_alt_norm)
    
  }
 
  return(ttList)
}

