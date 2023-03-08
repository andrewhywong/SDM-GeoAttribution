####################  GBM model (normalized, weighted) 
# calculating the case weights
train_mod <- function(final_train,final_train_Xnorm,env_test,Xnorm_testdf,traindf_xy,
                      clim_alt_norm,modelNow,s,objective,errorSpecies){
  # # create dir base on absence mode:
  # mdir <- paste0('./SPmodels_',objective,'_',absenseMode,'/model_outputs/',modelNow)
  mdir <- paste0('./model_outputs/',modelNow)
  
  if(!file.exists(mdir)){
    dir.create(file.path(mdir))
    print("The directory is created")
  }
  
  if (modelNow=='brt'){
    myseed <- sum(final_train$occ) + 300 + n
    set.seed(myseed)
    prNum <- as.numeric(table(final_train$occ)["1"]) # number of presences
    bgNum <- as.numeric(table(final_train$occ)["0"]) # number of backgrounds
    wt <- ifelse(final_train$occ == 1, 1, prNum / bgNum)
    
    ptm <- proc.time()
    paste0(mdir,'/',s,'_',modelNow,".png")
    png(filename = paste0(mdir,'/',s,'_',modelNow,".png"),width = 3000, height = 3000,res = 300)
    
    if(inherits(try(
      modeled <- gbm.step(data = final_train,
                      gbm.x = 2:ncol(final_train), # column indices for covariates
                      gbm.y = 1, # column index for response
                      family = "bernoulli",
                      tree.complexity = ifelse(prNum < 50, 1, 5),
                      learning.rate = 0.001,
                      bag.fraction = 0.75,
                      max.trees = 10000, # taking too long, may change 10000 to 5000
                      n.trees = 50,
                      n.folds = 5, # 5-fold cross-validation
                      site.weights = wt,
                      silent = TRUE) # avoid printing the cv results
    ), "try-error")){errorSpecies = append(errorSpecies,paste0("Error for species ", s," for model ", modelNow))}
    
    dev.off()
    t <- proc.time() - ptm
    

  }
  
  if (modelNow=='rf'){
    myseed <- sum(final_train$occ) + 300 + n
    set.seed(myseed)
    ptm <- proc.time()
    # convert the response to factor for producing class relative likelihood
    final_train$occ <- as.factor(final_train$occ)
    if(inherits(try(
      modeled <- randomForest(formula = occ ~.,
                         data = final_train,
                         ntree = 500) # the default number of trees
    ), "try-error")){errorSpecies = append(errorSpecies,paste0("Error for species ", s," for model ", modelNow))}
    
    png(filename = paste0(mdir,'/',s,'_',modelNow,".png"),width = 3000, height = 3000,res = 300)
    plot(modeled, main = paste0("RF_",s))
    dev.off()
    t <- proc.time() - ptm

  }
  
  if (modelNow=='rfdown'){
    myseed <- sum(final_train$occ) + 300 + n
    set.seed(myseed)
    ptm <- proc.time()
    # convert the response to factor for producing class relative likelihood
    final_train$occ <- as.factor(final_train$occ)
    prNum <- as.numeric(table(final_train$occ)["1"]) # number of presences
    bgNum <- as.numeric(table(final_train$occ)["0"]) # number of backgrounds
    # the sample size in each class; the same as presence number
    smpsize <- c("0" = prNum, "1" = prNum)
    if(inherits(try(
      modeled <- randomForest(formula = occ ~.,
                                    data = final_train,
                                    ntree = 1000,
                                    sampsize = smpsize,
                                    replace = TRUE)
    ), "try-error")){errorSpecies = append(errorSpecies,paste0("Error for species ", s," for model ", modelNow))}
    png(filename = paste0(mdir,'/',s,'_',modelNow,".png"),width = 3000, height = 3000,res = 300)
    plot(modeled, main = paste0("RFdownsampled_",s))
    dev.off()
    t <- proc.time() - ptm
    
    
  }

  if (modelNow=='biomod2'){  # BIOMOD EXPORT MAPS ITSELF
    myseed <- sum(final_train$occ) + 300 + n
    set.seed(myseed)
    ptm <- proc.time()

    myRespName <- paste0('biomod2_',s)
    colnames(final_train)[1] <- myRespName
    myResp <- as.numeric(final_train[, myRespName])
    
    myExpl <- data.frame(final_train[, 2:ncol(final_train)])
    myRespXY <- traindf_xy[, c("x", "y")]
    
    if (absenseMode == 'background'){
      myResp[which(myResp == 0)] <- NA
      myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                           expl.var = myExpl,
                                           resp.name = myRespName,
                                           resp.xy = myRespXY,
                                           PA.nb.absences = 50000,
                                           PA.strategy = 'random',
                                           na.rm = TRUE)
      mymodels <- c("GLM","GBM","GAM","CTA","ANN","FDA","MARS","RF","MAXENT.Phillips.2")
      
      
    }else if(absenseMode == 'pseudoAbsence'){
      #### 12142022 UPDATE: CHANGE PA SELECTIONS ARGUMENTS: 
      # if both presence and absence data are available, and there is enough absences : 
      #   set PA.nb.rep = 0 and no pseudo-absence will be selected.
      myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                           expl.var = myExpl,
                                           resp.name = myRespName,
                                           resp.xy = myRespXY,
                                           PA.nb.rep = 0,
                                           # PA.nb.absences = 50000, 
                                           # PA.strategy = 'random',
                                           na.rm = TRUE)
      mymodels <- c("GLM","GBM","GAM","CTA","ANN","FDA","MARS","RF") 
      
    }

    myBiomodOption <- BIOMOD_ModelingOptions()
    
    myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                        models = mymodels,
                                        models.options = myBiomodOption,
                                        NbRunEval = 1,
                                        DataSplit = 100, # use all since splitted data 80/20
                                        models.eval.meth = c("ROC"),
                                        SaveObj = TRUE,
                                        rescal.all.models = FALSE,
                                        do.full.models = TRUE,
                                        modeling.id = paste0(myRespName,'_biomod2'))
    # ensemble modeling using mean probability
    myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                          chosen.models = 'all',
                                          em.by = 'all',
                                          eval.metric = c("ROC"),
                                          eval.metric.quality.threshold = NULL, # since some species's auc are naturally low
                                          prob.mean = TRUE,
                                          prob.cv = FALSE,
                                          prob.ci = FALSE,
                                          prob.median = FALSE,
                                          committee.averaging = FALSE, 
                                          prob.mean.weight = FALSE)
    
    t <- proc.time() - ptm
    
  }
  
  if (modelNow=='maxent'){ #NOTE NOT NORMALIZED!
    maxentPath <- paste0('./tunedMaxent')
    myseed <- sum(final_train_Xnorm$occ) + 300 + n
    set.seed(myseed)
    
    ptm <- proc.time()
    # call function to tune maxent
    # number of folds
    nfolds <- ifelse(sum(final_train_Xnorm$occ) < 10, 2, 5)
    # tune maxent parameters
    param_optim <- maxent_param(data = final_train_Xnorm,
                                k = nfolds,
                                filepath = maxentPath)
    
    # fit a maxent model with the tuned parameters
    # library(rJava)
    # .jinit()
    if(inherits(try(
      modeled <- dismo::maxent(x = final_train_Xnorm[, names(final_train_Xnorm)[which(names(final_train_Xnorm) != 'occ')]],
                              p = final_train_Xnorm$occ, 
                              removeDuplicates = FALSE,
                              path = maxentPath,
                              args = param_optim)
    ), "try-error")){errorSpecies = append(errorSpecies,paste0("Error for species ", s," for model ", modelNow))}

    t <- proc.time() - ptm
  }
  
  # make prediction map
  # the rasters should be normalized with the mean and sd of the training data
  # predicting RF down-sampled on rasters with raster package
  # use index = 2 for the likelihood of presence
  
  #### add condition for full prediction
  if (objective=='evaluate'){
    out_file <- cbind(env_test, trueLabel=paraList$trueLabel)
    
    if (modelNow=='maxent'){
      # maxent use Xnorm
      out_file <- cbind(Xnorm_testdf, trueLabel=paraList$trueLabel)
      predNow <- predict(modeled, Xnorm_testdf, args = c("outputformat=cloglog"))
    } else if (modelNow=='brt') {
      predNow <- predict(modeled, env_test, n.trees = modeled$gbm.call$best.trees,type = "response")
    } else if (modelNow=='rf') {
      predNow <- predict(modeled, env_test, type = "prob")[,"1"]
    } else if (modelNow=='rfdown') {
      predNow <- predict(modeled, env_test, type = "prob")[,"1"]
    } else if (modelNow=='biomod2') {
      myBiomodProj_pred <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                             new.env = env_test, 
                                             proj.name = paste0(myRespName,'_biomod2_proj'),
                                             selected.models = "all",
                                             binary.meth = "ROC",
                                             compress = TRUE,
                                             clamping.mask = TRUE)
      myBiomodEnProj_pred <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj_pred,
                                                        EM.output = myBiomodEM,
                                                        selected.models = "all")
      # extracting the values for ensemble prediction
      myEnProj_DF <- as.data.frame(get_predictions(myBiomodEnProj_pred))
      predNow <- myEnProj_DF[,1]/1000 # biomod scale 0-1000
    }


    out_file$prediction <- predNow
    out_file$spid <- s
    out_file$model <- modelNow
    out_file$time <- t[3]
    
    outList <- list('model'=modelNow,"out_file" = out_file)
    
  } else if (objective=='full'){ # USE FULL TRAIN SET TO GET OUR OBJECTIVES, MAKE PRED HERE TO REDUCE TIME
    if (modelNow=='maxent'){ # Xnorm
      # cloglog gives an estimate between 0 and 1 of probability of presence
      pred_map = predict(modeled, clim_alt, args=c("outputformat=cloglog"))
    }
    if (modelNow=='brt'){
      pred_map <- predict(object = clim_alt_norm,model = modeled,type = "response",index = 2)
    }
    if (modelNow=='rf'){
      pred_map <- predict(object = clim_alt_norm,model = modeled,type = "prob",index = 2)
    }
    if (modelNow=='rfdown'){
      pred_map <- predict(object = clim_alt_norm,model = modeled,type = "prob",index = 2)
    }
    if (modelNow=='biomod2'){
      # project single models
      myBiomodProj_maps <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                             new.env = clim_alt_norm, # change to env_test for pred
                                             proj.name = paste0(myRespName,'_biomod2_proj'),
                                             selected.models = "all",
                                             binary.meth = "ROC",
                                             compress = TRUE,
                                             clamping.mask = TRUE)
      # project ensemble of all models
      # load to check .ird
      # tempget1 <- get(load(paste0(getwd(),"/biomod2.CITE2/","proj_biomod2_CITE2_biomod2_proj/",
      #                            "biomod2.CITE2.biomod2_CITE2_biomod2_proj.projection.out")))
      # tempget2 <- get(load(paste0(getwd(),"/biomod2.CITE2/",
      #                             "biomod2.CITE2.biomod2_CITE2_biomod2ensemble.models.out")))
      
      myBiomodEnProj_maps <- BIOMOD_EnsembleForecasting(projection.output =myBiomodProj_maps, #myBiomodProj_maps,tempget1,
                                                        EM.output =myBiomodEM,#myBiomodEM,tempget2,
                                                        output.format = '.img', # must get a formatted output .img/.tif, or can't open connection
                                                        selected.models = "all")
      myEnProj_map <- get_predictions(myBiomodEnProj_maps)
      pred_map = myEnProj_map[[1]]/1000 # biomod scale 0-1000
    }
    
    outList <- list('model'=modelNow,"out_map" = pred_map)
    
  }
  
  return(outList)

}

