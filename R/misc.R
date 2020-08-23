dump_models <- function(dir=".",fold="FINAL") {
  files <- list.files(path=dir)
  files <- grep("gbm_obj",files,value=TRUE)
  if(length(fold)==1) {
    fold <- paste0("-",fold)
    files <- grep(fold,files,value=TRUE)
    } else {
    stop("fold parameter must be length 1!")
  }
  out <- lapply(files,get_model,dir=dir)
  names(out) <- sapply(out,"[[","variable")
  class(out) <- c("gossisimputemodels","list")
  return(out)
}


get_model <- function(f,dir) {
  load(paste0(dir,f))
  out <- list()
  if(exists("gbm.obj2")) { # Models with two imputation models
    gbm.obj2$trainingData <- NULL
    #gbm.obj2$results <- NULL
    attr(gbm.obj2$terms,".Environment") <- NULL
    out$typeII <- gbm.obj2
  }
  f <- gsub("gbm_obj-",replacement = "",f)
  f <- gsub("\\-.*$","",f)
  out$variable <- f
  gbm.obj$trainingData <- NULL
  attr(gbm.obj$terms,".Environment") <- NULL
  #gbm.obj$control <- NULL
  #gbm.obj$results <- NULL
  out$typeI <- gbm.obj
  gc()
  return(out)
}
