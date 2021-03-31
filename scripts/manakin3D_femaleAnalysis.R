# This script loads female data and analyses it.
# You must have run manakin3D_preProcessing_females_manualVersion first!!
#
# Additionally, as it also appends female data to male data, this
# script must be run after the manakin3D_fittingEtc.R script for males, i.e. after
# male data xls files have been saved in the results directory (data/results_methodsMS/).
#
# cliodhna.quigley@vetmeduni.ac.at, 2020-2021

DEBUG = 0 # if == 1, do extra plotting for sanity checks
          # note: plots won't be displayed if run inside of for loop, so you'll need to manually run per file
          # or add code to save the plot to a file in the for loop

dataDir_2D_F = 'data/processedData/2Dprojections_femaleMan/' # where to load data (Rds files)
dataDir_results = 'data/results_methodsMS/'

# parameters:
n_values_corr <- 50 # number of parabola values used in correlation analysis
fps <- 60 # frame rate assumed fixed
spf <- 1/fps # seconds per frame
gravity <- 9.81 # ms-2; used the typical estimate, although a better estimate for panama is 9.78 according to wolfram alpha!
maxFrames_params2d <- 3 # up to how many frames difference between first and second point used to estimate take-off vector from 2d points

library(plotly)
library(dplyr)
library(ggplot2)

source('scripts/3d_utils.R') # contains functions to calculate distance and speed

birdCodes_F <- c('CAN03','CAN04','HIL02','HIL07','JUR12') # these are the arenas we used for the methods paper

# check in the data directory for which displays per bird
firstRun <- TRUE
fcounter <- 1
for (b in 1:length(birdCodes_F)) {
  bfiles_F <- unique(getDisplayVidCodesFromCSVDir(dataDir_2D_F,birdcode=birdCodes_F[b])) 
  
  for (f in 1:length(bfiles_F)) {
    cat('processing ',bfiles_F[f],'\n')
    
    # load 2D projected segmented jumps
    projJumps_man <- readRDS(paste0(dataDir_2D_F,bfiles_F[f],'_man.Rds'))
    
    man_reduced <- projJumps_man %>% select('jumpIndex','sap1','sap2','jumpIndSinceGflip') %>% unique()
    jIndsReversed <- man_reduced$jumpIndSinceGflip
    jInds <- man_reduced$jumpIndex
    sap1s <- man_reduced$sap1
    sap2s <- man_reduced$sap2
    
    # fit models to all jumps in the dataframe
    models_man <- fitAllJumps(projJumps_man)
    
    # quantify g.o.f. [adjusted r2 to account for differing sample sizes]
    gof_man <- lapply(models_man,function(x) ifelse(length(x)>1,summary(x)$adj.r.squared,NA))
    
    # plot them (this plot is in a function because it's more complicated)
    if (DEBUG==1) {
      g <- jumpPlotter(projJumps_man,models_man,paste(bfiles_F[f],'Manually annotated data + fits'))
      g
    }

    tmpStrs <- strsplit(bfiles_F[f],'_')[[1]] # birdcode, date, time, displN, vidGroup
    tmp_bcode <- paste0('female',formatC(fcounter,width=3,flag='0'))
    tmp_ccode <- tmpStrs[1] # court ID
    fcounter <- fcounter+1
    xx <- jumpParamsFrom2dJumpDF(projJumps_man,spf,maxFrames_params2d) 
    xx$courtID <- tmp_ccode
    xx$birdCode <- tmp_bcode
    xx$date <- tmpStrs[2]
    xx$filecode <- bfiles_F[f]
    xx$displayNum <- as.numeric(regmatches(tmpStrs[4],regexpr('[0-9]+',tmpStrs[4])))
    xx$daytime <- ifelse(as.POSIXct(tmpStrs[3],format='%H%M%S')<as.POSIXct('120000',format='%H%M%S'),'AM','PM')
    xx$jumpSinceGFlipIndex <- NULL # doesn't make sense for females!    
    yy <- jumpParamsFromFittedModels(models_man,gravity)
    yy$courtID <- tmp_ccode
    yy$birdCode <- tmp_bcode
    yy$date <- tmpStrs[2]
    yy$filecode <- bfiles_F[f]
    yy$displayNum <- as.numeric(regmatches(tmpStrs[4],regexpr('[0-9]+',tmpStrs[4])))
    yy$jumpIndex <- xx$jumpIndex  # info is not in model list, so copy from 2d jump params
    yy$sap1 <- xx$sap1
    yy$sap2 <- xx$sap2
    yy$daytime <- ifelse(as.POSIXct(tmpStrs[3],format='%H%M%S')<as.POSIXct('120000',format='%H%M%S'),'AM','PM')
    yy$gof <- unlist(gof_man) 
    yy$totalSpeed <- yy$distanceCurved / (xx$jumpDuration_frames * spf) # this is the curved distance divided by the time taken in seconds
    
    if (DEBUG==1) {
      g <- jumpPlotter(projJumps_man,models_man,paste(bfiles_F[f],'Female data + fits'))
      g
    }
    
    # and concatenate
    if (firstRun) {
      paramsFrom2D_F <- xx
      paramsFromParabolas_F <- yy
      firstRun <- FALSE
    } else {
      paramsFrom2D_F <- rbind(paramsFrom2D_F,xx)
      paramsFromParabolas_F <- rbind(paramsFromParabolas_F,yy)
      
    }
    
  }
  
}

# load saved male data in order to append female data
paramsFromParabolas_M <- read.table(file = paste0(dataDir_results,"paramsFromParabolas.xls"), sep = ",", header=TRUE)
paramsFromParabolas_M$jumpSinceGFlipIndex <- NULL # don't have this info for females so remove
paramsFromParabolas_M$vidGroup <- NULL # don't have for F
paramsFrom2D_M <- read.table(file = paste0(dataDir_results,"paramsFrom2dPoints.xls"), sep = ",", header=TRUE)
paramsFrom2D_M$jumpSinceGFlipIndex <- NULL # not relevant
paramsFrom2D_M$vidGroup <- NULL # not relevant

# join male and female
paramsFromParabolas_F$sex <- 'F' # add sex to both
paramsFromParabolas_M$sex <- 'M'
paramsFromParabolas_M$courtID <- paramsFromParabolas_M$birdCode # in males, courtID is synonomous with birdCode
paramsFromParabolas_MandF <- rbind(paramsFromParabolas_F,paramsFromParabolas_M)

paramsFrom2D_F$sex <- 'F'
paramsFrom2D_M$sex <- 'M'
paramsFrom2D_M$courtID <- paramsFrom2D_M$birdCode # in males, courtID is synonomous with birdCode
paramsFrom2D_MandF <- rbind(paramsFrom2D_F,paramsFrom2D_M)

# save them
write.table(paramsFrom2D_MandF, file = paste0(dataDir_results,"paramsFrom2dPoints_maleAndFemale.xls"), sep = ",", qmethod = "double", row.names=FALSE)
write.table(paramsFromParabolas_MandF, file = paste0(dataDir_results,"paramsFromParabolas_maleAndFemale.xls"), sep = ",", qmethod = "double", row.names=FALSE)