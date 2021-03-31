# loads pre-processed data and does main (non-statistical!) analysis for the methods paper.
# In addition to comparing manual with automated tracking, it also processes a set of
# automated data to allow comparison of some jump parameters between males. A separate
# script (manakin3D_femaleAnalysis.R) must be run after this one if you want to later compare
# male and non-male jump parameters.
#
# cliodhna.quigley@vetmeduni.ac.at, 2020-2021

DEBUG = 0 # if == 1, do extra plotting for sanity checks
          # note: plots won't be displayed if run inside of for loop, so you'll need to manually run per file
          # or add code to save the plot to a file in the for loop

dataDir_2D_MA = 'data/processedData/2Dprojections_manVsaut/' # where to load data (Rds files)
dataDir_2D_aut = 'data/processedData/2Dprojections_aut/'
dataDir_results = 'data/results_methodsMS/'

# parameters:
n_values_corr <- 50 # number of parabola values used in correlation analysis
fps <- 60 # frame rate assumed fixed
spf <- 1/fps # seconds per frame
gravity <- 9.81 # ms-2; used the typical estimate, although better estimate for panama is 9.78 according to wolfram alpha!
maxFrames_params2d <- 3 # up to how many frames difference between first and second point used to estimate take-off vector from 2d points
frameTol_jumpMatcher <- 10 # how many frames difference to allow between manual and automated jumps

library(plotly)
library(dplyr)
library(ggplot2)

source('scripts/3d_utils.R') # contains functions to calculate distance and speed


###########################################
###########################################
###########################################
# FIRST: manual vs automatic:
birdCodes_MvsA <- c('CAN03','CAN04','JUR12') # these are the ones we're interested in for the methods paper

# check in the data directory for which displays per bird
firstRun <- TRUE
for (b in 1:length(birdCodes_MvsA)) {
  bfiles_MvsA <- unique(getDisplayVidCodesFromCSVDir(dataDir_2D_MA,birdcode=birdCodes_MvsA[b])) # unique needed as they are doubled up (man, aut)

  for (f in 1:length(bfiles_MvsA)) {
    cat('processing ',bfiles_MvsA[f],'\n')
    
    # load 2D projected segmented jumps
    projJumps_aut <- readRDS(paste0(dataDir_2D_MA,bfiles_MvsA[f],'_aut.Rds'))
    projJumps_man <- readRDS(paste0(dataDir_2D_MA,bfiles_MvsA[f],'_man.Rds'))
    projJumps_man_full <- projJumps_man
    
    # exclude any G from man if there are no G in aut. This deals with the (currently 
    # only 1) case where we had to remove G from the sapling list for an automated file
    # as the jump was too sparse to be clustered correctly
    if (length(which(projJumps_aut$sap1=='G' | projJumps_aut$sap2=='G'))==0) { # no G rows in aut
      cat('Removing G from man as there is no G in aut\n')
      grows <- which(projJumps_man$sap1=='G' | projJumps_man$sap2=='G') # G rows in man
      projJumps_man <- projJumps_man[-grows,] # G rows are removed
      projJumps_man <- droplevels(projJumps_man) # and factor levels are re-computed
    }
    
    # check that everything matches
    aut_reduced <- projJumps_aut %>% select('jumpIndex','sap1','sap2','jumpIndSinceGflip') %>% unique()
    man_reduced <- projJumps_man %>% select('jumpIndex','sap1','sap2','jumpIndSinceGflip') %>% unique()
    # grab frame index of jump start for troubleshooting
    aut_reduced$frameIndex <- projJumps_aut$frameIndex[as.numeric(rownames(aut_reduced))] 
    man_reduced$frameIndex <- projJumps_man_full$frameIndex[as.numeric(rownames(man_reduced))] 
    
    # now find all jumps to and from the ground 
    gjumps_man <- which(man_reduced$sap1=='G' | man_reduced$sap2=='G') # G rows in man 
    gjumps_aut <- which(aut_reduced$sap1=='G' | aut_reduced$sap2=='G') # G rows in aut
  
    if ((nrow(aut_reduced)-length(gjumps_aut))!=(nrow(man_reduced)-length(gjumps_man))) {
      cat('\tman and aut have mis-match in number of jumps!!! Checking... \n')
      
      # this part is essentially to deal with the CAN04 display, in which
      # the perching threshold isn't optimal, leading to two jumps being 'missed' by
      # automated tracking. The transition between those two jumps involved a single frame
      # dropping under the threshold (the surrounding ones were close!) so was
      # excluded from perch classification. Although the two jumps were then counted
      # as a single jump, because it was a jump from 4->2 and 2->4, this is seen
      # as a big long jump from 4->4, so is excluded from our jump definition and lost.
      
      # go through each jump of each reduced table and confirm matches:
      match_aut_frametol <- sapply(aut_reduced$frameIndex, function(x) any(abs(man_reduced$frameIndex-x)<frameTol_jumpMatcher))
      match_man_frametol <- sapply(man_reduced$frameIndex, function(x) any(abs(aut_reduced$frameIndex-x)<frameTol_jumpMatcher))
      
      # note the ind of the problem ones
      bad_aut_inds <- which(match_aut_frametol == FALSE)
      bad_man_inds <- which(match_man_frametol == FALSE) 
      
      # kick out the rows that don't have a match in each reduced table:
      aut_reduced <- aut_reduced[match_aut_frametol,]
      man_reduced <- man_reduced[match_man_frametol,]
      
      # but note that the jump index might not match, and we leave it that way
      
      # now remove the offending jumps from the full data table
      if (length(bad_aut_inds)>0) {
        cat('\tRemoving non-matched jumps from automated tracking data... (jumps ',bad_aut_inds, ')\n')
        badrows_aut <- which(projJumps_aut$jumpIndex %in% bad_aut_inds) # find rows belonging to bad jumpinds in aut
        projJumps_aut <- projJumps_aut[-badrows_aut,] # those rows are removed
        projJumps_aut <- droplevels(projJumps_aut) # and factor levels are re-computed
        # re-calculate!
        aut_reduced <- projJumps_aut %>% select('jumpIndex','sap1','sap2','jumpIndSinceGflip') %>% unique()
      }
      if (length(bad_man_inds)>0) {
        cat('\tRemoving non-matched jumps from manually tracking data... (jumps ',bad_man_inds, ')\n')
        badrows_man <- which(projJumps_man$jumpIndex %in% bad_man_inds) 
        projJumps_man <- projJumps_man[-badrows_man,] # G rows are removed
        projJumps_man <- droplevels(projJumps_man) # and factor levels are re-computed
        man_reduced <- projJumps_man %>% select('jumpIndex','sap1','sap2','jumpIndSinceGflip') %>% unique()
      }
      
    }
    
    
    # now remove G rows
    # find all jumps to and from the ground 
    grows_man <- which(projJumps_man$sap1=='G' | projJumps_man$sap2=='G') # G rows in man 
    grows_aut <- which(projJumps_aut$sap1=='G' | projJumps_aut$sap2=='G') # G rows in aut
    if (length(grows_man)>0) {
      projJumps_man <- projJumps_man[-grows_man,] # G rows are removed
      projJumps_man <- droplevels(projJumps_man) # and factor levels are re-computed
      man_reduced <- projJumps_man %>% select('jumpIndex','sap1','sap2','jumpIndSinceGflip') %>% unique()
    }
    if (length(grows_aut)>0) {
      projJumps_aut <- projJumps_aut[-grows_aut,] # G rows are removed
      projJumps_aut <- droplevels(projJumps_aut) # and factor levels are re-computed
      # re-calculate!
      aut_reduced <- projJumps_aut %>% select('jumpIndex','sap1','sap2','jumpIndSinceGflip') %>% unique()
    }
    

    if (all(aut_reduced[,c('sap1','sap2')]==man_reduced[,c('sap1','sap2')])) { # compare all cols except frameIndex, jumpIndSinceGflip
      # it's ok!
      cat('\tman and aut match now\n')
      jIndsReversed <- man_reduced$jumpIndSinceGflip
      jInds <- man_reduced$jumpIndex
      sap1s <- man_reduced$sap1
      sap2s <- man_reduced$sap2
    } else {
      cat('\tman and aut dont match after removing G jumps; skipping file!!! \n')
      # browser()
      next
    }
    
    # fit models to all jumps in the dataframe
    models_man <- fitAllJumps(projJumps_man)
    models_aut <- fitAllJumps(projJumps_aut)
    # quantify g.o.f. [adjusted r2 to account for differing sample sizes]
    gof_man <- lapply(models_man,function(x) ifelse(length(x)>1,summary(x)$adj.r.squared,NA))
    gof_aut <- lapply(models_aut,function(x) ifelse(length(x)>1,summary(x)$adj.r.squared,NA))
    
    # plot them (this plot is in a function because it's more complicated)
    if (DEBUG==1) {
      g <- jumpPlotter(projJumps_man,models_man,paste(bfiles_MvsA[f],'Manually annotated data + fits')) + coord_fixed()
      g
    }
    if (DEBUG==1) {
      g <- jumpPlotter(projJumps_aut,models_aut,paste(bfiles_MvsA[f],'Automatically annotated data + fits')) + coord_fixed()
      g
    }
    
    # and compare them
    xmin <- round(sapply(unique(projJumps_man$jumpIndex),function(i) min(projJumps_man$horizontal[which(projJumps_man$jumpIndex==i)])),digits=2)
    xmax <- round(sapply(unique(projJumps_man$jumpIndex),function(i) max(projJumps_man$horizontal[which(projJumps_man$jumpIndex==i)])),digits=2)
    x_limits <- data.frame(xmin=xmin,xmax=xmax) # min and max taken from the manual file for that jump
    corrCoeffs_man_aut <- parabolaCorrelator(models_man,models_aut,x_limits,n_values_corr)
    
    # prepare for table entry
    nTing <- length(corrCoeffs_man_aut)
    tmpStrs <- strsplit(bfiles_MvsA[f],'_')[[1]] # birdcode, date, time, displN
    babyTable <- data.frame(birdCode=rep(tmpStrs[1],nTing), filename=rep(bfiles_MvsA[f],nTing), 
                            vidDate=rep(tmpStrs[2],nTing), vidTime=rep(tmpStrs[3],nTing), 
                            displayNum=rep(as.numeric(gsub("\\D", "", tmpStrs[4])),nTing), 
                            jumpInd=jInds, jumpIndSingGflip=jIndsReversed,sap1=sap1s, sap2=sap2s, 
                            pearsonR_manVSaut=corrCoeffs_man_aut,
                            gof_man=unlist(gof_man), gof_aut=unlist(gof_aut))
    
    if (firstRun) {
      # create table to be filled
      corrCoeffsAllFiles_MvsA <- babyTable
      firstRun = FALSE
    } else{
      corrCoeffsAllFiles_MvsA <- rbind(corrCoeffsAllFiles_MvsA,babyTable)
    }
    
  }
  
  
  
}
# save the data for later analysis/plotting:
fname <- 'manVsAutResult.Rds'
saveRDS(corrCoeffsAllFiles_MvsA, file = paste0(dataDir_results,fname))





###########################################
###########################################
###########################################
# THEN: automatic 
birdCodes_A <- c('CAN03','CAN04','JUR12') # these are the ones we're interested in for the methods paper
birdFirstRun_params <- TRUE
birdFirstRun <- TRUE
# check in the data directory for which displays per bird
for (b in 1:length(birdCodes_A)) {
  cat('processing bird ',birdCodes_A[b],'\n')
  bfiles_A <- getDisplayVidCodesFromCSVDir(dataDir_2D_aut,birdcode=birdCodes_A[b],vidgroup = TRUE) 
  # sapling naming is consistent across files by design.
  # so we can concatenate all videos into one, if they come from the same bird 

  # process files within bird, concatenating into single jump and single model table
  firstRun <- TRUE
  for (f in 1:length(bfiles_A)) {
    cat('processing ',bfiles_A[f],'\n')
    
    # load 2D projected segmented jumps
    projJumps_aut <- readRDS(paste0(dataDir_2D_aut,bfiles_A[f],'_aut.Rds'))
    # remove jumps to/from G
    grows <- which(projJumps_aut$sap1=='G' | projJumps_aut$sap2=='G') # G rows in man 
    if (length(grows)>0) {
      projJumps_aut <- projJumps_aut[-grows,] # G rows are removed
      projJumps_aut <- droplevels(projJumps_aut) # and factor levels are re-computed
    }

    # add a column to keep track of vidgroup (indicates those videos come from the same recording session) and fileind
    tmpStrs <- strsplit(bfiles_A[f],'_')[[1]] # birdcode, date, time, displN, vidGroup
    projJumps_aut$vidGroup <- tmpStrs[5]
    projJumps_aut$fileind <- f
    projJumps_aut$displayNum <- tmpStrs[4]
    # fit models to all jumps in the dataframe
    models_aut <- fitAllJumps(projJumps_aut)
    gof_aut <- lapply(models_aut,function(x) ifelse(length(x)>1,summary(x)$adj.r.squared,NA))
    
    # derive param for jumps:
    # 1. from 2D: jumpParamsFrom2dJumpDF(projectedJumpDF,spf,maxFramesUsed=5)
    xx <- jumpParamsFrom2dJumpDF(projJumps_aut,spf,maxFrames_params2d) 
    xx$birdCode <- birdCodes_A[b]
    xx$date <- tmpStrs[2]
    xx$filecode <- bfiles_A[f]
    xx$displayNum <- as.numeric(regmatches(tmpStrs[4],regexpr('[0-9]+',tmpStrs[4])))
    xx$vidGroup <- tmpStrs[5]
    xx$daytime <- ifelse(as.POSIXct(tmpStrs[3],format='%H%M%S')<as.POSIXct('120000',format='%H%M%S'),'AM','PM')
    # 2. from fits: jumpParamsFromFittedModels(modelled_jumps,gravity)
    yy <- jumpParamsFromFittedModels(models_aut,gravity)
    yy$birdCode <- birdCodes_A[b]
    yy$date <- tmpStrs[2]
    yy$filecode <- bfiles_A[f]
    yy$displayNum <- as.numeric(regmatches(tmpStrs[4],regexpr('[0-9]+',tmpStrs[4])))
    yy$vidGroup <- tmpStrs[5]
    yy$jumpIndex <- xx$jumpIndex  # info is not in model list, so copy from 2d jump params
    yy$jumpSinceGFlipIndex <- xx$jumpSinceGFlipIndex
    yy$sap1 <- xx$sap1
    yy$sap2 <- xx$sap2
    yy$daytime <- ifelse(as.POSIXct(tmpStrs[3],format='%H%M%S')<as.POSIXct('120000',format='%H%M%S'),'AM','PM')
    yy$gof <- unlist(gof_aut)
    yy$totalSpeed <- yy$distanceCurved / (xx$jumpDuration_frames * spf) # this is the curved distance divided by the time taken in seconds
      
    if (DEBUG==1) {
      g <- jumpPlotter(projJumps_aut,models_aut,paste(bfiles_A[f],'Automatically annotated data + fits'))
      g
    }
      
    # and concatenate
    if (firstRun) {
      projJumps_bird <- projJumps_aut
      models_bird <- models_aut
      jumpParamsFrom2d_bird <- xx
      jumpParamsFromParabolas_bird <- yy
      firstRun <- FALSE
    } else {
      projJumps_bird <- rbind(projJumps_bird,projJumps_aut)
      models_bird <- append(models_bird,models_aut)
      jumpParamsFrom2d_bird <- rbind(jumpParamsFrom2d_bird,xx)
      jumpParamsFromParabolas_bird <- rbind(jumpParamsFromParabolas_bird,yy)
    }
      
      
        
  } # for each file
    
  # combine parameter tables over birds:
  if (birdFirstRun_params) {
    params_2d <- jumpParamsFrom2d_bird
    params_fit <- jumpParamsFromParabolas_bird
    birdFirstRun_params <- FALSE
  } else {
    params_2d <- rbind(params_2d,jumpParamsFrom2d_bird)
    params_fit <- rbind(params_fit,jumpParamsFromParabolas_bird)
  }
      
  # now do the correlation analysis for this bird:
  # evaluate match between jumps within and between displays, for those jumps with the same sapling origin and destination
  jumpComparisonResult <- jumpComparisonWithinBetweenDisplays(models_bird,projJumps_bird,n_values_corr)
  if (length(jumpComparisonResult)==1 && is.na(jumpComparisonResult)) {
      cat('\tNo repeated jumps between same sapling pair for this bird, so nothing to compare. Not entering in comparison table\n')
  } else {
      # prepare for table entry
      nTing <- nrow(jumpComparisonResult)
      # it's bird focused now.... add birdCode to large table
      # jInd is jumpIndex within display in case it's a meaningful number
      babyTable <- data.frame(birdCode=rep(birdCodes_A[b],nTing), 
                              saplingCombo=jumpComparisonResult$saplingCombo,
                              dispInd1=jumpComparisonResult$dispIndex1,
                              dispInd2=jumpComparisonResult$dispIndex2,
                              jumpInd1=jumpComparisonResult$jumpIndex1,
                              jumpInd2=jumpComparisonResult$jumpIndex2,
                              vidGroup1=jumpComparisonResult$vidGroup1,
                              vidGroup2=jumpComparisonResult$vidGroup2,
                              jumpSinceFlip1=jumpComparisonResult$jumpSinceFlipIndex1,
                              jumpSinceFlip2=jumpComparisonResult$jumpSinceFlipIndex2,
                              pearsonR=jumpComparisonResult$corrCoeffs)

      if (birdFirstRun) {
        # create table to be filled
        corrCoeffsAllFiles_withinBird <- babyTable
        birdFirstRun <- FALSE
      } else{
        corrCoeffsAllFiles_withinBird <- rbind(corrCoeffsAllFiles_withinBird,babyTable)
      }
      
  }
    
} # birds

# save as excel for judith's stat functions:
write.table(params_fit, file = paste0(dataDir_results,"paramsFromParabolas.xls"), sep = ",", qmethod = "double", row.names=FALSE)
write.table(params_2d, file = paste0(dataDir_results,"paramsFrom2dPoints.xls"), sep = ",", qmethod = "double", row.names=FALSE)
write.table(corrCoeffsAllFiles_withinBird, file = paste0(dataDir_results,"corrCoeffWithinBird_automataticMales.xls"), sep = ",", qmethod = "double", row.names=FALSE)