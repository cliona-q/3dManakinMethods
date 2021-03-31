# general utilities for Judith's 3D data
# cliodhna.quigley@vetmeduni.ac.at, 2020-2021

# Remove rows with distance of zero from previous frame as they will result in NA for speed.
# Here we keep first frame (distance NA) to allow second frame distance / speed calculations later
cleanXYZdata <- function(xyzDF) {
  # exclude frames with repeated coding of location (output of loopy sometimes fills in 'blank' frames)
  tmpDF <- EuclidDistanceFromXYZ(xyzDF) # new column with distance is added
  # filter only the rows with non-zero distances between frames; keep NA entries too (first frame can't have a speed)
  tmpDF <-filter(tmpDF, tmpDF$EuclidDist3D>0 | is.na(tmpDF$EuclidDist3D)) # >0 is an ok assumption as distance can't be negative
  tmpDF <- select(tmpDF,-'EuclidDist3D') # remove column with distance
  return(tmpDF) # return row-reduced data-frame
}

# imports loopy tracking file, keeping only the necessary columns
# and removes useless frame estimates
parseJJcsv_1bird <- function(csvFilename) {
  # ASSUMES ONLY ONE BIRD IS TRACKED
  fulldata <-read.csv(csvFilename)
  # REMOVE ROWS THAT HAVE 0 CAMERAS -> X,Y,Z values are copied from previous frames so useless
  # and select only the columns we want
  usefuldata <- fulldata %>% filter(ncameras>0) %>% select('frame_number','x','y','z','err','oid','ncameras')
  # err might be used later; oid increments whenever the object is lost for >= 1 frame
  # ncameras: number of cameras used to estimated 3d point
  
  # give user some insight into number of frames etc.
  f1 <- usefuldata$frame_number[1]
  fend <-usefuldata$frame_number[nrow(usefuldata)]
  cat('parseJJcsv_1bird: file contained',nrow(usefuldata),'rows spanning',fend-f1+1,'frames [',f1,'-',fend,']',' \n')
  
  return(usefuldata)
}

# cleans data by checking for outlier speeds frame-by-frame
rejectHighSpeeds <- function(xyzDF,speed_threshold,spf) {
  # if there are single outlier points, blindly applying a threshold will remove outlier +
  # the subsequent 'return' position, as travel to and from outlier has high speed.
  # This should be able to deal with >1 frame of outlier - each arrival at outlier will
  # be rejected in subsequent runs of while loop.
  # For this reason, we do it row by row. It's still cheapest to recalculate speed
  # for entire table, I think, so need to make it recursive
  ok <- FALSE
  workingTable <- xyzDF
  rcounter <- 2 # start at the second row -> the first one has no defined speed (NA)
  while (!ok) {
    if (nrow(workingTable)==rcounter) {
      ok <- TRUE # we're done after checking this last row
    }
    tmp_speed <- workingTable$speed3d[rcounter]
    if (tmp_speed>speed_threshold) { # 'arrival speed' is too fast, so kick out this frame
        workingTable <- workingTable[-rcounter,] # removes the row
        # and DO not increment the row counter, because we need to check the next frame, which has jumped into this row
        
        # remove the distance and speed columns, and recaluclate speed
        workingTable <- select(workingTable,-'EuclidDist3D',-'speed3d')
        workingTable <- speedFromXYZ(workingTable,spf)
    } else { # 'arrival speed' is ok
        # do nothing to the table, but increment the row counter
        rcounter <- rcounter + 1
    }
  }
  okspeedtable <- workingTable
  return(okspeedtable)
}

getDisplayVidCodesFromCSVDir <- function(dirname,birdcode=NULL,vidgroup=FALSE) {
  # look into the directory DIRNAME and create a character vector of the date_time_displN strings
  # which summarise each csv file's contents. If optional argument birdcode is specified, limits
  # search to files prepended by birdcode
  fullfiles <- dir(dirname)
  # then grab the unique ID portion (8 digit date _ 6 digit time _ displ 1 (or 2) digit display index)
  if (is.null(birdcode)) {
    pattern <- '(\\d{8}[_]\\d{6}[_displ]{6}\\d*)'
  } else {
    if (vidgroup==FALSE) { # man vs aut
      pattern <- paste0('(',birdcode,'_\\d{8}[_]\\d{6}[_displ]{6}\\d*)') # this is the chunk we want to pull out of the full name
    } else { # aut are saved with a vidgroup code as well
      pattern <- paste0('(',birdcode,'_\\d{8}[_]\\d{6}[_displ]{6}\\d*[_]\\d*)') # this is the chunk we want to pull out of the full name
    }
  }
  m <- regexpr(pattern,fullfiles)
  fchunks <- regmatches(fullfiles,m)
  return(fchunks)
}

findLoadcsvJJ <- function(dirname,filecode) {
  # look into directory DIRNAME, find and load file containing
  # code FILECODE. Return data table if file is present, or
  # NA if file not found
  fullfiles <- dir(dirname) # all file names
  # searhc for code we want
  matchInd <- grep(filecode,fullfiles)
  if (length(matchInd)>1) {
    cat('ERROR: more than one file found:',filecode,'\n')
    data3d <-NA
  } else{
    if (length(matchInd)==0) {
      cat('ERROR: no file found:',filecode,'\n')
      data3d <- NA
    } else {
      data3d <- parseJJcsv_1bird(paste0(dirname,fullfiles[matchInd]))    
    }
  }
  return(data3d)
}

getSaplingGroundXYLocations_code <- function(filecode,saplingfname) {
  # uses a data_time_displayN code as input
  # figure out which sapling and ground locations are appropriate.
  # saplingfname should point to the csv file with all of this info for all files.
  sapdf <- read.csv(saplingfname)
  # now extract rows of table with filename matching fchunk
  goodrows <- which(sapdf$filename==filecode)
  tmpTable <- sapdf[goodrows,]
  smallDF <- select(tmpTable,-'filename') # return everything except the file code
  return(smallDF)
}

EuclidDistanceFromXYZ <- function(xyzdata) {
  # calculates Euclidean distance between xyz coordinates, adding it as a new column to dataframe
  xyz <- data.matrix(select(xyzdata,'x','y','z'))
  diffs <- diff(xyz,lag=1,differences = 1)
  sumsqdiffs <- rowSums(diffs^2)
  dist3d <- NULL
  dist3d[1] <- NA
  dist3d[2:(length(sumsqdiffs)+1)] <- sqrt(sumsqdiffs)
  
  xyzdata['EuclidDist3D'] <-dist3d
  return(xyzdata)
}

speedFromXYZ <- function(xyzdata,spf) {  # m/s
  # calculates speed from xyz coordinates
  # first need distance
  xyzdist <- EuclidDistanceFromXYZ(xyzdata)
  # then time between samples (cannot assume no missing frames)
  frdiff <- diff(xyzdist$frame_number,lag = 1,differences = 1)
  dt <- frdiff*spf
  if (any(dt==0)) { # sanity check for dt of 0 which produces inf speeds
    cat('Problem: 0 time difference found between frames (repeated frame), cannot compute speed\n')
    xyzdist <- NA #  otherwise we'd have an infinite speed
  } else {
    # speed = dist / time
    xyzdist$speed3d <- NA # first frame can't have a speed
    xyzdist$speed3d[2:(nrow(xyzdist))] <- (xyzdist$EuclidDist3D[2:(nrow(xyzdist))])/dt
  }
  return(xyzdist)
}

# NEAREST NEIGHBOUR VERSION, conservatively excluding brief (1 frame) dips below speed threshold
assignLowSpeedsToSaplings <- function(speedDF,sitting_threshold,saplings_xy) {
  # mark and select frames with low speed, i.e. bird has moved little since last frame
  slowrows <- which(speedDF$speed3d<sitting_threshold) 
  
  # exclude any occurrences of fly, perch, fly for adjacent frames. This is likely to be a technical issue with the recordings (e.g. JUR12)
  mask <- sapply(slowrows, function(x) speedDF$frame_number[x+1]-speedDF$frame_number[x-1]==2 & speedDF$speed3d[x-1]>sitting_threshold & speedDF$speed3d[x+1]>sitting_threshold)
  # the mask marks slowrows that have available adjacent frames, for which both are above the speed threshold
  suspect_slowrows <- slowrows[mask] # grab the row indices
  slowrows <- setdiff(slowrows,suspect_slowrows) # and remove them from the vector
  
  tmp_data <- speedDF[slowrows,]
  # use the sapling + ground xys to cluster the slow xy data into separate sapling / ground
  xy_slow <- matrix(c(tmp_data$x,tmp_data$y),ncol=2,byrow=FALSE)
  
  # calculate Euclidean distance between each perch point (xy_slow row) and each sapling xy, by sapling
  eucdist <- matrix(nrow=nrow(tmp_data),ncol = nrow(saplings_xy))
  for (s in 1:nrow(saplings_xy)) { # for each sapling
    sap <- c(saplings_xy$x[s],saplings_xy$y[s])
    eucdist[,s] <- sqrt((xy_slow[,1]-sap[1])^2 + (xy_slow[,2]-sap[2])^2) 
  }
  mininds <- apply(eucdist, 1, function(x) which.min(x)) # indexes of closest saplings (NOT THEIR IDS)
  
  # add new column to table containing NA for probable jump frames and sampling ID for slow frames
  speedDF$sapling_num <- NA
  speedDF$sapling_num[slowrows] <- as.character(saplings_xy$sapID[mininds]) # put ID of sapling in to new column as we have it from the Judith excel table
  
  return(speedDF)
}

# segment full tracking into jumps
extractTrajectories <- function(speedDF) {
  JUMPLENGTHTHRESH = 3 # minimum duration of a jump in frames
  # add new columns jumpIndex, startSap, endSap
  nrows <- nrow(speedDF)
  # vectors we'll use
  saps <- speedDF$sapling_num
  
  # speed attached to a given frame is between this and the previous frame.
  # so a below-threshold speed has not moved from the previous position
  jumping <- is.na(saps) # T means bird is in the air or just was (might be on sapling now), F means bird is most likely on a sapling (from sapling to sapling)
  jumping <- paste(as.character(as.integer(jumping)),collapse = "") # convert to a string of 1s and 0s
  
  takeoffs <- gregexpr("01",jumping)[[1]] # this returns the indexes of the frame of last sit on sapling before jump
  landings <- gregexpr("10",jumping)[[1]] # this returns the indexes of the frame of first sit on sapling
                                          # first index of 01 is likely the landing on the sapling, as the 1 indicates that previous frame was v close to current
  if (landings[1] <= takeoffs[1]) {   # if the first landing is earlier than the first take-off, exclude it as it's the bird entering the arena
    landings <- landings[-1] # remove the first one
  }
  
  # vectors we'll fill
  jumpIndex <- integer(nrows) # intialised to have all zeros
  startSap <- integer(nrows)
  endSap <- integer(nrows)
  
  # need a for loop :(
  jumpNum <- 1
  for (j in 1:length(landings)) {
    # for each jump
    # which frames are included?
    frinds <- takeoffs[j]:landings[j]
    if ((speedDF$frame_number_readable[frinds[length(frinds)]]-speedDF$frame_number_readable[frinds[1]])<=JUMPLENGTHTHRESH) { # if the jump is too short, ignore
      # ignore
    } else if (saps[takeoffs[j]]==saps[landings[j]+1]) { # if the jump is from and to the same sapling, ignore
      # note: need to look at the frame after the landing to figure out which sapling we're on!
      # ignore
    } else {
      # each frame-in-table including and between the takeoff and landing index gets the following values:
      jumpIndex[frinds] <- jumpNum # index of jump (0=no jump happening right now)
      startSap[frinds] <- speedDF$sapling_num[takeoffs[j]] # which sapling he took off from
      endSap[frinds] <- speedDF$sapling_num[landings[j]+1] # which sapling he will land on
    
      jumpNum <- jumpNum+1
    }
  }
  
  speedDF$jumpIndex <- jumpIndex
  speedDF$startSapling <- startSap
  speedDF$stopSapling <- endSap
  return(speedDF)
}

# adds a column to log reverse index of the jump wrt jump to ground (G)
addJumpIndexSinceG <- function(trajDF) {
  G_jInds <- unique(trajDF$jumpIndex[which(trajDF$startSapling=='G')])
  if (length(G_jInds)==0) { # no jumps from G 
    cat('Warning: addJumpIndexSinceG: Did not find any jumps from G in this file\n')
    trajDF$jumpIndSinceGflip <- NA
  } else {
    if (length(G_jInds)==1) {
      subtractor <- G_jInds
      trajDF$jumpIndSinceGflip <- trajDF$jumpIndex - subtractor
    } else {
      cat('Warning: addJumpIndexSinceG: multiple jumps from G in this file\n')
      trajDF$jumpIndSinceGflip <- NA
    }
  }
  return(trajDF)
}
 
flattenJumps <- function(trajDF) {
  # takes 3D jumps and projects them onto a 2D plane which is defined in its vertical aspect by the vertical coord from original 3D,
  # horizontal is a slice through the space containing start and end points of each jump
  njumps <- max(trajDF$jumpIndex)
  njumpedframes <- sum(trajDF$jumpIndex>0) # jumpIndex is 0 when bird is sitting [except for last and first sit before and after jump! They are in the jump]
  # create empty vectors (zeros)
  jumpIndex <- integer(njumpedframes)
  jumpIndSinceGflip <- integer(njumpedframes)
  sap1 <- integer(njumpedframes)
  sap2 <- integer(njumpedframes)
  horizontal <- integer(njumpedframes)
  vertical <- integer(njumpedframes)
  frameIndex <- integer(njumpedframes)
  # and a counter
  jcounter <- 1
  # and fill them
  for (n in 1:njumps) {
    jinds <- which(trajDF$jumpIndex==n)
    xyz <- matrix(c(trajDF$x[jinds],trajDF$y[jinds],trajDF$z[jinds]),nrow=3,byrow=T) # 3 x n matrix of 3d coords, xyz
   
    # set it so that x,y of first frame is 0,0 [height must remain unchanged]
    xyz <- xyz - c(xyz[1,1],xyz[2,1],0) # subtract first x from all x, first y from all y, 0 from all z
    
    blah <- projectToPlane(xyz,c(xyz[,1]),c(xyz[,ncol(xyz)]))
    projd <- blah$projected
    # we also have blah$projectionMat if needed later!
    
    jtinds <- jcounter:(jcounter+ncol(projd)-1) # target indices into new vectors
    jumpIndex[jtinds] <- n
    jumpIndSinceGflip[jtinds] <- trajDF$jumpIndSinceGflip[jinds[1]]
    sap1[jtinds] <- trajDF$startSapling[jinds[1]]
    sap2[jtinds] <- trajDF$stopSapling[jinds[1]]
    horizontal[jtinds] <- projd[1,]
    vertical[jtinds] <- projd[2,]
    frameIndex[jtinds] <- trajDF$frame_number_readable[jinds]
    jcounter <- jtinds[length(jtinds)]+1
  }
  jumpDF <- data.frame(jumpIndex,jumpIndSinceGflip,sap1,sap2,horizontal,vertical,frameIndex)
  return(jumpDF)
}

# called by flattenJumps; does simple linear algebra 
projectToPlane <- function(xyz,startxyz,endxyz) {
  # take a matrix of xyz points (each column is 1 point) and project them onto a 2d plane
  # calculated from start and end points, and orthogonal to z axis
  se_diff <- endxyz - startxyz  
  # gram-schmidt method: v2 = x2 - (x2.v1/v1.v1)v1
  v1 <- c(0,0,1) # first vector is z unit vector
  v2 <- se_diff - ((sum(se_diff*v1)/sum(v1*v1))*v1) # this just eliminates the z part in order to orthogonalize this dimension
  # and make them orthonormal by normalising:
  u1 <- v1 # I made it a unit vector from the beginning
  u2 <- v2 / (sqrt(sum(v2*v2))) # divide it by its norm
  
  # now make it into a projection matrix (already transposed for the upcoming multiplication)
  Pt <- rbind(u2,u1) # u2 goes first because it's the horizontal component
  rownames(Pt) <- c('xy','z')
  # and project it:
  projected <- Pt %*% xyz # 2x3 proj mat times 3xn coordinate mat
  outputs <- list('projected' = projected,'projectionMat' =Pt) # also return projection matrix for re-projecting it if necessary
  return(outputs)
}


fitAllJumps <- function(jumpsXYdf,minN=5) {
  # takes a table of xy jump coordinates and returns a list of fitted models.
  # default minimum number of points is 5; any fewer and we don't bother fitting.
  # Equation is poly(x,2)
  jInds <- unique(jumpsXYdf$jumpIndex)
  njumps <- length(jInds)
  modelled_jumps <- vector("list",njumps)
  for (j in 1:njumps) {
    # minitable of those rows
    tmp_df <- jumpsXYdf[which(jumpsXYdf$jumpIndex==jInds[j]),]
    # check sufficiently many data-points present
    if (nrow(tmp_df)<minN) {
      # can't fit it so don't try
      modelled_jumps[[j]] <- NA
    } else {
      # EXTREMELY IMPORTANT: if we want to interpret the coefficients later, use raw and not orthogonal poly [uncorrelated]
      tmp_model <- lm(vertical ~ poly(horizontal, 2, raw=T), tmp_df, x=TRUE) # save x values - need them later
      modelled_jumps[[j]] <- tmp_model      
    }
    
  }
  return(modelled_jumps)
}

jumpPlotter <- function(jumpsXYdf,modelled_jumps,titleText) {
  # takes a table of xy jump coords and the list of fitted models
  # and plots them
  jInds <- unique(jumpsXYdf$jumpIndex)
  njumps <- length(jInds)
  jumpSaps <- unique(select(jumpsXYdf,c('sap1','sap2'))) # table containing unique combos of sap1 and sap2
  sapLevels <- paste0(jumpSaps$sap1,'->',jumpSaps$sap2) # to translate index into a meaningful label later
  
  # initialise stuff to be filled in during for loop
  jumpsXYdf$sapIndex <- NA # new column of table
  # future data frame
  fittedXYs <- data.frame( jumpIndex = numeric(0), x = numeric(0), fit = numeric(0), lwr = numeric(0),  upr = numeric(0),  sapIndex = as.factor(numeric(0)))
  
  for (j in 1:njumps) {
    # identify rows for this jump:
    jrows <- which(jumpsXYdf$jumpIndex==jInds[j])
    # which sapling combo is it
    tmp_sapIndex <- which((jumpSaps$sap1==jumpsXYdf$sap1[jrows[1]]) & (jumpSaps$sap2==jumpsXYdf$sap2[jrows[1]])) # index into jumpSaps table
    # add it to the table in the new column
    jumpsXYdf$sapIndex[jrows] <- tmp_sapIndex
    
    # and add to the fitted data table
    tmp_model <- modelled_jumps[[j]]
    if (length(tmp_model)>1) { 
      smooth_x <- data.frame(horizontal=seq(from=0,by=0.01,to=max(jumpsXYdf$horizontal[jrows])))
      tmp_fit <- predict(tmp_model, smooth_x, interval="confidence")
      tmp_fitDF <- data.frame(tmp_fit)
      tmp_fitDF$x <- smooth_x$horizontal
      tmp_fitDF$sapIndex <- tmp_sapIndex
      tmp_fitDF$jumpIndex <- jInds[j]
      fittedXYs <- merge(fittedXYs,tmp_fitDF,all=T)
    } else {
      # do nothing      
    }
  }
  
  jumpsXYdf$sapIndex <- as.factor(sapLevels[jumpsXYdf$sapIndex])
  fittedXYs$sapIndex <- as.factor(sapLevels[fittedXYs$sapIndex])
  
  # facet_wrap
  # first the fits
  g <- ggplot(fittedXYs) + 
    geom_ribbon(aes(x=x,ymin=lwr,ymax=upr,fill=as.factor(jumpIndex),color=NULL),alpha=0.15) + guides(fill=F) +
    geom_line(aes(x=x,y=fit,color=as.factor(jumpIndex)),size=0.5)
  g <- g + facet_wrap(~ sapIndex)  + 
    coord_fixed() # keep x units same as y units to allow visualisation of angle
  # and then the data on top, wrapping it and dealing with labels
  g <- g + geom_point(data=jumpsXYdf,aes(x=horizontal,y=vertical,color=as.factor(jumpIndex))) + 
    xlab('Horizontal distance from jump origin (m)') +
    ylab('Height above ground (m)') +
    labs(color='Jump Index') +
    ggtitle(titleText)
  g
  return(g)
}

parabolaCorrelator <- function(modelled_jumps_1,modelled_jumps_2,x_limits,n_values) {
  # takes two lists of jump models which must have the same length and correlates the predicted y values
  # for the given x_limits (table containing min and max for each model pair), using n_values.
  # Ignore the variable naming; it was originally written for manual and automated correlation but
  # order of input or type of data doesn't really matter at all
  
  nmodels <- length(modelled_jumps_1)
  corCoeffs <- numeric(nmodels)*NA
  
  for (m in 1:nmodels) {
    # first check whether the manual model really worked. Some of them didn't fit because of too
    # few points so there won't be model coefficients to use later
    tmp_mod_man <- modelled_jumps_1[[m]]
    tmp_mod_aut <- modelled_jumps_2[[m]]
    
    if (length(tmp_mod_man)==1 | length(tmp_mod_aut)==1) { # too few data points in one or other
      next # do nothing, corCoeffs[m] will be a NA; skip to next model
    }
    
    if (all(tmp_mod_man$residuals==0)  | all(tmp_mod_aut$residuals==0)) {
      # do nothing, corCoeffs[m] will be a NA; skip to next model
    } else {
    
    toCor <- data.frame(manY=numeric(n_values)*NA,autY=numeric(n_values)*NA)
    
    # create x values for this model:
    xvals <- data.frame(horizontal=seq(from=x_limits$xmin[m],to=x_limits$xmax[m],length.out=n_values))
    # use each model to create y values:
    tmp_fit_aut <- predict(tmp_mod_aut, newdata=xvals)
    toCor$autY <- tmp_fit_aut
    
    tmp_fit_man <- predict(tmp_mod_man,newdata=xvals)
    toCor$manY <- tmp_fit_man
    
    # and correlate
    r <- cor(toCor,method='pearson')
    corCoeffs[m] <- r[lower.tri(r)] # take relevent value from full corr matrix
    }
  }
  return(corCoeffs)
}

jumpComparisonWithinDisplay <- function(modelled_jumps,jumpsXYdf, n_values) {
# this function compares jumps within a display which have the same sapling start and end.
# If there are no matching jumps, it returns an NA instead of the data frame summarisng the result.
# Third input is n_values used to determine how many xvalues to compare in the correlation of fits
  
  jIDdID <- unique(select(jumpsXYdf,c('jumpIndex','displayNum')))
  njumps <- nrow(jIDdID)
  jumpSaps <- unique(select(jumpsXYdf,c('sap1','sap2'))) # table containing unique combos of sap1 and sap2
  sapLevels <- paste0(jumpSaps$sap1,'->',jumpSaps$sap2) # to translate index into a meaningful label later
  
  # initialise stuff to be filled in during for loop
  jumpsXYdf$sapIndex <- NA # new column of table

  sapIndex <- numeric(njumps)
  jumpIndex <- numeric(njumps)
  dispIndex <- numeric(njumps)
  jumpIndexSinceGFlip <-numeric(njumps)
  for (j in 1:njumps) {
    # identify rows for this jump:
    jrows <- which(jumpsXYdf$jumpIndex==jIDdID[j,1] & jumpsXYdf$displayNum==jIDdID[j,2])
    # which sapling combo is it
    tmp_sapIndex <- which((jumpSaps$sap1==jumpsXYdf$sap1[jrows[1]]) & (jumpSaps$sap2==jumpsXYdf$sap2[jrows[1]])) # index into jumpSaps table
    # which jump since G flip is it
    tmp_jSinceInd <- unique(jumpsXYdf$jumpIndSinceGflip[jrows])
    # add it to the table in the new column
    jumpsXYdf$sapIndex[jrows] <- tmp_sapIndex
    sapIndex[j] <- tmp_sapIndex
    jumpIndex[j] <- jIDdID[j,1]
    dispIndex[j] <- jIDdID[j,2]
    jumpIndexSinceGFlip[j] <- tmp_jSinceInd
  }
  repeatedSapIDs <- unique(sapIndex[duplicated(sapIndex)])
  
  if (length(repeatedSapIDs)==0) {
    cat('WARNING jumpComparisonWithinDisplay: no repeated saplings found, returning NA\n')
    jumpComparisonResult <- NA
    return(jumpComparisonResult)
  }
  
  result_sapIDs <- numeric(0)
  result_Rs <- numeric(0)*NA
  result_jumpIDs_1 <- numeric(0)
  result_jumpIDs_2 <- numeric(0)
  result_dispIDs_1 <- numeric(0)
  result_dispIDs_2 <- numeric(0)
  result_jumpSinceFlipIDs_1 <- numeric(0)*NA
  result_jumpSinceFlipIDs_2 <- numeric(0)*NA
  
  for (s in 1:length(repeatedSapIDs)) {
    # for each repeated sapIndex
    tmp_sapIndex <- repeatedSapIDs[s]
    # which jumps was it?
    jumpInds <- which(tmp_sapIndex==sapIndex) # index
    jumpIDs <- jumpIndex[jumpInds]            # value
    dispIDs <- dispIndex[jumpInds]
    
    # how many repetitions?
    nreps <- length(jumpInds)
    combos <- combn(jumpInds,2) # generates all possible pairs; index
    npairs <- ncol(combos)
    model_list_1 <- list()
    model_list_2 <- list()
    xmin <- 0
    xmax <- 0
    for (p in 1:npairs) {
      # grab the models and the limits and get the corr coeff from the function
      model_list_1[[p]] <- modelled_jumps[[combos[1,p]]]
      model_list_2[[p]] <- modelled_jumps[[combos[2,p]]]
      
      # identify rows for this jump pair:
      jrows1 <- which(jumpsXYdf$jumpIndex==jIDdID[combos[1,p],1] & jumpsXYdf$displayNum==jIDdID[combos[1,p],2])
      jrows2 <- which(jumpsXYdf$jumpIndex==jIDdID[combos[2,p],1] & jumpsXYdf$displayNum==jIDdID[combos[2,p],2])
      
      # here use unique limits for each paired comparison
      xmin[p] <- min(jumpsXYdf$horizontal[c(jrows1,jrows2)])
      xmax[p] <- max(jumpsXYdf$horizontal[c(jrows1,jrows2)])
    }
    x_limits <- data.frame(xmin=xmin,xmax=xmax)
    tmp_Rs <- parabolaCorrelator(model_list_1,model_list_2,x_limits,n_values) 
    
    result_sapIDs <- c(result_sapIDs,rep(tmp_sapIndex,npairs))
    result_jumpIDs_1 <- c(result_jumpIDs_1,jumpIndex[combos[1,]])
    result_jumpIDs_2 <- c(result_jumpIDs_2,jumpIndex[combos[2,]])
    result_dispIDs_1 <- c(result_dispIDs_1,dispIndex[combos[1,]])
    result_dispIDs_2 <- c(result_dispIDs_2,dispIndex[combos[2,]])
    result_jumpSinceFlipIDs_1 <- c(result_jumpSinceFlipIDs_1,jumpIndexSinceGFlip[combos[1,]])
    result_jumpSinceFlipIDs_2 <- c(result_jumpSinceFlipIDs_2,jumpIndexSinceGFlip[combos[2,]])
    
    result_Rs <- c(result_Rs,tmp_Rs)
    
  }
  jumpComparisonResult <- data.frame('saplingCombo'=as.factor(result_sapIDs),
                                     'jumpIndex1'=result_jumpIDs_1,'jumpIndex2'=result_jumpIDs_2,
                                     'jumpSinceFlipIndex1'=result_jumpSinceFlipIDs_1,'jumpSinceFlipIndex2'=result_jumpSinceFlipIDs_2,
                                     'dispIndex1'=result_dispIDs_1,'dispIndex2'=result_dispIDs_2,
                                     'corrCoeffs'=result_Rs)
  levels(jumpComparisonResult$saplingCombo) <- sapLevels
  return(jumpComparisonResult)
}

# new version that compares within and across recordings and records
# whether it was within/across for later analysis. Note that you need
# to specify 'videoGroup' in the sapling table for it to work. Use a dummy
# value for it but don't leave it empty, as that will default to NA and the 
# matching procedure will fail, causing a cryptic error!
jumpComparisonWithinBetweenDisplays <- function(modelled_jumps,jumpsXYdf, n_values) {
  # this function compares jumps within a display which have the same sapling start and end.
  # If there are no matching jumps, it returns an NA instead of the data frame summarisng the result.
  # Third input is n_values used to determine how many xvalues to compare in the correlation of fits
  
  jIDdID <- unique(select(jumpsXYdf,c('jumpIndex','displayNum','vidGroup'))) # added vidgroup as well
  njumps <- nrow(jIDdID)
  jumpSaps <- unique(select(jumpsXYdf,c('sap1','sap2'))) # table containing unique combos of sap1 and sap2
  sapLevels <- paste0(jumpSaps$sap1,'->',jumpSaps$sap2) # to translate index into a meaningful label later
  
  # initialise stuff to be filled in during for loop
  jumpsXYdf$sapIndex <- NA # new column of table
  
  sapIndex <- numeric(njumps)
  jumpIndex <- numeric(njumps)
  dispIndex <- numeric(njumps)
  jumpIndexSinceGFlip <-numeric(njumps)
  vidGroup <- numeric(njumps)
  for (j in 1:njumps) {
    # identify rows for this jump:
    jrows <- which(jumpsXYdf$jumpIndex==jIDdID[j,1] & jumpsXYdf$displayNum==jIDdID[j,2] & jumpsXYdf$vidGroup==jIDdID[j,3])
    # which sapling combo is it (uses only the first row as it's the same for all jump entries)
    tmp_sapIndex <- which((jumpSaps$sap1==jumpsXYdf$sap1[jrows[1]]) & (jumpSaps$sap2==jumpsXYdf$sap2[jrows[1]])) # index into jumpSaps table
    # which jump since G flip is it
    tmp_jSinceInd <- unique(jumpsXYdf$jumpIndSinceGflip[jrows])
    # add it to the table in the new column
    jumpsXYdf$sapIndex[jrows] <- tmp_sapIndex
    sapIndex[j] <- tmp_sapIndex
    jumpIndex[j] <- jIDdID[j,1]
    dispIndex[j] <- jIDdID[j,2]
    vidGroup[j] <- jIDdID[j,3]
    jumpIndexSinceGFlip[j] <- tmp_jSinceInd
  }
  repeatedSapIDs <- unique(sapIndex[duplicated(sapIndex)])
  
  if (length(repeatedSapIDs)==0) {
    cat('WARNING jumpComparisonWithinDisplay: no repeated saplings found, returning NA\n')
    jumpComparisonResult <- NA
    return(jumpComparisonResult)
  }
  
  result_sapIDs <- numeric(0)
  result_Rs <- numeric(0)*NA
  result_jumpIDs_1 <- numeric(0)
  result_jumpIDs_2 <- numeric(0)
  result_dispIDs_1 <- numeric(0)
  result_dispIDs_2 <- numeric(0)
  result_vidGroups_1 <- numeric(0)
  result_vidGroups_2 <- numeric(0)
  result_jumpSinceFlipIDs_1 <- numeric(0)*NA
  result_jumpSinceFlipIDs_2 <- numeric(0)*NA
  
  for (s in 1:length(repeatedSapIDs)) {
    # for each repeated sapIndex
    tmp_sapIndex <- repeatedSapIDs[s]
    # which jumps was it?
    jumpInds <- which(tmp_sapIndex==sapIndex) # index
    jumpIDs <- jumpIndex[jumpInds]            # value, not used, kept for debugging
    dispIDs <- dispIndex[jumpInds]            # value, not used, kept for debugging
    vidGroups <- vidGroup[jumpInds]            # value, not used, kept for debugging
    
    # how many repetitions?
    nreps <- length(jumpInds)
    combos <- combn(jumpInds,2) # generates all possible pairs; index
    npairs <- ncol(combos)
    model_list_1 <- list()
    model_list_2 <- list()
    xmin <- 0
    xmax <- 0
    for (p in 1:npairs) {
      # grab the models and the limits and get the corr coeff from the function
      model_list_1[[p]] <- modelled_jumps[[combos[1,p]]]
      model_list_2[[p]] <- modelled_jumps[[combos[2,p]]]
      
      # identify rows for this jump pair:
      jrows1 <- which(jumpsXYdf$jumpIndex==jIDdID[combos[1,p],1] & jumpsXYdf$displayNum==jIDdID[combos[1,p],2] & jumpsXYdf$vidGroup==jIDdID[combos[1,p],3])
      jrows2 <- which(jumpsXYdf$jumpIndex==jIDdID[combos[2,p],1] & jumpsXYdf$displayNum==jIDdID[combos[2,p],2] & jumpsXYdf$vidGroup==jIDdID[combos[2,p],3])
      
      # here use unique limits for each paired comparison
      xmin[p] <- min(jumpsXYdf$horizontal[c(jrows1,jrows2)])
      xmax[p] <- max(jumpsXYdf$horizontal[c(jrows1,jrows2)])
    }
    x_limits <- data.frame(xmin=xmin,xmax=xmax)
    tmp_Rs <- parabolaCorrelator(model_list_1,model_list_2,x_limits,n_values) 
    
    result_sapIDs <- c(result_sapIDs,rep(tmp_sapIndex,npairs))
    result_jumpIDs_1 <- c(result_jumpIDs_1,jumpIndex[combos[1,]])
    result_jumpIDs_2 <- c(result_jumpIDs_2,jumpIndex[combos[2,]])
    result_dispIDs_1 <- c(result_dispIDs_1,dispIndex[combos[1,]])
    result_dispIDs_2 <- c(result_dispIDs_2,dispIndex[combos[2,]])
    result_vidGroups_1 <- c(result_vidGroups_1,vidGroup[combos[1,]])
    result_vidGroups_2 <- c(result_vidGroups_2,vidGroup[combos[2,]])
    
    result_jumpSinceFlipIDs_1 <- c(result_jumpSinceFlipIDs_1,jumpIndexSinceGFlip[combos[1,]])
    result_jumpSinceFlipIDs_2 <- c(result_jumpSinceFlipIDs_2,jumpIndexSinceGFlip[combos[2,]])
    
    result_Rs <- c(result_Rs,tmp_Rs)
    
  }
  jumpComparisonResult <- data.frame('saplingCombo'=as.factor(result_sapIDs),
                                     'jumpIndex1'=result_jumpIDs_1,'jumpIndex2'=result_jumpIDs_2,
                                     'jumpSinceFlipIndex1'=result_jumpSinceFlipIDs_1,'jumpSinceFlipIndex2'=result_jumpSinceFlipIDs_2,
                                     'dispIndex1'=result_dispIDs_1,'dispIndex2'=result_dispIDs_2,
                                     'vidGroup1'=result_vidGroups_1,'vidGroup2'=result_vidGroups_2,
                                     'corrCoeffs'=result_Rs)
  levels(jumpComparisonResult$saplingCombo) <- sapLevels
  return(jumpComparisonResult)
}



jumpParamsFrom2dJumpDF <- function(projectedJumpDF,spf,maxFramesUsed=5) {
  # takes a DF with 2D positions and frame numbers for all jumps of a display and returns 
  # parameterised jumps. Required 2nd input spf (secs per frame).
  # Optional input defines the maximum number of frames used in the calculation (default: 5)
  
  jumpInds <- unique(projectedJumpDF$jumpIndex)
  njumps <- length(jumpInds)
  
  # initialise things:
  jumpIndex <-integer(njumps)*NA
  jumpSinceGFlipIndex <- integer(njumps)*NA
  sap1 <- integer(njumps)*NA
  sap2 <- integer(njumps)*NA
  frame0 <- integer(njumps)*NA
  frame1 <- integer(njumps)*NA
  TOangles <- integer(njumps)*NA
  TOvelX <- integer(njumps)*NA
  TOvelY <- integer(njumps)*NA
  TOvelXY <- integer(njumps)*NA
  TOheight <- integer(njumps)*NA
  Lheight <- integer(njumps)*NA
  maxTrackedHeight <- integer(njumps)*NA
  distanceStraightLine <- integer(njumps)*NA
  distanceCurved <- integer(njumps)*NA
  jumpDuration_frames <- integer(njumps)*NA
  for (j in 1:njumps) {
    current_jumpInd <- jumpInds[j]
    jumpIndex[j] <- current_jumpInd
    jDF <- filter(projectedJumpDF,jumpIndex==current_jumpInd) # extract rows of interest
    sap1[j] <- unique(jDF$sap1)
    sap2[j] <- unique(jDF$sap2)
    jumpDuration_frames[j] <- jDF$frameIndex[nrow(jDF)]-jDF$frameIndex[1] # don't correct by adding 1 as the first frame the bird is still perched
    jumpSinceGFlipIndex[j] <- unique(jDF$jumpIndSinceGflip)
    TOheight[j] <- jDF$vertical[1]        # take-off height is first frame (maybe bad assumption!)
    Lheight[j] <- jDF$vertical[nrow(jDF)] # landing height is last frame (maybe bad assumption!)
    maxTrackedHeight[j] <- max(jDF$vertical) # max height is simply a max (biased assumption) 
    distanceStraightLine[j] <- sqrt((jDF$vertical[nrow(jDF)]-jDF$vertical[1])^2+(jDF$horizontal[nrow(jDF)]-jDF$horizontal[1])^2) # distance between first and last points (euclidean)
    distanceCurved[j] <- sum(sqrt(diff(jDF$horizontal)^2+diff(jDF$vertical)^2)) # sum of paired distances (euclidean)
    f0 <- jDF$frameIndex[2] # use the first tracked point in the air, not the one at rest
    FOI <- which((jDF$frameIndex-f0)<maxFramesUsed & (jDF$frameIndex-f0)>0) # use the next possible point to calculate take-off angle+vel, unless it's too far away in time
    if (length(FOI)==0) { # no next inital frame within the maximum time window
      next
    }
    f1 <- jDF$frameIndex[min(FOI)] # use next possible point
    dt <- (f1-f0)*spf
    if (dt==0) { # not enough frames, skip it
      next
    }
    frame0[j] <- f0
    frame1[j] <- f1
    dy <- jDF$vertical[jDF$frameIndex==f1]-jDF$vertical[jDF$frameIndex==f0]
    dx <- jDF$horizontal[jDF$frameIndex==f1]-jDF$horizontal[jDF$frameIndex==f0]
    TOangles[j] <- atan2(dy,dx)*180/pi    # take-off angle in degrees
    TOvelX[j] <- dx/dt                    # horizontal take-off velocity
    TOvelY[j] <- dy/dt                    # vertical take-off velocity
    TOvelXY[j] <- sqrt(TOvelX[j]^2+TOvelY[j]^2)      # take-off velocity using Pythagorus
  }
  jumpParamsDF <- data.frame(jumpIndex,jumpSinceGFlipIndex,sap1,sap2,frame0,frame1,
                             TOangles,TOvelX,TOvelY,TOvelXY,
                             TOheight,Lheight,maxTrackedHeight,
                             distanceStraightLine,distanceCurved,jumpDuration_frames) 
  return(jumpParamsDF)
}

# Most common mistake while writing this was to convert angles to degrees and forget to convert back if using trig functions again. Double and triple-checked!
jumpParamsFromFittedModels <- function(modelled_jumps,gravity) {
  # takes a DF with fitted models for all jumps of a display and returns 
  # parameterised jumps. Required 2nd input gravity (value of force of gravity).
  nmodels <- length(modelled_jumps)
  # initialise things:
  TOangles <- integer(nmodels)*NA
  TOvelXY <- integer(nmodels)*NA
  TOheight <- integer(nmodels)*NA
  Lheight <- integer(nmodels)*NA
  maxHeightFit <- integer(nmodels)*NA
  distanceStraightLine <- integer(nmodels)*NA
  distanceCurved <- integer(nmodels)*NA
  
  for (j in 1:nmodels) {
    # first check whether the model really worked
    tmp_mod <- modelled_jumps[[j]]
    if (length(tmp_mod)==1 && is.na(tmp_mod)) { # probably there were too few data points so fitting wasn't attempted
      next # do nothing, all entries will be NA
    }
    TOheight[j] <- tmp_mod$fitted.values[1]   # take-off height is first fitted y value
    Lheight[j] <- tmp_mod$fitted.values[length(tmp_mod$fitted.values)] # landing height is last fitted y value

    coeffs <- as.numeric(tmp_mod$coefficients) # 1:b0 intercept, 1:b1 linear, 3:b2 squared term
    # coefficient of linear x term is tan of take-off angle
    TOangles[j] <- atan(coeffs[2])*180/pi     
    #  note: single argument arctan is fine here as we projected all jumps to be from left to right
    
    # x-coord of vertex (highest point of jump) is -b1/2b2 [linear term / sq term]
    xVertx <- -coeffs[2]/(2*coeffs[3])
    maxHeightFit[j] <- predict(tmp_mod,newdata = data.frame(horizontal=xVertx),raw=T) # same as coeffs[1]+coeffs[2]*xVertx+coeffs[3]*(xVertx^2)
    
    tmp_x <- tmp_mod$x[,2]
    tmp_y <- tmp_mod$fitted.values
    distanceStraightLine[j] <- sqrt((tmp_x[length(tmp_x)]-tmp_x[1])^2+(tmp_y[length(tmp_y)]-tmp_y[1])^2) # distance between first and last points (euclidean)
    
    # coefficient of x2 term is -g/(2 * v^2 cos^2angle) => v = sqrt(-g*(1+b1^2)/2*b2) using either trig equality of sec and tan, or cos(atanx) substitution
    TOvelXY[j] <- sqrt( -(gravity*(1+coeffs[2]^2)) / (2*coeffs[3]) )
    # can also use sqrt (-gravity   / 2 * b2 * cos sq theta REMEMBER TO RADIANISE IT AGAIN!), equivalent anyway, so leave it out
    # i.e.  sqrt(-gravity/(2 * coeffs[3] * ((cos(TOangles[j]*pi/180))^2)))
    
    
    # Distance is the integral of the square root of 1 + the square of the derivative of the parabola! 
    # Parabola is ax2 + bx + c , therefore 2ax + b is the derivative 
    # e.g. https://www.math24.net/arc-length/ https://www.intmath.com/applications-integration/11-arc-length-curve.php
    
    integralResult <- integrate(function(x,a,b,c) {sqrt(1+(2*a*x+b)^2)},lower=tmp_x[1],upper=tmp_x[length(tmp_x)],a=tmp_mod$coefficients[[3]], b=tmp_mod$coefficients[[2]])
    distanceCurved[j] <- integralResult[[1]]
    }
  
  jumpParamsModel <- data.frame(TOangles,TOvelXY,  
                             TOheight,Lheight,maxHeightFit,
                             distanceStraightLine,distanceCurved) 
  return(jumpParamsModel)
}