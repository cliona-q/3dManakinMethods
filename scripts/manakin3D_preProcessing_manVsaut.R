# This script performs preprocessing for later comparison of manual and automated 3d tracking 
#
# cliodhna.quigley@vetmeduni.ac.at, 2020-2021

autManStartDiffThresh = 20 # threshold for allowed difference between manual and automated start frames; more then we cut

FDIR_man = 'data/manualdata/' # where to find the csv files
FDIR_aut = 'data/automaticdata/'

dataDir_2D = 'data/processedData/2Dprojections_manVsaut/' # where to save data
dataDir_3D = 'data/processedData/3Ddata_manVsaut/'

library(plotly)
library(dplyr)
library(ggplot2)

source('scripts/3d_utils.R') # contains functions to calculate distance and speed

DEBUG = 0 # make extra plots to check what's going on; 1 is minimum plots; 2 is everything
          # note: plots won't be displayed if run inside of for loop, so you'll need to manually run per file
          # or add code to save the plot to a file in the for loop

# Format: each row is an xy coordinate with corresponding filename and sapling index (number or G for ground)
mapTable <- 'data/methodsPaper_videoMetadata.csv'

# parameters:
fps <- 60 # frame rate assumed fixed
spf <- 1/fps # seconds per frame
speed_threshold <- 7 #  this is m/s, to exclude outliers; chosen using visual inspection of histogram
sitting_threshold <- 0.6 # this is m/s, to classify position as sitting / jumping; chosen using visual inspection of histogram

n_values_corr <- 50 # number of x values for which we correlate y values of fitted parabolas when comparing annotations

# FIRST STEP: which files do we have? The limiting factor is the manual annotations, we
# can always get new automated ones by running the model
man_files <- getDisplayVidCodesFromCSVDir(FDIR_man) # this is a character vector

# for each manual file:
for (m in 1:length(man_files)) {
  current_filecode <- man_files[m]
  cat('Processing',current_filecode,'\n')
  
  # # #
  # 1: load data
  # check there is an automated file by trying to load it
  manakin3D_aut <- findLoadcsvJJ(FDIR_aut,current_filecode)
  if (any(is.na(manakin3D_aut))) {
    cat('PROBLEM WITH AUTOMATIC ANNOTATION FILE; SKIPPING\n')
    next # skip this file for now, there's no automatic annotation 
  }
  # and load the manual one too
  manakin3D_man <- findLoadcsvJJ(FDIR_man,current_filecode)
  if (any(is.na(manakin3D_man))) {
    cat('PROBLEM WITH MANUAL ANNOTATION FILE; SKIPPING\n')
    next # skip this file for now, there's maybe more than one manual annotation 
  }
  # restrict the analysis to the frame limits given in the manual data table
  min_frame <- manakin3D_man$frame_number[1] # first manual frame is first interesting one
  max_frame <- manakin3D_man$frame_number[nrow(manakin3D_man)] # last manual frame is last interesting one
  manakin3D_aut <- manakin3D_aut[which(manakin3D_aut$frame_number>=min_frame & manakin3D_aut$frame_number<=max_frame),]
  
  # for one of the files, we have a later start in the automatic than the manual
  if (abs(manakin3D_aut$frame_number[1]-manakin3D_man$frame_number[1])>autManStartDiffThresh) {
    # ignore the frames of manual before the start of the automated tracking in these cases
    cat('AUT starts after MAN so cutting the start of AUT\n')
    min_frame <- manakin3D_aut$frame_number[1] # first automated frame is first interesting one
    manakin3D_man <- manakin3D_man[which(manakin3D_man$frame_number>=min_frame),]
  }
  
  # load sapling data
  saplings_xy <- getSaplingGroundXYLocations_code(current_filecode,mapTable)
  if (nrow(saplings_xy)==0) {
    cat('ERROR, NO SAPLINGS FOUND IN TABLE, skipping file\n')
    next
  }
  # and grab bird code for data saving
  birdcode <- unique(as.character(saplings_xy$ID))
  
  # # #
  # 2: plot if debugging
  if (DEBUG==2) {
    p <- plot_ly(manakin3D_man, x=~x, y=~y, z=~z, type = 'scatter3d', mode='lines+markers',
                 opacity=0.5)
    p %>% layout(title=paste('Manual annotation',current_filecode))
  }
  if (DEBUG==2) {
    p2 <- plot_ly(manakin3D_aut, x=~x, y=~y, z=~z, type = 'scatter3d', mode='lines+markers',
                  opacity=0.5)
    p2 %>% layout(title=paste('Automatic annotation',current_filecode))
  }

  # Exclude position estimates with distances that are 0 
  manakin3D_aut <- cleanXYZdata(manakin3D_aut)
  manakin3D_man <- cleanXYZdata(manakin3D_man)
  # Most of these frames are excluded on import by ignoring frames with ncameras == 0  
  # There are sometimes still subsequent frames with distance of 0, and if they are not
  # rejected they will lead to false positive classifications of landings during a jump
  
  # # #
  # 3: deal with outliers
  # Uses a speed threshold
  # Get speed of bird per tracked frame, using assumption that it's 60 fps video
  manakin3D_speed_man <- speedFromXYZ(manakin3D_man,spf) # adds new columns to dataframe for distance and speed
  if (all(is.na(manakin3D_speed_man))) {
    cat('problem calculating speed for manual file; skipping\n')
    next
  }
  manakin3D_speed_aut <- speedFromXYZ(manakin3D_aut,spf) # adds new columns to dataframe for distance and speed
  if (DEBUG==2) {
    # CHECK SPEED DISTRIBUTION, pre-cleaning
    man_speeds <- manakin3D_speed_man$speed3d
    aut_speeds <-  manakin3D_speed_aut$speed3d
    hist(aut_speeds,seq(from=0,to=ceiling(max(c(aut_speeds,man_speeds),na.rm=T)),by=0.01))
    hist(man_speeds,seq(from=0,to=ceiling(max(c(aut_speeds,man_speeds),na.rm=T)),by=0.01))
  }
  # Threshold speed, rejecting any frames containing speeds above the threshold (and recomputing speed)
  manakin3D_speed_aut <- rejectHighSpeeds(manakin3D_speed_aut,speed_threshold,spf)
  if (DEBUG==2) {
    # CHECK SPEED DISTRIBUTION, post-cleaning
    man_speeds <- manakin3D_speed_man$speed3d
    aut_speeds <-  manakin3D_speed_aut$speed3d
    hist(aut_speeds,seq(from=0,to=ceiling(max(c(aut_speeds,man_speeds),na.rm=T)),by=0.01))
    hist(man_speeds,seq(from=0,to=ceiling(max(c(aut_speeds,man_speeds),na.rm=T)),by=0.01))
  }
  if (DEBUG==1) {
    # CHECK SPEED DISTRIBUTION, post-cleaning, zoom in to region where perching is assigned
    man_speeds <- manakin3D_speed_man$speed3d
    aut_speeds <-  manakin3D_speed_aut$speed3d
    
    mytitle = paste(birdcode,current_filecode,'automated')
    hist(aut_speeds,seq(from=0,to=ceiling(max(c(aut_speeds,man_speeds),na.rm=T)),by=0.01),xlim=c(0,1),main=mytitle)
    
    mytitle = paste(birdcode,current_filecode,'manual')
    hist(man_speeds,seq(from=0,to=ceiling(max(c(aut_speeds,man_speeds),na.rm=T)),by=0.01),xlim=c(0,1),main=mytitle)
    
  }
  
  
  
  # Create new frame numbering column to start with 1, for easier readibility
  # BUT, ensure that frame numbering remains matched across manual / annotated files
  min_frame <- min(manakin3D_speed_aut$frame_number[1],manakin3D_speed_man$frame_number[1])
  manakin3D_speed_aut$frame_number_readable <- manakin3D_speed_aut$frame_number-min_frame+1
  manakin3D_speed_man$frame_number_readable <- manakin3D_speed_man$frame_number-min_frame+1
  
  # # #
  # 4: mark resting frames and assign them to a sapling
  #
  manakin3D_speed_aut_cl <- assignLowSpeedsToSaplings(manakin3D_speed_aut,sitting_threshold,saplings_xy)
  if (all(is.na(manakin3D_speed_aut_cl))) { # clustering failed 
    cat('clustering of resting frames to saplings failed for automated file\n')
    cat('trying without G sapling\n')
    Gind <- which(saplings_xy$sapID=='G')
    saplings_xy_noG <- saplings_xy[-Gind,]
    
    manakin3D_speed_aut_cl <- assignLowSpeedsToSaplings(manakin3D_speed_aut,sitting_threshold,saplings_xy_noG)
    if (all(is.na(manakin3D_speed_aut_cl))) { # clustering failed 
      cat('clustering failed again; skipping\n')
      next
    }
  }
  manakin3D_speed_aut <- manakin3D_speed_aut_cl
  manakin3D_speed_man <- assignLowSpeedsToSaplings(manakin3D_speed_man,sitting_threshold,saplings_xy)
  if (all(is.na(manakin3D_speed_man))) { # clustering failed 
    cat('clustering of resting frames to saplings failed for manual file; skipping\n')
    next
  }
  
  if (DEBUG==2) {
    p <- plot_ly(manakin3D_speed_man, x=~x, y=~y, z=~z, type = 'scatter3d', mode='lines+markers',
                 opacity=0.5,
                 marker=list(size = ~speed3d, sizeref = 0.02, sizemode = 'area'),
                 text = ~frame_number_readable,
                 hovertemplate = paste('<i>Speed</i>: %{marker.size:,.2f}m/s',
                                       '<br><b>Framenum</b>: %{text}<br>')
    )
    p %>% layout(title=paste('Manual, clean, speed',current_filecode))
  }
  if (DEBUG==2) {
    p2 <- plot_ly(manakin3D_speed_aut, x=~x, y=~y, z=~z, type = 'scatter3d', mode='lines+markers',
                  opacity=0.5,
                  marker=list(size = ~speed3d, sizeref = 0.02, sizemode = 'area'),
                  text = ~frame_number_readable,
                  hovertemplate = paste('<i>Speed</i>: %{marker.size:,.2f}m/s',
                                        '<br><b>Framenum</b>: %{text}<br>'))
    p2 %>% layout(title=paste('Automatic, clean, speed',current_filecode))
  }
  
  # # #
  # 5: derive trajectories
  #
  # extract 3D trajectories
  trajDF_man <- extractTrajectories(manakin3D_speed_man)  # data frame now includes jumpIndex, startSapling, stopSapling columns
  trajDF_man <- addJumpIndexSinceG(trajDF_man)
  trajDF_aut <- extractTrajectories(manakin3D_speed_aut)  # data frame now includes jumpIndex, startSapling, stopSapling columns
  trajDF_aut <- addJumpIndexSinceG(trajDF_aut)
  # turn into 2D trajectories: extract plane described by start and end points
  # of jump, project all jump locations onto this plane
  projectedJumpDF_man <- flattenJumps(trajDF_man)
  projectedJumpDF_aut <- flattenJumps(trajDF_aut)
  if (DEBUG==1) { # plot them 
    g <- ggplot(projectedJumpDF_aut,aes(horizontal,vertical,color=as.factor(jumpIndex))) + geom_point() + coord_fixed()
    g <- g + facet_grid(rows = vars(sap1), cols = vars(sap2), drop=F) +
      ggtitle(paste(current_filecode,'automatic (rows: origin; cols: dest)')) + 
      xlab('Horizontal distance from jump origin (m)') +
      ylab('Height above ground (m)') +
      labs(color='Jump Index')
    g
  }
  if (DEBUG==1) { # plot them
    g <- ggplot(projectedJumpDF_man,aes(horizontal,vertical,color=as.factor(jumpIndex))) + geom_point() + coord_fixed()
    g <- g + facet_grid(rows = vars(sap1), cols = vars(sap2), drop=F) +
      ggtitle(paste(current_filecode,'manual (rows: origin; cols: dest)')) + 
      xlab('Horizontal distance from jump origin (m)') +
      ylab('Height above ground (m)') +
      labs(color='Jump Index')
    g
  }
 
  # # #
  # 6: save the data in appropriate dataDir
  fname_chunk <- paste0(birdcode,'_',current_filecode,'_')
  saveRDS(projectedJumpDF_aut, file = paste0(dataDir_2D,fname_chunk,'aut.Rds'))
  saveRDS(projectedJumpDF_man, file = paste0(dataDir_2D,fname_chunk,'man.Rds'))
  saveRDS(manakin3D_speed_aut, file = paste0(dataDir_3D,fname_chunk,'clean_aut.Rds'))
  saveRDS(manakin3D_speed_man, file = paste0(dataDir_3D,fname_chunk,'clean_man.Rds'))
  saveRDS(manakin3D_aut, file = paste0(dataDir_3D,fname_chunk,'noSpeedExclusion_aut.Rds'))
  cat('data saved successfully\n')
}

