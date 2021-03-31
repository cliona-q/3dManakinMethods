# This script performs preprocessing for all available automated 3d tracking files 
#
# cliodhna.quigley@vetmeduni.ac.at, 2020-2021

FDIR = 'data/automaticdata/' # where to look for the data

dataDir_2D = 'data/processedData/2Dprojections_aut/' # where to save data
dataDir_3D = 'data/processedData/3Ddata_aut/'

library(plotly)
library(dplyr)
library(ggplot2)

source('scripts/3d_utils.R') # contains functions to calculate distance and speed

# Format: each row is an xy coordinate with corresponding filename and sapling index (number or G for ground)
mapTable <- 'data/methodsPaper_videoMetadata.csv'

# parameters:
fps <- 60 # frame rate assumed fixed
spf <- 1/fps # seconds per frame
speed_threshold <- 7 # this is m/s, to exclude outliers
sitting_threshold <- 0.6 # to decide when the bird is not moving

# FIRST STEP: which files do we have?
aut_files <- getDisplayVidCodesFromCSVDir(FDIR) 

DEBUG = 0 # make some plots if 1, make way more plots but not saved if 2
          # note: plots won't be displayed if run inside of for loop, so you'll need to manually run per file
          # or add code to save the plot to a file in the for loop

# for each file:
for (f in 1:length(aut_files)) {
  current_filecode <- aut_files[f]
  cat('Processing',current_filecode,'\n')
  
  # # #
  # 1: load data
  manakin3D_aut <- findLoadcsvJJ(FDIR,current_filecode)
  if (any(is.na(manakin3D_aut))) {
    cat('PROBLEM LOADING AUTOMATIC ANNOTATION FILE; SKIPPING\n')
    next # skip this file for now, there's no automatic annotation 
  }
  # load sapling data
  saplings_xy <- getSaplingGroundXYLocations_code(current_filecode,mapTable)
  if (nrow(saplings_xy)==0) {
    cat('ERROR, NO SAPLINGS FOUND, skipping file\n')
    next
  }
  # and grab bird code and vidGroup for data saving
  birdcode <- unique(as.character(saplings_xy$ID))
  vidGroup <- unique(as.integer(saplings_xy$vidGroup))
  
  # # #
  # 3: clean data 
  # Exclude position estimates with distances that are 0 
  manakin3D_aut <- cleanXYZdata(manakin3D_aut)
  
  # Get speed of bird per tracked frame, using assumption that it's 60 fps video
  manakin3D_speed_aut <- speedFromXYZ(manakin3D_aut,spf) # adds new columns to dataframe for distance and speed
  # Threshold speed, rejecting any frames containing speeds above the threshold (and recomputing speed)
  manakin3D_speed_aut <- rejectHighSpeeds(manakin3D_speed_aut,speed_threshold,spf)
  # Reset frame numbering to start with 1, for easier readibility
  manakin3D_speed_aut$frame_number_readable <- manakin3D_speed_aut$frame_number-manakin3D_speed_aut$frame_number[1]+1

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

  # # #
  # 5: derive trajectories
  #
  # extract 3D trajectories
  trajDF_aut <- extractTrajectories(manakin3D_speed_aut)  # data frame now includes jumpIndex, startSapling, stopSapling columns
  trajDF_aut <- addJumpIndexSinceG(trajDF_aut)
  # turn into 2D trajectories: extract plane described by start and end points
  # of jump, project all jump locations onto this plane
  projectedJumpDF_aut <- flattenJumps(trajDF_aut)

  if (DEBUG==2) {
    p2 <- plot_ly(manakin3D_speed_aut, x=~x, y=~y, z=~z, type = 'scatter3d', mode='lines+markers',
                  opacity=0.5,
                  marker=list(size = ~speed3d, sizeref = 0.02, sizemode = 'area'),
                  text = ~frame_number_readable,
                  hovertemplate = paste('<i>Speed</i>: %{marker.size:,.2f}m/s',
                                        '<br><b>Framenum</b>: %{text}<br>'))
    p2 %>% layout(title=paste('Automatic, clean, speed',current_filecode))
  }
  if (DEBUG==1) { # plot them
    g <- ggplot(projectedJumpDF_aut,aes(horizontal,vertical,color=as.factor(jumpIndex))) + geom_point() + coord_fixed(ratio=1)
    g <- g + facet_grid(rows = vars(sap1), cols = vars(sap2), drop=F) +
      ggtitle(paste(current_filecode,'automatic (rows: origin; cols: dest)')) + 
      xlab('Horizontal distance from jump origin (m)') +
      ylab('Height above ground (m)') +
      labs(color='Jump Index')
    g
  }
  
  # # #
  # 6: save the data in appropriate dataDir
  fname_chunk <- paste0(birdcode,'_',current_filecode,'_',vidGroup,'_')
  saveRDS(projectedJumpDF_aut, file = paste0(dataDir_2D,fname_chunk,'aut.Rds'))
  saveRDS(manakin3D_speed_aut, file = paste0(dataDir_3D,fname_chunk,'clean_aut.Rds'))
  saveRDS(manakin3D_aut, file = paste0(dataDir_3D,fname_chunk,'noSpeedExclusion_aut.Rds'))
  cat('data saved successfully\n\n\n')
}

