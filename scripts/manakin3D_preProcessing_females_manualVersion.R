# This script performs preprocessing for all available manual 3d tracking files
# from non-males (females or maybe juveniles; cannot visually distinguish them!)
#
# cliodhna.quigley@vetmeduni.ac.at, 2020-2021

FDIR_females = 'data/manualdata/females/' # where to find the csv files

dataDir_2D = 'data/processedData/2Dprojections_femaleMan/' # where to save data
dataDir_3D = 'data/processedData/3Ddata_femaleMan/'

library(plotly)
library(dplyr)
library(ggplot2)

source('scripts/3d_utils.R') # contains functions to calculate distance and speed

DEBUG = 0 # 1=make extra plots to check what's going on.  2= additionally plot extra stuff
          # note: plots won't be displayed if run inside of for loop, so you'll need to manually run per file
          # or add code to save the plot to a file in the for loop

# Format: each row is an xy coordinate with corresponding filename and sapling index (number or G for ground)
mapTable <- 'data/female_saplings.csv'

# parameters:
fps <- 60 # frame rate assumed fixed
spf <- 1/fps # seconds per frame
speed_threshold <- 7 # this is m/s, to exclude outliers
sitting_threshold_man <- 0.6 # to decide when the bird is not moving
n_values_corr <- 50 # number of x values for which we correlate y values of fitted parabola  when comparing annotations

# FIRST STEP: which files do we have? 
man_files <- getDisplayVidCodesFromCSVDir(FDIR_females) # this is a character vector

# for each file:
for (m in 1:length(man_files)) {
  current_filecode <- man_files[m]
  cat('Processing',current_filecode,'\n')
  
  # # #
  # 1: load data
  manakin3D_man <- findLoadcsvJJ(FDIR_females,current_filecode)
  if (any(is.na(manakin3D_man))) {
    cat('PROBLEM WITH MANUAL ANNOTATION FILE; SKIPPING\n')
    next # skip this file for now, there's maybe more than one manual annotation 
  }

  # load sapling data
  saplings_xy <- getSaplingGroundXYLocations_code(current_filecode,mapTable)
  if (nrow(saplings_xy)==0) {
    cat('ERROR, NO SAPLINGS FOUND in table, skipping\n')
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


  # Get speed of bird per tracked frame, using assumption that it's 60 fps video
  manakin3D_speed_man <- speedFromXYZ(manakin3D_man,spf) # adds new columns to dataframe for distance and speed
  if (all(is.na(manakin3D_speed_man))) {
    cat('problem calculating speed; skipping\n')
    next
  }

  if (DEBUG==1) {
    # CHECK SPEED DISTRIBUTION, manual annotation so no cleaning done
    man_speeds <- manakin3D_speed_man$speed3d
    mytitle = paste(birdcode,current_filecode,'manual, female')
    hist(man_speeds,seq(from=0,to=ceiling(max(man_speeds,na.rm=T)),by=0.01),xlim=c(0,1),main=mytitle)
  }
  
  # Create new frame numbering column to start with 1, for easier readibility
  manakin3D_speed_man$frame_number_readable <- manakin3D_speed_man$frame_number-manakin3D_speed_man$frame_number[1]+1
  
  # # #
  # 4: mark resting frames and assign them to a sapling
  #
  manakin3D_speed_man <- assignLowSpeedsToSaplings(manakin3D_speed_man,sitting_threshold_man,saplings_xy)
  if (all(is.na(manakin3D_speed_man))) { # clustering failed 
    cat('clustering of resting frames to saplings failed; skipping\n')
    next
  } # can't try and rescue it by removing G as females/juveniles don't 'jump' to ground
  
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

  # # #
  # 5: derive trajectories
  #
  # extract 3D trajectories
  trajDF_man <- extractTrajectories(manakin3D_speed_man)  # data frame now includes jumpIndex, startSapling, stopSapling columns
  trajDF_man <- addJumpIndexSinceG(trajDF_man)
  # turn into 2D trajectories: extract plane described by start and end points
  # of jump, project all jump locations onto this plane
  projectedJumpDF_man <- flattenJumps(trajDF_man)
  if (DEBUG==1) { # plot them
    g <- ggplot(projectedJumpDF_man,aes(horizontal,vertical,color=as.factor(jumpIndex))) + geom_point()  + coord_fixed()
    g <- g + facet_grid(rows = vars(sap1), cols = vars(sap2), drop=F) +
      ggtitle(paste(current_filecode,'manual (rows: origin; cols: dest), Female')) + 
      xlab('Horizontal distance from jump origin (m)') +
      ylab('Height above ground (m)') +
      labs(color='Jump Index')
    g
  }
  
  # # #
  # 6: save the data in appropriate dataDir
  fname_chunk <- paste0(birdcode,'_',current_filecode,'_')
  saveRDS(projectedJumpDF_man, file = paste0(dataDir_2D,fname_chunk,'man.Rds'))
  saveRDS(manakin3D_speed_man, file = paste0(dataDir_3D,fname_chunk,'clean_man.Rds'))
  cat('data saved successfully\n')
}