# script to plot manakin 3d jumps, in order to check automated tracking of display
# Note:  in rstudio, some plots will be in viewer tab and some in plots tab
#
# cliodhna.quigley@vetmeduni.ac.at, 2020-2021

# ENTER THE NAME OF THE FILE YOU WANT TO PLOT HERE:
# it needs to be BIRDCODE_DATE_TIME_displN
current_file = 'CAN04_20180220_160809_displ13_5'


dataDir_2D = 'data/processedData/2Dprojections_aut/' # where to find data
dataDir_3D = 'data/processedData/3Ddata_aut/'

sitting_threshold <- 0.6 # 0.6 is final threshold

library(plotly)
library(dplyr)
library(ggplot2)


# load the data
projectedJumpDF_aut<-readRDS(file = paste0(dataDir_2D,current_file,'_aut.Rds'))
manakin3D_speed_aut<-readRDS(file = paste0(dataDir_3D,current_file,'_clean_aut.Rds'))
manakin3D_aut<-readRDS(file = paste0(dataDir_3D,current_file,'_noSpeedExclusion_aut.Rds'))

# plots:
# 3d plot of automated data before speed threshold is implemented:
p2 <- plot_ly(manakin3D_aut, x=~x, y=~y, z=~z, type = 'scatter3d', mode='lines+markers',
                opacity=0.5)
p2 %>% layout(title=paste('Automatic annotation',current_file,'before speed threshold'))

# histogram of remaining speeds, post-cleaning
aut_speeds <-  manakin3D_speed_aut$speed3d
hist(aut_speeds,seq(from=0,to=ceiling(max(aut_speeds,na.rm=T)),by=0.01))

# 3d plot of automated data after speed threshold is implemented:
p2 <- plot_ly(manakin3D_speed_aut, x=~x, y=~y, z=~z, type = 'scatter3d', mode='lines+markers',
              opacity=0.5)
p2 %>% layout(title=paste('Automatic annotation',current_file,'after speed threshold'))

# and again with speed indicated by marker size
p2 <- plot_ly(manakin3D_speed_aut, x=~x, y=~y, z=~z, type = 'scatter3d', mode='lines+markers',
                opacity=0.5,
                marker=list(size = ~speed3d, sizeref = 0.02, sizemode = 'area'),
                text = ~frame_number_readable,
                hovertemplate = paste('<i>Speed</i>: %{marker.size:,.2f}m/s',
                                      '<br><b>Framenum</b>: %{text}<br>'))
p2 %>% layout(title=paste('Automatic, clean, speed is dot size',current_file))

# and again with different colours for jump and sit
colinds_tmp <- manakin3D_speed_aut$speed3d<sitting_threshold
colinds_tmp <- colinds_tmp + 1
cols <- c("red", "black")
p2 <- plot_ly(manakin3D_speed_aut, x=~x, y=~y, z=~z, type = 'scatter3d', mode='lines+markers',
              marker=list(color = colinds_tmp, colours = cols, opacity=0.5),
              text = ~frame_number_readable)
p2 %>% layout(title=paste('Automatic, clean, colour is jump/sit',current_file))


# projected jumps
g <- ggplot(projectedJumpDF_aut,aes(horizontal,vertical,color=as.factor(jumpIndex))) + geom_point() + coord_fixed()
g <- g + facet_grid(rows = vars(sap1), cols = vars(sap2), drop=F) +
    ggtitle(paste(current_file,'automatic (rows: origin; cols: dest)')) + 
    xlab('Horizontal distance from jump origin (m)') +
    ylab('Height above ground (m)') +
    labs(color='Jump Index')
g
