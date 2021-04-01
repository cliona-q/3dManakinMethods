# 3dManakinMethods

Code and some sample data to derive jump trajectory parameters 
from 3d tracking data of courting golden-collared manakins.

This project is licensed under the terms of the GNU GPLv3 license.

Project is published in the journal Ethology as Janisch, Perinot, 
Fusani & Quigley (2021) "Deciphering choreographies of elaborate 
courtship displays of golden-collared manakins using markerless motion capture."

To use the code, it's easiest to open the R project in Rstudio. You 
must run the preprocessing files before running fitting/analysis 
scripts. 

You must run the male analysis file (manakin3D_fittingEtc.R) before
running the female analysis file (manakin3D_femaleAnalysis.R).

Sample data is provided in the data/ directory; empty sub 
directories are intentional and will be used to store processed data 
after you run the preprocessing scripts.

Warnings during execution are normal (and informative!); errors should 
not occur. All scripts are thoroughly documented; see the publication 
for more info on the method. A plotting function is also provided to give
some hands-on insight into the method.