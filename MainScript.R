### Overarchign script, tying everythign together

### Let's get all the libraries in

library(dplyr)
library(pool)
library(postGIStools)
library(sf)
library(getPass)
library(rgbif)
library(dplyr) # for data-wrangling
library(mapedit)
library(rio)
library(raster)
library(FNN)
library(readr)
library(ggplot2)
library(greta)
library(stringr)
source("ExtraFunctions.R")


### Script 1 downloads everything for you. 
species_name <- 'Esox lucius'

# All the data will be transferred to one file. Create that file if it hasn't already been made.
if (dir.exists(paste0("./Data/",gsub(' ','_',species_name))) == FALSE
    ) {dir.create(paste0("./Data/",gsub(' ','_',species_name)))}


#   If the download has already been initiated, set the initiate_download value to FALSE. 
#   If you have already downloaded lake data and it is sitting in your Data folder, set 
#   download_lakes to FALSE. If you have already downloaded all lake location data in 
#   order to create absences, set all_lakes to FALSE. If you want to get rid of all lakes
#   from Nordland, Troms and Finnmark, set delete_north to FALSE.
#   If the species does not have a native range in Norway, select native_range to FALSE.

initiate_download <- FALSE
download_lakes <- FALSE
all_lakes <-  FALSE
delete_north <-  TRUE
download_native_range <- TRUE

#   dist_threshold sets the distance between a point and its designated lake which
#   is acceptable.
dist_threshold <- 50

source("./R/1_GetIntroductions.R")

# Next step is to add in all the covariates.

source("./R/2_BioticDataAddition.R")

# Now we run the preliminary model. Only thing we need to choose is what our size limit 
# on lakes will be.

size_threshold <- 0.02

# This script will take a while, as it's running a model on up to 250,000 lakes. Grab a 
# coffee. Teach it to do algebra.
source("./R/3_FullModelConstruct.R")

# Now we run a second model, using buffered data
source("./R/4_BufferedModelConstruct.R")

# The output from this is the same as the second, except that you get a map as well
# showing what the buffered region looks like.
plot(comb_bif)

# Next one gives you uncertainty from each lake, convergence diagnostic,
# beta intervals, and model deviance.
# 
use_buffered_model <- FALSE
source("./R/5_ModelAnalysis.R")

# Now we can run our simulations. need to define which model we want to use.
# n_loops simply tells us how many iterations we want to run. Be careful,
# because the run time can be enormous.
n_loops <- 100

source("./R/6_Looping.R")

# Last step is simply doing some data visualisation.

prob_threshold <- 0.5
axis_limit <- 150

source("./R/7_Visualisations.R")




maps$predicted_appearances$layers



