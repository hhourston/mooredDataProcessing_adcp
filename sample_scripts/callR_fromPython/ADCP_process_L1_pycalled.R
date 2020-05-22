# This script was adapted by Hana Hourston (IOS) from one written by Emily Chisholm (BIO)

# Inputs: a raw ADCP file (.000, .664, or .pd0) and a csv file containing metadata
# Corrects for magnetic declination
# Flags bad leading and trailing ensembles
# Creates a netCDF file

# WARNING: Use caution when working with ADCPs missing a pressure sensor. If the serial number is known
# for the instrument, the safest choice is to add it to line 82 so that pressure values may be calculated

# Adaptations to this code to make it executable from Python come from https://www.kdnuggets.com/2015/10/integrating-python-r-executing-part2.html

library(ADCP) 
library(ncdf4)
library(oce)
library(tools)
library(english)

################################################

mean_orientation <- function(orientation){
  upward = 0
  downward = 0
  for (i in 1:length(orientation)){
    if (orientation[i] == 'upward'){
      upward = upward + 1
    } else if (orientation[i] == 'downward'){
      downward = downward + 1
    } else {
      stop(paste('Invalid orientation value for ensemble', i))
    }
  }
  if (upward > downward){
    return('upward')
  } else {
    return('downward')
  }
}


# Use to calculate mode of pressure and determine if pressure sensor missing
getmode <- function(v){
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


########################################

# Make the R script executable from Python:

# Fetch command line arguments
myArgs <- commandArgs(trailingOnly = TRUE)

# Assign paths to the raw ADCP file and its corresponding csv metadata file
wd <- myArgs[1]
f <- myArgs[2]
meta <- myArgs[3]

#print(wd); print(f); print(meta)

setwd(wd)


# Create standard name from raw file name (WITH data type and processing level add-on)
ncname_L <- paste0(substr(basename(f), 1, nchar(basename(f))-4), '.adcp.L1')
# ncname below is for use in navigating to the directory of the file name
ncname <- substr(basename(f), 1, nchar(basename(f))-4)


# Read in the data and metadata
adp <- read.adp.easy(f, meta)


# Take average orientation since it is now vector with an orientation for every ensemble
adp <- oceSetMetadata(adp, 'orientation', mean_orientation(adp[['orientation']]))


# Convert numeric metadata items from string to numeric
adp <- oceSetMetadata(adp, 'country_institute_code', as.numeric(adp[['country_institute_code']]))
adp <- oceSetMetadata(adp, 'instrument_depth', as.numeric(adp[['instrument_depth']]))
adp <- oceSetMetadata(adp, 'latitude', as.numeric(adp[['latitude']]))
adp <- oceSetMetadata(adp, 'longitude', as.numeric(adp[['longitude']]))
adp <- oceSetMetadata(adp, 'processing_level', as.numeric(adp[['processing_level']]))
adp <- oceSetMetadata(adp, 'serialNumber', toString(adp[['serialNumber']])) # since some numbers have a leading 0
adp <- oceSetMetadata(adp, 'station_number', as.numeric(adp[['station_number']]))
adp <- oceSetMetadata(adp, 'water_depth', as.numeric(adp[['water_depth']]))


#For serial numbers starting with 0 so they have leading zero added back in
if (nchar(adp[['serialNumber']]) == 3){
  adp <- oceSetMetadata(adp, 'serialNumber', paste0('0', adp[['serialNumber']]))
}


# If missing pressure sensor, calculate pressure based on static (instrument) depth
# In these cases pressure may be all zeros, or it may have some NAs, or one very large negative or positive value with the rest zeros
# Other anomalour value combinations may occur with missing pressure sensors
# Use sum() to find the number of occurrences of zero values
if (sum(adp[['pressure']] == 0) > length(adp[['pressure']])/2 | adp[['instrumentSubtype']] == 'BroadBand' | adp[['serialNumber']] == '1588' | adp[['serialNumber']] == '3694' | adp[['serialNumber']] == '5588'){
  z <- -adp[['instrument_depth']] #since z positive is up and depth positive is down
  p <- round(gsw_p_from_z(z, adp[['latitude']]), digits = 0) #to match number of significant figures of instrument_depth
  
  adp[['pressure']] <- rep_len(p, length(adp[['pressure']])) # array of pres value in shape of pressure vector
  adp@processingLog <- processingLogAppend(adp@processingLog, value = sprintf('Pressure values calculated from static instrument depth (%s m) using the TEOS-10 75-term expression for specific volume and rounded to %s significant digits', adp[['instrument_depth']], as.english(nchar(p))))
  
  #print(adp[['pressure']][1:50])
}


# Create a subdirectory for pre- and post-processing R plots
dir.create('R_plots')
plotDir <- './R_plots'
setwd(plotDir)


# Plots
startPlots(adp, path = getwd())

# To save the map plot
plotMap(adp)

# Make bin plots for u, v, w, and er
for (i in 1:4){
  binPlot(adp, path = getwd(), x = adp[['v']][,,i])
}


#uv scatter
pdf('UV Scatter PreProcessing.pdf')
plot(adp, which = 28, main = 'UV Scatter: PreProcessing') #uv scatter plot
dev.off()


# To save progressive vector plot
pdf('Progressive_vec.pdf')
plot(adp, which = 23, col = 'red')
par(new = TRUE)
mtext('Progressive Vector (U)', side = 3) #progressive vector for u


# Remove leading and trailing ensembles from deployment and recovery
# based off cut_lead_ensembles and cut_trail_ensembles
# adp[['cut_lead_ensembles']] and trail are both of type 'character' so as.numeric() is needed to carry out calculations
st_index <- as.numeric(adp[['cut_lead_ensembles']])+1 #+1 is the first time entry after the cut ones
et_index <- length(adp[['time']])-as.numeric(adp[['cut_trail_ensembles']])

# Adjust time format if needed for mag decl function, limit time function
st <- format(adp[['time']][st_index], '%Y-%m-%d %H:%M:%S', usetz = T)
et <- format(adp[['time']][et_index], '%Y-%m-%d %H:%M:%S', usetz = T)

#st and et values weren't showing up in processing log, just value=st or value=et
adp <- oceSetMetadata(adp, 'time_coverage_start', st)
adp <- oceSetMetadata(adp, 'time_coverage_end', et)


# Convert coordinate system, oceCoordinate, if not in enu
# since enu is the only coordinate system accepted by applyMagneticDeclinationAdp()
coord <- adp[['oceCoordinate']]
if (coord == 'enu'){
} else {
  adp <- toEnuAdp(adp)
}


# apply magnetic declination
adp<- applyMagneticDeclinationAdp(adp, lat=adp[['latitude']], lon=adp[['longitude']], st=adp[['time_coverage_start']], et=adp[['time_coverage_end']], tz='UTC', type='average')


#check plots
# looking for rotation compared to plots before magnetic declination applied
par(new = TRUE)
plot(adp, which = 23, axes = FALSE) #progressive vector for u
legend('topright', legend = list('PreProcessing', 'PostProcessing'), col = c('red', 'black'), lty = 1, cex = 0.8)
dev.off()


#mtext('Progressive Vector (U): PostProcessing', side = 3)
# par(new = TRUE)
# plot(adp, which = 28) #uv scatter plot


# Set flags to 0 to indicate no processing had been done
adp <- adpFlag(adp, adp[['percentgd_threshold']], adp[['error_threshold']])
adpClean <- handleFlags(adp, flags = list(all=c(0)), actions = list('NA')) #flags=4


# Check data visually post processing
# Dataset still maintains complete integrity
#endPlots(adpClean, path = plotDir) 


# Make bin plots for u, v, w, er
for (qc in c('u', 'v', 'w', 'er')){
  qcPlots(adp, QC = qc, path = getwd())
}


# exit R_Plots subdirectory for netCDF file output
setwd(wd)


#these set values to NA and don't remove them
adp <- limit_depthbytime(adp, tz = "UTC")


# limit time, v, pressure, salinity, temp, pitch, roll, and heading based on deployment/recovery times
adp <- limit_time(adp)


# Fix orientation name so it can be used by oceNc_create() to calculate geospatial_vertical_min and max
# and so it aligns with previously made .adcp files (why the function code for oceNc_create wasn't changed instead)
if (adp[['orientation']] == "upward"){
  adp <- oceSetMetadata(adp, "orientation", "up")
} else if (adp[['orientation']] == "downward"){
  adp <- oceSetMetadata(adp, "orientation", "down")
}


# processingLog export
adp <-exportPL(adp)

# export to netCDF file
oceNc_create(adp, ncname_L)

# Return completion message to Python
cat('Done')
# Remove everything from R workspace
# rm(list=ls(all=TRUE))
