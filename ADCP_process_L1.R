# Inputs: a raw ADCP file (.000, .664, or .pd0) and a csv file containing metadata
# Corrects for magnetic declination
# Flags bad leading and trailing ensembles
# Creates a netCDF file

library(ADCP) 
library(ncdf4)
library(oce)
library(tools)


# Calculates average orientation
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


# Set working directory; a subfolder for plots will be created in this directory
wd <- 'Your/path/here'
setwd(wd)

# Write paths to the raw ADCP file and its corresponding csv metadata file
f <- './sample_data/a1_20050503_20050504_0221m.000'
meta = './sample_data/a1_20050503_20050504_0221m_meta_L1.csv'

# Create standard name from raw file name (WITH data type and processing level add-on)
ncname_L <- paste0(substr(basename(meta), 1, nchar(basename(meta))-12), '.adcp.L1')
# ncname below is for use in navigating to the directory of the file name
ncname <- substr(basename(meta), 1, nchar(basename(meta))-12)

# Set wd
wd <- sprintf('/home/hourstonh/Documents/Hana_D_drive/ADCP_processing/ADCP/%s/P01', ncname)
setwd(wd)

# Read in the data and metadata
adp <- read.adp.easy(f, meta)


# Take average orientation since it is now vector with an orientation for every ensemble
adp <- oceSetMetadata(adp, 'orientation', mean_orientation(adp[['orientation']]))

# Convert numeric metadata items a  s needed
adp <- oceSetMetadata(adp, 'country_institute_code', as.numeric(adp[['country_institute_code']]))
adp <- oceSetMetadata(adp, 'instrument_depth', as.numeric(adp[['instrument_depth']]))
adp <- oceSetMetadata(adp, 'latitude', as.numeric(adp[['latitude']]))
adp <- oceSetMetadata(adp, 'longitude', as.numeric(adp[['longitude']]))
adp <- oceSetMetadata(adp, 'processing_level', as.numeric(adp[['processing_level']]))
adp <- oceSetMetadata(adp, 'serialNumber', toString(adp[['serialNumber']])) # since some numbers have a leading 0
adp <- oceSetMetadata(adp, 'station_number', as.numeric(adp[['station_number']]))
adp <- oceSetMetadata(adp, 'water_depth', as.numeric(adp[['water_depth']]))

# Update naming_authority from CF v52 to CF v72
adp <- oceSetMetadata(adp, 'naming_authority', 'BODC, MEDS, CF v72')
# Correct flag_meanings attribute to follow BODC SeaDataNet
adp <- oceSetMetadata(adp, 'flag_meaning', 'no_quality_control, good_value, probably_good_value, probably_bad_value, bad_value, changed_value, value_below_detection, value_in_excess, interpolated_value, missing_value')

#For numbers starting with 0 so they have leading zero added back in
if (nchar(adp[['serialNumber']]) == 3){
  adp <- oceSetMetadata(adp, 'serialNumber', paste0('0', adp[['serialNumber']]))
}


# If missing pressure sensor, calculate pressure based on static (instrument) depth
# In these cases pressure may be all zeros, or it may have some NAs, or one very large negative or positive value with the rest zeros
if (getmode(adp[['pressure']]) == 0 | adp[['instrumentSubtype']] == 'BroadBand'){
  z <- -adp[['instrument_depth']] #since z positive is up and depth positive is down
  p <- round(gsw_p_from_z(z, adp[['latitude']]), digits = 0) #to match number of significant figures of instrument_depth
  
  adp[['pressure']] <- rep_len(p, length(adp[['pressure']])) # array of pres value in shape of pressure vector
  adp@processingLog <- processingLogAppend(adp@processingLog, value = sprintf('Pressure values calculated from static instrument depth (%s m) using the TEOS-10 75-term expression for specific volume and rounded to %s significant digits', adp[['instrument_depth']], toString(nchar(p))))
  
  print(sprintf('Pressure values calculated from static instrument depth (%s m)', adp[['instrument_depth']]))
  print(adp[['pressure']][1:50])
}

# Check instrument_depth from metadata csv file: compare with pressure values
depths_check <- mean(adp[['pressure']], na.rm = T) - adp[['distance']]
inst_depth_check <- (depths_check[1] + adp[['distance']][1])[1]
abs_difference <- abs(inst_depth_check - adp[['instrument_depth']])
# Calculate percent difference in relation to total water depth
if ((abs_difference / adp[['water_depth']] * 100) > 0.05){
  warning("Difference between calculated instrument depth and metadata instrument_depth 
          exceeds 0.05% of the total water depth")
}

# Create a subdirectory for pre- and post-processing R plots
dir.create('R_plots')
plotDir <- './R_plots'
setwd(plotDir)


# Plots
startPlots(adp, path = getwd())

# map plot
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
# adp[['cut_lead_ensembles']] and trail are both of type 'character'
# so as.numeric() is needed to carry out calculations
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

# Flag pressure separately
# Based on start/end times
PRESPR01_QC <- rep(0, length(adp@data$pressure))
PRESPR01_QC[1:st_index-1] = 4
if (et_index != length(PRESPR01_QC)){
  PRESPR01_QC[et_index+1:length(PRESPR01_QC)] = 4
}

# And negative pressure values
for (i in 1:length(adp@data$pressure)){
  if (adp@data$pressure[i] < 0){
    PRESPR01_QC[i] = 4
  }
}

adp@data$pressure[PRESPR01_QC == 4] <- NA
adp <- oceSetData(adp, 'PRESPR01_QC', PRESPR01_QC)
processingLogAppend(adp@processingLog, 'Negative pressure values flagged as \"bad_data\" and set to nan\'s.')

# Check data visually post processing
# Dataset still maintains complete integrity
endPlots(adpClean, path = getwd()) 


# Make bin plots for u, v, w, er
for (qc in c('u', 'v', 'w', 'er')){
  qcPlots(adp, QC = qc, path = getwd())
}


# exit R_Plots subdirectory for netCDF file output
setwd(wd)


#these set values to NA and don't remove them.
#should not be after qcPlot()
#Limit depth by time should be before adpFlag()
adp <- limit_depthbytime(adp, tz = "UTC")

# limit time, v, pressure, salinity, temp, pitch, roll, and heading based on
# deployment/recovery times
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

  
#Remove everything from R workspace
#rm(list=ls(all=TRUE))
