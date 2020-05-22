'''
author: Hana Hourston
date: March 2, 2020

The following code was adapted from
https://www.kdnuggets.com/2015/10/integrating-python-r-executing-part2.html

'''

import subprocess

# Define command
command ='Rscript'

# Define arguments
# Path to the R ADCP processing script on your computer that you want to execute from Python
path2script ='/home/hourstonh/Documents/Hana_D_drive/ADCP_processing/R_scripts/ADCP_lvl1_process_pycalled.R'
# Working directory
wd = '/home/hourstonh/Documents/Hana_D_drive/ADCP_processing/callR_fromPython/'
# Path to the raw ADCP file to be processed
path2raw = '/home/hourstonh/Documents/Hana_D_drive/ADCP_processing/callR_fromPython/20568-A1-56.000'
# Path to the csv metadata file for the raw ADCP file
path2meta = '/home/hourstonh/Documents/Hana_D_drive/ADCP_processing/ADCP/a1_20160713_20170513_0480m/P01/a1_20160713_20170513_0480m_meta_L1.csv'

# Create a list of arguments: the working directory, raw ADCP file and associated metadata file paths
args = [wd, path2raw, path2meta]

# Build subprocess command
cmd = [command, path2script] + args

# check_output will run the command and store to result
x = subprocess.check_output(cmd, universal_newlines=True)

print(x)
