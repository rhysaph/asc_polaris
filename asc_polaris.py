#!/usr/bin/env python

# usage python asc_polaris.py image_name.fits

import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy import coordinates as coord
from photutils.detection import find_peaks
from photutils.detection import find_peaks
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils.centroids import (centroid_1dg, centroid_2dg,
                                 centroid_com, centroid_quadratic)
import math
from astropy.io import fits

# Check input parameters
if ( len(sys.argv) != 2):
    print("There is a problem with your arguments.")
    print("Argument 1 should be an existing fits file.")
    sys.exit('Parameters are wrong')

# If all OK, then first argument is the file to open
imagefile = sys.argv[1]
print("Opening ",imagefile)

# Open image
hdu = fits.open(imagefile)[0]


# We don't need a WCS for this, just the pixel coordinates
# enough, and they need to be approximate to detect the movement.

# Approximate x, y coordinates of Polaris
polaris_x = 1467
polaris_y = 1732

print("Approximate x,y for Polaris is ",polaris_x, polaris_y)

# Set size of search box
starbox_size = 40

# Define a box
# For some reason, xdirection is point[0][1] and y direction is point[0][0]

startarrx = int(polaris_y - starbox_size/2)
endarrx = int(polaris_y + starbox_size/2)
startarry = int(polaris_x - starbox_size/2)
endarry = int(polaris_x + starbox_size/2)

# Get the data from the image array
skybox = hdu.data[startarrx:endarrx,startarry:endarry]

# Check for saturated values
#saturated = 65535
saturated = 65534

if np.any(x  >= saturated for x in skybox):
    print("WARNING! Saturated value found\n")

mean, median, std = sigma_clipped_stats(skybox, sigma=3.0)  


# Create mask array wth all values set to True
mask = np.ones(hdu.data.shape, dtype=bool)

# Set the subarray to false (no masking)
mask[startarrx:endarrx,startarry:endarry] = False

# set up the daofind parameters
daofind = DAOStarFinder(fwhm=3.0, threshold=5. * std)

# Run the DAO starfinder software
sources = daofind(hdu.data - median, mask=mask)

for col in sources.colnames:  
    sources[col].info.format = '%.8g'  # for consistent table output

print(sources)
#help(sources)
print(" ")
# Sort sources by flux
sources.sort('flux', reverse=True)
print(sources)

print("Polaris x, y is ",sources['xcentroid'][0],", ",sources['ycentroid'][0])

#print("\nmean, median, std are \n",mean, median, std)

print("Remember FITS images start at (1,1) while python is (0,0)")
