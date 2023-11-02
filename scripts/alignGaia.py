# Copied/adapted from Gaia_alignment.ipynb 
# from gaia_alignment repository on spacetelescope github
# https://github.com/spacetelescope/gaia_alignment/blob/master/Gaia_alignment.ipynb
# needs to be run in the "driz" environment

import astropy.units as u
import glob
import numpy as np
import os

from astropy.io import fits
from astropy.table import Table
from astropy.units import Quantity
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astropy.wcs import WCS

from drizzlepac import tweakreg
from stwcs import updatewcs
from stsci.tools import teal

import math

def get_footprints(im_name):
    """Calculates positions of the corners of the science extensions of some image 'im_name' in sky space"""
    footprints = []
    hdu = fits.open(im_name)
    
    flt_flag = 'flt.fits' in im_name or 'flc.fits' in im_name
    
    # Loop ensures that each science extension in a file is accounted for. This is important for 
    # multichip imagers like WFC3/UVIS and ACS/WFC
    for ext in hdu:
        if 'SCI' in ext.name:
            hdr = ext.header
            wcs = WCS(hdr, hdu)
            footprint = wcs.calc_footprint(hdr, undistort=flt_flag)
            footprints.append(footprint)
    
    hdu.close()
    
    return footprints

# ----------------------------------------------------------------------------------------------------------
def bounds(footprint_list):
    """Calculate RA/Dec bounding box properties from multiple RA/Dec points"""
    
    # Flatten list of extensions into numpy array of all corner positions
    merged = [ext for image in footprint_list for ext in image]
    merged = np.vstack(merged)
    ras, decs = merged.T
    
    # Compute width/height
    delta_ra = (max(ras)-min(ras))
    delta_dec = max(decs)-min(decs)

    # Compute midpoints
    ra_midpt = (max(ras)+min(ras))/2.
    dec_midpt = (max(decs)+min(decs))/2.
    
    return ra_midpt, dec_midpt, delta_ra, delta_dec
# ----------------------------------------------------------------------------------------------------------

def get_error_mask(catalog, max_error):
    """Returns a mask for rows in catalog where RA and Dec error are less than max_error"""
    ra_mask = catalog['ra_error']< max_error
    dec_mask = catalog['dec_error'] < max_error
    mask = ra_mask & dec_mask
#     print('Cutting sources with error higher than {}'.format(max_error))
#     print('Number of sources befor filtering: {}\nAfter filtering: {}\n'.format(len(mask),sum(mask)))
    
    return mask

# ----------------------------------------------------------------------------------------------------------

def alignGaia(targname,filt,mainDir='../'):

    tDir = os.path.join(mainDir,'data',f'{targname}_f{filt}w')

    str_in = os.path.join(tDir,'*fl?.fits')
    images = glob.glob(str_in)
    print(images)
    footprint_list = list(map(get_footprints, images))

    ra_midpt, dec_midpt, delta_ra, delta_dec = bounds(footprint_list)

    coord = SkyCoord(ra=ra_midpt, dec=dec_midpt, unit=u.deg)
    print(coord)

    width = Quantity(delta_ra, u.deg)
    height = Quantity(delta_dec, u.deg)

    r = Gaia.query_object_async(coordinate=coord, width=width, height=height)

    mask = get_error_mask(r, 2.)

    ras = r['ra']
    decs = r['dec']

    tbl = Table([ras[mask], decs[mask]]) # Make a temporary table of just the positions
    outname = os.path.join(tDir,'gaia.cat')
    tbl.write(outname, format='ascii.fast_commented_header') # Save the table to a file
    
    str_in = os.path.join(tDir,'*flc.fits')
    input_images = sorted(glob.glob(str_in)) 

    # This would update the WCS, but it tends to break things and should be unnecessary 
    # if the files are newly downloaded from MAST
    # if first:
    #     derp = list(map(updatewcs.updatewcs, input_images))

    cat = os.path.join(tDir,'gaia.cat')
    wcsname ='GAIA'
    teal.unlearn('tweakreg')
    teal.unlearn('imagefindpars')

    cw = 3.5 # psf width measurement (2*FWHM).  Use 3.5 for WFC3/UVIS and ACS/WFC and 2.5 for WFC3/IR

    tweakreg.TweakReg(input_images, # Pass input images
                  updatehdr=True, # update header with new WCS solution
                  imagefindcfg={'threshold':500.,'conv_width':cw},# Detection parameters, threshold varies for different data
                  separation=0.0, # Allow for very small shifts
                  refcat=cat, # Use user supplied catalog (Gaia)
                  clean=True, # Get rid of intermediate files
                  interactive=False,
                  see2dplot=False,
                  writecat=False,
                  residplot='No plot',
                  #shiftfile=False, # Save out shift file (so we can look at shifts later)
                  wcsname=wcsname, # Give our WCS a new name
                  reusename=True,
                  fitgeometry='general', # Use the 6 parameter fit
                  configobj=None,
                  searchrad=10, tolerance=2.0) 

    return None


if __name__=='__main__':
    targname = input("Enter targname (ex. HOROLOGIUM I): ")
    filt_arr = ['606','814']
    for ff in filt_arr:
        alignGaia(targname,ff,mainDir='../')
    # os.remove('./residuals*.png')
    # os.remove('./vector*.png')
    
