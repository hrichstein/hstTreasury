import argparse as ap
import os
import shutil

import numpy as np
import pandas as pd

from astropy.io import fits
from astropy.convolution import convolve
from photutils.segmentation import make_2dgaussian_kernel, detect_sources
from photutils.segmentation import SourceCatalog
from photutils.utils import circular_footprint

import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

def maskACS(targname,filt,dataDir,qualDir):
    fname = f'{targname}_f{filt}w.fits'
    print(f'Masking {fname}...')
    with fits.open(os.path.join(dataDir,fname)) as hdu:
        try:
            data = hdu[0].data
        except: # May throw an Attribute Error?
            data = hdu.data
            
    kernel = make_2dgaussian_kernel(3.0,size=5)
    convolved_data = convolve(data,kernel)
    segment_map = detect_sources(convolved_data,100,npixels=200)
    
    # Make plot for checking
    fig, ax = plt.subplots(figsize=(6,6))
    ax.imshow(segment_map,origin='lower',cmap=segment_map.cmap,interpolation='nearest')
    
    ax.set_title('Segmentation Image')
    ax.set_rasterization_zorder(1)
    
    plt.savefig(os.path.join(qualDir,f'segm_f{filt}w.pdf'),bbox_inches='tight',dpi=150)
    plt.close()
    
    cat = SourceCatalog(data,segment_map,convolved_data=convolved_data)
    
    print('Converting catalog to table...')
    tbl = cat.to_table()
    arr = np.asarray(tbl['area'])
    idx = np.argwhere(arr>=700)
    idx = idx.flatten()
    
    segment_map.keep_labels(labels=tbl['label'][idx])

    footprint = circular_footprint(radius=3)
    mask = segment_map.make_source_mask(footprint=footprint)
    
    tmp2 = data*[~mask]
    tmp3 = np.where(tmp2==0,1e11,tmp2)
    
    norm = ImageNormalize(stretch=SqrtStretch(),vmax=8000)

    # Making plot of fits image after the masking
    fig, ax = plt.subplots(figsize=(15, 12.5))
    ax.imshow(tmp3[0],origin='lower',norm=norm,cmap='Greys_r')

    ax.set_title('Masked Image')
    ax.set_rasterization_zorder(1)
    
    plt.savefig(os.path.join(qualDir,f'masked_f{filt}w.pdf'),bbox_inches='tight',dpi=150)
    plt.close()
    
    print('Writing masked FITS...')
    with fits.open(os.path.join(dataDir,fname)) as hdu:
        hdu[0].data = tmp3
        hdu.writeto(os.path.join(dataDir,f'{targname}_mf{filt}w.fits'),overwrite=True)
    
    return None


def mvFits(targname,filt,dataDir,resDir):
    mname = os.path.join(dataDir,f'{targname}_mf{filt}w.fits')  # Masked file
    if not os.path.exists(mname):
        print('The masked file doesn\'t seem to exist. Check some things.')
        return None
    
    oname = os.path.join(dataDir,f'{targname}_f{filt}w.fits')  # Old, unmasked file
    nname = os.path.join(dataDir,f'{targname}_f{filt}w_nomask.fits')  # New, unmasked file
    rname = os.path.join(resDir,f'{targname}_f{filt}w.fits') # Masked file for results dir.
    
    c0 = os.path.exists(nname)
    if c0:
        print(f'{nname} already exists. Check yourself.')
        return None
    
    os.rename(oname,nname)
    print('Moved original file to nomask fits.')
    
    c1 = os.path.exists(oname)
    if c1:
        print('Error in moving files.')
        return None
    
    shutil.copy(mname,rname)  # Copying masked file to results dir.
    os.rename(mname,oname)  # Renaming masked file in data dir.
    print('Moved files successfully.')
    
    return None


def main(args):
    config = pd.read_json(args.config)

    targname = config.main.targname    
    filt_arr = [
        f'{config.main.filt1}',
        f'{config.main.filt2}'
        ]
    
    resDir = os.path.join(config.script.resDir,targname)
    qualDir = os.path.join(config.script.qualDir,targname)
    
    for ff in filt_arr:
        dataDir = os.path.join(config.script.dataDir,f'{targname}_f{ff}w')
        # Checking if this code was run before, which would mean the file was already masked
        if not os.path.exists(
            os.path.join(dataDir,f'{targname}_f{ff}w_nomask.fits') ):
            maskACS(targname,ff,dataDir,qualDir)
        mvFits(targname,ff,dataDir,resDir)
        
    return None


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Drizzle the FLC files'
    )
    _ = parser.add_argument(
        '-c',
        '--config',
        help='Name of the config json file.\
        (Default: config.json)',
        default='../config.json',
        type=str,
    )
    
    args = parser.parse_args()
    
    raise SystemExit(main(args))

