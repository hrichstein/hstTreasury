import argparse as ap
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from astropy.convolution import convolve
from astropy.io import fits
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize

from photutils.segmentation import make_2dgaussian_kernel, detect_sources
from photutils.segmentation import SourceCatalog
from photutils.utils import circular_footprint


def mask_acs(targname, filt, data_dir, qual_dir):
    fname = f'{targname}_f{filt}w.fits'
    print(f'Masking {fname}...')
    with fits.open(os.path.join(data_dir, fname)) as hdu:
        try:
            data = hdu[0].data
        except:  # May throw an Attribute Error?
            data = hdu.data

    kernel = make_2dgaussian_kernel(3.0, size=5)
    convolved_data = convolve(data, kernel)
    segment_map = detect_sources(convolved_data, 100, npixels=200)

    # Make plot for checking
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(segment_map, origin='lower',
              cmap=segment_map.cmap, interpolation='nearest')

    ax.set_title('Segmentation Image')
    ax.set_rasterization_zorder(1)

    plt.savefig(os.path.join(qual_dir, f'segm_f{filt}w.pdf'),
                bbox_inches='tight', dpi=150)
    plt.close()

    cat = SourceCatalog(data, segment_map,
                        convolved_data=convolved_data)

    print('Converting catalog to table...')
    tbl = cat.to_table()
    arr = np.asarray(tbl['area'])
    idx = np.argwhere(arr >= 700)
    idx = idx.flatten()

    segment_map.keep_labels(labels=tbl['label'][idx])

    footprint = circular_footprint(radius=3)
    mask = segment_map.make_source_mask(footprint=footprint)

    tmp2 = data*[~mask]
    tmp3 = np.where(tmp2 == 0, 1e11, tmp2)

    norm = ImageNormalize(stretch=SqrtStretch(), vmax=8000)

    # Making plot of fits image after the masking
    fig, ax = plt.subplots(figsize=(15, 12.5))
    ax.imshow(tmp3[0], origin='lower', norm=norm, cmap='Greys_r')

    ax.set_title('Masked Image')
    ax.set_rasterization_zorder(1)

    plt.savefig(os.path.join(qual_dir, f'masked_f{filt}w.pdf'),
                bbox_inches='tight', dpi=150)
    plt.close()

    print('Writing masked FITS...')
    with fits.open(os.path.join(data_dir, fname)) as hdu:
        hdu[0].data = tmp3
        hdu.writeto(os.path.join(data_dir,
                                 f'{targname}_mf{filt}w.fits'),
                    overwrite=True)

    return None


def mv_fits(targname, filt, data_dir, res_dir):
    # Masked file
    mname = os.path.join(data_dir, f'{targname}_mf{filt}w.fits')
    if not os.path.exists(mname):
        print('The masked file doesn\'t seem to exist. ')
        print('Check some things.')
        return None

    # Old, unmasked file
    oname = os.path.join(data_dir, f'{targname}_f{filt}w.fits')
    # New, unmasked file
    nname = os.path.join(data_dir, f'{targname}_f{filt}w_nomask.fits')
    # Masked file for results dir.
    rname = os.path.join(res_dir, f'{targname}_f{filt}w.fits')

    c0 = os.path.exists(nname)
    if c0:
        print(f'{nname} already exists.')
        return None

    os.rename(oname, nname)
    print('Moved original file to nomask fits.')

    c1 = os.path.exists(oname)
    if c1:
        print('Error in moving files.')
        return None

    # Copying masked file to results dir.
    shutil.copy(mname, rname)
    # Renaming masked file in data dir.
    os.rename(mname, oname)
    print('Moved files successfully.')

    return None


def main(args):
    config = pd.read_json(args.config)

    targname = config.main.targname
    filt_arr = [
        f'{config.main.filt1}',
        f'{config.main.filt2}'
    ]

    res_dir = os.path.join(config.script.res_dir, targname)
    qual_dir = os.path.join(config.script.qual_dir, targname)

    for ff in filt_arr:
        data_dir = os.path.join(config.script.data_dir,
                                f'{targname}_f{ff}w')
        # Checking if this code was run before,
        # which would mean the file was already masked
        if not os.path.exists(
                os.path.join(data_dir, f'{targname}_f{ff}w_nomask.fits')):
            mask_acs(targname, ff, data_dir, qual_dir)
        mv_fits(targname, ff, data_dir, res_dir)

    return None


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Mask the DRC files.'
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
