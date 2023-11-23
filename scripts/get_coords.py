import argparse as ap
import os
import shutil

from astropy.io import fits, ascii
from astropy import wcs
import numpy as np
import pandas as pd


def get_coords(targname, cat_name, filt_arr,
               res_dir, xpos='xcenter', ypos='ycenter', drc_loc=None):
    """
    cat_name: string that has entire path of the filename 
    that has pixel positions to be converted to coordinates
    drc_loc: path to where the drizzled image is located if
    it is not in the res_dir
    """

    # Code to get the filename with everything except .dat
    root_name = None  # fill in with proper code
    orig_file = os.path.join(f'{root_name}_orig.dat')

    # Making a copy of the original file for safety
    shutil.copy(cat_name, orig_file)

    dat = ascii.read(cat_name)

    for ff in filt_arr:
        xkey = f'{xpos}_{ff}'
        ykey = f'{ypos}_{ff}'
        # Getting pixel values
        xdat = dat[f'{xkey}']
        ydat = dat[f'{ykey}']

        pix_in = np.vstack((xdat, ydat))

        # Getting drc file
        if drc_loc is not None:
            drc_file = os.path.join(drc_loc, f'{targname}_f{ff}w.fits')
        else:
            drc_file = os.path.join(res_dir, f'{targname}_f{ff}w.fits')

        # Grabbing WCS info
        with fits.open(drc_file, names=True) as dd:
            try:
                wf = wcs.WCS(dd[0].header, naxis=2)
            except:
                wf = wcs.WCS(dd.header, naxis=2)

        coords = wf.all_pix2world(pix_in.T, 1)

        dat[f'RA_{ff}'] = coords[:, 0]
        dat[f'DEC_{ff}'] = coords[:, 1]

    ascii.write(dat, cat_name, overwrite=True, format='commented_header')

    return None


def main(args):
    config = pd.read_json(args.config)
    psf = eval(args.psf)

    targname = config.main.targname
    filt_arr = [
        f'{config.main.filt1}',
        f'{config.main.filt2}'
    ]

    res_dir = os.path.join(config.script.res_dir, targname)

    if psf:
        cat_name = os.path.join(res_dir, f'{targname}_dp_matched.dat')
        get_coords(targname, cat_name, filt_arr, res_dir, xpos='X',
                   ypos='Y')

    else:
        cat_name = os.path.join(res_dir, f'{targname}_matched_drc.dat')
        get_coords(targname, cat_name, filt_arr, res_dir)

    return None


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='TBD.'
    )
    _ = parser.add_argument(
        '-c',
        '--config',
        help='Name of the config json file.\
        (Default: config.json)',
        default='../config.json',
        type=str,
    )
    _ = parser.add_argument(
        '-p',
        '--psf',
        help='Flag for denoting PSF input catalog.',
        default=False,
        type=str,
    )
    args = parser.parse_args()

    raise SystemExit(main(args))
