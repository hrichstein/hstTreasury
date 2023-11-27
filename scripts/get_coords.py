import argparse as ap
from datetime import datetime
from glob import glob
import os
import shutil

from astropy.io import fits, ascii
from astropy import wcs
import numpy as np
import pandas as pd

def get_coords(targname, cat_name, filt_arr,
               res_dir, xpos='xcenter', ypos='ycenter', suffix='.dat', drc_loc=None, force=False):
    """
    cat_name: string that has entire path of the filename 
    that has pixel positions to be converted to coordinates
    drc_loc: path to where the drizzled image is located if
    it is not in the res_dir
    """

    dat = ascii.read(cat_name)
    if not force:
        try:
            _ = dat[f'RA_{filt_arr[0]}']
            print('There already seem to be coordinates \
                  in this catalog.\n')
            print('If you would like to still run this, \
                  use the "force" flag option.')
            return None
        except:
            # Code to get the filename with everything except suffix
            root_name, _ = os.path.splitext(cat_name)
            orig_file = os.path.join(f'{root_name}_orig{suffix}')
            # Making a copy of the original file for safety
            shutil.copy(cat_name, orig_file)

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
            if os.path.exists(os.path.join(res_dir, f'{targname}_f{ff}w.fits')):
                drc_file = os.path.join(res_dir, f'{targname}_f{ff}w.fits')
            else:
                drc_file = os.path.join(
                    os.path.dirname(res_dir), f'{targname}_f{ff}w.fits')

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
    suffix = config.aper.suffix
    filt_arr = [
        f'{config.main.filt1}',
        f'{config.main.filt2}'
    ]

    res_dir = os.path.join(config.script.res_dir, targname)
    if args.date is not None:
        phot_dir = os.path.join(res_dir, f'drc_phot{args.date}')
    else:
        phot_dir_list = glob(os.path.join(res_dir, f'drc_phot*'))
        if len(phot_dir_list) == 1:
            phot_dir = phot_dir_list[0]
        else:
            phot_dir = max(phot_dir_list, key=lambda p:
                           datetime.strptime(p, os.path.join(res_dir, f'drc_phot%d%b%y')))

    if psf:
        cat_name = os.path.join(res_dir, f'{targname}_dp_matched.dat')
        get_coords(targname, cat_name, filt_arr, res_dir,
                   xpos='X', ypos='Y', suffix=suffix, force=eval(args.force))

    else:
        cat_name = os.path.join(phot_dir,
                                f'{targname}_matched_drc.dat')
        get_coords(targname, cat_name, filt_arr, phot_dir,
                   suffix=suffix, force=eval(args.force))

    return None


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Gets RA and DEC coordinates.'
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
        '-d',
        '--date',
        help='Date of aperture photometry in format \
        DDMonYY (01Jan24).',
        type=str,
    )
    _ = parser.add_argument(
        '-f',
        '--force',
        help='Flag for forcing coordinate \
        calculation, even if it has previously been done.',
        default='False',
        type=str,
    )
    _ = parser.add_argument(
        '-p',
        '--psf',
        help='Flag for denoting PSF input catalog.',
        default='False',
        type=str,
    )
    args = parser.parse_args()

    raise SystemExit(main(args))
