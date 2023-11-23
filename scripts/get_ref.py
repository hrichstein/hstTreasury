import argparse as ap
from datetime import date as dd
import os

from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np


def get_ref(aper_file, psf_file, res_dir, qual_dir):
    """
    Outputs a list of reference stars for matching between
    PSF and aperture photometry catalogs
    """
    today = dd.today()
    date = today.strftime("%d%b%y")

    outname = os.path.join(res_dir, f'ref_aper_psf_{date}')
    reg_606 = os.path.join(qual_dir)

    aper = ascii.read(aper_file)
    psf = ascii.read(psf_file)

    return None


def main(args):
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
    args = parser.parse_args()

    raise SystemExit(main(args))
