import argparse as ap
from datetime import datetime
from glob import glob
import os

from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def make_cmd(targname, filt1, filt2, xmin, xmax, ymin, ymax, phot_dir):
    drc_file = ascii.read(os.path.join(phot_dir,
                                       f'{targname}_matched_drc.dat'))

    mask = np.logical_or(drc_file[f'six_4_flag_{filt1}'] == 1,
                         drc_file[f'six_4_flag_{filt2}'] == 1)
    drc = drc_file[mask]

    fig, ax = plt.subplots(figsize=(4, 6.5))
    ax.scatter(drc_file[f'magr_{filt1}']-drc_file[f'magr_{filt2}'],
               drc_file[f'magr_{filt1}'],
               s=10,
               label='All',
               color='gray',
               alpha=0.5)

    ax.scatter(drc[f'magr_{filt1}']-drc[f'magr_{filt2}'],
               drc[f'magr_{filt1}'],
               s=10,
               label='Likely Stars',
               color='k',
               marker='*')

    ax.set_ylim(ymax, ymin)
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel(f'F{filt1}W - F{filt2}W')
    ax.set_ylabel(f'F{filt1}W')
    ax.set_title(f'{targname} DRC APERTURE CMD')

    ax.legend()
    ax.set_rasterization_zorder(1)
    plt.savefig(os.path.join(phot_dir, f'{targname}_drc.pdf'),
                dpi=150, bbox_inches='tight')
    plt.close()

    return None


def main(args):
    config = pd.read_json(args.config)
    targname = config.main.targname
    filt_arr = [
        f'{config.main.filt1}',
        f'{config.main.filt2}'
    ]

    res_dir = os.path.join(config.aper.res_dir, targname)
    if args.date is not None:
        phot_dir = os.path.join(res_dir, f'drc_phot{args.date}')
    else:
        phot_dir_list = glob(os.path.join(res_dir, f'drc_phot*'))
        if len(phot_dir_list) == 1:
            phot_dir = phot_dir_list[0]
        else:
            phot_dir = max(phot_dir_list,
                           key=lambda p: datetime.strptime(p,
                                                           os.path.join(res_dir, f'drc_phot%d%b%y')))

    make_cmd(targname, filt_arr[0], filt_arr[1],
             config.cmd.xmin, config.cmd.xmax,
             config.cmd.ymin, config.cmd.ymax,
             phot_dir)

    return None


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Creates a CMD using the \
        matched DRC photometry.'
    )
    _ = parser.add_argument(
        '-c',
        '--config',
        help='Name of the config json file.\
        (Default: config.json)',
        default='../../config.json',
        type=str,
    )
    _ = parser.add_argument(
        '-d',
        '--date',
        help='Date of aperture photometry in format \
        DDMonYY (01Jan24).',
        type=str,
    )
    args = parser.parse_args()

    raise SystemExit(main(args))
