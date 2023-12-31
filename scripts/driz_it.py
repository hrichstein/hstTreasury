# needs to be run in the "driz" environment
import argparse as ap
from datetime import date
from glob import glob
import os
import shutil

from drizzlepac import astrodrizzle
import pandas as pd


def driz_it(targname, filt, data_dir, cfg, scale, fpf):
    print(f'Going to drizzle images for {targname} in F{filt}W')

    today = date.today()
    cDate = today.strftime("%d%b")

    flc_list = glob(os.path.join(data_dir, '*_flc.fits'))
    print(flc_list)

    astrodrizzle.AstroDrizzle(flc_list,
                              output=os.path.join(data_dir,
                                                  f'F{filt}W_{cDate}'),
                              configobj=cfg,
                              final_scale=scale,
                              final_pixfrac=fpf)

    shutil.copy(os.path.join(data_dir, f'F{filt}W_{cDate}_drc_sci.fits'),
                os.path.join(data_dir, f'{targname}_f{filt}w.fits'))

    return None


def main(args):
    config = pd.read_json(args.config)

    targname = config.main.targname
    cfg = config.driz.cfg
    scale = config.driz.scale  # final_scale
    fpf = config.driz.fpf  # final_pixfrac
    filt_arr = [
        f'{config.main.filt1}',
        f'{config.main.filt2}'
    ]

    for ff in filt_arr:
        data_dir = os.path.join(config.script.data_dir, f'{targname}_f{ff}w')
        driz_it(targname, ff, data_dir, cfg, scale, fpf)

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
