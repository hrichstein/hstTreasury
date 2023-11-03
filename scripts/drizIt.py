# needs to be run in the "driz" environment
import argparse as ap
from datetime import date
from glob import glob
import os
import shutil

import pandas as pd

from drizzlepac import astrodrizzle

def drizIt(targname,filt,dataDir,cfg,scale,fpf):
    print(f'Going to drizzle images for {targname} in F{filt}W')

    today = date.today()
    cDate = today.strftime("%d%b")

    flc_list = glob(os.path.join(dataDir,'j*_flc.fits'))
    print(flc_list)

    astrodrizzle.AstroDrizzle(flc_list,
                              output=os.path.join(dataDir,f'F{filt}W_{cDate}'),
                              configobj=cfg,final_scale=scale,
                              final_pixfrac=fpf)

    shutil.copy(os.path.join(dataDir,f'F{filt}W_{cDate}_drc_sci.fits'),
                os.path.join(dataDir,f'{targname}_f{filt}w.fits'))
    
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
        dataDir = os.path.join(config.script.dataDir,f'{targname}_f{ff}w')
        drizIt(targname,ff,dataDir,cfg,scale,fpf)

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

    # for running just from the code itself... (if you want to run both filters in one go)?

    # args = parser.parse_args('')
    # args.filter = '606'
    # main(args)

    # args.filter = '814'
    # main(args)
