# needs to be run in the "driz" environment
import argparse as ap

import os
import shutil
import sys

from datetime import date
from glob import glob

import numpy as np
import pandas as pd

from drizzlepac import astrodrizzle

def main(args):

    config = pd.read_json(args.config)
    filt = args.filter
    cfg = config.driz.cfg
    scale = config.driz.scale

    print(filt)

    targname = config.main.targname
    dataDir = config.script.dataDir
    resDir = config.script.resDir

    tDir = os.path.join(dataDir,f'{targname}_f{filt}w')
    oDir = os.path.join(resDir,targname)

    today = date.today()
    cDate = today.strftime("%d%b")

    flc_list = glob(os.path.join(tDir,'j*_flc.fits'))
    print(flc_list)

    astrodrizzle.AstroDrizzle(flc_list,
                              output=os.path.join(tDir,f'F{filt}W_{cDate}'),
                              configobj=cfg,final_scale=scale)

    shutil.copy(os.path.join(tDir,f'F{filt}W_{cDate}_drc_sci.fits'),os.path.join(oDir,f'{targname}_f{filt}w.fits'))

    return None


# if __name__=='__main__': 
#     filt_arr = ['606','814']
#     for ff in filt_arr:
#         doDriz(sys.argv[1],filt=ff)

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

    _ = parser.add_argument(
        '-f',
        '--filter',
        help='String integer of filter (e.g., \'606\' for F606W)',
        default = '606',
        type=str,
    )
    args = parser.parse_args()

    # for running just from the code itself... (if you want to run both filters in one go)?

    # args = parser.parse_args('')
    # args.filter = '606'
    # main(args)

    # args.filter = '814'
    # main(args)
