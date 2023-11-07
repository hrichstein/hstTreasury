import argparse as ap
from datetime import datetime
from glob import glob
import os

from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def makeCMD(targname,filt1,filt2,xmin,xmax,ymin,ymax,photDir):
    drcFile = ascii.read(os.path.join(photDir,
                                      f'{targname}_matchedDRC.dat'))
    
    mask = np.logical_or(drcFile[f'six_4_flag_f{filt1}w']==1,
                         drcFile[f'six_4_flag_f{filt2}w']==1)
    drc = drcFile[mask]
    
    fig,ax = plt.subplots(figsize=(4,6.5))
    ax.scatter(drcFile[f'magr_f{filt1}w']-drcFile[f'magr_f{filt2}w'],
               drcFile[f'magr_f{filt1}w'],
               s=10,
               label='All',
               color='gray',
               alpha=0.5)
    
    ax.scatter(drc[f'magr_f{filt1}w']-drc[f'magr_f{filt2}w'],
               drc[f'magr_f{filt1}w'],
               s=20,
               label='Likely Stars',
               color='k')
    
    ax.set_ylim(ymax,ymin)
    ax.set_xlim(xmin,xmax)
    ax.set_xlabel(f'F{filt1}W-F{filt2}W')
    ax.set_ylabel(f'F{filt1}W')
    ax.set_title(f'{targname} DRC APERTURE CMD')
    
    ax.legend()
    ax.set_rasterization_zorder(1)
    plt.savefig(os.path.join(photDir,f'{targname}_drc.pdf'),
                dpi=150,bbox_inches='tight')
    plt.close()

    return None 


def main(args):
    config = pd.read_json(args.config)
    targname = config.main.targname
    filt_arr = [
        f'{config.main.filt1}',
        f'{config.main.filt2}'
        ]
    
    resDir = os.path.join(config.aper.resDir,targname)
    if args.date is not None:
        photDir = os.path.join(resDir,f'drcPhot{args.date}')
    else:    
        photDirList = glob(os.path.join(resDir,f'drcPhot*'))
        if len(photDirList) == 1:
            photDir = photDirList[0]
        else:
            photDir = max(photDirList, 
                          key=lambda p: datetime.strptime(p, 
                                        os.path.join(resDir,f'drcPhot%d%b%y')))
    
    makeCMD(targname,filt_arr[0],filt_arr[1],
            config.cmd.xmin,config.cmd.xmax,
            config.cmd.ymin,config.cmd.ymax,
            photDir)
    
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