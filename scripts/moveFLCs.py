import argparse as ap
from glob import glob

from astropy.io import fits
from astropy.table import Table

from pathlib import Path
import os
import pandas as pd

# Add functionality to delete empty MAST folders

def main(args):
    config = pd.read_json(args.config)

    targname = config.main.targname
    dataDir = config.script.dataDir
    resDir = config.script.resDir

    filt_arr = [f'{config.main.filt1}',f'{config.main.filt2}']

    filelist = glob(f'{dataDir}/mast*/**/*flc.fits',recursive=True)

    # If the files were downloaded via curl method
    if len(filelist) == 0:
        filelist = glob(f'{dataDir}/MAST*/**/*flc.fits',recursive=True)

    if len(filelist) == 0:
        print('No FLCs to be moved!')
        return None

    if not os.path.exists(os.path.join(resDir,targname)):
        os.mkdir(os.path.join(resDir,targname))      

    for ff in filt_arr:
        if not os.path.exists(os.path.join(dataDir,f'{targname}_f{ff}w')):
            os.mkdir(os.path.join(dataDir,f'{targname}_f{ff}w'))

    for file in filelist:
        with fits.open(file,names=True) as hdu:
            if hdu[0].header['TARGNAME'] == targname:
                if hdu[0].header['EXPTIME'] < 200:
                    os.remove(file)
                    print('Removed short exposure')
                if hdu[0].header['FILTER1'][0] == 'F':
                    filt= hdu[0].header['FILTER1']
                else:
                    filt = hdu[0].header['FILTER2']
                fileroot = os.path.basename(file)
                ff = filt[1:4] 
                os.rename(file,os.path.join(dataDir,f'{targname}_f{ff}w',fileroot))
            
            else:
                print(f'Skipping {file}...')
                continue
    
    print('Finished moving the FLCs.')    

    return None


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Move downloaded FLC files'
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

    df = main(args)