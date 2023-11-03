# Must be in base environment
# acszpt is particular
import os
import sys

sys.path.insert(0,'../../')

import argparse as ap
import pandas as pd

from acstools import acszpt
from astropy.io import fits
from astropy.time import Time
from datetime import date

from src.photred import sky
from src.photred import io2 as io

# path = '/Volumes/Samsung_T5/hstTreasury'
# dataDir = os.path.join(path,'data')
# resDir = os.path.join(path,'results')
# currDir = os.path.relpath(dataDir,'.')
# outDir = os.path.relpath(resDir,'.')

def getZPT(date,filt,detector='WFC'):
    """ Get date-dependent STmag ZPT """
    q = acszpt.Query(date=date,detector=detector,filt=filt)
    zpt_table = q.fetch()
    try:
        return zpt_table['STmag'].value[0],zpt_table['VEGAmag'].value[0]
    except:
        return zpt_table['STmag'],zpt_table['VEGAmag']


def main(args):

    config = pd.read_json(args.config)

    targname = config.main.targname
    dataDir = config.aper.dataDir
    resDir = config.aper.resDir

    filt = args.filter

    # Directory where the target data is stored
    tDir = os.path.join(dataDir,f'{targname}_f{filt}w')
    
    # Directory that we will output the file in
    oDir = os.path.join(resDir,targname)

    with fits.open(os.path.join(tDir,f'{targname}_f{filt}w.fits')) as hdu:
        expend = Time(hdu[0].header['EXPEND'],format='mjd')
        dateI = expend.to_value('iso',subfmt='date')
        full_name = hdu[0].header['TARGNAME']
        
    im,head = io.readfile(os.path.join(tDir,f'{targname}_f{filt}w.fits'))
    skymode, skysig = sky.getsky(im)
        
    stmag,vegamag = getZPT(dateI,f'F{filt}W')

    today = date.today()

    d1 = today.strftime("%Y-%m-%d")
    
    if not os.path.exists(os.path.join(oDir,f'zeropoints.dat')):
        with open(os.path.join(oDir,f'zeropoints.dat'),'w') as zz:
            zz.write('# TARGNAME FILTER OBS-DATE STMAG VEGAMAG RETR-DATE SKYSIG\n')
        
    with open(os.path.join(oDir,f'zeropoints.dat'),'a') as zz:
        zz.write(f'{full_name} F{filt}W {dateI} {stmag:.3f} {vegamag:.3f} {d1} {skysig:.2f}\n')

    return None


# if __name__=='__main__':
#     targname = sys.argv[1]
# #     zptOut(sys.argv[1],filt='606',dir=os.path.join('../data','sys.argv[1]_f606w'),first=True)
# #     zptOut(sys.argv[1],filt='814',dir=os.path.join('../data',sys.argv[1]_f814w'))

#     zptOut(targname,filt='606',dir=currDir,first=True)
#     zptOut(targname,filt='814',dir=currDir)


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Creates a file with zeropoints and sky sigma values'
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

    df = main(args)

    # args = parser.parse_args('')
    # args.config = '../../config.json'
    # args.filter = '606'
    # main(args)

    # args.filter = '814'
    # main(args)