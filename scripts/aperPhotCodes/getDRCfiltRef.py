import argparse as ap
from datetime import datetime
from glob import glob
import os

from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.matchlistUNQ import matchlistID

def getDRCfiltRef(targname,filt1,filt2,photDir,matchtol,slice,suffix='.dat',plot=True):
    """
    Getting the reference stars between the two filters so that I can do a 6D
    linear transform.
    """
    xx = os.listdir(photDir)
    for ii in xx:
        if (ii.startswith("._")) and (ii.endswith(suffix)):
            s0 = os.path.join(photDir,ii)
            os.remove(s0)
            print('Deleted ._ files')
        elif ii.endswith(f'aper{filt1}' + suffix):
            filt1_file = ii
        elif ii.endswith(f'aper{filt2}' + suffix):
            filt2_file = ii
            
    filt1f = ascii.read(os.path.join(photDir,filt1_file))
    filt2f = ascii.read(os.path.join(photDir,filt2_file))
    
    # Taking the brightest stars, to some 'count' number (~500 usually good)
    cut1 = np.argsort(filt1f['magr'])[:slice]
    f1_cut = filt1f[cut1]
    
    cut2 = np.argsort(filt2f['magr'])[:slice]
    f2_cut = filt2f[cut2]
    f2_cut['id'] = np.arange(0,len(f2_cut),1,dtype=int).flatten()
    
    prime_in = Table()
    prime_in[f'x_f{filt1}w'] = f1_cut['xcenter'].flatten()
    prime_in[f'y_f{filt1}w'] = f1_cut['ycenter'].flatten()
    
    p_in = Table(prime_in,copy=True)
    cat = f2_cut
    
    nF_out = True  # not False
    while nF_out:
        prime, matchids = matchlistID(p_in,cat,matchtol,f'x_f{filt1}w',
                                      f'y_f{filt1}w','xcenter','ycenter',
                                      'id',tag='DRCref',saveDir=photDir,
                                      plot=plot)
        
        if len(prime)>=int(12):  # Because it's going to be a 6D transformation, let's get at least double
            nF_out = False
            print(f'Minimum Number Reached for {targname}: {len(prime):d}')
        else:
            print(f'Need more stars, pixel tolerance at {matchtol:d}')
            p_in = Table(prime_in,copy=True)
            matchtol += 0.5
            if matchtol >= 4:
                nF_out = False
                print('Not enough matches with pixel tolerance at {matchtol:d}.')
                print('Look for other issues.')
                return None
                
    idxCol = np.asarray(matchids,int)
    pts2 = f2_cut[idxCol]
    
    prime[f'x_f{filt2}w'] = pts2['xcenter'].flatten()
    prime[f'y_f{filt2}w'] = pts2['ycenter'].flatten()
    
    for name in prime.colnames:
        prime[name] = prime[name].flatten()
        
    outName = os.path.join(photDir,f'drcFiltRef_{targname}')
    ascii.write(prime,f'{outName}.dat',overwrite=True,format='commented_header')
    
    # Reference quality check plot
    fig, ax = plt.subplots(figsize=(6,6))

    ax.scatter(prime[f'x_f{filt1}w'],prime[f'y_f{filt1}w'],label=f'F{filt1}W',s=50,color='black')
    ax.scatter(prime[f'x_f{filt2}w'],prime[f'y_f{filt2}w'],label=f'F{filt2}W',s=20,color='magenta')

    ax.legend()
    ax.set_title(targname)
    ax.set_rasterization_zorder(1)
    plt.savefig(f'{outName}.pdf',dpi=150,bbox_inches='tight')
    plt.close()

    return None


def main(args):
    plot = args.plot
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
            photDir = max(photDirList, key=lambda p: datetime.strptime(p, os.path.join(resDir,f'drcPhot%d%b%y')))

    print(f'Getting DRC references for files in {photDir}')    
    getDRCfiltRef(targname,filt_arr[0],filt_arr[1],photDir,matchtol=config.aper.drctol,
                  slice=config.aper.slice,suffix=config.aper.suffix,plot=plot)

    return None


if __name__ == '__main__':
    parser = ap.ArgumentParser(description='Find reference stars in DRC images.')
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
        help='Date of aperture photometry in format DDMonYY (01Jan24).',
        type=str,
    )
    _ = parser.add_argument(
        '-p',
        '--plot',
        help='Plot intermediate/end XY of reference stars.',
        default=True,
        type=str,
    )
    args = parser.parse_args()

    raise SystemExit(main(args))