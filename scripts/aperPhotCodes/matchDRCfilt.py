import argparse as ap
from datetime import datetime
from glob import glob
import os

from astropy.io import ascii
from astropy.table import Table
import numpy as np
import pandas as pd

from src.matchlistUNQ import matchlistID

def matchDRCfilt(targname,filt1,filt2,photDir,
                 matchtol,maxtol,suffix='.dat',plot=False):
    """
    Matching sources between DRC filters using 6D-transformed positions
    In the previous code, I transformed sources from filter2 into 
    filter1 space. I will take the original filter1 catalog and the 
    new filter2 catalog.
    """
    xx = os.listdir(photDir)
    for ii in xx:
        if (ii.startswith("._")) and (ii.endswith(suffix)):
            s0 = os.path.join(photDir,ii)
            os.remove(s0)
            print("Deleted ._ files")
        elif ii.endswith(f'aper{filt1}' + suffix):
            filt1_file = ii
            
    filt1f = ascii.read(os.path.join(photDir,filt1_file))
    filt2f = ascii.read(os.path.join(photDir,
                                     f'drcTrans_{targname}_F{filt2}W.dat'))
    filt2f['id'] = np.arange(0,len(filt2f),1,dtype=int)
    
    prime_in = Table()
    prime_in[f'id_{filt1}'] = np.arange(0,len(filt1f),1,dtype=int)
    prime_in[f'x_f{filt1}w'] = filt1f['xcenter']
    prime_in[f'y_f{filt1}w'] = filt1f['ycenter']
    
    p_in = Table(prime_in,copy=True)
    cat = filt2f
    
    lenFilt1 = len(filt1f)
    lenFilt2 = len(filt2f)
    minLen = np.min([lenFilt1,lenFilt2])
    
    nF_out = True
    while nF_out:
        prime, matchids = matchlistID(p_in,cat,matchtol,
                                      f'x_f{filt1}w',f'y_f{filt1}w',
                                      f'x_f{filt1}wTrans',
                                      f'y_f{filt1}wTrans',
                                      'id',tag='matchDRC',
                                      saveDir=photDir,
                                      plot=plot)
        
        if len(prime)>=int(0.2*minLen):
            nF_out = False
            print(f'Minimum Number Reached for {targname}: {len(prime):d}')
        else:
            if matchtol < maxtol+1:
                print('Looking for more matches')    
                print(f'Pixel tolerance: {matchtol:.1f},')
                print(f'Sources: {len(prime):d}')
                matchtol += 1
                p_in = Table(prime_in,copy=True)
            else:
                print(f'Maximum tolerance ({maxtol:.1f} pixels) reached.')
                print(f'Final number of sources: {len(prime):d}')
                nF_out = False
    
    # The percentage of sources matched
    print(f'{len(prime)/minLen*100:.2f}% of sources matched.')  
    print(f'{len(prime):d} matched sources.')
    
    idx1 = np.asarray(prime[f'id_{filt1}'],int)
    reg1 = filt1f[idx1]
    
    idx2 = np.asarray(matchids,int)
    reg2 = filt2f[idx2]
    
    outTable = Table()
    for name in reg1.colnames:
        outTable[f'{name}_f{filt1}w'] = reg1[name].flatten()
    for name in reg2.colnames:
        outTable[f'{name}_f{filt2}w'] = reg2[name].flatten()
        
    ascii.write(outTable,os.path.join(photDir,
                                      f'{targname}_matchedDRC.dat'),
                                      overwrite=True,
                                      format='commented_header')
    
    return None


def main(args):
    plot = args.plot
    config = pd.read_json(args.config)
    if plot != eval(config.aper.plot):
        print(f'Using argument plot value {plot} rather than \
              config file value {eval(config.aper.plot)}.')
        
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
            photDir = max(photDirList, key=lambda p: 
                          datetime.strptime(p, os.path.join(resDir,f'drcPhot%d%b%y')))
            
    matchDRCfilt(targname,filt_arr[0],filt_arr[1],
                 photDir,matchtol=config.aper.filtol,
                 maxtol=config.aper.maxtol,suffix=config.aper.suffix,
                 plot=plot)
    
    return None


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Matches sources between the two filters.'
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
    _ = parser.add_argument(
        '-p',
        '--plot',
        help='Flag for plotting intermediate matching XY figures.',
        default=False,
        type=str,
    )
    args = parser.parse_args()

    raise SystemExit(main(args))