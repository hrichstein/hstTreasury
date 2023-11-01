from astropy.io import ascii,fits
from astropy.table import Table
from glob import glob
from pathlib import Path
import os

filelist = glob('../data/MAST*/HST/*/*flc.fits',recursive=True)
filt_arr = ['606','814']

def mvFLCs(targname,mainDir='../'):

    for ff in filt_arr:
        if not os.path.exists(os.path.join(mainDir,'data',targname,'_f{ff}w')):
             os.mkdir(os.path.join(mainDir,'data',targname,f'_f{ff}w'))
        if not os.path.exists(os.path.join(mainDir,'results',targname)):
            os.mkdir(os.path.join(mainDir,'results',targname))
    for file in filelist:
        fullfile = os.path.relpath(file)
        with fits.open(file,names=True) as hdu:
            if hdu[0].header['TARGNAME'] == targname:
                if hdu[0].header['EXPTIME'] < 200:
                    os.remove(fullfile)
                if hdu[0].header['FILTER1'][0] == 'F':
                    filt= hdu[0].header['FILTER1']
                else:
                    filt = hdu[0].header['FILTER2']
                fileroot = os.path.basename(file)
                for ff in filt_arr:
                    if filt == f'F{ff}W':
                        print(fullfile)
                        os.rename(fullfile,f'../data/{targname}_f{ff}w/{fileroot}')
            else: 
                print(f'Skipping {file}...')
                continue
        
    print('Finished moving the FLCs.')    

    return None
            

if __name__ =='__main__':
    targ = input('Enter the target name (ex. HOROLOGIUM-I): ')
    mvFLCs(targ)