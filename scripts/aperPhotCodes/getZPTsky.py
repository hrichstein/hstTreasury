from datetime import date
import os
import sys

from acstools import acszpt
from astropy.io import fits
from astropy.time import Time
from photred import sky
from photred import io2 as io

path = '/Volumes/Samsung_T5/hstTreasury'
dataDir = os.path.join(path,'data')
resDir = os.path.join(path,'results')
currDir = os.path.relpath(dataDir,'.')
outDir = os.path.relpath(resDir,'.')

def getZPT(date,filt,detector='WFC'):
    """ Get date-dependent STmag ZPT """
    q = acszpt.Query(date=date,detector=detector,filt=filt)
    zpt_table = q.fetch()

    return zpt_table['STmag'].value[0],zpt_table['VEGAmag'].value[0]


def zptOut(targname, filt, dir=currDir, first=False, outDir = outDir):
    # Directory where the target data is stored
    tdir = os.path.join(dir,f'{targname}_f{filt}w')
    
    # Directory that we will output the file in
    odir = os.path.join(outDir,f'{targname}')
    
    with fits.open(os.path.join(tdir,f'{targname}_f{filt}w.fits')) as hdu:
        expend = Time(hdu[0].header['EXPEND'],format='mjd')
        dateI = expend.to_value('iso',subfmt='date')
        full_name = hdu[0].header['TARGNAME']
        
    im,head = io.readfile(os.path.join(tdir,f'{targname}_f{filt}w.fits'))
    skymode, skysig = sky.getsky(im)
        
    stmag,vegamag = getZPT(dateI,f'F{filt}W')

    today = date.today()

    d1 = today.strftime("%Y-%m-%d")
    
    if first:
        with open(os.path.join(odir,f'zeropoints.dat'),'w') as zz:
            zz.write('# TARGNAME FILTER OBS-DATE STMAG VEGAMAG RETR-DATE SKYSIG\n')
        
    with open(os.path.join(odir,f'zeropoints.dat'),'a') as zz:
        zz.write(f'{full_name} F{filt}W {dateI} {stmag:.3f} {vegamag:.3f} {d1} {skysig:.2f}\n')

    return None


if __name__=='__main__':
    targname = sys.argv[1]
#     zptOut(sys.argv[1],filt='606',dir=os.path.join('../data','sys.argv[1]_f606w'),first=True)
#     zptOut(sys.argv[1],filt='814',dir=os.path.join('../data',sys.argv[1]_f814w'))

    zptOut(targname,filt='606',dir=currDir,first=True)
    zptOut(targname,filt='814',dir=currDir)