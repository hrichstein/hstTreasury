# needs to be run in the "driz" environment
from drizzlepac import astrodrizzle
import numpy as np
import os
import shutil
import sys
from datetime import date
from glob import glob

def doDriz(targname,filt='606',mainDir='../',dataDir='data'):

    today = date.today()
    cDate = today.strftime("%d%b")

    dDir = os.path.join(mainDir,dataDir,f'{targname}_f{filt}w')


    flc_list = glob(os.path.join(dDir,'j*_flc.fits'))
    print(flc_list)

    astrodrizzle.AstroDrizzle(flc_list,
                              output=os.path.join(dDir,f'F{filt}W_{cDate}'),
                              configobj='./refFiles/astrodrizzle_new08.cfg',final_scale=0.035,)

    shutil.copy(os.path.join(dDir,f'F{filt}W_{cDate}_drc_sci.fits'),os.path.join(dDir,f'{targname}_f{filt}w.fits'))

    return None


if __name__=='__main__': 
    filt_arr = ['606','814']
    for ff in filt_arr:
        doDriz(sys.argv[1],filt=ff)