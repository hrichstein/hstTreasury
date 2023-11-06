import argparse as ap
from datetime import date
import os
import warnings

from astropy.io import fits,ascii
from astropy.stats import sigma_clipped_stats
import numpy as np
import pandas as pd
from photutils.aperture import aperture_photometry,CircularAnnulus,CircularAperture
from photutils.detection import DAOStarFinder
from photutils.utils import calc_total_error

warnings.filterwarnings('ignore',category=RuntimeWarning)

def runPhotUtils(targname,filt,zpt,sky,eeBand,resDir,
                 fwhm=3.5,aperture=0.2,scale=0.035,native=0.05):
    today = date.today()
    cDate = today.strftime("%d%b%y")
    photDir = os.path.join(resDir,f'drcPhot{cDate}')
    if not os.path.exists(photDir):
        print(f'Creating {photDir}')
        os.makedirs(photDir)
        
    with fits.open(os.path.join(resDir,f'{targname}_f{filt}w.fits')) as hdu:
        data = hdu[0].data
        exptime = hdu[0].header['EXPTIME']
        readA = hdu[0].header['READNSEA']
        readB = hdu[0].header['READNSEB']
        readC = hdu[0].header['READNSEC']
        readD = hdu[0].header['READNSED']
        
    radius = aperture/scale
    rnMean = np.mean(np.array([readA,readB,readC,readD]))
    rnIn = np.sqrt(sky**2 + (rnMean*(scale/native))**2)
    effGain = 1.
    
    if len(data) < 10:  # If the data has dimensions like (1,x,y), this will get it to (x,y)
        data = data[0]
        
    # For getting a flux error based on readnoise
    bkg_arr = np.zeros((data.shape[0],data.shape[1]))    
    bkg_arr.fill(rnIn)
    
    error = calc_total_error(data,bkg_arr,effGain)
    
    _, median, std = sigma_clipped_stats(data, sigma=3.0,
                                         maxiters=10)

    daofind = DAOStarFinder(fwhm=fwhm, threshold=4.0*std, peakmax=9e4)
    sources = daofind(data - median)
    loc = np.array([sources['xcentroid'], sources['ycentroid']])
    positions = np.transpose(loc)
    apertures_rad = CircularAperture(positions, r=radius)
    rawflux_rad = aperture_photometry(data, apertures_rad, error=error)
    rawflux_rad['roundness1'] = sources['roundness1']
    rawflux_rad['roundness2'] = sources['roundness2']
    rawflux_rad['sharpness'] = sources['sharpness']
    
    annulus_apertures = CircularAnnulus(positions, r_in=9.,
                                        r_out=12.)
    annulus_masks = annulus_apertures.to_mask(method='center')
    
    bkg_median = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(data)
        annulus_data_1d = annulus_data[mask.data > 0]
        _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)
    bkg_median = np.array(bkg_median)
    
    rawflux_rad['annulus_median'] = bkg_median
    rawflux_rad['aper_bkg'] = bkg_median * apertures_rad.area
    rawflux_rad['final_phot'] = rawflux_rad['aperture_sum'] \
                                - rawflux_rad['aper_bkg']

    # Looking at the excess flux (r+) and trying to use it as a star/galaxy indicator.
    apertures_rp = CircularAperture(positions, r=radius+(0.1/scale))
    rawflux_rp = aperture_photometry(data, apertures_rp)
    rawflux_rp['aper_bkg'] = rawflux_rad['annulus_median'] \
                             * apertures_rp.area
    rawflux_rp['final_phot'] = rawflux_rp['aperture_sum'] \
                               - rawflux_rp['aper_bkg']
    
    # Looking at concentrated flux (r-)
    apertures_rm = CircularAperture(positions, r=radius-(0.1/scale))
    rawflux_rm = aperture_photometry(data, apertures_rm)
    rawflux_rm['aper_bkg'] = rawflux_rad['annulus_median'] \
                             * apertures_rm.area
    rawflux_rm['final_phot'] = rawflux_rm['aperture_sum'] \
                               - rawflux_rm['aper_bkg']
    
    # Getting items with positive mag differences between normal and r+ radius
    mask_negative = (rawflux_rad['final_phot'] > 0) \
                     & (rawflux_rp['final_phot'] > 0)
    rawflux_pos_rad = rawflux_rad[mask_negative]
    rawflux_pos_rp = rawflux_rp[mask_negative]
    rawflux_pos_rm = rawflux_rm[mask_negative]
    
    # Converting flux to magnitude
    mag_rad = -2.5 * np.log10(rawflux_pos_rad['final_phot'])
    mag_rp = -2.5 * np.log10(rawflux_pos_rp['final_phot'])
    mag_rm = -2.5 * np.log10(rawflux_pos_rm['final_phot'])
    
    # Finding differences between magnitude in the fiducial radius and the larger radius (r+)
    delta_mag = mag_rad - mag_rp
    # Differences must be positive and less than 0.1 mag
    mask_1 = (delta_mag > 0) & (delta_mag < 0.1)
    _, median, std = sigma_clipped_stats(delta_mag[mask_1],
                                         sigma=3.0, maxiters=5)
    # Checking that any one source's difference is within 1.5 std of the median of all sources
    mask_rp = (delta_mag > 0) & (delta_mag < median + 1.5 * std)
    flagCol = np.zeros((len(mask_rp),1),dtype=int)
    for ii in range(len(mask_rp)):
        if mask_rp[ii]:
            flagCol[ii]=int(1)
        else:
            flagCol[ii]=int(0)
    
    rawflux_pos_rad['six_4_flag'] = flagCol.flatten()
    rawflux_pos_rad['four_pix_imag'] = mag_rad.flatten()
    rawflux_pos_rad['six_pix_imag'] = mag_rp.flatten()
    rawflux_pos_rad['two_pix_imag'] = mag_rm.flatten()
    
    # Converting magnitudes from instrumental to observational STMAG
    final_phot = -2.5 * np.log10(rawflux_pos_rad['final_phot']/eeBand) \
                 + zpt + 2.5 * np.log10(exptime)
    flux_err = rawflux_pos_rad['aperture_sum_err']/rawflux_pos_rad['aperture_sum']
    mag_err = 2.5/np.log(10) * flux_err
    rawflux_pos_rad['magr'] = final_phot.flatten()
    rawflux_pos_rad['magr_err'] = mag_err.flatten()
    rawflux_pos_rad['id'] = np.arange(0,len(rawflux_pos_rad),1,dtype=int).flatten()
    
    form_dict = {}
    for name in rawflux_pos_rad.dtype.names:
        form_dict[f'{name}'] = '%.5f'
    
    form_dict['id'] = '%d'
    form_dict['six_4_flag'] = '%d'
    
    ascii.write(rawflux_pos_rad,os.path.join(photDir,f'{targname}_aper{filt}.dat'),
                format='commented_header',formats=form_dict,overwrite=True)

    return None


def main(args):
    config = pd.read_json(args.config)
    
    targname = config.main.targname
    filt_arr = [
        f'{config.main.filt1}',
        f'{config.main.filt2}'
        ]
    ee_arr = [
        config.aper.ee606,
        config.aper.ee814
        ]
    
    resDir = os.path.join(config.aper.resDir,targname)
    paramFile = os.path.join(resDir,config.aper.param)
    
    param = ascii.read(paramFile)
    paraDF = param.to_pandas()
    df = paraDF.set_index('FILTER',drop=False)
    
    for ff,filt in enumerate(filt_arr):
        zpt = df.loc[f"F{filt}W","STMAG"]
        sky = df.loc[f"F{filt}W","SKYSIG"]
        try:
            if len(zpt) > 1:  # If the zeropoint was added to the file multiple times
                zpt = zpt[-1]
                sky = sky[-1]
        except TypeError:
            zpt = zpt
            sky = sky
        eeBand = ee_arr[ff]
        runPhotUtils(targname,filt,zpt,sky,eeBand,resDir,
                     fwhm=config.aper.fwhm,aperture=config.aper.aperrad,
                     scale=config.driz.scale,native=config.driz.native)
    
    return None


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Creates a file with zeropoints and sky sigma values'
    )
    _ = parser.add_argument(
        '-c',
        '--config',
        help='Name of the config json file.\
        (Default: config.json)',
        default='../../config.json',
        type=str,
    )

    args = parser.parse_args()

    raise SystemExit(main(args))
