# Must be in base (treas) environment
# acszpt is particular
import argparse as ap
from datetime import date
import os

from acstools import acszpt
from astropy.io import fits
from astropy.time import Time
import pandas as pd

from src.photred import sky
from src.photred import io2 as io


def get_zpt(date, filt, detector='WFC'):
    """ Get date-dependent STmag ZPT """
    q = acszpt.Query(date=date, detector=detector, filt=filt)
    zpt_table = q.fetch()
    try:
        return zpt_table['STmag'].value[0], zpt_table['VEGAmag'].value[0]
    except:
        return zpt_table['STmag'], zpt_table['VEGAmag']


def out_zpt_sky(targname, filt, data_dir, res_dir):
    with fits.open(os.path.join(data_dir,
                                f'{targname}_f{filt}w.fits')) as hdu:
        expend = Time(hdu[0].header['EXPEND'], format='mjd')
        dateI = expend.to_value('iso', subfmt='date')
        full_name = hdu[0].header['TARGNAME']

    im, _ = io.readfile(os.path.join(data_dir,
                                     f'{targname}_f{filt}w.fits'))
    try:
        _, skysig = sky.getsky(im)
    except ValueError:
        # The im tends to have 3 dimensions, for some reason
        _, skysig = sky.getsky(im[0])

    today = date.today()
    d1 = today.strftime("%Y-%m-%d")
    stmag, vegamag = get_zpt(dateI, f'F{filt}W')

    if not os.path.exists(os.path.join(res_dir, f'zeropoints.dat')):
        with open(os.path.join(res_dir, f'zeropoints.dat'), 'w') as zz:
            zz.write('# TARGNAME FILTER OBS-DATE STMAG ')
            zz.write('VEGAMAG RETR-DATE SKYSIG\n')

    with open(os.path.join(res_dir, f'zeropoints.dat'), 'a') as zz:
        zz.write(f'{full_name} F{filt}W {dateI} {stmag:.3f} ')
        zz.write(f'{vegamag:.3f} {d1} {skysig:.2f}\n')

    return None


def main(args):
    config = pd.read_json(args.config)

    targname = config.main.targname
    filt_arr = [
        f'{config.main.filt1}',
        f'{config.main.filt2}'
    ]

    res_dir = os.path.join(config.aper.res_dir, targname)
    for ff in filt_arr:
        data_dir = os.path.join(config.aper.data_dir, f'{targname}_f{ff}w')
        out_zpt_sky(targname, ff, data_dir, res_dir)

    return None


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Creates a file with zeropoints \
        and sky sigma values'
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
