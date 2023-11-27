import argparse as ap
from datetime import datetime
from glob import glob
import os

from astropy.io import ascii
from astropy.table import Table
import pandas as pd


def output_cat(targname, filt_arr, phot_dir, suffix='.dat'):
    match_file = os.path.join(phot_dir,
                              f'{targname}_matched_drc{suffix}')

    cat = ascii.read(match_file)
    try:
        _ = cat[f'RA_{filt_arr[0]}']
    except:
        print('No coordinates found. Please run get_coords.py routine.')
        return None

    out_table = Table()
    form_dict = {}
    for ff in filt_arr:
        out_table[f'id_{ff}'] = cat[f'id_{ff}']
        form_dict[f'id_{ff}'] = '%d'

        out_table[f'RA_{ff}'] = cat[f'RA_{ff}']
        form_dict[f'RA_{ff}'] = '%1.7f'
        out_table[f'DEC_{ff}'] = cat[f'DEC_{ff}']
        form_dict[f'DEC_{ff}'] = '%1.7f'

        out_table[f'x_{ff}'] = cat[f'xcenter_{ff}']
        form_dict[f'x_{ff}'] = '%1.4f'
        out_table[f'y_{ff}'] = cat[f'ycenter_{ff}']
        form_dict[f'y_{ff}'] = '%1.4f'

        out_table[f'magr_{ff}'] = cat[f'magr_{ff}']
        form_dict[f'magr_{ff}'] = '%1.4f'
        out_table[f'magr_err_{ff}'] = cat[f'magr_err_{ff}']
        form_dict[f'magr_err_{ff}'] = '%1.5f'

        out_table[f'flag_{ff}'] = cat[f'six_4_flag_{ff}']
        form_dict[f'flag_{ff}'] = '%d'

    ascii.write(out_table, os.path.join(phot_dir,
                                        f'{targname}_aper_cat{suffix}'),
                overwrite=True,
                format='commented_header',
                formats=form_dict)

    return None


def main(args):
    config = pd.read_json(args.config)

    targname = config.main.targname
    suffix = config.aper.suffix
    filt_arr = [
        f'{config.main.filt1}',
        f'{config.main.filt2}'
    ]

    res_dir = os.path.join(config.script.res_dir, targname)
    if args.date is not None:
        phot_dir = os.path.join(res_dir, f'drc_phot{args.date}')
    else:
        phot_dir_list = glob(os.path.join(res_dir, f'drc_phot*'))
        if len(phot_dir_list) == 1:
            phot_dir = phot_dir_list[0]
        else:
            phot_dir = max(phot_dir_list, key=lambda p:
                           datetime.strptime(p, os.path.join(res_dir, f'drc_phot%d%b%y')))

    output_cat(targname, filt_arr, phot_dir, suffix=suffix)

    return None


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Outputs a catalog with the \
        most relevant information.'
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
        '-d',
        '--date',
        help='Date of aperture photometry in format \
        DDMonYY (01Jan24).',
        type=str,
    )
    _ = parser.add_argument(
        '-p',
        '--psf',
        help='Flag for denoting PSF input catalog.',
        default='False',
        type=str,
    )
    args = parser.parse_args()

    raise SystemExit(main(args))
