""" Code to do 6D linear transform """
import argparse as ap
from datetime import datetime
from glob import glob
import os

from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.linear6d import test_linear


def mk_plot(targname, x1, y1, x2, y2, x3, y3, label_1,
            label_2, label_3, outname=None):

    fig, ax = plt.subplots(figsize=(6, 6))

    ax.scatter(x1, y1, label=label_1, s=60)
    ax.scatter(x2, y2, label=label_2, s=25)
    ax.scatter(x3, y3, label=label_3, s=10)

    ax.legend()
    ax.set_title(targname)
    ax.set_rasterization_zorder(1)

    plt.savefig(f'{outname}.pdf', dpi=150, bbox_inches='tight')
    plt.close()

    return None


def drc_filt_trans(targname, filt1, filt2, phot_dir, suffix='.dat'):
    xx = os.listdir(phot_dir)
    for ii in xx:
        if (ii.startswith("._")) and (ii.endswith(suffix)):
            s0 = os.path.join(dir, ii)
            os.remove(s0)
            print("Deleted ._ files")
        elif ii.endswith(f'aper{filt2}' + suffix):
            filt2_file = ii

    # Going from filter 2 to filter 1
    ref = ascii.read(os.path.join(phot_dir,
                                  f'drc_filt_ref_{targname}.dat'))
    all = ascii.read(os.path.join(phot_dir, filt2_file))

    # Putting the filter 2 positions into the match array
    # to be used as refs
    match_arr = np.zeros((len(ref), 2))
    match_arr[:, 0] = ref[f'x_f{filt2}w'].flatten()
    match_arr[:, 1] = ref[f'y_f{filt2}w'].flatten()

    # Putting the filter 1 positions into the prime array
    prime_arr = np.zeros((len(ref), 2))
    prime_arr[:, 0] = ref[f'x_f{filt1}w'].flatten()
    prime_arr[:, 1] = ref[f'y_f{filt1}w'].flatten()

    # Weights required for the 6D transformation function
    # We haven't assigned useful values to them.
    weights = np.ones(len(prime_arr))

    # Creating an array that will hold the original
    # and transformed values
    all_arr = np.zeros((len(all), 2))
    all_arr[:, 0] = all['xcenter'].flatten()
    all_arr[:, 1] = all['ycenter'].flatten()

    print(f'Transforming {targname}')
    # This takes the xy positions from the points you have matched
    # to the prime, the xy positions of the corresponding points in
    # the prime, weights, and the xy points in the match_arr frame
    # that will be transformed into the prime_arr frame.
    new_match, new_all = test_linear(match_arr[:, 0], match_arr[:, 1],
                                     prime_arr[:, 0], prime_arr[:, 1],
                                     weights, weights,
                                     all_arr[:, 0], all_arr[:, 1])

    all[f'x_f{filt1}w_trans'] = np.round(new_all[:, 0], 5).flatten()
    all[f'y_f{filt1}w_trans'] = np.round(new_all[:, 1], 5).flatten()

    outname = os.path.join(phot_dir, f'drc_trans_{targname}_F{filt2}W')
    ascii.write(all, f'{outname}.dat', overwrite=True,
                format='commented_header')

    mk_plot(targname, match_arr[:, 0], match_arr[:, 1],
            prime_arr[:, 0], prime_arr[:, 1],
            new_match[:, 0], new_match[:, 1],
            label_1=f'Original in F{filt2}W',
            label_2=f'Original in F{filt1}W',
            label_3=f'New in F{filt2}W to F{filt1}W',
            outname=f'{outname}_match_check')

    return None


def main(args):
    config = pd.read_json(args.config)

    targname = config.main.targname
    filt_arr = [
        f'{config.main.filt1}',
        f'{config.main.filt2}'
    ]

    res_dir = os.path.join(config.aper.res_dir, targname)
    if args.date is not None:
        phot_dir = os.path.join(res_dir, f'drc_phot{args.date}')
    else:
        phot_dir_list = glob(os.path.join(res_dir, f'drc_phot*'))
        if len(phot_dir_list) == 1:
            phot_dir = phot_dir_list[0]
        else:
            phot_dir = max(phot_dir_list,
                           key=lambda p: datetime.strptime(p,
                                                           os.path.join(res_dir,
                                                                        f'drc_phot%d%b%y')))

    print(f'Linearly transforming DRC catalogs in {phot_dir}')
    drc_filt_trans(targname, filt_arr[0], filt_arr[1],
                   phot_dir, suffix=config.aper.suffix)

    return None


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Calculates and outputs transformations of \
        sources in filter 2 to filter 1.'
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
    args = parser.parse_args()

    raise SystemExit(main(args))
