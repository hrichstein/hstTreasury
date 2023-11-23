import argparse as ap
from glob import glob
import os
import shutil

from astropy.io import fits
import pandas as pd

# To do: Add functionality to delete empty MAST folders


def move_flc(targname, filt_arr, exp_cut, data_dir, res_dir):
    filelist = glob(f'{data_dir}/mast*/**/*flc.fits', recursive=True)

    # If the files were downloaded via curl method:
    if len(filelist) == 0:
        filelist = glob(f'{data_dir}/MAST*/**/*flc.fits', recursive=True)

    if len(filelist) == 0:
        print('No FLCs to be moved!')
        return None

    if not os.path.exists(os.path.join(res_dir, targname)):
        os.mkdir(os.path.join(res_dir, targname))
        os.mkdir(os.path.join(res_dir, 'qual_figs', targname))

    for filt in filt_arr:
        if not os.path.exists(os.path.join(data_dir, f'{targname}_f{filt}w')):
            os.mkdir(os.path.join(data_dir, f'{targname}_f{filt}w'))

    for file in filelist:
        with fits.open(file, names=True) as hdu:
            if hdu[0].header['TARGNAME'] == targname:
                if hdu[0].header['EXPTIME'] < exp_cut:
                    os.remove(file)
                    print('Removed short exposure')
                else:
                    if hdu[0].header['FILTER1'][0] == 'F':
                        filter = hdu[0].header['FILTER1']
                    else:
                        filter = hdu[0].header['FILTER2']
                    fileroot = os.path.basename(file)
                    ff = filter[1:4]
                    os.rename(file, os.path.join(data_dir,
                                                 f'{targname}_f{ff}w',
                                                 fileroot))

            else:
                print(f'Skipping {file}...')
                continue

    # shutil.rmtree('mastDownload')

    print('Finished moving the FLCs.')

    return None


def main(args):
    config = pd.read_json(args.config)

    targname = config.main.targname
    exp_cut = config.obs.exp_cut
    data_dir = config.script.data_dir
    res_dir = config.script.res_dir
    filt_arr = [
        f'{config.main.filt1}',
        f'{config.main.filt2}'
    ]

    # for ff in filt_arr:
    move_flc(targname, filt_arr, exp_cut, data_dir, res_dir)

    return None


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Move downloaded FLC files'
    )
    _ = parser.add_argument(
        '-c',
        '--config',
        help='Name of the config json file.\
        (Default: config.json)',
        default='../config.json',
        type=str,
    )
    args = parser.parse_args()

    raise SystemExit(main(args))
