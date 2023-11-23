import argparse as ap
import os

import pandas as pd

from src.photred import mkopt


def main(args):
    config = pd.read_json(args.config)

    targname = config.main.targname
    hilim = config.obs.hlim
    varpsf = config.obs.varpsf
    scale = config.driz.scale

    filt_arr = [
        f'{config.main.filt1}',
        f'{config.main.filt2}'
    ]

    res_dir = os.path.join(config.script.res_dir, targname)

    for ff in filt_arr:
        fname = os.path.join(res_dir, f'{targname}_f{ff}w.fits')
        mkopt(fname, hilimit=hilim, va=varpsf, scale=scale)

    return None


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Make the .opt files for DAOPHOT.'
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
