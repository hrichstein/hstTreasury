import argparse as ap
import os

from astroquery.mast import Observations
import pandas as pd


def getObs(targname, proposal_id, instrument, data_dir):
    obs_table = Observations.query_criteria(proposal_id=proposal_id,
                                            obs_collection='HST',
                                            dataproduct_type='IMAGE',
                                            instrument_name=instrument,
                                            target_name=targname)

    prod_list = Observations.get_product_list(obs_table)

    cut_list = Observations.filter_products(prod_list,
                                            productType='SCIENCE',
                                            productSubGroupDescription='FLC',
                                            project='CALACS')

    manifest = Observations.download_products(cut_list,
                                              extension='flc.fits',
                                              download_dir=data_dir)

    print(manifest['Local Path'])

    return None


def main(args):
    config = pd.read_json(args.config)

    targname = config.main.targname
    prop_id = config.obs.prop_id
    instrument = config.obs.instrument
    data_dir = config.script.data_dir

    getObs(targname, prop_id, instrument, data_dir)

    return None


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Download the FLC files\
        from MAST.'
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
