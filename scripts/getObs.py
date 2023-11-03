import argparse as ap
import os

from astroquery.mast import Observations
import pandas as pd

def getObs(targname,dataDir):
    obs_table = Observations.query_criteria(proposal_id='14734',
                                        obs_collection='HST',
                                        dataproduct_type='IMAGE',
                                        instrument_name='ACS/WFC',
                                        target_name=targname)
    
    prod_list = Observations.get_product_list(obs_table)
    
    cut_list = Observations.filter_products(prod_list, productType='SCIENCE',
                                       productSubGroupDescription='FLC',
                                       project='CALACS')
    
    manifest = Observations.download_products(
        cut_list,extension='flc.fits',download_dir=dataDir)
    
    print(manifest['Local Path'])
    
    return None


def main(args):
    config = pd.read_json(args.config)
    
    targname = config.main.targname
    dataDir = config.script.dataDir
    
    getObs(targname,dataDir)
    
    return None


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Drizzle the FLC files'
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
