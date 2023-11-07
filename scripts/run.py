import argparse as ap
import subprocess as sp

import pandas as pd

def switch_conda_environment(env_name):
    sp.run(
        f"source activate {env_name}",
        shell=True
    )
    print(f"Switched to conda environment: {env_name}")
    
    return None

def main(args):
    config = pd.read_json(args.config)
    mainEnv = config.main.mainEnv
    drizEnv = config.main.drizEnv
    # sfdEnv = config.main.sfdEnv
    
    switch_conda_environment(mainEnv)
    sp.run(
        'python getObs.py',
        shell=True
    )
    sp.run(
        'python moveFLCs.py',
        shell=True
    )
    
    switch_conda_environment(drizEnv)
    sp.run(
        'python alignGaia.py',
        shell=True
    )
    sp.run(
        'python drizIt.py',
        shell=True
    )
    
    switch_conda_environment(mainEnv)
    sp.run(
        'python maskACS.py',
        shell=True
    )
    
    sp.run(
        'python getZPTsky.py',
        shell=True, 
        cwd='./aperPhotCodes/'
    )
    sp.run(
        'python runPhotUtils.py',
        shell=True, 
        cwd='./aperPhotCodes/'
    )
    sp.run(
        'python getDRCfiltRef.py',
        shell=True, 
        cwd='./aperPhotCodes/'
    )
    sp.run(
        'python drcFiltLinTrans.py',
        shell=True, 
        cwd='./aperPhotCodes/'
    )
    sp.run(
        'python matchDRCfilt.py',
        shell=True, 
        cwd='./aperPhotCodes/'
    )
    sp.run(
        'python makeCMD.py',
        shell=True, 
        cwd='./aperPhotCodes/'
    )
    
    
if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Run hstTreasury Pipeline.')
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