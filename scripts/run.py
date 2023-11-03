import argparse as ap
import subprocess as sp

import pandas as pd

def switch_conda_environment(env_name):
    sp.run(
        f"conda activate {env_name}",
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
        'cd aperPhotCodes/'    
    )
