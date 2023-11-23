import argparse as ap
import subprocess as sp

import pandas as pd


# def run_env(env_path, script):
#     # sp.run(
#     #     f"source activate {env_name}",
#     #     shell=True
#     # )
#     # sp.Popen(
#     #     [f'{env_path}/bin/{env_comm}', f'{script}'])
#     sp.run(
#         f"{env_path}/bin/{env_comm} {script}")
#     # os.system(f"source activate {env_name}")
#     # print(f"Switched to conda environment: {env_name}")

#     return None


def main(args):
    config = pd.read_json(args.config)
    main_env = config.main.main_env
    driz_env = config.main.driz_env
    # sfdEnv = config.main.sfd_env

    # env_comm = f"{main_env}/bin/python"
    # # # sp.run(
    # # #     f'{env_comm} get_obs.py',
    # # #     shell=True
    # # # )
    # sp.run(
    #     f'{env_comm} move_flc.py',
    #     shell=True
    # )

    env_comm = f"{driz_env}/bin/python"
    # sp.run(
    #     f'{env_comm} align_gaia.py',
    #     shell=True
    # )
    sp.run(
        f'{env_comm} driz_it.py',
        shell=True
    )

    env_comm = f"{main_env}/bin/python"
    sp.run(
        f'{env_comm} mask_acs.py',
        shell=True
    )
    sp.run(
        f'{env_comm} get_zpt_sky.py',
        shell=True,
        cwd='./aperphot/'
    )

    # sp.run(
    #     f'{env_comm} run_photutils.py',
    #     shell=True,
    #     cwd='./aperphot/'
    # )
    # sp.run(
    #     f'{env_comm} drc_filt_ref.py',
    #     shell=True,
    #     cwd='./aperphot/'
    # )
    # sp.run(
    #     f'{env_comm} drc_filt_trans.py',
    #     shell=True,
    #     cwd='./aperphot/'
    # )
    # sp.run(
    #     f'{env_comm} match_drc_filt.py',
    #     shell=True,
    #     cwd='./aperphot/'
    # )
    # sp.run(
    #     f'{env_comm} make_cmd.py',
    #     shell=True,
    #     cwd='./aperphot/'
    # )


if __name__ == '__main__':
    parser = ap.ArgumentParser(
        description='Run hst_treasury Pipeline.')
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
