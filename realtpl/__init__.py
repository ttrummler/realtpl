import argparse
import pandas as pd
import time
import warnings

from realtpl import config
from realtpl import fluid_properties
from realtpl import nasa
from realtpl.fluid_properties \
    import fluid_properties_from_coolprop_and_data_base
from realtpl.ref_data_from_coolprop import ref_data_from_coolprop
from realtpl.calc_all import calc_eos_data
from realtpl.visualization import vis_data, vis_deviation
from realtpl.write_data_to_files import write_csv

# do not provide anything for * imports
__all__ = []

_description = ('Computes thermodynamic quantities using a thermodynamic model'
                ' based on cubic equations of state (SRK, PR, RKPR).'
                ' Additionally reference values are extracted from CoolProp.')
_parser = argparse.ArgumentParser(_description)
_parser.add_argument(
    '--config-file',
    help='Path to configuration file.',
    default='config.yaml'
)


def main():
    args = vars(_parser.parse_args())

    time_start = time.process_time()

    # read and check config file
    cfg = config.load_config(args)
    # write config data to file
    config.write_config(cfg)

    # nasa data for cp calculation
    data_nasa = nasa.NasaCoefficients.from_name_and_coeff(
        cfg['fluid_name'], cfg['n_nasa_coeff'])

    _check_temp_range(data_nasa, cfg)

    # get fluid properties
    fp = fluid_properties_from_coolprop_and_data_base(cfg['fluid_name'],
                                                      data_nasa)
    fluid_properties.save_fp_to_file(fp, cfg['output_dir'])

    # df is main data frame
    df = pd.DataFrame(
        columns=['kind', 'press_Pa', 'temp_K', 'rho_kg/m3', 'cp_J/(kgK)',
                 'sound_m/s', 'visc_Pas', 'cond_W/(mK)'])

    time_after_setup = time.process_time()

    # ref data
    if cfg['include_ref_data']:
        df_ref = ref_data_from_coolprop(cfg['fluid_name'], cfg['temp_array'],
                                        cfg['pressure_array'])
        df = pd.concat([df, df_ref])

    time_after_ref = time.process_time()

    # eos data
    for eos in cfg['eos_list']:
        df_add = calc_eos_data(eos, fp, cfg['temp_array'],
                               cfg['pressure_array'])
        df = pd.concat([df, df_add])

    time_after_eos = time.process_time()

    # plot and optionally save fig
    if cfg['save_plots'] or cfg['show_plots']:
        vis_data(df, fp, cfg['save_plots'], cfg['show_plots'],
                 cfg['output_dir'])

    # optionally: show and/or save deviation
    if cfg['show_deviation'] or cfg['save_deviation']:
        vis_deviation(df, cfg['output_dir'], cfg['fluid_name'],
                      cfg['show_deviation'], cfg['save_deviation'])

    time_after_figs = time.process_time()

    # save data to csv
    if cfg['save_data_to_csv']:
        write_csv(df, fp, cfg['output_dir'])

    time_after_save = time.process_time()

    if cfg['performance_tracking']:
        filename = 'performance.out'
        print(f'Performance evaluation saved to {filename}')
        with open(filename, 'w') as f:
            performance_data = {
                'num_eval': len(cfg['temp_array']),
                'time_total [s]': time_after_save - time_start,
                'time_setup [s]': time_after_setup - time_start,
                'time_ref_data [s]': time_after_ref - time_after_setup,
                'time_eos_data [s]': time_after_eos - time_after_ref,
                'time_figs [s]': time_after_figs - time_after_eos,
                'time_write_csv [s]': time_after_save - time_after_figs
            }
            for key, value in performance_data.items():
                f.write('%s: %s\n' % (key, value))

    print('...successfully finished')


def _check_temp_range(data_nasa, cfg):
    data_nasa_temp_range = data_nasa.get_temp_range()
    if data_nasa_temp_range[0] > cfg['temperature_start_K'] \
            or data_nasa_temp_range[1] < cfg['temperature_end_K']:
        warnings.warn(f'Nasa coefficients not valid for entire temperature '
                      f'range! \n'
                      f'Coefficients provided for {data_nasa_temp_range[0]} '
                      f'to {data_nasa_temp_range[1]} K.')
