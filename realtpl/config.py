import numpy as np
import os
import yaml
import warnings

_CFG_DEFAULT = {'eos_list': ['SRK', 'PR', 'RKPR'],
                'include_ref_data': True,
                'temperature_step_K': 1,
                'pressure_step_Pa': 1e5,
                'n_nasa_coeff': 7,
                'output_dir': 'results',
                'save_data_to_csv': True,
                'show_plots': True,
                'save_plots': False,
                'show_deviation': False,
                'save_deviation': False,
                'performance_tracking': False}


def load_config(args):
    cfg = _CFG_DEFAULT.copy()
    file = args['config_file']
    with open(file, 'r') as f:
        cfg_user = yaml.safe_load(f)
    cfg.update(cfg_user)

    # check config for mandatory input
    for key in ['fluid_name',
                'temperature_start_K',
                'temperature_end_K']:
        if not cfg.get(key):
            raise RuntimeError(f'{key} is mandatory in config file. \n'
                               f'Revise the config file {file}.')

    if not cfg.get('pressure_Pa') and not cfg.get('pressure_start_Pa'):
        raise RuntimeError(f'Either pressure_Pa or pressure_start_Pa has to be'
                           f' provided in the config file.\n'
                           f'Revise the config file {file}.')

    if not cfg.get('pressure_start_Pa'):
        cfg['pressure_start_Pa'] = cfg['pressure_Pa']
    if not cfg.get('pressure_Pa'):
        cfg['pressure_Pa'] = cfg['pressure_start_Pa']
    if not cfg.get('pressure_end_Pa'):
        cfg['pressure_end_Pa'] = cfg['pressure_Pa']

    # make sure expo numbers (e.g. 1e5) are interpreted as floats
    for key in ['temperature_start_K',
                'temperature_end_K',
                'temperature_step_K',
                'pressure_Pa',
                'pressure_start_Pa',
                'pressure_end_Pa',
                'pressure_step_Pa']:
        cfg[key] = float(cfg[key])

    # check consistency of config input
    if cfg['temperature_start_K'] > cfg['temperature_end_K']:
        raise RuntimeError(f'wrong input: temperature_start_K >'
                           ' temperature_end_K.\n'
                           f'Revise the config file {file}.')

    if cfg['pressure_start_Pa'] > cfg['pressure_end_Pa']:
        raise RuntimeError(f'wrong input: pressure_start_Pa >'
                           ' pressure_end_Pa.\n'
                           f'Revise the config file {file}.')

    if (cfg['pressure_end_Pa'] > cfg['pressure_start_Pa'] and
            (cfg['show_plots'] or cfg['save_plots'] or
             cfg['show_deviation'] or cfg['save_deviation'] or
             not cfg['save_data_to_csv'])):
        raise RuntimeError(f'wrong input: pressure array (e.g. pressure_end_Pa'
                           f' > pressure_Pa or pressure_end_Pa > '
                           f'pressure_start_Pa) does not work with show/save '
                           f'plots and deviation, but requires save data to '
                           f'csv. \n'
                           f'Revise the config file {file}. ')

    if ((cfg['show_deviation'] or cfg['save_deviation']) and
            not cfg['include_ref_data']):
        raise RuntimeError(f'Deviation can only be evaluated with '
                           f'include_ref_data. Revise the config file {file}.')

    cfg['temp_array'] = np.arange(
        cfg['temperature_start_K'],
        cfg['temperature_end_K'] + cfg['temperature_step_K'],
        cfg['temperature_step_K']
    )
    cfg['pressure_array'] = np.arange(
        cfg['pressure_start_Pa'],
        cfg['pressure_end_Pa'] + cfg['pressure_step_Pa'],
        cfg['pressure_step_Pa']
    )

    if len(cfg['temp_array'])*len(cfg['pressure_array']) > 1e9:
        warnings.warn('More than 1e9 pressure/temperature data points!\n'
                      'Calculation will take more than about 30 minutes!')

    return cfg


def write_config(args):
    outdir = os.path.join(args['output_dir'], args['fluid_name'])
    os.makedirs(outdir, exist_ok=True)

    with open(os.path.join(outdir, 'config_data.out'), 'w') as f:
        for key, value in args.items():
            f.write('%s: %s\n' % (key, value))
