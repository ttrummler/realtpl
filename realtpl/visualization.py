from dataclasses import dataclass

import pandas as pd
import matplotlib.pyplot as plt

import os

# color cycle for plots
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=['r', 'g', 'b', 'k'])


def vis_data(df: pd.DataFrame, fp: dataclass, save_fig: bool, show_fig: bool,
             output_dir: str):
    pressure_MPa = df['press_Pa'].max() / 1e6
    pressure_rel = round(df['press_Pa'].max() / fp.p_c, 2)
    title = (fp.name + ' at p =' + str(pressure_MPa) + ' MPa '
             + '(p/p_c = ' + str(pressure_rel) + ')')

    for item in ['rho_kg/m3', 'cp_J/(kgK)', 'sound_m/s', 'visc_Pas',
                 'cond_W/(mK)']:
        fig, ax = plt.subplots(figsize=(8, 5))

        for kind, dff in df.groupby('kind'):
            ax.plot(dff['temp_K'], dff[item], label=kind)
        plt.legend()
        plt.xlabel('temperature [K]')
        plt.ylabel(str(item.split("_")[0]) + ' [' + str(item.split("_")[1]) +
                   ']')
        plt.xlim([df['temp_K'].min(), df['temp_K'].max()])
        plt.grid(True)
        plt.title(title)

        if save_fig:
            path_graphics = os.path.join(output_dir, fp.name, 'graphics')
            os.makedirs(path_graphics, exist_ok=True)
            fig.savefig(os.path.join(path_graphics,
                                     str(item.split("_")[0]) + '.png'),
                        dpi=200)
    if show_fig:
        plt.show()


def vis_deviation(df: pd.DataFrame, output_dir: str, name: str,
                  flag_show: bool, flag_save: bool):
    ref = df[df['kind'] == 'ref_data']

    for item in ['rho_kg/m3', 'cp_J/(kgK)', 'sound_m/s', 'visc_Pas',
                 'cond_W/(mK)']:
        fig, ax = plt.subplots(figsize=(8, 5))

        for kind, dff in df.groupby('kind'):
            ax.plot(dff['temp_K'], (dff[item] - ref[item]) / ref[item] * 100,
                    label=kind)
        plt.legend()
        plt.xlabel('temperature [K]')
        plt.ylabel('deviation ' + str(item.split("_")[0]) + ' [%]')
        plt.xlim([df['temp_K'].min(), df['temp_K'].max()])
        plt.grid(True)
        if flag_save:
            path_graphics = os.path.join(output_dir, name, 'graphics')
            os.makedirs(path_graphics, exist_ok=True)
            fig.savefig(os.path.join(path_graphics, 'deviation_'
                                     + str(item.split("_")[0]) + '.png'),
                        dpi=200)

    if flag_show:
        plt.show()
