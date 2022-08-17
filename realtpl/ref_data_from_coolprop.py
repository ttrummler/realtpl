import numpy as np
import pandas as pd
import CoolProp as CP


def ref_data_from_coolprop(name: str, temp: np.ndarray, press: np.ndarray):
    """
    Imports reference data from CoolProp

    References:
    -----------
    ..[1] Bell, Wronski, Quoilin, Lemort, “Pure and pseudo-pure fluid
    thermophysical property evaluation and the open-source thermophysical
    property library coolprop,” Industrial & engineering chemistry research
    53, 2498–2508 (2014).
    """
    df = pd.DataFrame(
        columns=['kind', 'press_Pa', 'temp_K', 'rho_kg/m3', 'cp_J/(kgK)',
                 'sound_m/s', 'visc_Pas', 'cond_W/(mK)'])

    kind = 'ref_data'
    heos = CP.AbstractState('HEOS', name)
    for j, press_step in enumerate(press):
        for i, temp_step in enumerate(temp):
            heos.update(CP.PT_INPUTS, press_step, temp_step)
            location = j*len(temp) + i
            df.loc[location] = (
                [kind, press_step, temp_step,
                 heos.rhomass(),
                 heos.cpmass(),
                 heos.speed_sound(),
                 heos.viscosity(),
                 heos.conductivity()]
            )

    return df
