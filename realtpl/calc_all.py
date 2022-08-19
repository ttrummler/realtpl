import numpy as np
import pandas as pd

from realtpl.fluid_properties import FluidProperties
from realtpl.thermophysical_constants import R_UNIV
from realtpl.eos_data import eos_parameter_from_eos_name
from realtpl.eos_data import alpha_functions_from_eos_name
from realtpl.calc_compressibility import calc_compressibility
from realtpl.calc_cv_cp_sound import calc_cv_cp_sound
from realtpl.calc_visc_cond_chung import calc_visc_cond_chung


def calc_eos_data(eos: str, fp: FluidProperties, temp_array: np.ndarray,
                  pressure_array: np.ndarray):
    """
    calc_eos_data - calculates all thermodynamic quantities base on the eos

    First, current eos data and alpha functions are evaluated based on the eos-
    name and the fluid properties. Then the compressibility z is computed and
    with it the molar volume v and the density rho. Afterwards, the caloric
    properties cp and cv are computed using departure function formalism. In
    this step, also the speed of sound is evaluated. The last step is the
    evaluation of the transport properties viscosity and heat conductivity.

    Parameters:
    -----------
    eos: str
        eos-name
    fp:  FluidProperties
        dataclass with all fluid properties
    temp_array: numpy array
        temperature range in Kelvin
    pressure_array: numpy array
        pressure in Pascal where data is evaluated

    Returns:
    --------
    df: Pandas DataFrame
        dataframe with all the relevant data

    Authors
    ----------
    Trummler, Glatzle

    References:
    -----------
    ..[1] Trummler, Glatzle, Doehring, Urban, Klein (2022), Thermodynamic
    modeling for numerical simulations based on the generalized cubic equation
    of state

    ..[2] Matheis (2018), PhD. Thesis, Numerical Simulation of Fuel
    Injection and Turbulent Mixing Under High-Pressure Conditions,
    http://mediatum.ub.tum.de/?id=1363601

    ...[3] Matheis and Hickel (2017), Multi-component vapor-liquid equilibrium
    model for LES of high-pressure fuel injection and application to ECN
    Spray A. Int. J. Multiph. Flow,
    https://doi.org/10.1016/j.ijmultiphaseflow.2017.11.001
    """

    current_eos_data = eos_parameter_from_eos_name(eos, fp)
    alpha_funcs = alpha_functions_from_eos_name(eos, fp)

    df_temp = pd.DataFrame(columns=['kind', 'press_Pa', 'temp_K', 'rho_kg/m3',
                                    'cp_J/(kgK)', 'sound_m/s', 'visc_Pas',
                                    'cond_W/(mK)'])

    for pressure in pressure_array:
        z = calc_compressibility(current_eos_data, alpha_funcs.alpha,
                                 temp_array, pressure)
        vol = z * R_UNIV * temp_array/pressure
        rho = fp.mass/vol

        cv_cp_sound = calc_cv_cp_sound(fp, temp_array, current_eos_data,
                                       alpha_funcs, vol)

        visc_and_cond = calc_visc_cond_chung(fp, temp_array, rho,
                                             cv_cp_sound[0])

        df_temp_ = pd.DataFrame({
                'kind': eos,
                'press_Pa': pressure,
                'temp_K': temp_array,
                'rho_kg/m3': rho,
                'cp_J/(kgK)': cv_cp_sound[1],
                'sound_m/s': cv_cp_sound[2],
                'visc_Pas': visc_and_cond[0],
                'cond_W/(mK)': visc_and_cond[1]
            })
        df_temp = pd.concat([df_temp, df_temp_])

    return df_temp
