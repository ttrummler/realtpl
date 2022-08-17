import numpy as np

from realtpl.thermophysical_constants import R_UNIV
from realtpl.fluid_properties import FluidProperties
from realtpl.eos_data import EosParameter
from realtpl.eos_data import AlphaFunctions
from realtpl.calc_cp_ref_nasa import calc_cp_ref_nasa


def calc_cv_cp_sound(fp: FluidProperties, temp: np.array, ed: EosParameter,
                     alpha_funcs: AlphaFunctions, vol: np.array):
    """
    Calculation of the heat capacities  cv, cp, and the speed of sound with
    departure functions

    Parameters:
    -----------
    fp: FluidProperties
    temp: np.array
    ed: EosParameter
    alpha_funcs: AlphaFunctions
    vol: np.array

    Returns:
    --------
    cv, cp, sound: np.array
       cv, cp, and sound values for temperature array

    Authors
    ----------
    Trummler, Glatzle

    References:
    -----------
    ..[1] Trummler, Glatzle, Doehring, Urban, Klein (2022), Thermodynamic
    modeling for numerical simulations based on the generalized cubic equation
    of state

    .. [2] Matheis (2018), PhD. Thesis, Numerical Simulation of Fuel
     Injection and Turbulent Mixing Under High-Pressure Conditions,
     http://mediatum.ub.tum.de/?id=1363601

    .. [3] Matheis and Hickel (2017), Multi-component vapor-liquid equilibrium
     model for LES of high-pressure fuel injection and application to ECN
     Spray A
     Int. J. Multiph. Flow, doi: 10.1016/j.ijmultiphaseflow.2017.11.001
    """

    # a*alpha(temp) and derivatives
    a_alpha = ed.a*alpha_funcs.alpha(temp)
    d_a_alpha = ed.a*alpha_funcs.d_alpha_d_temp(temp)
    dd_a_alpha = ed.a*alpha_funcs.d2_alpha_d2_temp(temp)

    # denominator of the cubic eos
    denom = (vol**2 + ed.d_1_p_d_2*ed.b*vol + ed.d_1_t_d_2*ed.b**2)

    # dp/dT_c_V dp/dv_c_T
    d_p_d_temp_c_v = (R_UNIV/(vol - ed.b) - d_a_alpha/denom)
    d_p_d_v_c_temp = -(R_UNIV*temp/(vol - ed.b)**2
                       - a_alpha*(2*vol + ed.d_1_p_d_2*ed.b)/denom**2)

    cp_ref = calc_cp_ref_nasa(fp.data_nasa, temp)
    cv_ref = cp_ref - R_UNIV

    right = np.log((vol + ed.b*ed.d_2)/(vol + ed.b*ed.d_1))
    dcv = -temp*dd_a_alpha*right/(ed.b*ed.d_1_m_d_2)

    cv = cv_ref + dcv
    cp = (cv - temp*d_p_d_temp_c_v**2/d_p_d_v_c_temp)/fp.mass
    sound = vol*(-cp/cv*d_p_d_v_c_temp)**0.5

    return cv, cp, sound
