import numpy as np

from realtpl.thermophysical_constants import R_UNIV


def calc_cp_ref_nasa(data_nasa, temp: np.array):
    """
    calculates ideal reference state from NASA polynomials for specific heat
    at constant pressure.

    Parameters:
    -----------
    data_nasa: dict
        data with the nasa coefficients
    temp: np.array
        temperature array in Kelvin

    Returns:
    --------
    cp_ref: np.array
       reference cp values for temperature array

    Authors:
    ---------
    Trummler, Glatzle

    References:
    ------------
    .. [1] E. Goos, A. Burcat, and B. Ruscic. Third Millennium Ideal Gas
          and Condensed Phase Thermochemical Database for Combustion, 2009.
          URL http://burcat. technion.ac.il/dir/.
    .. [2] http://combustion.berkeley.edu/gri-mech/data/nasa_plnm.html
    .. [3] Matheis (2018), PhD. Thesis, Numerical Simulation of Fuel
     Injection and Turbulent Mixing Under High-Pressure Conditions,
     http://mediatum.ub.tum.de/?id=1363601
    .. [4] Matheis and Hickel (2017), Multi-component vapor-liquid equilibrium
     model for LES of high-pressure fuel injection and application to ECN
     Spray A
     Int. J. Multiph. Flow, doi: 10.1016/j.ijmultiphaseflow.2017.11.001
    """

    temp_2 = temp*temp
    temp_3 = temp_2*temp
    temp_4 = temp_3*temp

    if data_nasa.n_coeff == 7:
        cp = (data_nasa.get_coeff(0, temp)
              + data_nasa.get_coeff(1, temp)*temp
              + data_nasa.get_coeff(2, temp)*temp_2
              + data_nasa.get_coeff(3, temp)*temp_3
              + data_nasa.get_coeff(4, temp)*temp_4)
    elif data_nasa.n_coeff == 9:
        temp_inv = 1/temp
        temp_inv_2 = temp_inv/temp
        cp = (data_nasa.get_coeff(0, temp)*temp_inv_2
              + data_nasa.get_coeff(1, temp)*temp_inv
              + data_nasa.get_coeff(2, temp)
              + data_nasa.get_coeff(3, temp)*temp
              + data_nasa.get_coeff(4, temp)*temp_2
              + data_nasa.get_coeff(5, temp)*temp_3
              + data_nasa.get_coeff(6, temp)*temp_4)
    else:
        raise ValueError(f"Unknown n_coeff: {data_nasa.n_coeff}.")
    return cp*R_UNIV
