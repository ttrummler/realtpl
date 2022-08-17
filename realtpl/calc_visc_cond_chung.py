import numpy as np
from realtpl.fluid_properties import FluidProperties
from realtpl.thermophysical_constants import J_PER_CAL, R_MOL


def calc_visc_cond_chung(fp: FluidProperties,
                         temp: np.array,
                         rho_kg_p_m3: np.array,
                         cv_joule_p_kmol_p_kelvin: np.array):
    """
    visc_and_cond_chung - calculates viscosity and heat conductivity for
    dense fluids based on paper by Chung et al. 1988.

    Parameters:
    -----------
    fp: FluidProperties
    temp: np.array
    rho_kg_p_m3: np.array
    cv_joule_p_kmol_p_kelvin: np.array

    Returns:
    --------
    visc, cond: np.array
       visc and cond values for temperature array

    Authors
    ----------
    Trummler, Glatzle

    References
    ----------
    ..[1] Chung, Ajlan, Lee, Starling. (1988).
    Generalized Multiparameter Correlation for Nonpolar and Polar Fluid
    Transport Properties. Ind. Eng. Chem. Res., 27, 671?679.
    http://doi.org/10.1021/ie00076a024
    """

    # Convert input data to correct units
    v_c_cm3_p_mol = fp.v_c*1e3  # cm3/mol
    rho_mol_p_cm3 = rho_kg_p_m3/fp.mass*1e-3  # mol/cm3
    cv_cal_p_mol_p_kelvin = (cv_joule_p_kmol_p_kelvin
                             / (1000*J_PER_CAL))  # cal/mol K

    # Flag for extended calculation (hydrogen bounding, dipole)
    extended_calc = False
    if fp.dipole_moment != 0 or fp.association_parameter != 0:
        extended_calc = True

    # get reduced dipole moment
    mu_r = 131.3*fp.dipole_moment/(v_c_cm3_p_mol*fp.temp_c)**0.5

    # Constants for A and B
    a0 = np.array([6.32402, 0.12102e-2, 5.28346, 6.62263, 19.74540,
                   -1.89992, 24.27450, 0.79716, -0.23816, 0.68629e-1])
    a1 = np.array([50.41190, -0.11536e-2, 254.20900, 38.09570, 7.63034,
                   -12.53670, 3.44945, 1.11764, 0.67695e-1, 0.34793])
    a2 = np.array([-51.68010, -0.62571e-2, -168.481, -8.46414, -14.35440,
                   4.98529, -11.29130, 0.12348e-1, -0.81630, 0.59256])
    a3 = np.array([1189.020, 0.37283e-1, 3898.27, 31.4178, 31.5267,
                   -18.15070, 69.3466, -4.11661, 4.02528, -0.72663])
    b0 = np.array([2.41657, -0.50924, 6.61069, 14.54250, 0.79274, -5.86340,
                   81.17100])
    b1 = np.array([0.74824, -1.50936, 5.62073, -8.91387, 0.82019, 12.80050,
                   114.15800])
    b2 = np.array([-0.91858, -49.9912, 64.7599, -5.63794, -0.69369, 9.58926,
                   -60.841])
    b3 = np.array([121.721, 69.9834, 27.0389, 74.3435, 6.31734, -65.52920,
                   466.775])

    # Calculation of A and B
    a_vec = a0 + a1*fp.omega
    if extended_calc:
        a_vec += a2*mu_r**4 + a3*fp.association_parameter

    b_vec = b0 + b1*fp.omega
    if extended_calc:
        b_vec += b2*mu_r**4 + b3*fp.association_parameter

    # Collision Integral(ci, original paper Omega*)
    aa = 1.16145
    bb = 0.14874
    cc = 0.52487
    dd = 0.77320
    ee = 2.16178
    ff = 2.43787
    gg = -6.435e-4
    hh = 7.27371
    ss = 18.0323
    ww = -0.76830

    temp_star = 1.2593*temp/fp.temp_c
    ci = ((aa/temp_star**bb) +
          cc/np.exp(dd*temp_star) +
          ee/np.exp(ff*temp_star) +
          gg*temp_star**bb*np.sin(ss*temp_star**ww - hh))
    fc = 1 - 0.2756*fp.omega + 0.059035*mu_r**4 + fp.association_parameter

    # Calculation visc_ref
    visc_ref = 4.0785e-5*(fp.mass*temp)**0.5/(
            v_c_cm3_p_mol**(2/3)*ci)*fc

    # Calculation visc
    y = rho_mol_p_cm3 * v_c_cm3_p_mol / 6
    g1 = (1 - 0.5*y)/(1 - y)**3
    g2 = ((a_vec[0]*(1 - np.exp(-a_vec[3]*y))/y
           + a_vec[1]*g1*np.exp(a_vec[4]*y) + a_vec[2]*g1)
          / (a_vec[0]*a_vec[3] + a_vec[1] + a_vec[2]))
    visc_k = visc_ref*(1/g2 + a_vec[5]*y)
    visc_p = ((36.344e-6*(fp.mass*fp.temp_c)**0.5/v_c_cm3_p_mol**(2/3))
              * a_vec[6]*y**2*g2
              * np.exp(a_vec[7] + a_vec[8]/temp_star + a_vec[9]/temp_star**2))
    visc = visc_k + visc_p  # P

    # Calculation cond_ref
    alpha = (cv_cal_p_mol_p_kelvin/R_MOL) - (3/2)
    beta = 0.7862 - 0.7109*fp.omega + 1.3168*fp.omega**2
    temp_r = temp/fp.temp_c
    zeta = 2 + 10.5*temp_r**2
    psi = (1 + alpha*((0.215 + 0.28288*alpha - 1.061*beta + 0.26665*zeta)
                      / (0.6366 + beta*zeta + 1.061*alpha*beta)))
    cond_ref = 7.452*(visc_ref/fp.mass)*psi

    # Calculation cond
    h2 = ((b_vec[0]*(1 - np.exp(-b_vec[3]*y))/y
           + b_vec[1]*g1*np.exp(b_vec[4]*y) + b_vec[2]*g1)
          / (b_vec[0]*b_vec[3] + b_vec[1] + b_vec[2]))
    cond_k = cond_ref*(1/h2 + b_vec[5]*y)
    cond_p = ((3.039e-4*(fp.temp_c/fp.mass)**0.5/v_c_cm3_p_mol**(2/3))
              * b_vec[6]*y**2*h2*temp_r**0.5)
    cond = cond_k + cond_p  # cal/cm s K

    # Conversion to correct unit for main program
    visc_pascal_s = visc/10
    cond_watt_p_meter_p_kelvin = cond*J_PER_CAL*100

    return visc_pascal_s, cond_watt_p_meter_p_kelvin
