import numpy as np

from realtpl.eos_data import EosParameter
from realtpl.thermophysical_constants import R_UNIV


def calc_compressibility(ed: EosParameter, alpha: callable,
                         temp: np.ndarray, press: float):
    """
        compressibility - solves cubic equation for compressibility factor

        Computes compressibility factor Z (in code z) from cubic EOS
                    Z^3 + c_2*Z^2 + c_1*Z + c_0 = 0
        Parameters c_2, c_1, c_0 depend on the eos parameter d_1 and d_2
        and A and B (in code aa and bb).

        This implementation is based on the work by Matheis [1,2]. We
        reformulated the parameters of the cubic eos to hold for the
        generalized formulation and added some tweaks and tricks to conduct
        the calculations in vectorized form. For the general approach to solve
        a cubic eos, see also [3].

        First Q and R (qq and rr) are calculated and then based on the sign of
        the discriminant D = R² - Q³ either one or three real roots exist.

        If only one real root exists, this the solution, which is calculated as
        described in [1] and in [3] at the bottom. Note that in OpenFOAM the
        implementation is slightly different (as in [3] in the top), however,
        yielding the same solution.

        If three real roots may exist, these first have to calculated. Then the
        maximum root is the solution for the vapor/gas phase and the minimum for
        the liquid phase. (The center root is thermodynamically meaningless.)
        Afterwards it has to be checked if the liquid root (z_l) makes sense from
        a thermodynamical point of view, i.e. z_l >= bb, since it can not be
        smaller than the co-volume bb.  If a thermodynamically meaningful
        solutions for z_v and z_l exists, the correct root is determined by
        evaluating, where the Gibbs energy is smaller, see [1,2].

        Parameters
        -----------
        ed: EosParameter
        alpha: callable
            function that returns the temperature dependent alpha values and
            their derivation
        temp: np.ndarray
            temperature in Kelvin
        press: np.ndarray
            pressure in Pascal

        Returns
        ----------
        z: np.ndarray
            compressibility factor

        Authors
        ----------
        Trummler, Glatzle

        References
        ----------
        .. [1] Matheis (2018), PhD. Thesis, Numerical Simulation of Fuel
        Injection and Turbulent Mixing Under High-Pressure Conditions,
        http://mediatum.ub.tum.de/?id=1363601

        .. [2] Matheis and Hickel (2017), Multi-component vapor-liquid
        equilibrium model for LES of high-pressure fuel injection and
        application to ECN Spray A Int. J. Multiph. Flow,
        https://doi.org/10.1016/j.ijmultiphaseflow.2017.11.001

        .. [3] https://www.e-education.psu.edu/png520/m11_p6.html

        .. [5] Elliot & Lira (2012), Introductory Chemical Engineering
        Thermodynamics

        .. [6] Michelsen & Mollerup (2007), Thermodynamic Models: Fundamentals
        \& Computational Aspects, 87-989961-1-8

        .. [7] Trummler, Glatzle, Doehring, Urban, Klein (2022), Thermodynamic
         modeling for numerical simulations based on the generalized cubic
         equation of state.
    """
    aa = (ed.a*alpha(temp)*press)/(R_UNIV*temp)**2
    bb = (ed.b*press)/(R_UNIV*temp)

    c_2 = bb*(ed.d_1_p_d_2 - 1) - 1
    c_1 = aa + bb*(ed.d_1_t_d_2*bb - ed.d_1_p_d_2*(bb + 1))
    c_0 = -bb*(ed.d_1_t_d_2*(bb**2 + bb) + aa)

    qq = (c_2**2 - 3*c_1)/9
    rr = (2*c_2**3 - 9*c_2*c_1 + 27*c_0)/54
    dd = rr**2 - qq**3

    # flags for the number of roots
    one_real_root = (dd >= 0.)
    three_real_roots = (dd < 0)

    # one real root, two imaginary roots exist
    sqrt_dd = np.abs(dd)**0.5
    ee = -np.sign(rr)*(abs(rr) + sqrt_dd)**(1/3)
    ff = qq/ee
    z_one_real = ee + ff - c_2/3

    # check if solution of three real roots is required
    if np.count_nonzero(three_real_roots) == 0:
        z_three_real = np.zeros_like(three_real_roots)

    else:
        sqrt_qq = np.abs(qq)**0.5
        phi = np.arccos((rr/(sqrt_qq*qq))*three_real_roots)

        x1 = (-2*sqrt_qq*np.cos(phi/3) - c_2/3)*three_real_roots
        x2 = (-2*sqrt_qq*np.cos(
            (phi + 2*np.pi)/3) - c_2/3)*three_real_roots
        x3 = (-2*sqrt_qq*np.cos(
            (phi - 2*np.pi)/3) - c_2/3)*three_real_roots

        z_solutions = np.matrix([x1, x2, x3])

        # min root: liquid, max: vapor, center: thermodynamically meaningless
        z_l = z_solutions.min(0)
        z_v = z_solutions.max(0)
        z_l = np.squeeze(np.asarray(z_l))
        z_v = np.squeeze(np.asarray(z_v))

        # volume cannot be smaller than the co-volume
        is_z_v = (z_l < bb)*three_real_roots

        # check if all solutions ar already found
        if sum(is_z_v) == sum(three_real_roots):
            z_three_real = z_v

        # otherwise: determine correct root based on Gibbs
        # see Ref [1] Eq. 2.58
        else:
            eval_gibbs = three_real_roots & ~is_z_v

            # some clipping and mods to avoid error messages
            z_l_minus_bb = np.clip((z_l - bb), 1e-16, np.inf)
            z_v_minus_bb = np.clip((z_v - bb), 1e-16, np.inf)
            dd_l_1 = (z_l + ed.d_1*bb)*eval_gibbs + 1e-16*~eval_gibbs
            dd_l_2 = (z_l + ed.d_2*bb)*eval_gibbs + 1e-16*~eval_gibbs
            dd_v_1 = (z_v + ed.d_1*bb)*eval_gibbs + 1e-16*~eval_gibbs
            dd_v_2 = (z_v + ed.d_2*bb)*eval_gibbs + 1e-16*~eval_gibbs

            dg = (np.log(z_l_minus_bb/z_v_minus_bb)
                  + aa/(bb*(ed.d_1 - ed.d_2))
                  * np.log(dd_l_1/dd_l_2*dd_v_2/dd_v_1)
                  - (z_l - z_v))

            is_z_l = (dg >= 0)*eval_gibbs
            is_z_v = is_z_v + (dg < 0)*eval_gibbs

            z_three_real = is_z_l*z_l + is_z_v*z_v

    return one_real_root*z_one_real + three_real_roots*z_three_real
