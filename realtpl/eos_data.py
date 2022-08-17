from dataclasses import dataclass, field
import numpy as np

from realtpl.fluid_properties import FluidProperties
from realtpl.thermophysical_constants import R_UNIV


@dataclass(frozen=True)
class EosParameter:
    """"
        Dataclass EosParameters - contains Eos Parameter

        Function eos_parameter_from_eos_name  - sets Eos Parameters according
        to the name

        References:
        -----------
        ..[1] Kim, S. K., Choi, H. S., & Kim, Y. (2012). Thermodynamic modeling
         based on a generalized cubic equation of state for kerosene/LOx rocket
        combustion. Combustion and Flame, 159(3), 1351-1365.
    """
    name: str
    d_1: float
    d_2: float = field(init=False)
    d_1_p_d_2: float = field(init=False)
    d_1_t_d_2: float = field(init=False)
    d_1_m_d_2: float = field(init=False)
    a: float
    b: float

    def __post_init__(self):
        object.__setattr__(self, 'd_2', (1 - self.d_1)/(1 + self.d_1))
        object.__setattr__(self, 'd_1_p_d_2', self.d_1 + self.d_2)
        object.__setattr__(self, 'd_1_m_d_2', self.d_1 - self.d_2)
        object.__setattr__(self, 'd_1_t_d_2', self.d_1*self.d_2)


def eos_parameter_from_eos_name(name: str,
                                props: FluidProperties) -> EosParameter:
    if name == 'SRK':
        d_1 = 1
        a_coeff = 0.42747
        b_coeff = 0.08664

    elif name == 'PR':
        d_1 = 1 + np.sqrt(2)
        a_coeff = 0.45724
        b_coeff = 0.07780

    elif name == 'RKPR':
        c_z = 1.168
        d1 = 0.428363
        d2 = 18.496215
        d3 = 0.338426
        d4 = 0.660000
        d5 = 789.723105
        d6 = 2.512392

        d_1 = (d1 + d2*(d3 - c_z*props.Z_c)**d4
               + d5*(d3 - c_z*props.Z_c)**d6)

        d_rkpr = (1 + d_1**2)/(1 + d_1)

        y_rkpr = 1 + (2*(1 + d_1))**(1/3) + (4/(1 + d_1))**(1/3)

        a_coeff = ((3*y_rkpr**2 + 3*y_rkpr*d_rkpr
                    + d_rkpr**2 + d_rkpr - 1)/(3*y_rkpr + d_rkpr - 1)**2)

        b_coeff = 1/(3*y_rkpr + d_rkpr - 1)

    else:
        raise ValueError(f'Unknown EOS: {name}')

    return EosParameter(
        name,
        d_1,
        a_coeff*R_UNIV**2*props.temp_c**2/props.p_c,
        b_coeff*R_UNIV*props.temp_c/props.p_c
    )


class AlphaFunctions:

    def __init__(
            self,
            alpha: callable,
            d_alpha_d_temp: callable,
            d2_alpha_d2_temp: callable):
        self.alpha = alpha
        self.d_alpha_d_temp = d_alpha_d_temp
        self.d2_alpha_d2_temp = d2_alpha_d2_temp


def alpha_functions_from_eos_name(name: str,
                                  props: FluidProperties) -> AlphaFunctions:
    omega = props.omega
    temp_c = props.temp_c
    z_c = props.Z_c

    if name == 'SRK':
        c_alpha = 0.48508 + 1.55171*omega - 0.15613*omega**2

    elif name == 'PR':
        c_alpha = 0.37464 + 1.54226*omega - 0.26992*omega**2

    elif name == 'RKPR':
        c_z = 1.168

        a1 = 66.125
        a0 = -23.359
        b1 = -40.594
        b0 = 16.855
        c1 = 5.27345
        c0 = -0.25826

        c_alpha = ((c_z*z_c*a1 + a0)*omega**2
                   + (c_z*z_c*b1 + b0)*omega
                   + (c_z*z_c*c1 + c0))

    else:
        raise ValueError(f'Unknown EOS: {name}')

    if name == 'SRK' or name == 'PR':
        def alpha(temp: np.ndarray):
            return (1 + c_alpha*(1 - (temp/temp_c)**0.5))**2

        def d_alpha_d_temp(temp: np.ndarray):
            return ((1 + c_alpha*(1 - (temp/temp_c)**0.5))
                    * (-c_alpha/(temp_c*temp)**0.5))

        def d2_alpha_d2_temp(temp: np.ndarray):
            return ((1 + c_alpha*(1 - (temp/temp_c)**0.5))
                    * c_alpha/(2*(temp**3*temp_c)**0.5)
                    + c_alpha**2/(2*temp_c*temp))

    elif name == 'RKPR':
        def alpha(temp: np.ndarray):
            return (3/(2 + temp/temp_c))**c_alpha

        def d_alpha_d_temp(temp: np.ndarray):
            return (-3**c_alpha*c_alpha
                    / (temp_c*(2 + temp/temp_c)**(c_alpha + 1)))

        def d2_alpha_d2_temp(temp: np.ndarray):
            return (3**c_alpha*c_alpha*(c_alpha + 1)
                    / (temp_c**2*(2 + temp/temp_c)**(c_alpha + 2)))

    else:
        raise ValueError(f'Unknown EOS: {name}')

    return AlphaFunctions(
        alpha,
        d_alpha_d_temp,
        d2_alpha_d2_temp
    )
