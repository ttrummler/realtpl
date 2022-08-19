import os
from dataclasses import dataclass
from CoolProp.CoolProp import PropsSI

from realtpl.thermophysical_constants import R_UNIV

from realtpl.data_association_parameter import data_association_parameter
from realtpl.data_dipole_moment import data_dipole_moment
from realtpl import nasa


@dataclass
class FluidProperties:
    name: str
    mass: float
    omega: float
    p_c: float
    temp_c: float
    rho_c: float
    data_nasa: nasa.NasaCoefficients
    association_parameter: float = 0.0
    dipole_moment: float = 0.0

    def __post_init__(self):
        self.v_c = 1/self.rho_c  # m^3/kmol
        self.Z_c = self.p_c*self.v_c/(R_UNIV*self.temp_c)  # -

    def __str__(self):
        return (self.name + '\n'
                + 'mass: ' + str(self.mass) + ' kg/kmol\n'
                + 'acentric factor: ' + str(self.omega) + '\n'
                + 'critical pressure: ' + str(self.p_c) + ' Pa\n'
                + 'critical temperature: ' + str(self.temp_c) + ' K\n'
                + 'critical density: ' + str(self.rho_c) + ' kmol/m3\n'
                + 'critical volume: ' + str(self.v_c) + ' m3/kmol\n'
                + 'critical compressibility: ' + str(self.Z_c) + '  ')


def fluid_properties_from_coolprop_and_data_base(name: str,
                                                 data_nasa: nasa.NasaCoefficients
                                                 ) -> FluidProperties:
    return FluidProperties(
        name,
        mass=PropsSI('M', name)*1e3,  # kg/kmol
        omega=PropsSI('acentric', name),  # -
        p_c=PropsSI('pcrit', name),  # Pa
        temp_c=PropsSI('Tcrit', name),  # K
        rho_c=PropsSI('rhomolar_critical', name)*1e-3,  # kmol/m^3
        data_nasa=data_nasa,
        association_parameter=data_association_parameter[name],
        dipole_moment=data_dipole_moment[name]
    )


def save_fp_to_file(fp: FluidProperties, output_dir: str):
    path = os.path.join(os.getcwd(), output_dir, fp.name)
    os.makedirs(path, exist_ok=True)
    print('Starting evaluation for ' + fp.name + '\n' +
          'output is written to ' + path)
    with open(os.path.join(path, fp.name + '.out'), 'w') as f:
        f.writelines(str(fp))
