# What is _realtpl_?
_realtpl_ stands for real gas thermodynamic python library and computes the density (rho), heat capacity (cp), speed of sound (sound), viscosity (visc) and heat conductivity (lambda) using a thermodynamic model based on cubic equations of state (EoS) and compares it with reference data from the open source library CoolProp. The cubic EoS are the well-known Soave-Redlich-Kwong (SRK) and Peng-Robinson (PR) EoS and the more sophisticated three-parameter Redlich-Kwong-Peng-Robinson (RKPR) EoS, which uses the critical compressibility factor for the EoS parameters, for details see reference _Trummler et al._. The current implementation is designed to evaluate results over a temperature range (with defined number of temperature steps) for a given pressure level. The data is displayed graphically and can also be exported to a _csv_ files for further processing. The only mandatory input is the fluid name, the temperature range and the pressure. Additionally, also evaluation over a desired temperature and pressure range is possible. This feature aims to support for table generations. 

Fluid properties, such as mass, critical pressure etc., are directly extracted from CoolProp. This data is also saved to the output directory (specified in the configuration file) as _Fluidname_.out. For the calculation of the heat capacity a reference value has to be evaluated, which is done using the NASA 7 or 9 coefficient polynomials. Therefore the coefficients have to be provided and for selected fluids we have already generated a database in the corresponding file - nasa_7.yaml or nasa_9.yam. More data can be found in the references provided in the header of these files.

A special feature of _realtpl_ is that all calculations with the suggested thermodynamic model are vectorized, making them very fast.

# Motivation and targeted audience
_realtpl_ is open-source software (under GNU GPL License) for physicists, engineers, scientists, technicians and anyone interested in real gas thermodynamics. It runs on all operating systems which support Python, is quick to install, free of charge and designed to be easy to use. The indented use is to run the executable and modify the configuration files only.

_realtpl_ has been specifically designed to be used for CFD simulations. Checking the accuracy of a thermodynamic model in advance is a central step before conducting CFD simulations. For different fluids as well as different pressure and temperature ranges, such an evaluation can be complicated and especially time consuming. To this end, we have written _realtpl_ to easily compare the results obtained with a thermodynamic model based on cubic EoS, for details see reference _Trummler et al._.

# Installation

_realtpl_ is available as python package. Using pip you can simply install it with:
````bash
 pip install realtpl
 ````

 Additionally, _realtpl_ is also available on github https://github.com/ttrummler/realtpl. Clone the latest version of _realtpl_ with
````bash
git clone https://github.com/ttrummler/realtpl
````

You will now have an executable named _realtpl_, which you can run anywhere.

# Testing
Testing is configured using `tox`.
If not already installed, install `tox` using pip
````bash
pip install tox
````
Then, you can simply enter the `tox` and the tests will automatically be executed
````bash
tox
````
**Note** that it might not run for all tested python versions (py38, py39, py310) on your OS. 

# Running _realtpl_
To run the test example use
````bash
realtpl --config /path/to/tests/example/config.yaml
```` 
In the configuration file everything for the computation is specified. 

# Configuration file

See `tests/example/` for example configuration files.
A minimum working example for the configuration file reads:
````yaml
fluid_name: nHexane
pressure_Pa: 6.0e+06
temperature_start_K: 300
temperature_end_K: 600
````
Note that pressure and temperature input has to be in the SI units Kelvin and Pascal, as also indicated in the parameter names. 

Additional optional configuration parameters are:
````yaml
fluid_name: nHexane
eos_list: [SRK, PR, RKPR] # optional; default: SRK, PR, RKPR;
include_ref_data: true # optional; default: true
pressure_Pa: 6.0e+06
temperature_start_K: 300
temperature_end_K: 600
temperature_step_K: 1 # optional; default: 1 
n_nasa_coeff: 7 # optional; default: 7
output_dir: results # optional; default: results
save_data_to_csv: true # optional;
show_plots: false # optional; default: true
save_plots: true # optional; default: false
show_deviation: true # optional; default: false
save_deviation: true # optional; default: false
performance_tracking: true # optional; default: false
````
Apart from that also evaluations for a pressure and temperature range are possible. However, then no graphical output can be activated. 
````yaml
fluid_name: nHexane
eos_list: [SRK, PR, RKPR]
include_ref_data: true
pressure_start_Pa: 4.0e+06
pressure_end_Pa: 8.0e+06
pressure_step_Pa: 1.0e+06
temperature_start_K: 300
temperature_end_K: 600
temperature_step_K: 100
n_nasa_coeff: 7
output_dir: results
save_data_to_csv: true
show_plots: false
save_plots: false
````
The used configurations are written out to the output directory.

# Latest source code

The latest development version of _realtpl_ can be obtained at

    https://github.com/ttrummler/realtpl

# Bug reports

To report bugs, please use _realtpl_’s Bug Tracker at:

    https://github.com/ttrummler/realtpl

# License information

See LICENSE for information on the terms & conditions for usage of this software, and a DISCLAIMER OF ALL WARRANTIES.

Cite _realtpl_ if used in your work or to generate data required for your work. Cite as: _Trummler, T., Glatzle, M., Doehring, A., Urban, N., Klein, M., Thermodynamic modeling for numerical simulations based on the generalized cubic equation of state._

**NOTE** the license for the included data base of the 7 coefficient NASA polynomials used for the calculation of the caloric properties. Using this data requires proper citation to be included in the pertinent publications. Cite: _Goos, E., Burcat, A., Ruscic, B.. New NASA Thermodynamic Polynomials Database With Active Thermochemical Tables updates. Report ANL 05/20 TAE 960._

# Citation

To cite _realtpl_ in publications use:

_Trummler, T., Glatzle, M., Doehring, A., Urban, N., Klein, M., Thermodynamic modeling for numerical simulations based on the generalized cubic equation of state._