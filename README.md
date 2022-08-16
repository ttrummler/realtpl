# What is `realtpl`?
`realtpl` stands for real gas thermodynamic python library. It computes the density (rho), heat capacity (cp), speed of sound (sound), viscosity (visc) and heat conductivity (lambda) using a thermodynamic model based on cubic equations of state (EoS) and compares it with reference data from the open source library CoolProp. The cubic EoS are the well-known Soave-Redlich-Kwong (SRK) and Peng-Robinson (PR) EoS and the more sophisticated three-parameter Redlich-Kwong-Peng-Robinson (RKPR) EoS, which uses the critical compressibility factor for the EoS parameters, for details see reference _Trummler et al._. The current implementation is designed to evaluate results over a temperature range (with defined number of temperature steps) for a given pressure level. The data is displayed graphically and can also be exported to a _csv_ files for further processing. The only mandatory input is the fluid name, the temperature range and the pressure. Fluid properties, such as mass, critical pressure etc., are directly extracted from CoolProp. This data is also saved to the output directory (specified in the configuration file) as _Fluidname_.out. For the calculation of the heat capacity a reference value has to be evaluated, which is done using the NASA polynomials. Therefore the coefficients have to be provided and for selected fluids we have already generated a database in the corresponding nasa_X_ yaml files. More data can be found in the references provided in the header of these files.

A special feature of `realtpl` is that all calculations with the suggested thermodynamic model are vectorized, making them very fast. 

# Motivation and targeted audience
`realtpl` is open-source software (under GNU GPL License) for physicists, engineers, scientists, technicians and anyone interested in real gas thermodynamics. It runs on all operating systems which support Python, is quick to install, free of charge and designed to be easy to use. The indented use is to run the executable and modify the configuration files only. 

`realtpl` has been specifically designed to be used for CFD simulations. Checking the accuracy of a thermodynamic model in advance is a central step before conducting CFD simulations. For different fluids as well as different pressure and temperature ranges, such an evaluation can be complicated and especially time consuming. To this end, we have written `realtpl` to easily compare the results obtained with a thermodynamic model based on cubic EoS, for details see reference _Trummler et al._.

With some small extensions it could also be used as a table generator. Furthermore, `realtpl` can serve as an inspiration for implementing the present model into an internal flow solver. 

# Installation

`realtpl` is available as python package. Using pip you can simply install it with:
````bash
 pip install realtpl
 ````

 Additionally, `realtpl` is also available on github https://github.com/ttrummler/realtpl. Clone the latest version of `realtpl` with
````bash
git clone https://github.com/ttrummler/realtpl
````

You will now have an executable named `real-tpl`, which you can run anywhere.

# Running `realtpl`
To run the test example use 
````bash
real-tpl --config /path/to/tests/example/config.yaml
```` 
In the configuration file everything for the computation is specified. 

# Configuration file

See `tests/example/` for example configuration files.
A minimum working example for the configuration file reads:
````yaml

````
Additional optional configuration parameters are:


The used configurations are written out to the output directory.

# Latest source code

The latest development version of `realtpl` can be obtained at

    https://github.com/ttrummler/realtpl

# Bug reports

To report bugs, please use `realtpl`’s Bug Tracker at:

    https://github.com/ttrummler/realtpl

# License information

See LICENSE for information on the terms & conditions for usage of this software, and a DISCLAIMER OF ALL WARRANTIES.

Cite `realtpl` if used in your work or to generate data required for your work. Cite as: _Trummler, T., Glatzle, M., Doehring, A., Urban, N., Klein, M., Thermodynamic modeling for numerical simulations based on the generalized cubic equation of state._

**NOTE** the license for the included data base of the 7 coefficient NASA polynomials used for the calculation of the caloric properties. Using this data requires proper citation to be included in the pertinent publications. Cite: _Goos, E., Burcat, A., Ruscic, B.. New NASA Thermodynamic Polynomials Database With Active Thermochemical Tables updates. Report ANL 05/20 TAE 960._ 

# Citation

To cite `realtpl` in publications use:

_Trummler, T., Glatzle, M., Doehring, A., Urban, N., Klein, M., Thermodynamic modeling for numerical simulations based on the generalized cubic equation of state._