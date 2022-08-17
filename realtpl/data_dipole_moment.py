# Database to calculate viscosity and heat conductivity for
# dense fluids based on paper by Chung et al. 1988.

# Macro - https://macro.lsu.edu/HowTo/solvents/Dipole%20Moment.htm

# Dipole moment [D]

from realtpl.ZeroDict import ZeroDict

data_dipole_moment = ZeroDict({
    'Methanol': 2.87,
    'Ethanol':  1.66,
    'Water':    1.87
})
