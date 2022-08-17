# Database to calculate viscosity and heat conductivity for
# dense fluids based on paper by Chung et al. 1988.

# REF: Chung, T.-H., Ajlan, M., Lee, L. L., & Starling, K. E. (1988).
# Generalized Multiparameter Correlation for Nonpolar and Polar Fluid
# Transport Properties. Ind. Eng. Chem. Res., 27, 671?679.
# http://doi.org/10.1021/ie00076a024

# Association Parameters Kappa

from realtpl.ZeroDict import ZeroDict

data_association_parameter = ZeroDict({
    'Methanol': 0.215175,
    'Ethanol':  0.174823,
    'Propanol': 0.143453,
    'Butanol':  0.131671,
    'Pentanol': 0.121555,
    'Hexanol':  0.114230,
    'Heptanol': 0.108674,
    'Acid':     0.091549,
    'Water':    0.075908
})
