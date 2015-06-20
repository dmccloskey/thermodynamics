Thermodynamics
============
A thermodynamics module for cobrapy
============
Douglas McCloskey
-----------------

thermodynamics module for quantitative metabolomics data

To install, please follow the [instructions](INSTALL.md).

Methods
-----------------
Methods are based off of the following work:

Henry CS, Broadbelt LJ, Hatzimanikatis V (2007) Thermodynamics-based metabolic flux analysis. Biophysical journal 92: 1792–805.

Noor E, Bar-Even A, Flamholz A, Lubling Y, Davidi D, et al. (2012) An integrated open framework for thermodynamics of reactions that combines accuracy and coverage. Bioinformatics (Oxford, England) 28: 2037–2044.

Noor E, Haraldsdo´ ttir HS, Milo R, Fleming RMT (2013) Consistent Estimation of Gibbs Energy Using Component Contributions. PLoS Comput Biol 9(7): e1003098. doi:10.1371/journal.pcbi.1003098

Notes:
-----------------
A modified version of the component-contribution method https://github.com/dmccloskey/component-contribution/ by Elad Noor (https://github.com/eladnoor/component-contribution)

the modified version does not require openbabel (currently, openbabel only supports 32bit python)

the modified version does not require oct2py

several of the function calls to calculate dG0_r and dG_r have been modified to return dG0_f and dG_f instead

all compounds covered by the reactant contribution method are provided as a flat file for faster calculations
