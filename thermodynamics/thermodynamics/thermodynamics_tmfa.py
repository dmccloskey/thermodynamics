from __future__ import with_statement
from math import floor,ceil,log,sqrt,pow,exp,fabs
from copy import deepcopy
from cobra.core.Metabolite import Metabolite
from cobra.core.Reaction import Reaction
from collections import Counter

# Other dependencies
import csv,json,sys

class thermodynamics_tfba():    
    """1. Runs thermodynamic flux balance analysis analysis on a cobra.Model object
    2. Runs thermodynamic flux variabiity balance analysis analysis on a cobra.Model object
    3. Runs thermodynamic dG_r variability analysis analysis on a cobra.Model object
    4. Runs thermodynamic metabolite variability analysis analysis on a cobra.Model object

    #1 Calculate the optimal solution with additional constraints from thermodynamics

    measured_concentrations: measured concentration values
                             metabolite.id: {'concentration': float,
                                             'concentration_lb': float,
                                             'concentration_ub': float,
                                             'concentration_var': float,
                                             'concentration_units': string}
                     NOTE: all measured concentrations and variances are computed in ln space (i.e. geometric mean and variance)

    estimated_concentrations: estimated concentrations for those metabolites that were not measured
                     metabolite.id: {'concentration': float,
                                     'concentration_lb': float,
                                     'concentration_ub': float,
                                     'concentration_var': float,
                                     'concentration_units': string}

    calculated reaction free energies dG_r = {reaction.id: {'dG_r': float,
                          'dG_r_var': float,
                          'dG_r_lb': float,
                          'dG_r_ub': float,
                          'Keq': float,
                          'ri': float,
                          'ri_estimate': float,
                          'Q': float,
                          'Q_estimate': float,
                          'dG_r_units': string}}

    #2 Calculate the flux ranges with additional constraints from thermodynamics

        flux_bounds: {reaction.id: {'maximum': float, 'minimum': float, 'flux_units': 'mmol*gDW-1*hr-1'}}
        
    #3 Calculate the dG_r ranges with additional constraints from thermodynamics

        dG_r_bounds: {reaction.id: {'dG_r_lb': float, 'dG_r_ub': float, 'flux_units': string}}

    #4 Calculate the metabolite concentration ranges with additional constraints from thermodynamics

        concentration_bounds: {metabolite.id: {'concentration_lb': float, 'concentration_ub': float, 'concentration_units': string}}
    """

    def __init__(self):
        return;
    def tfba(self, cobra_model, measured_concentration, estimated_concentration, dG_r):
        '''performs thermodynamic flux balance analysis'''
        return;
    def tfva(self, cobra_model, measured_concentration, estimated_concentration, dG_r):
        '''performs thermodynamic flux variability analysis'''
        return;
    def tfva_dG_r(self, cobra_model, measured_concentration, estimated_concentration, dG_r):
        '''performs thermodynamic dG_r variability analysis'''
        return;
    def tfva_concentrations(self, cobra_model, measured_concentration, estimated_concentration, dG_r):
        '''performs thermodynamic metabolite concentration variability analysis'''
        return;