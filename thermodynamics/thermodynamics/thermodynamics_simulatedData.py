# Dependencies from cobra
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis import single_deletion
from cobra.flux_analysis.objective import update_objective

import json, csv
from math import sqrt,exp,pow
from numpy import average, var, log

from thermodynamics_io import thermodynamics_io

class thermodynamics_simulatedData(thermodynamics_io):
    """Class to hand input of simulated Data"""

    def __init__(self):
        self.fva_data = {}
        self.sra_data = {}

    def check_data(self):
        '''check data integrity'''
        return

    def generate_sra_data(self, cobra_model, element_list=None,
                            method='fba', the_problem='return',
                            element_type='reaction', solver='gurobi',
                            error_reporting=None):

        print 'Single Reaction Deletion...'

        # single reaction deletion
        single_reaction_deletions = single_deletion(cobra_model, element_list=None,
                            method='fba', the_problem='return',
                            element_type='reaction', solver='gurobi',
                            error_reporting=None);

        # FBA
        cobra_model.optimize(solver='gurobi');

        for k,v in single_reaction_deletions[0].iteritems():
            self.sra_data[k.id] = {'gr':None,'gr_ratio':None};
            if v:
                self.sra_data[k.id] = {'gr':v,'gr_ratio':v/cobra_model.solution.f};

    def export_sra_data(self, filename):
        '''export sra data'''
        self.export_values_json(filename, self.sra_data);

    def generate_fva_data(self, cobra_model, fraction_of_optimum=0.9,
                                      objective_sense='maximize', the_reactions=None,
                                      allow_loops=True, solver='gurobi',
                                      the_problem='return', tolerance_optimality=1e-6,
                                      tolerance_feasibility=1e-6, tolerance_barrier=1e-8,
                                      lp_method=1, lp_parallel=0, new_objective=None,
                                      relax_b=None, error_reporting=None,
                                      number_of_processes=1, copy_model=True):

        print 'FVA...'
        # calculate the reaction bounds using FVA
        self.fva_data = flux_variability_analysis(cobra_model, fraction_of_optimum=0.9,
                                      objective_sense='maximize', the_reactions=None,
                                      allow_loops=True, solver='gurobi',
                                      the_problem='return', tolerance_optimality=1e-6,
                                      tolerance_feasibility=1e-6, tolerance_barrier=1e-8,
                                      lp_method=1, lp_parallel=0, new_objective=None,
                                      relax_b=None, error_reporting=None,
                                      number_of_processes=1, copy_model=True);

    def export_fva_data(self, filename):
        '''export fva data'''
        self.export_values_json(filename, self.fva_data);

    def import_sra_data(self,filename):
        self.sra_data = json.load(open(filename));

    def import_fva_data(self,filename):
        fva = json.load(open(filename));
        self.fva_data = self._convert_fluxBounds2var(fva);
        
    def _convert_fluxBounds2var(self, flux_bounds):
        """
        convert flux bounds from FVA to median and variance
    
        variance = (max - median)^2

        flux_bounds: {reaction.id: {'maximum': float, 'minimum': float}}

        returns a dictionary: {reaction.id: {'flux': float, 'flux_var': float, 'flux_units': 'mmol*gDW-1*hr-1'}

        """

        flux_bounds_O = {};
        for k,v in flux_bounds.iteritems():
            median = (v['maximum'] - v['minimum'])/2;
            variance = (v['maximum'] - median)*(v['maximum'] - median);
            flux_bounds_O[k] = {'flux': median, 'flux_var': variance, 'flux_units': 'mmol*gDW-1*hr-1',
                                'flux_lb': v['minimum'], 'flux_ub': v['maximum']};

        return flux_bounds_O