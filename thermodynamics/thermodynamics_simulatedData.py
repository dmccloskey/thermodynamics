# -*- coding: utf-8 -*-
# Dependencies from cobra
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis.single_deletion import single_reaction_deletion,single_gene_deletion
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
from cobra.flux_analysis.loopless import construct_loopless_model

import json, csv
from math import sqrt,exp,pow
from numpy import average, var, log

from .thermodynamics_io import thermodynamics_io

class thermodynamics_simulatedData(thermodynamics_io):
    """Class to generate and handle COBRA simulated data"""

    def __init__(self,fva_data_I={},
                 sra_data_I={},
                 sga_data_I={},
                 fba_primal_data_I={},
                 fba_dual_data_I={}):
        if fva_data_I:
            self.fva_data = self._convert_fluxBounds2var(fva_data_I)
        else:
            self.fva_data = {}
        if sra_data_I:
            self.sra_data = sra_data_I
        else:
            self.sra_data = {}
        if fba_primal_data_I:
            self.fba_primal_data = fba_primal_data_I
        else:
            self.fba_primal_data = {}
        if fba_dual_data_I:
            self.fba_dual_data = fba_dual_data_I
        else:
            self.fba_dual_data = {}
        if sga_data_I:
            self.sga_data = sra_data_I
        else:
            self.sga_data = {}

    def check_data(self):
        """check data integrity"""
        return

    def generate_sra_data(self, cobra_model, reaction_list=None,
                            method_I='fba', solver='glpk', verbose_I=True):
        """Single reaction deletion analysis

        Args:
            cobra_model (cobra.Model): cobra model object
            reaction_list (list(cobra.Reaction)): list of cobra model reactions to use with SRA
            method_I (string): 'fba', 'moma'
            solver (string): default = 'glpk'
            verbose_I (boolean): print messages to the console
        """
        if verbose_I:
            print('Single Reaction Deletion...')

        # single reaction deletion
        cobra_model.solver = solver
        single_reaction_deletions = single_reaction_deletion(cobra_model,
                        #reaction_list=reaction_list,
                        method=method_I
                        )

        # FBA
        cobra_model.optimize()

        # for k,v in single_reaction_deletions[0].items():
        #     self.sra_data[k] = {'gr':None,'gr_ratio':None,'method':method_I}
        #     if v:
        #         self.sra_data[k] = {'gr':v,'gr_ratio':v/cobra_model.objective.value}
        for i in range(len(single_reaction_deletions.index)):
            self.sra_data[single_reaction_deletions.index[i]] = {'gr':None,'gr_ratio':None,'method':method_I}
            if single_reaction_deletions.flux[i]:
                self.sra_data[single_reaction_deletions.index[i]] = {
                    'gr':single_reaction_deletions.flux[i],
                    'gr_ratio':single_reaction_deletions.flux[i]/cobra_model.objective.value}


    def export_sra_data(self, filename):
        """export sra data"""
        with open(filename, 'w') as outfile:
            json.dump(self.sra_data, outfile, indent=4);

    def import_sra_data(self,filename):
        """import sra data"""
        self.sra_data = json.load(open(filename))

    def generate_fva_data(self, cobra_model, fraction_of_optimum=0.9,
            objective_sense='maximize', reaction_list=None,
            allow_loops=True, solver='glpk', verbose_I=True):

        if verbose_I:
            print('FVA...')
        #add in loop law constrain
        if not allow_loops: cobra_model=construct_loopless_model(cobra_model)
        # calculate the reaction bounds using FVA
        cobra_model.solver = solver
        fva_data = flux_variability_analysis(cobra_model, fraction_of_optimum=0.9,
                                      objective_sense='maximize',
                                      reaction_list=reaction_list,
                                      )        
        self.fva_data = dict(zip(list(fva_data.index),fva_data.to_dict('records')))

    def export_fva_data(self, filename):
        """export fva data"""
        with open(filename, 'w') as outfile:
            json.dump(self.fva_data, outfile, indent=4);


    def import_fva_data(self,filename):
        """import fva data"""
        fva = json.load(open(filename))
        self.fva_data = self._convert_fluxBounds2var(fva)
        
    def _convert_fluxBounds2var(self, flux_bounds):
        """convert flux bounds from FVA to median and variance
    
        variance = (max - median)^2

        Args:
            flux_bounds (dict): {reaction.id: {'maximum': float, 'minimum': float}}

        Returns:             
            dict: output:  {reaction.id: {'flux': float, 'flux_var': float, 'flux_units': 'mmol*gDW-1*hr-1'}

        """

        flux_bounds_O = {}
        for k,v in flux_bounds.items():
            median = (v['maximum'] - v['minimum'])/2
            variance = (v['maximum'] - median)*(v['maximum'] - median)
            flux_bounds_O[k] = {'flux': median, 'flux_var': variance, 'flux_units': 'mmol*gDW-1*hr-1',
                                'flux_lb': v['minimum'], 'flux_ub': v['maximum']}

        return flux_bounds_O

    def generate_fba_data(self,cobra_model,allow_loops=True, method_I='fba',solver='glpk'):
        """perform FBA simulation on the model
        """
        #add in loop law constrain
        if not allow_loops: cobra_model=construct_loopless_model(cobra_model)
        #check for the optimization method:
        cobra_model.solver = solver
        if method_I=='fba' or method_I=='loopless-fba':
            sol = cobra_model.optimize()
        elif method_I =='pfba' or method_I=='loopless-pfba':
            sol = optimize_minimal_flux(model=cobra_model,solver=solver)
        else:
            print('method not recognized.')
            return
        
        self.fba_primal_data={}
        for k,v in sol.x_dict.items():
            self.fba_primal_data[k] = v
        self.fba_dual_data={}
        for k,v in sol.y_dict.items():
            self.fba_dual_data[k] = v

    def generate_sga_data(self, cobra_model, gene_list=None,
                            method_I='fba', solver='glpk'):
        """Single gene deletion analysis
        
        Args:
            gene_list (list): list of genes
            method_I (str): 'fba', 'moma'
        """

        print('Single Reaction Deletion...')

        # single gene deletion
        single_gene_deletions = single_gene_deletion(cobra_model,
                        #gene_list=gene_list,
                        method=method_I,
                        solver=solver
                        )

        # FBA
        cobra_model.solver = solver
        cobra_model.optimize()

        for k,v in single_gene_deletions[0].items():
            self.sga_data[k] = {'gr':None,'gr_ratio':None,'method':method_I}
            if v:
                self.sga_data[k] = {'gr':v,'gr_ratio':v/cobra_model.objective.value}

    def export_sga_data(self, filename):
        """export sga data"""
        with open(filename, 'w') as outfile:
            json.dump(self.sga_data, outfile, indent=4);

    def import_sga_data(self,filename):
        """import sga data"""
        self.sga_data = json.load(open(filename))

    #TODO:
    def reduce_model(self,cobra_model,method_I='fba',solver='glpk'):
        """reduce the model"""
        pass
    def generate_fluxSum_data(self,cobra_model,method_I='fba',solver='glpk'):
        """
        perform a fluxSum analysis        
        """
        pass