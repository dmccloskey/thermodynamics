from __future__ import with_statement
from math import floor,ceil,log,sqrt,pow,exp,fabs
from copy import deepcopy
from cobra.core.Metabolite import Metabolite
from cobra.core.Reaction import Reaction
from cobra.io import save_matlab_model, write_sbml_model
from collections import Counter
from math import log, exp
from warnings import warn
import numpy
import scipy

from six import iteritems, string_types
from cobra.solvers import solver_dict, get_solver_name
from thermodynamics_utility import find_transportRxns, null
from thermodynamics_sampling import thermodynamics_sampling

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
        # set min/max for dG_r
        self.dG_r_min = -1e4;
        self.dG_r_max = 1e4;
        # set min/max for dG0_r
        self.dG0_r_min = -1e4;
        self.dG0_r_max = 1e4;
        # set min/max for metabolite activity
        y = 1000.0
        conc = 5.0e-5
        self.conc_min = 1/sqrt(y)*conc;
        self.conc_max = sqrt(y)*conc;
        self.conc_min_ln = log(self.conc_min);
        self.conc_max_ln = log(self.conc_max);
        # initialize constants
        self.R = 8.314e-3; # gas constant 8.3144621 [kJ/K/mol]
        self.K=10*self.dG_r_max; # an arbitrary constant larger than dG_r_max
        # initialize data structures
        self.tfba_data = {};
        self.tfva_data = {};
        self.tfva_dG_r_data = {};
        self.tfva_concentration_data = {};
        self.tfva_analysis = {};
        self.tsampling_dG_r_data = {};
        
    def export_tfva_data(self, filename):
        '''export tfva data'''
        self.export_values_json(filename, self.tfva_data);
    def export_tfva_dG_r_data(self, filename):
        '''export tfva_dG_r data'''
        self.export_values_json(filename, self.tfva_dG_r_data);
    def export_tfva_concentrations_data(self, filename):
        '''export tfva_concentrations data'''
        self.export_values_json(filename, self.tfva_concentrations_data);
    def export_tfva_analysis(self, filename):
        '''export tfva_analysis'''
        self.export_values_json(filename, self.tfva_analysis);
    def _scale_dG_r(self,dG_r):
        '''scale dG_r lb/ub to be within pre-defined bounds'''
        # scale down the magnitude of dG_r
        for k,v in dG_r.iteritems():
            if v['dG_r_lb']<self.dG_r_min:
                dG_r[k]['dG_r_lb']=self.dG_r_min;
            elif v['dG_r_lb']>self.dG_r_max:
                dG_r[k]['dG_r_lb']=1.0;
            if v['dG_r_ub']>self.dG_r_max:
                dG_r[k]['dG_r_ub']=self.dG_r_max;  
            if v['dG_r_ub']<self.dG_r_min:
                dG_r[k]['dG_r_ub']=-1.0;
        return dG_r;
    def _scale_conc(self,concentrations):
        '''scale conc lb/ub to be within pre-defined bounds'''
        # scale down the magnitude of conc
        for k,v in concentrations.iteritems():
            if v['concentration_lb']<self.conc_min:
                concentrations[k]['concentration_lb']=self.conc_min;
            #elif v['concentration_lb']>self.conc_max:
            #    concentrations[k]['concentration_lb']=self.conc;
            if v['concentration_ub']>self.conc_max:
                concentrations[k]['concentration_ub']=self.conc_max;  
            #elif v['concentration_ub']<self.conc_min:
            #    concentrations[k]['concentration_ub']=self.conc;
        return concentrations;
    def _add_dG_r_constraints(self, cobra_model_irreversible, dG_r, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check,
                              measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99, use_measured_dG_r=True, return_dG_r_variables=False):
        '''add constraints for dG_r to the model'''
        # pre-process the data
        dG_r = self._scale_dG_r(dG_r);
        # bounds
        dG_r_indicator_constraint = 1-1e-6;
        # find system boundaries and the objective reaction
        system_boundaries = [x.id for x in cobra_model_irreversible.reactions if x.boundary == 'system_boundary'];
        objectives = [x.id for x in cobra_model_irreversible.reactions if x.objective_coefficient == 1];
        transporters = find_transportRxns(cobra_model_irreversible);
        # add variables and constraints to model for tfba
        reactions = [r for r in cobra_model_irreversible.reactions];
        dG_r_variables = {};
        variables_break = []
        for i,r in enumerate(reactions):
            # ignore system_boundary, objective, and transport reactions
            if r.id in system_boundaries or r.id in objectives or r.id in transporters:
                continue;
            # make a boolean indicator variable
            indicator = Reaction('indicator_' + r.id);
            indicator.lower_bound = 0;
            indicator.upper_bound = 1;
            indicator.variable_kind = 'integer';
            # make a continuous variable for dG_r
            dG_rv = Reaction('dG_rv_' + r.id);
            if use_measured_dG_r and thermodynamic_consistency_check[r.id]: # ignore inconsistent reactions
                #if dG_r_coverage[r.id]>measured_dG_f_coverage_criteria:
                #if metabolomics_coverage[r.id] > measured_concentration_coverage_criteria and \
                #dG_r_coverage[r.id]>measured_dG_f_coverage_criteria:
                    dG_rv.lower_bound = dG_r[r.id]['dG_r_lb'];
                    dG_rv.upper_bound = dG_r[r.id]['dG_r_ub'];
                    dG_rv.objective_coefficient = 0;
                #else: 
                #    dG_rv.lower_bound = self.dG_r_min;
                #    dG_rv.upper_bound = self.dG_r_max;
                #    dG_rv.objective_coefficient = 0;
            else:
                dG_rv.lower_bound = self.dG_r_min;
                dG_rv.upper_bound = self.dG_r_max;
                dG_rv.objective_coefficient = 0;
            dG_rv.variable_kind = 'continuous';
            # create a constraint for vi-zi*vmax<=0
            indicator_plus = Metabolite(r.id + '_plus');
            indicator_plus._constraint_sense = 'L';
            ## create a constraint for vi+zi*vmax>=0
            #indicator_minus = Metabolite(r.id + '_minus');
            #indicator_minus._constraint_sense = 'G';
            # create additional constraint for dG_ri/K+zi<=1-1e-4
            dG_r_constraint = Metabolite(r.id + '_dG_r');
            dG_r_constraint._constraint_sense = 'L';
            dG_r_constraint._bound = self.K - dG_r_indicator_constraint
            # add constraints to the variables
            indicator.add_metabolites({indicator_plus: -r.upper_bound});#,indicator_minus: -r.lower_bound});
            indicator.add_metabolites({dG_r_constraint: self.K});#,indicator_minus: -r.lower_bound});
            dG_rv.add_metabolites({dG_r_constraint: 1.0})
            #indicator.add_metabolites({dG_r_constraint: 1});#,indicator_minus: -r.lower_bound});
            #dG_rv.add_metabolites({dG_r_constraint: 1.0/self.K})
            cobra_model_irreversible.reactions.get_by_id(r.id).add_metabolites({indicator_plus: 1});#,indicator_minus: 1});
            # add indicator reactions to the model
            cobra_model_irreversible.add_reaction(indicator);
            cobra_model_irreversible.add_reaction(dG_rv);
            # check to see if the model broke
            cobra_model_irreversible.optimize(solver='gurobi');
            if not cobra_model_irreversible.solution.f:
                print dG_rv.id + ' broke the model!';
                variables_break.append(dG_rv.id);
                #cobra_model_irreversible.remove_reactions(indicator)
                cobra_model_irreversible.remove_reactions(dG_rv)
                #cobra_model_irreversible.reactions.get_by_id(r.id).subtract_metabolites({indicator_plus:1})
                #indicator.subtract_metabolites({dG_r_constraint: self.K})
                dG_rv.lower_bound = self.dG_r_min;
                dG_rv.upper_bound = self.dG_r_max;
                #indicator.add_metabolites({dG_r_constraint: self.K})
                dG_rv.add_metabolites({dG_r_constraint: 1.0})
                cobra_model_irreversible.add_reaction(dG_rv);
            else:
                # record dG_rv variables
                dG_r_variables[r.id] = dG_rv;
        if return_dG_r_variables:
            return dG_r_variables;
    def _add_conc_ln_constraints_v1(self,cobra_model_irreversible, measured_concentration, estimated_concentration, dG0_r, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check,
                              measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,use_measured_concentrations=True,use_measured_dG0_r=True,return_concentration_variables=False,return_dG0_r_variables=False):
        # pre-process the data
        dG0_r = self._scale_dG_r(dG0_r);
        #measured_concentration = self._scale_conc(measured_concentration);
        #estimated_concentration = self._scale_conc(estimated_concentration);
        # initialize hydrogens:
        hydrogens = [];
        compartments = list(set(cobra_model_irreversible.metabolites.list_attr('compartment')));
        for compart in compartments:
             hydrogens.append('h_' + compart);
        # find system boundaries and the objective reaction
        system_boundaries = [x.id for x in cobra_model_irreversible.reactions if x.boundary == 'system_boundary'];
        objectives = [x.id for x in cobra_model_irreversible.reactions if x.objective_coefficient == 1];
        transporters = find_transportRxns(cobra_model_irreversible);
        # bounds
        dG_r_indicator_constraint = 1-1e-6;
        # add variables and constraints to model for tfba
        reactions = [r for r in cobra_model_irreversible.reactions];
        conc_lnv_dict = {}; # dictionary to record conc_ln variables
        dG_r_variables = {};
        dG0_r_dict = {};
        variables_break = [];
        for i,r in enumerate(reactions):
            if r.id in system_boundaries or r.id in objectives or r.id in transporters:
                continue;
            # create a constraint for vi-zi*vmax<=0
            indicator_plus = Metabolite(r.id + '_plus');
            indicator_plus._constraint_sense = 'L';
            ## create a constraint for vi+zi*vmax>=0
            #indicator_minus = Metabolite(r.id + '_minus');
            #indicator_minus._constraint_sense = 'G';
            # create additional constraint for dG_ri+K*zi<=K-1e-4
            dG_r_constraint = Metabolite(r.id + '_dG_r');
            dG_r_constraint._constraint_sense = 'L';
            dG_r_constraint._bound = self.K - dG_r_indicator_constraint
            # create additional constraint for dG0_ri + RT*SUM[sij*ln(xj)] = dG_ri
            conc_ln_constraint = Metabolite(r.id + '_conc');
            conc_ln_constraint._constraint_sense = 'E';
            conc_ln_constraint._bound = 0
            # make continuous variables for conc
            products = [p for p in r.products];
            # Method 2:
            for p in products:
                if not(p.id in hydrogens): # exclude hydrogen because it has already been accounted for when adjusting for the pH
                    if use_measured_concentrations:# and thermodynamic_consistency_check[r.id]: # ignore inconsistent reactions:
                        if p.id in measured_concentration.keys():
                            # check that if the p_id has already been created:
                            if p.id in conc_lnv_dict.keys():
                                # add constraints
                                conc_lnv_dict[p.id].add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)}); #reaction is also updated in model automatically
                                #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[p.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + p.id);
                                conc_lnv.lower_bound = log(measured_concentration[p.id]['concentration_lb']);
                                conc_lnv.upper_bound = log(measured_concentration[p.id]['concentration_ub']);
                                conc_lnv.variable_kind = 'continuous';
                                conc_lnv.objective_coefficient = 0;
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                                # record the new variable
                                conc_lnv_dict[p.id] = conc_lnv;
                                # add the new variable to the model
                                cobra_model_irreversible.add_reaction(conc_lnv);
                        elif p.id in estimated_concentration.keys():
                            # check that if the p_id has already been created:
                            if p.id in conc_lnv_dict.keys():
                                # add constraints
                                conc_lnv_dict[p.id].add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                                #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[p.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + r.id + '_' + p.id);
                                conc_lnv.lower_bound = log(estimated_concentration[p.id]['concentration_lb']);
                                conc_lnv.upper_bound = log(estimated_concentration[p.id]['concentration_ub']);
                                conc_lnv.variable_kind = 'continuous';
                                conc_lnv.objective_coefficient = 0;
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                                # record the new variable
                                conc_lnv_dict[p.id] = conc_lnv;
                                # add the new variable to the model
                                cobra_model_irreversible.add_reaction(conc_lnv);
                    else:
                        # check that if the p_id has already been created:
                        if p.id in conc_lnv_dict.keys():
                            # add constraints
                            conc_lnv_dict[p.id].add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                            #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[p.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                        else:
                            conc_lnv = Reaction('conc_lnv_' + p.id);
                            conc_lnv.lower_bound = self.conc_min_ln;
                            conc_lnv.upper_bound = self.conc_max_ln;
                            conc_lnv.variable_kind = 'continuous';
                            conc_lnv.objective_coefficient = 0;
                            # add constraints
                            conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                            # record the new variable
                            conc_lnv_dict[p.id] = conc_lnv;
                            # add the new variable to the model
                            cobra_model_irreversible.add_reaction(conc_lnv);
            reactants = [react for react in r.reactants];
            # Method 2:
            for react in reactants:
                if not(react.id in hydrogens): # exclude hydrogen because it has already been accounted for when adjusting for the pH
                    if use_measured_concentrations:# and thermodynamic_consistency_check[r.id]: # ignore inconsistent reactions:
                        if react.id in measured_concentration.keys():
                            # check that if the react_id has already been created:
                            if react.id in conc_lnv_dict.keys():
                                # add constraints
                                conc_lnv_dict[react.id].add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                                #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[react.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + react.id);
                                conc_lnv.lower_bound = log(measured_concentration[react.id]['concentration_ub']);
                                conc_lnv.upper_bound = log(measured_concentration[react.id]['concentration_lb']);
                                conc_lnv.variable_kind = 'continuous';
                                conc_lnv.objective_coefficient = 0;
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                                # record the new variable
                                conc_lnv_dict[react.id] = conc_lnv;
                                # add the new variable to the model
                                cobra_model_irreversible.add_reaction(conc_lnv);
                        elif react.id in estimated_concentration.keys():
                            # check that if the react_id has already been created:
                            if react.id in conc_lnv_dict.keys():
                                # add constraints
                                conc_lnv_dict[react.id].add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                                #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[react.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + r.id + '_' + react.id);
                                conc_lnv.lower_bound = log(estimated_concentration[react.id]['concentration_ub']);
                                conc_lnv.upper_bound = log(estimated_concentration[react.id]['concentration_lb']);
                                conc_lnv.variable_kind = 'continuous';
                                conc_lnv.objective_coefficient = 0;
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                                # record the new variable
                                conc_lnv_dict[react.id] = conc_lnv;
                                # add the new variable to the model
                                cobra_model_irreversible.add_reaction(conc_lnv);
                    else:
                        # check that if the react_id has already been created:
                        if react.id in conc_lnv_dict.keys():
                            # add constraints
                            conc_lnv_dict[react.id].add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                            #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[react.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                        else:
                            conc_lnv = Reaction('conc_lnv_' + react.id);
                            conc_lnv.lower_bound = self.conc_min_ln;
                            conc_lnv.upper_bound = self.conc_max_ln;
                            conc_lnv.variable_kind = 'continuous';
                            conc_lnv.objective_coefficient = 0;
                            # add constraints
                            conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                            # record the new variable
                            conc_lnv_dict[react.id] = conc_lnv;
                            # add the new variable to the model
                            cobra_model_irreversible.add_reaction(conc_lnv);
            ## Method 1:
            #metabolites = [p for p in r.products] + [react for react in r.reactants];
            #for met in metabolites:
            #    if not(met.id in hydrogens): # exclude hydrogen because it has already been accounted for when adjusting for the pH
            #        if use_measured_concentrations:# and thermodynamic_consistency_check[r.id]: # ignore inconsistent reactions:
            #            if met.id in measured_concentration.keys():
            #                # check that if the met_id has already been created:
            #                if met.id in conc_lnv_dict.keys():
            #                    # add constraints
            #                    conc_lnv_dict[met.id].add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)});
            #                else:
            #                    conc_lnv = Reaction('conc_lnv_' + met.id);
            #                    conc_lnv.lower_bound = log(measured_concentration[met.id]['concentration_lb']);
            #                    conc_lnv.upper_bound = log(measured_concentration[met.id]['concentration_ub']);
            #                    conc_lnv.variable_kind = 'continuous';
            #                    conc_lnv.objective_coefficient = 0;
            #                    # add constraints
            #                    conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)});
            #                    # record the new variable
            #                    conc_lnv_dict[met.id] = conc_lnv;
            #            elif met.id in estimated_concentration.keys():
            #                # check that if the met_id has already been created:
            #                if met.id in conc_lnv_dict.keys():
            #                    # add constraints
            #                    conc_lnv_dict[met.id].add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)});
            #                else:
            #                    conc_lnv = Reaction('conc_lnv_' + r.id + '_' + met.id);
            #                    conc_lnv.lower_bound = log(estimated_concentration[met.id]['concentration_lb']);
            #                    conc_lnv.upper_bound = log(estimated_concentration[met.id]['concentration_ub']);
            #                    conc_lnv.variable_kind = 'continuous';
            #                    conc_lnv.objective_coefficient = 0;
            #                    # add constraints
            #                    conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)});
            #                    # record the new variable
            #                    conc_lnv_dict[met.id] = conc_lnv;
            #        else:
            #            # check that if the met_id has already been created:
            #            if met.id in conc_lnv_dict.keys():
            #                # add constraints
            #                conc_lnv_dict[met.id].add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)});
            #            else:
            #                conc_lnv = Reaction('conc_lnv_' + met.id);
            #                conc_lnv.lower_bound = self.conc_min_ln;
            #                conc_lnv.upper_bound = self.conc_max_ln;
            #                conc_lnv.variable_kind = 'continuous';
            #                conc_lnv.objective_coefficient = 0;
            #                # add constraints
            #                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)});
            #                # record the new variable
            #                conc_lnv_dict[met.id] = conc_lnv;
            # make a boolean indicator variable
            indicator = Reaction('indicator_' + r.id);
            indicator.lower_bound = 0;
            indicator.upper_bound = 1;
            indicator.variable_kind = 'integer';
            # make a continuous variable for dG_r
            dG_rv = Reaction('dG_rv_' + r.id);
            #if use_measured_dG_r and thermodynamic_consistency_check[r.id]: # ignore inconsistent reactions
            if False and thermodynamic_consistency_check[r.id]: # ignore inconsistent reactions
                #if dG_r_coverage[r.id]>measured_dG_f_coverage_criteria:
                #if metabolomics_coverage[r.id] > measured_concentration_coverage_criteria and \
                #dG_r_coverage[r.id]>measured_dG_f_coverage_criteria:
                    dG_rv.lower_bound = dG_r[r.id]['dG_r_lb'];
                    dG_rv.upper_bound = dG_r[r.id]['dG_r_ub'];
                    dG_rv.objective_coefficient = 0;
                #else: 
                #    dG_rv.lower_bound = self.dG_r_min;
                #    dG_rv.upper_bound = self.dG_r_max;
                #    dG_rv.objective_coefficient = 0;
            else:
                dG_rv.lower_bound = self.dG_r_min
                dG_rv.upper_bound = self.dG_r_max;
            # make a continuous variable for dG0_r
            dG0_rv = Reaction('dG0_rv_' + r.id);
            if use_measured_dG0_r and thermodynamic_consistency_check[r.id] and r.id != 'NTD4': # ignore inconsistent reactions:
                if dG_r_coverage[r.id]>measured_dG_f_coverage_criteria:
                    dG0_rv.lower_bound = dG0_r[r.id]['dG_r_lb'];
                    dG0_rv.upper_bound = dG0_r[r.id]['dG_r_ub'];
                    dG0_rv.objective_coefficient = 0;
                else:
                    dG0_rv.lower_bound = self.dG0_r_min;
                    dG0_rv.upper_bound = self.dG0_r_max;
                    dG0_rv.objective_coefficient = 0;
            else:
                dG0_rv.lower_bound = self.dG0_r_min;
                dG0_rv.upper_bound = self.dG0_r_max;
                dG0_rv.objective_coefficient = 0;
            dG0_rv.variable_kind = 'continuous';
            # add constraints to the variables
            indicator.add_metabolites({indicator_plus: -r.upper_bound});#,indicator_minus: -r.lower_bound});
            indicator.add_metabolites({dG_r_constraint: self.K});#,indicator_minus: -r.lower_bound});
            dG_rv.add_metabolites({dG_r_constraint: 1.0})
            dG_rv.add_metabolites({conc_ln_constraint: -1.0})
            dG0_rv.add_metabolites({conc_ln_constraint: 1.0})
            cobra_model_irreversible.reactions.get_by_id(r.id).add_metabolites({indicator_plus: 1});#,indicator_minus: 1});
            # add indicator reactions to the model
            cobra_model_irreversible.add_reaction(indicator);
            cobra_model_irreversible.add_reaction(dG_rv);
            cobra_model_irreversible.add_reaction(dG0_rv);
            # record dG_rv variables
            dG0_r_dict[r.id] = dG0_rv
            # check to see if the model broke
            cobra_model_irreversible.optimize(solver='gurobi');
            if not cobra_model_irreversible.solution.f:
                print dG_rv.id + ' broke the model!';
                variables_break.append(dG_rv.id);
                #cobra_model_irreversible.remove_reactions(indicator)
                cobra_model_irreversible.remove_reactions(dG_rv)
                #cobra_model_irreversible.reactions.get_by_id(r.id).subtract_metabolites({indicator_plus:1})
                #indicator.subtract_metabolites({dG_r_constraint: self.K})
                dG_rv.lower_bound = self.dG_r_min;
                dG_rv.upper_bound = self.dG_r_max;
                #indicator.add_metabolites({dG_r_constraint: self.K})
                dG_rv.add_metabolites({dG_r_constraint: 1.0})
                cobra_model_irreversible.add_reaction(dG_rv);
            else:
                # record dG_rv variables
                dG_r_variables[r.id] = dG_rv;
            ## check to see if the model broke
            #cobra_model_irreversible.optimize(solver='gurobi');
            #if not cobra_model_irreversible.solution.f:
            #    print dG0_rv.id + ' broke the model!';
            #    variables_break.append(dG0_rv.id);
            #    #cobra_model_irreversible.remove_reactions(indicator)
            #    cobra_model_irreversible.remove_reactions(dG0_rv)
            #    #cobra_model_irreversible.reactions.get_by_id(r.id).subtract_metabolites({indicator_plus:1})
            #    dG0_rv.lower_bound = self.dG0_r_min;
            #    dG0_rv.upper_bound = self.dG0_r_max;
            #    dG0_rv.add_metabolites({conc_ln_constraint: 1.0})
            #    cobra_model_irreversible.add_reaction(dG0_rv);
        #for dG0_rv in dG0_r_dict.values():
        #    # add indicator reactions to the model
        #    cobra_model_irreversible.add_reaction(dG0_rv);
        #    # check to see if the model broke
        #    cobra_model_irreversible.optimize(solver='gurobi');
        #    if not cobra_model_irreversible.solution.f:
        #        print dG0_rv.id + ' broke the model!';
        #        variables_break.append(dG0_rv.id);
        #        #cobra_model_irreversible.remove_reactions(indicator)
        #        cobra_model_irreversible.remove_reactions(dG0_rv)
        #        #cobra_model_irreversible.reactions.get_by_id(r.id).subtract_metabolites({indicator_plus:1})
        #        dG0_rv.lower_bound = self.dG0_r_min;
        #        dG0_rv.upper_bound = self.dG0_r_max;
        #        dG0_rv.add_metabolites({conc_ln_constraint: 1.0})
        #        cobra_model_irreversible.add_reaction(dG0_rv);
        if return_concentration_variables and not return_dG0_r_variables:
            return conc_lnv_dict;
        if return_dG0_r_variables and not return_concentration_variables:
            return dG0_r_dict;
        if return_concentration_variables and return_dG0_r_variables:
            return conc_lnv_dict,dG0_r_dict;
    def _copy_solution(self,cobra_model_irreversible,cobra_model_copy):
        '''copy out the model solution'''
        x_dict_copy = {};
        #y_dict_copy = {};
        solution_copy = cobra_model_copy.solution.f;
        status_copy = cobra_model_copy.solution.status;
        cobra_model_irreversible.solution.f = solution_copy;
        cobra_model_irreversible.solution.status = status_copy;
        if solution_copy:
            for rxn in cobra_model_irreversible.reactions:
                x_dict_copy[rxn.id]=cobra_model_copy.solution.x_dict[rxn.id];
            #for met in cobra_model_irreversible.metabolites:
            #    y_dict_copy[met.id]=cobra_model_copy.solution.y_dict[met.id];
            cobra_model_irreversible.solution.x_dict = x_dict_copy;
            #cobra_model_irreversible.solution.y_dict = y_dict_copy;

    def tfba(self, cobra_model_irreversible, dG_r, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99, use_measured_dG_r=True, solver=None):
        '''performs thermodynamic flux balance analysis'''

        """based on the method described in 10.1529/biophysj.106.093138
        max Z=c'v
        Sv=0
        0<=vi<=zi*vmax, {i=1,...,r}
        dG_ri-K+K*zi<0, {i=1,...,r}
        dG0_ri+RT*SUM[sij*ln(xj)]=dG_ri, {i=1,...,r}

        simplified:
        vi>=0, {i=1,...,r}
        vi-zi*vmax<=0, {i=1,...,r}
        dG_ri-K+K*zi<0, {i=1,...,r}
        dG_ri/K-1+zi<0, {i=1,...,r}
        dG_ri/K+zi<1, {i=1,...,r}
        dG_ri/K+zi<=1-1e-4, {i=1,...,r}
        
        where:
        zi is a binary variable, zi {0,1}
        K is always large enough such that dG_ri-K<0 or K>dG_ri
        """    
        # copy the model:
        cobra_model_copy = cobra_model_irreversible.copy();
        # add constraints
        self._add_dG_r_constraints(cobra_model_copy,dG_r, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria, measured_dG_f_coverage_criteria, use_measured_dG_r);
        # optimize
        cobra_model_copy.optimize(solver='gurobi');
        # record the results
        self._copy_solution(cobra_model_irreversible,cobra_model_copy);
    def tfba_conc_ln(self,cobra_model_irreversible, measured_concentration, estimated_concentration, dG0_r, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check,
                              measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,
                              use_measured_concentrations=True,use_measured_dG0_r=True, solver=None):
        '''performs thermodynamic flux balance analysis with bounds on metabolite activity insteady of dG_r'''

        """based on the method described in 10.1529/biophysj.106.093138
        max Z=c'v
        Sv=0
        0<=vi<=zi*vmax, {i=1,...,r}
        dG_ri-K+K*zi<0, {i=1,...,r}
        dG0_ri+RT*SUM[sij*ln(xj)]=dG_ri, {i=1,...,r}

        simplified:
        vi>=0, {i=1,...,r}
        vi-zi*vmax<=0, {i=1,...,r}

        dG_ri-K+K*zi<0, {i=1,...,r}
        dG_ri/K-1+zi<0, {i=1,...,r}
        dG_ri/K+zi<1, {i=1,...,r}
        dG_ri/K+zi<=1-1e-6, {i=1,...,r}
        dG0_ri/K+RT*SUM[sij*ln(xj)]/K+zi<=1-1e-6, {i=1,...,r}
        
        dG0_ri+RT*SUM[sij*ln(xj)]=dG_ri, {i=1,...,r}
        -dG_ri+RT*SUM[sij*ln(xj)]=-dG0_ri, {i=1,...,r}
        
        where:
        zi is a binary variable, zi {0,1}
        K is always large enough such that dG_ri-K<0 or K>dG_ri
        """    
        # copy the model:
        cobra_model_copy = cobra_model_irreversible.copy();
        # add constraints
        self._add_conc_ln_constraints(cobra_model_copy,measured_concentration, estimated_concentration, dG0_r, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria, measured_dG_f_coverage_criteria,
                     use_measured_concentrations,use_measured_dG0_r);
        # optimize
        cobra_model_copy.optimize(solver='gurobi');
        # record the results
        self._copy_solution(cobra_model_irreversible,cobra_model_copy);
    def tfva(self, cobra_model_irreversible, dG_r, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99, use_measured_dG_r=True,
             reaction_list=None,fraction_of_optimum=1.0, solver=None,
             objective_sense="maximize", **solver_args):
        """performs thermodynamic flux variability analysis to find max/min flux values

        cobra_model : :class:`~cobra.core.Model`:

        reaction_list : list of :class:`~cobra.core.Reaction`: or their id's
            The id's for which FVA should be run. If this is None, the bounds
            will be comptued for all reactions in the model.

        fraction_of_optimum : fraction of optimum which must be maintained.
            The original objective reaction is constrained to be greater than
            maximal_value * fraction_of_optimum

        solver : string of solver name
            If None is given, the default solver will be used.

        """
        # copy the model:
        cobra_model_copy = cobra_model_irreversible.copy();
        reactions_copy = [r for r in cobra_model_copy.reactions];
        if reaction_list is None and "the_reactions" in solver_args:
            reaction_list = solver_args.pop("the_reactions")
            from warnings import warn
            warn("the_reactions is deprecated. Please use reaction_list=")
        if reaction_list is None:
            reaction_list = [r for r in cobra_model_copy.reactions];
        else:
            reaction_list = [cobra_model_copy.reactions.get_by_id(i) if isinstance(i, string_types) else i for i in reaction_list]
        # add dG_r constraints: # adding constraints here is slower!
        self._add_dG_r_constraints(cobra_model_copy,dG_r, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria, measured_dG_f_coverage_criteria,use_measured_dG_r);

        solver = solver_dict[get_solver_name() if solver is None else solver]
        lp = solver.create_problem(cobra_model_copy)
        solver.solve_problem(lp, objective_sense=objective_sense)
        solution = solver.format_solution(lp, cobra_model_copy)
        if solution.status != "optimal":
            raise ValueError("TFVA requires the solution status to be optimal, not "
                             + solution.status)
        # set all objective coefficients to 0
        for i, r in enumerate(cobra_model_copy.reactions):
            if r.objective_coefficient != 0 and r in reactions_copy: # check that we are not messing with added variables from dG_r constraints
                f = solution.x_dict[r.id]
                new_bounds = (f * fraction_of_optimum, f)
                solver.change_variable_bounds(lp, i, min(new_bounds), max(new_bounds))
                solver.change_variable_objective(lp, i, 0.)
        ## add dG_r constraints: # adding constraints here is faster!
        #self._add_dG_r_constraints(cobra_model_copy,dG_r, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria, measured_dG_f_coverage_criteria,use_measured_dG_r);
        # perform fva
        for r in reaction_list:
            ## print reaction for debugging
            #print 'TFVA for rxn ' + r.id
            i = cobra_model_copy.reactions.index(r)
            self.tfva_data[r.id] = {}
            solver.change_variable_objective(lp, i, 1.)
            solver.solve_problem(lp, objective_sense="maximize", **solver_args)
            self.tfva_data[r.id]["flux_ub"] = solver.get_objective_value(lp)
            solver.solve_problem(lp, objective_sense="minimize", **solver_args)
            self.tfva_data[r.id]["flux_lb"] = solver.get_objective_value(lp)
            self.tfva_data[r.id]['flux_units']= 'mmol*gDW-1*hr-1';
            # revert the problem to how it was before
            solver.change_variable_objective(lp, i, 0.)

    def tfva_dG_r(self, cobra_model_irreversible, dG_r, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99, use_measured_dG_r=True,
             reaction_list=None,fraction_of_optimum=1.0, solver=None,
             objective_sense="maximize", **solver_args):
        """performs thermodynamic dG_r variability analysis to find max/min dG_r values

        cobra_model : :class:`~cobra.core.Model`:

        reaction_list : list of :class:`~cobra.core.Reaction`: or their id's
            The id's for which FVA should be run. If this is None, the bounds
            will be comptued for all reactions in the model.

        fraction_of_optimum : fraction of optimum which must be maintained.
            The original objective reaction is constrained to be greater than
            maximal_value * fraction_of_optimum

        solver : string of solver name
            If None is given, the default solver will be used.

        """
        # copy the model:
        cobra_model_copy = cobra_model_irreversible.copy();
        reactions_copy = [r for r in cobra_model_copy.reactions];
        if reaction_list is None and "the_reactions" in solver_args:
            reaction_list = solver_args.pop("the_reactions")
            from warnings import warn
            warn("the_reactions is deprecated. Please use reaction_list=")
        if reaction_list is None:
            reaction_list = [r for r in cobra_model_copy.reactions];
        else:
            reaction_list = [cobra_model_copy.reactions.get_by_id(i) if isinstance(i, string_types) else i for i in reaction_list]
        # add dG_r constraints: # adding constraints here is slower!
        dG_r_variables = self._add_dG_r_constraints(cobra_model_copy,dG_r, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria, measured_dG_f_coverage_criteria,use_measured_dG_r,True);

        solver = solver_dict[get_solver_name() if solver is None else solver]
        lp = solver.create_problem(cobra_model_copy)
        solver.solve_problem(lp, objective_sense=objective_sense)
        solution = solver.format_solution(lp, cobra_model_copy)
        if solution.status != "optimal":
            raise ValueError("TFVA requires the solution status to be optimal, not "
                             + solution.status)
        # set all objective coefficients to 0
        for i, r in enumerate(cobra_model_copy.reactions):
            if r.objective_coefficient != 0:
                f = solution.x_dict[r.id]
                new_bounds = (f * fraction_of_optimum, f)
                solver.change_variable_bounds(lp, i, min(new_bounds), max(new_bounds))
                solver.change_variable_objective(lp, i, 0.)
        ## add dG_r constraints: # adding constraints here is faster!
        #dG_r_variables = self._add_dG_r_constraints(cobra_model_copy,dG_r, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria, measured_dG_f_coverage_criteria,use_measured_dG_r,True);
        # perform tfva on dG_r
        for r in dG_r_variables.values():
            # print reaction for debugging
            print 'TFVA_dG_r for rxn ' + r.id
            i = cobra_model_copy.reactions.index(r)
            self.tfva_dG_r_data[r.id] = {}
            solver.change_variable_objective(lp, i, 1.)
            solver.solve_problem(lp, objective_sense="maximize", **solver_args)
            self.tfva_dG_r_data[r.id]["dG_r_ub"] = solver.get_objective_value(lp)
            solver.solve_problem(lp, objective_sense="minimize", **solver_args)
            self.tfva_dG_r_data[r.id]["dG_r_lb"] = solver.get_objective_value(lp)
            self.tfva_dG_r_data[r.id]['flux_units']= 'mmol*gDW-1*hr-1';
            # revert the problem to how it was before
            solver.change_variable_objective(lp, i, 0.)

    def tfva_concentrations(self, cobra_model_irreversible, measured_concentration, estimated_concentration, dG0_r, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check,
                              measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,
                             use_measured_concentrations=True,use_measured_dG0_r=True, reaction_list=None,fraction_of_optimum=1.0, solver=None,
                             objective_sense="maximize", **solver_args):
                '''performs thermodynamic metabolite concentration variability analysis'''
                # copy the model:
                cobra_model_copy = cobra_model_irreversible.copy();
                reactions_copy = [r for r in cobra_model_copy.reactions];
                if reaction_list is None and "the_reactions" in solver_args:
                    reaction_list = solver_args.pop("the_reactions")
                    from warnings import warn
                    warn("the_reactions is deprecated. Please use reaction_list=")
                if reaction_list is None:
                    reaction_list = [r for r in cobra_model_copy.reactions];
                else:
                    reaction_list = [cobra_model_copy.reactions.get_by_id(i) if isinstance(i, string_types) else i for i in reaction_list]
                # add constraints
                conc_ln_variables = self._add_conc_ln_constraints(cobra_model_copy,measured_concentration, estimated_concentration, dG0_r, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria, measured_dG_f_coverage_criteria,
                             use_measured_concentrations,use_measured_dG0_r,True,False);

                solver = solver_dict[get_solver_name() if solver is None else solver]
                lp = solver.create_problem(cobra_model_copy)
                solver.solve_problem(lp, objective_sense=objective_sense)
                solution = solver.format_solution(lp, cobra_model_copy)
                if solution.status != "optimal":
                    raise ValueError("TFVA requires the solution status to be optimal, not "
                                     + solution.status)
                # set all objective coefficients to 0
                for i, r in enumerate(cobra_model_copy.reactions):
                    if r.objective_coefficient != 0:
                        f = solution.x_dict[r.id]
                        new_bounds = (f * fraction_of_optimum, f)
                        solver.change_variable_bounds(lp, i, min(new_bounds), max(new_bounds))
                        solver.change_variable_objective(lp, i, 0.)
                ## add constraints
                #conc_ln_variables = self._add_conc_ln_constraints(cobra_model_copy,measured_concentration, estimated_concentration, dG0_r, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria, measured_dG_f_coverage_criteria,
                #             use_measured_concentrations,use_measured_dG0_r, True,False);
                # perform tfva on dG_r
                for r in conc_ln_variables.values():
                    # print reaction for debugging
                    print 'TFVA_dG_r for met ' + r.id
                    i = cobra_model_copy.reactions.index(r)
                    self.tfva_concentration_data[r.id] = {}
                    solver.change_variable_objective(lp, i, 1.)
                    solver.solve_problem(lp, objective_sense="maximize", **solver_args)
                    self.tfva_concentration_data[r.id]["concentration_ub"] = exp(solver.get_objective_value(lp)) #convert from ln
                    solver.solve_problem(lp, objective_sense="minimize", **solver_args)
                    self.tfva_concentration_data[r.id]["concentration_lb"] = exp(solver.get_objective_value(lp)) #convert from ln
                    self.tfva_concentration_data[r.id]['flux_units']= 'M';
                    # revert the problem to how it was before
                    solver.change_variable_objective(lp, i, 0.)

    def analyze_tfva_results(self,flux_threshold=1e-6):
        '''Determine what reactions are
        1. Blocked
        2. Essential
        3. Substitutable
        4. Constrained'''

        blocked_list = [];
        essential_list = [];
        substitutable_list = [];
        constrained_list = [];
        blocked_cnt = 0;
        essential_cnt = 0;
        substitutable_cnt = 0;
        constrained_cnt = 0;
        for k,v in self.tfva_data.iteritems():
            if v['flux_lb']<threshold and v['flux_ub']<threshold:
                self.tfva_analysis[k]['blocked'] = True;
                blocked_list.append(k);
                blocked_cnt+=1;
            else:
                self.tfva_analysis[k]['blocked'] = False;
            if v['flux_lb']>threshold and v['flux_ub']>threshold:
                self.tfva_analysis[k]['essential'] = True;
                essential_list.append(k);
                essential_cnt+=1;
            else:
                self.tfva_analysis[k]['essential'] = False;
            if v['flux_lb']<threshold and v['flux_ub']>threshold:
                self.tfva_analysis[k]['substitutable'] = True;
                substitutable_list.append(k);
                substitutable_cnt+=1;
            else:
                self.tfva_analysis[k]['substitutable'] = False;
            if v['flux_lb']-v['flux_ub']<threshold:
                self.tfva_analysis[k]['constrained'] = True;
            else:
                self.tfva_analysis[k]['constrained'] = False;
                constrained_list.append(k);
                constrained_cnt+=1;

        print "blocked reactions (" + str(blocked_cnt) + "): " + blocked_list;
        print "essential reactions (" + str(essential_cnt) + "): " + essential_list;
        print "substitutable reactions (" + str(substitutable_cnt) + "): " + substitutable_list;
        print "constrained reactions (" + str(constrained_cnt) + "): " + constrained_list;

    def tsampling_matlab_export(self, cobra_model_irreversible, dG_r, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99, use_measured_dG_r=True, solver=None,
                  fraction_optimal = 1.0):
        '''performs sampling with bounds on dG_r values'''

        # copy the model:
        cobra_model_copy = cobra_model_irreversible.copy();
        # add constraints
        self._add_dG_r_constraints(cobra_model_copy,dG_r, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria, measured_dG_f_coverage_criteria, use_measured_dG_r);
        # optimize
        cobra_model_copy.optimize(solver='gurobi');
        # confine the objective to a fraction of maximum optimal
        objective = [x.id for x in cobra_model_copy.reactions if x.objective_coefficient == 1]
        cobra_model_copy.reactions.get_by_id(objective[0]).upper_bound = fraction_optimal * cobra_model_copy.solution.f;
        # write model to mat
        save_matlab_model(cobra_model_copy,'data\\tsampling\\tsampler.mat');
        ## write model to xml
        #write_sbml_model(cobra_model_copy,'data\\tsampling\\tsampler.xml');
        # write the sampling script to file\
        mat_script = "% initialize with Tomlab_CPLEX\n"+\
                      "load('C:\\Users\\dmccloskey-sbrg\\Documents\\MATLAB\\tsampling\\tsampler.mat')\n"+\
                      "initCobraToolbox();\n"+\
                      "% sample\n"+\
                      "[tsampler_out, mixedFrac] = gpSampler(" + cobra_model_irreversible.description + ", [], [], [], [], [], true);\n"+\
                      "[tsampler_out, mixedFrac] = gpSampler(tsampler_out, [], [], [], 20000, [], true);";
        with open('data\\tsampling\\tsampler_run.m','w') as f:
            f.write(mat_script);

    def tsampling_matlab_import(self,cobra_model_irreversible, dG_r, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99, use_measured_dG_r=True, solver=None,
                  fraction_optimal = 1.0, matlab_data='data\\tsampling\\tsampler_out.mat',sampler_model_name='tsampler_out',plot_reactions=[]):
        '''performs sampling with bounds on dG_r values'''

        # copy the model:
        cobra_model_copy = cobra_model_irreversible.copy();
        # add constraints
        dG_r_variables = self._add_dG_r_constraints(cobra_model_copy,dG_r, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria, measured_dG_f_coverage_criteria, use_measured_dG_r, True);
        # optimize
        cobra_model_copy.optimize(solver='gurobi');
        # confine the objective to a fraction of maximum optimal
        objective = [x.id for x in cobra_model_copy.reactions if x.objective_coefficient == 1]
        cobra_model_copy.reactions.get_by_id(objective[0]).upper_bound = fraction_optimal * cobra_model_copy.solution.f;
        # process matlab sampling data
        tsampling = thermodynamics_sampling();
        tsampling.get_points_matlab(matlab_data,sampler_model_name);
        if tsampling.check_loops(cobra_model_copy):
            tsampling.simulate_loops(cobra_model_copy,'data\\tsampling\\tsampler_loops.json')
            tsampling.find_loops('data\\tsampling\\tsampler_loops.json');
            tsampling.remove_loops();
        # tsampling.export_points('data\\tsampling\\tsampler_points_loopless.json');
        dG_r_constraints = [];
        if plot_reactions:
            for rxn in plot_reactions:
                dG_r_constraints.append('dG_rv_'+rxn);
        else:
            for k,v in dG_r_variables.iteritems():
                dG_r_constraints.append(v.id);
        tsampling.plot_points(dG_r_constraints);
        # extract out dG_r data
        for k,v in dG_r_variables.iteritems():
            self.tsampling_dG_r_data[k] = {'dG_r_points':tsampling.points[v.id]['points'],
                                           'dG_r':tsampling.points[v.id]['mean'],
                                           'dG_r_var':pow(tsampling.points[v.id]['std'],2),
                                           'dG_r_lb':tsampling.points[v.id]['lb'],
                                           'dG_r_ub':tsampling.points[v.id]['ub']};

    def tsampling_conc_ln_matlab_export(self,cobra_model_irreversible, measured_concentration, estimated_concentration, dG0_r, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check,
                     measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,
                     use_measured_concentrations=True,use_measured_dG0_r=True, solver=None, fraction_optimal = 1.0):
        '''performs sampling with bounds on metabolite activity and dG0_r insteady of dG_r'''

       
        # copy the model:
        cobra_model_copy = cobra_model_irreversible.copy();
        # add constraints
        self._add_conc_ln_constraints_v1(cobra_model_copy,measured_concentration, estimated_concentration, dG0_r, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria, measured_dG_f_coverage_criteria,
                     use_measured_concentrations,use_measured_dG0_r);
        # optimize
        cobra_model_copy.optimize(solver='gurobi');
        # confine the objective to a fraction of maximum optimal
        objective = [x.id for x in cobra_model_copy.reactions if x.objective_coefficient == 1]
        cobra_model_copy.reactions.get_by_id(objective[0]).upper_bound = fraction_optimal * cobra_model_copy.solution.f;
        # write model to mat
        save_matlab_model(cobra_model_copy,'data\\tsampling\\tsampler_conc_ln.mat');
        ## write model to xml
        #write_sbml_model(cobra_model_copy,'data\\tsampling\\tsampler_conc_ln.xml');
        # write the sampling script to file\
        mat_script = "% initialize with Tomlab_CPLEX\n"+\
                      "load('C:\\Users\\dmccloskey-sbrg\\Documents\\MATLAB\\tsampling\\tsampler_conc_ln.mat')\n"+\
                      "initCobraToolbox();\n"+\
                      "% sample\n"+\
                      "[tsampler_out, mixedFrac] = gpSampler(" + cobra_model_irreversible.description + ", [], [], [], [], [], true);\n"+\
                      "[tsampler_out, mixedFrac] = gpSampler(tsampler_out, [], [], [], 20000, [], true);";
        with open('data\\tsampling\\tsampler_conc_ln_run.m','w') as f:
            f.write(mat_script);

    def tfba_looplaw(self, cobra_model_irreversible, solver=None):
        '''performs thermodynamic flux balance analysis using the looplaw'''

        #TODO: Constraint formulation is not yet correct

        """based on the method described in 0006-3495/11/02/0544/10
        """    
        # copy the model:
        cobra_model_copy = cobra_model_irreversible.copy();
        # add constraints
        self._add_looplaw_constraints(cobra_model_copy);
        # optimize
        cobra_model_copy.optimize(solver='gurobi');
        # record the results
        self._copy_solution(cobra_model_irreversible,cobra_model_copy);
    def _add_looplaw_constraints(self, cobra_model_irreversible, return_looplaw_variables=False):
        '''add constraints for looplaw to the model'''

        # TODO: Correct constraints

        # find system boundaries and the objective reaction
        system_boundaries = [x.id for x in cobra_model_irreversible.reactions if x.boundary == 'system_boundary'];
        transporters = find_transportRxns(cobra_model_irreversible);
        ## determine the compartments:
        #compartments = list(set(cobra_model_irreversible.metabolites.list_attr('compartment')));
        # add variables and constraints to model for tfba
        reactions = [r for r in cobra_model_irreversible.reactions];
        looplaw_variables = {};
        internal_dict = {};
        # calculate the null basis:
        cobra_model_array = cobra_model_irreversible.to_array_based_model();
        N = null(cobra_model_array.S.todense()) #convert S from sparse to full and compute the nullspace
        not_in_null = [];
        for i in range(N.shape[0]):
            met_j = N[i,:];
            met_sum = met_j.sum();
            if met_sum<1e-6: not_in_null.append(cobra_model_array.reactions[i].id)
        # make constraints for Nint*looplaw=0
        internal_constraint = Metabolite('internal');
        internal_constraint._constraint_sense = 'E';
        internal_constraint._bound = 0.0;
        for i,r in enumerate(reactions):
            # ignore system_boundary, objective, and transport reactions
            if r.id in system_boundaries or r.id in not_in_null: #r.id in transporters:
                continue;
            # create a constraint for vi-zi*vmax<=0
            indicator_plus = Metabolite(r.id + '_plus');
            indicator_plus._constraint_sense = 'L';
            # create a constraint for vi+zi*vmax>=0
            indicator_minus = Metabolite(r.id + '_minus');
            indicator_minus._constraint_sense = 'G';
            # create additional constraint for looplawi <= zi - 1000(1-zi) if zi = 1, gi <=1 
            # create additional constraint for looplawi <= 1001*zi - 1000
            # create additional constraint for looplawi -1001*zi <= -1000
            looplaw_constraint_plus = Metabolite(r.id + '_looplaw_plus');
            looplaw_constraint_plus._constraint_sense = 'L';
            looplaw_constraint_plus._bound = -1000.0;
            # create additional constraint for looplawi >= -1000*zi + (1-zi)
            # create additional constraint for looplawi >= -1001*zi + 1
            # create additional constraint for looplawi + 1001*zi >= 1
            looplaw_constraint_minus = Metabolite(r.id + '_looplaw_minus');
            looplaw_constraint_minus._constraint_sense = 'G';
            looplaw_constraint_minus._bound = 1;
            ## create additional constraint for looplaw >= 1
            #looplaw_constraint_gzero = Metabolite(r.id + '_looplaw_gzero');
            #looplaw_constraint_gzero._constraint_sense = 'G';
            #looplaw_constraint_gzero._bound = 1;
            ## create additional constraint for looplaw <= -1
            #looplaw_constraint_lzero = Metabolite(r.id + '_looplaw_lzero');
            #looplaw_constraint_lzero._constraint_sense = 'L';
            #looplaw_constraint_lzero._bound = -1;
            # make a boolean indicator variable
            indicator = Reaction('indicator_' + r.id);
            indicator.lower_bound = 0;
            indicator.upper_bound = 1;
            indicator.variable_kind = 'integer';
            # make a continuous variable for looplaw
            looplaw = Reaction('looplaw_' + r.id);
            looplaw.lower_bound = -999.0;
            looplaw.upper_bound = 999.0;
            looplaw.objective_coefficient = 0;
            looplaw.variable_kind = 'continuous';
            # add constraints to the variables
            indicator.add_metabolites({indicator_plus: -r.upper_bound});#,indicator_minus: -r.lower_bound});
            indicator.add_metabolites({looplaw_constraint_plus: -1001});
            indicator.add_metabolites({looplaw_constraint_minus: 1001});
            looplaw.add_metabolites({looplaw_constraint_plus: 1.0})
            looplaw.add_metabolites({looplaw_constraint_minus: 1.0})
            for j in range(N[i,:].shape[1]):
                met_stoich = N[i,j]
                if met_stoich>1e-6 or met_stoich<-1e-6:
                    looplaw.add_metabolites({internal_constraint: met_stoich})
            #looplaw.add_metabolites({looplaw_constraint_gzero: 1.0})
            #looplaw.add_metabolites({looplaw_constraint_lzero: 1.0})
            cobra_model_irreversible.reactions.get_by_id(r.id).add_metabolites({indicator_plus: 1.0},{indicator_minus: 1.0});
            # add indicator reactions to the model
            cobra_model_irreversible.add_reaction(indicator);
            cobra_model_irreversible.add_reaction(looplaw);
            # record looplaw_plus variables
            looplaw_variables[r.id] = looplaw;
        if return_looplaw_variables:
            return looplaw_variables;
    def _add_conc_ln_constraints(self,cobra_model_irreversible, measured_concentration, estimated_concentration, dG0_r, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check,
                              measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,use_measured_concentrations=True,use_measured_dG0_r=True,return_concentration_variables=False,return_dG0_r_variables=False):
        # pre-process the data
        dG0_r = self._scale_dG_r(dG0_r);
        #measured_concentration = self._scale_conc(measured_concentration);
        #estimated_concentration = self._scale_conc(estimated_concentration);
        # initialize hydrogens:
        hydrogens = [];
        compartments = list(set(cobra_model_irreversible.metabolites.list_attr('compartment')));
        for compart in compartments:
             hydrogens.append('h_' + compart);
        # find system boundaries and the objective reaction
        system_boundaries = [x.id for x in cobra_model_irreversible.reactions if x.boundary == 'system_boundary'];
        objectives = [x.id for x in cobra_model_irreversible.reactions if x.objective_coefficient == 1];
        transporters = find_transportRxns(cobra_model_irreversible);
        # bounds
        dG_r_indicator_constraint = 1-1e-6;
        # add variables and constraints to model for tfba
        reactions = [r for r in cobra_model_irreversible.reactions];
        conc_lnv_dict = {}; # dictionary to record conc_ln variables
        dG0_r_dict = {};
        variables_break = [];
        for i,r in enumerate(reactions):
            if r.id in system_boundaries or r.id in objectives or r.id in transporters:
                continue;
            # create a constraint for vi-zi*vmax<=0
            indicator_plus = Metabolite(r.id + '_plus');
            indicator_plus._constraint_sense = 'L';
            ## create a constraint for vi+zi*vmax>=0
            #indicator_minus = Metabolite(r.id + '_minus');
            #indicator_minus._constraint_sense = 'G';
            # create additional constraint for RT*SUM[sij*ln(xj)]/K+zi<=1-1e-6-dG0_ri/K
            conc_ln_constraint = Metabolite(r.id + '_conc');
            conc_ln_constraint._constraint_sense = 'L';
            conc_ln_constraint._bound = self.K - dG_r_indicator_constraint
            # make continuous variables for conc
            products = [p for p in r.products];
            # Method 2:
            for p in products:
                if not(p.id in hydrogens): # exclude hydrogen because it has already been accounted for when adjusting for the pH
                    if use_measured_concentrations:# and thermodynamic_consistency_check[r.id]: # ignore inconsistent reactions:
                        if p.id in measured_concentration.keys():
                            # check that if the p_id has already been created:
                            if p.id in conc_lnv_dict.keys():
                                # add constraints
                                conc_lnv_dict[p.id].add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)}); #reaction is also updated in model automatically
                                #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[p.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + p.id);
                                conc_lnv.lower_bound = log(measured_concentration[p.id]['concentration_lb']);
                                conc_lnv.upper_bound = log(measured_concentration[p.id]['concentration_ub']);
                                conc_lnv.variable_kind = 'continuous';
                                conc_lnv.objective_coefficient = 0;
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                                # record the new variable
                                conc_lnv_dict[p.id] = conc_lnv;
                                # add the new variable to the model
                                cobra_model_irreversible.add_reaction(conc_lnv);
                        elif p.id in estimated_concentration.keys():
                            # check that if the p_id has already been created:
                            if p.id in conc_lnv_dict.keys():
                                # add constraints
                                conc_lnv_dict[p.id].add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                                #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[p.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + r.id + '_' + p.id);
                                conc_lnv.lower_bound = log(estimated_concentration[p.id]['concentration_lb']);
                                conc_lnv.upper_bound = log(estimated_concentration[p.id]['concentration_ub']);
                                conc_lnv.variable_kind = 'continuous';
                                conc_lnv.objective_coefficient = 0;
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                                # record the new variable
                                conc_lnv_dict[p.id] = conc_lnv;
                                # add the new variable to the model
                                cobra_model_irreversible.add_reaction(conc_lnv);
                    else:
                        # check that if the p_id has already been created:
                        if p.id in conc_lnv_dict.keys():
                            # add constraints
                            conc_lnv_dict[p.id].add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                            #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[p.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                        else:
                            conc_lnv = Reaction('conc_lnv_' + p.id);
                            conc_lnv.lower_bound = self.conc_min_ln;
                            conc_lnv.upper_bound = self.conc_max_ln;
                            conc_lnv.variable_kind = 'continuous';
                            conc_lnv.objective_coefficient = 0;
                            # add constraints
                            conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                            # record the new variable
                            conc_lnv_dict[p.id] = conc_lnv;
                            # add the new variable to the model
                            cobra_model_irreversible.add_reaction(conc_lnv);
            reactants = [react for react in r.reactants];
            # Method 2:
            for react in reactants:
                if not(react.id in hydrogens): # exclude hydrogen because it has already been accounted for when adjusting for the pH
                    if use_measured_concentrations:# and thermodynamic_consistency_check[r.id]: # ignore inconsistent reactions:
                        if react.id in measured_concentration.keys():
                            # check that if the react_id has already been created:
                            if react.id in conc_lnv_dict.keys():
                                # add constraints
                                conc_lnv_dict[react.id].add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                                #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[react.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + react.id);
                                conc_lnv.lower_bound = log(measured_concentration[react.id]['concentration_ub']);
                                conc_lnv.upper_bound = log(measured_concentration[react.id]['concentration_lb']);
                                conc_lnv.variable_kind = 'continuous';
                                conc_lnv.objective_coefficient = 0;
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                                # record the new variable
                                conc_lnv_dict[react.id] = conc_lnv;
                                # add the new variable to the model
                                cobra_model_irreversible.add_reaction(conc_lnv);
                        elif react.id in estimated_concentration.keys():
                            # check that if the react_id has already been created:
                            if react.id in conc_lnv_dict.keys():
                                # add constraints
                                conc_lnv_dict[react.id].add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                                #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[react.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + r.id + '_' + react.id);
                                conc_lnv.lower_bound = log(estimated_concentration[react.id]['concentration_ub']);
                                conc_lnv.upper_bound = log(estimated_concentration[react.id]['concentration_lb']);
                                conc_lnv.variable_kind = 'continuous';
                                conc_lnv.objective_coefficient = 0;
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                                # record the new variable
                                conc_lnv_dict[react.id] = conc_lnv;
                                # add the new variable to the model
                                cobra_model_irreversible.add_reaction(conc_lnv);
                    else:
                        # check that if the react_id has already been created:
                        if react.id in conc_lnv_dict.keys():
                            # add constraints
                            conc_lnv_dict[react.id].add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                            #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[react.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                        else:
                            conc_lnv = Reaction('conc_lnv_' + react.id);
                            conc_lnv.lower_bound = self.conc_min_ln;
                            conc_lnv.upper_bound = self.conc_max_ln;
                            conc_lnv.variable_kind = 'continuous';
                            conc_lnv.objective_coefficient = 0;
                            # add constraints
                            conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                            # record the new variable
                            conc_lnv_dict[react.id] = conc_lnv;
                            # add the new variable to the model
                            cobra_model_irreversible.add_reaction(conc_lnv);
            ## Method 1:
            #metabolites = [p for p in r.products] + [react for react in r.reactants];
            #for met in metabolites:
            #    if not(met.id in hydrogens): # exclude hydrogen because it has already been accounted for when adjusting for the pH
            #        if use_measured_concentrations:# and thermodynamic_consistency_check[r.id]: # ignore inconsistent reactions:
            #            if met.id in measured_concentration.keys():
            #                # check that if the met_id has already been created:
            #                if met.id in conc_lnv_dict.keys():
            #                    # add constraints
            #                    conc_lnv_dict[met.id].add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)});
            #                else:
            #                    conc_lnv = Reaction('conc_lnv_' + met.id);
            #                    conc_lnv.lower_bound = log(measured_concentration[met.id]['concentration_lb']);
            #                    conc_lnv.upper_bound = log(measured_concentration[met.id]['concentration_ub']);
            #                    conc_lnv.variable_kind = 'continuous';
            #                    conc_lnv.objective_coefficient = 0;
            #                    # add constraints
            #                    conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)});
            #                    # record the new variable
            #                    conc_lnv_dict[met.id] = conc_lnv;
            #            elif met.id in estimated_concentration.keys():
            #                # check that if the met_id has already been created:
            #                if met.id in conc_lnv_dict.keys():
            #                    # add constraints
            #                    conc_lnv_dict[met.id].add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)});
            #                else:
            #                    conc_lnv = Reaction('conc_lnv_' + r.id + '_' + met.id);
            #                    conc_lnv.lower_bound = log(estimated_concentration[met.id]['concentration_lb']);
            #                    conc_lnv.upper_bound = log(estimated_concentration[met.id]['concentration_ub']);
            #                    conc_lnv.variable_kind = 'continuous';
            #                    conc_lnv.objective_coefficient = 0;
            #                    # add constraints
            #                    conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)});
            #                    # record the new variable
            #                    conc_lnv_dict[met.id] = conc_lnv;
            #        else:
            #            # check that if the met_id has already been created:
            #            if met.id in conc_lnv_dict.keys():
            #                # add constraints
            #                conc_lnv_dict[met.id].add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)});
            #            else:
            #                conc_lnv = Reaction('conc_lnv_' + met.id);
            #                conc_lnv.lower_bound = self.conc_min_ln;
            #                conc_lnv.upper_bound = self.conc_max_ln;
            #                conc_lnv.variable_kind = 'continuous';
            #                conc_lnv.objective_coefficient = 0;
            #                # add constraints
            #                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)});
            #                # record the new variable
            #                conc_lnv_dict[met.id] = conc_lnv;
            # make a boolean indicator variable
            indicator = Reaction('indicator_' + r.id);
            indicator.lower_bound = 0;
            indicator.upper_bound = 1;
            indicator.variable_kind = 'integer';
            # add constraints
            indicator.add_metabolites({indicator_plus: -r.upper_bound});#,indicator_minus: -r.lower_bound});
            indicator.add_metabolites({conc_ln_constraint: self.K});#,indicator_minus: -r.lower_bound});
            # add indicator reactions to the model
            cobra_model_irreversible.add_reaction(indicator);
            cobra_model_irreversible.reactions.get_by_id(r.id).add_metabolites({indicator_plus: 1});#,indicator_minus: 1});
            # add constraints to the boolean indicator variable
            # make a continuous variable for dG0_r
            dG0_rv = Reaction('dG0_rv_' + r.id);
            if use_measured_dG0_r and thermodynamic_consistency_check[r.id] and r.id != 'NTD4': # ignore inconsistent reactions:
                if dG_r_coverage[r.id]>measured_dG_f_coverage_criteria:
                    dG0_rv.lower_bound = dG0_r[r.id]['dG_r_lb'];
                    dG0_rv.upper_bound = dG0_r[r.id]['dG_r_ub'];
                    dG0_rv.objective_coefficient = 0;
                else:
                    dG0_rv.lower_bound = self.dG0_r_min;
                    dG0_rv.upper_bound = self.dG0_r_max;
                    dG0_rv.objective_coefficient = 0;
            else:
                dG0_rv.lower_bound = self.dG0_r_min;
                dG0_rv.upper_bound = self.dG0_r_max;
                dG0_rv.objective_coefficient = 0;
            dG0_rv.variable_kind = 'continuous';
            # add constraints
            dG0_rv.add_metabolites({conc_ln_constraint: 1.0})
            # record dG_rv variables
            dG0_r_dict[r.id] = dG0_rv
            # add indicator reactions to the model
            cobra_model_irreversible.add_reaction(dG0_rv);
            ## check to see if the model broke
            #cobra_model_irreversible.optimize(solver='gurobi');
            #if not cobra_model_irreversible.solution.f:
            #    print dG0_rv.id + ' broke the model!';
            #    variables_break.append(dG0_rv.id);
            #    #cobra_model_irreversible.remove_reactions(indicator)
            #    cobra_model_irreversible.remove_reactions(dG0_rv)
            #    #cobra_model_irreversible.reactions.get_by_id(r.id).subtract_metabolites({indicator_plus:1})
            #    dG0_rv.lower_bound = self.dG0_r_min;
            #    dG0_rv.upper_bound = self.dG0_r_max;
            #    dG0_rv.add_metabolites({conc_ln_constraint: 1.0})
            #    cobra_model_irreversible.add_reaction(dG0_rv);
        #for dG0_rv in dG0_r_dict.values():
        #    # add indicator reactions to the model
        #    cobra_model_irreversible.add_reaction(dG0_rv);
        #    # check to see if the model broke
        #    cobra_model_irreversible.optimize(solver='gurobi');
        #    if not cobra_model_irreversible.solution.f:
        #        print dG0_rv.id + ' broke the model!';
        #        variables_break.append(dG0_rv.id);
        #        #cobra_model_irreversible.remove_reactions(indicator)
        #        cobra_model_irreversible.remove_reactions(dG0_rv)
        #        #cobra_model_irreversible.reactions.get_by_id(r.id).subtract_metabolites({indicator_plus:1})
        #        dG0_rv.lower_bound = self.dG0_r_min;
        #        dG0_rv.upper_bound = self.dG0_r_max;
        #        dG0_rv.add_metabolites({conc_ln_constraint: 1.0})
        #        cobra_model_irreversible.add_reaction(dG0_rv);
        if return_concentration_variables and not return_dG0_r_variables:
            return conc_lnv_dict;
        if return_dG0_r_variables and not return_concentration_variables:
            return dG0_r_dict;
        if return_concentration_variables and return_dG0_r_variables:
            return conc_lnv_dict,dG0_r_dict;