
from math import floor,ceil,log,sqrt,pow,exp,fabs
from copy import deepcopy
from cobra.core.metabolite import Metabolite
from cobra.core.reaction import Reaction
from cobra.io import save_matlab_model, write_sbml_model
from collections import Counter
from math import log, exp
from warnings import warn
import numpy
import scipy

from six import iteritems, string_types
from cobra.solvers import solver_dict, get_solver_name
from thermodynamics.thermodynamics_utility import find_transportRxns, null, find_transportMetsAndRxns
from thermodynamics.thermodynamics_sampling import thermodynamics_sampling
from thermodynamics.thermodynamics_io import thermodynamics_io
from cobra_utilities.cobra_simulatedData import cobra_simulatedData

# Other dependencies
import csv,json,sys

class thermodynamics_tfba(thermodynamics_io):    
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
        #y = 1000.0
        y = 500000.0
        #conc = 5.0e-5
        conc = 1.0e-5
        self.conc_min = 1/sqrt(y)*conc;
        self.conc_max = sqrt(y)*conc;
        self.conc_max_e = 1.0; # 1 M
        self.conc_min_ln = log(self.conc_min);
        self.conc_max_ln = log(self.conc_max);
        self.conc_max_e_ln = log(self.conc_max_e);
        # initialize constants
        self.R = 8.314e-3; # gas constant 8.3144621 [kJ/K/mol]
        self.K=10*self.dG_r_max; # an arbitrary constant larger than dG_r_max
        # T = 298.15; # temperature [K] (Not constant b/w compartments)
        self.A = 0.5093; # parameter A of the extended Debye-Huckel equation
                    # [mol^.5*L^-.5]
        self.B = 1.6; # parameter B of the extended Debye-Huckel equation
                 # [mol^.5*L^-.5]
        self.F = 96.485; # faraday's constant [kJ/V/mol]
        self.Rkcal = 1.9858775e-3; # gas constant 8.3144621 [kcal/K/mol]
        self.Fkcal = 23.062e-3; # faraday's constant [kcal/mV/mol]
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
        for k,v in dG_r.items():
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
        for k,v in concentrations.items():
            if v['concentration_lb']<self.conc_min:
                concentrations[k]['concentration_lb']=self.conc_min;
            #elif v['concentration_lb']>self.conc_max:
            #    concentrations[k]['concentration_lb']=self.conc;
            if v['concentration_ub']>self.conc_max:
                concentrations[k]['concentration_ub']=self.conc_max;  
            #elif v['concentration_ub']<self.conc_min:
            #    concentrations[k]['concentration_ub']=self.conc;
        return concentrations;

    def _add_dG_r_constraints(self, cobra_model_irreversible, dG_r, dG_r_coverage, thermodynamic_consistency_check,
                              use_measured_dG_r=True, return_dG_r_variables=False):
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
            if use_measured_dG_r and r.id in thermodynamic_consistency_check.keys() and thermodynamic_consistency_check[r.id]: # ignore inconsistent reactions
                #if r.id in dG_r_coverage.keys() and dG_r_coverage[r.id]>measured_dG_f_coverage_criteria:
                #if metabolomics_coverage[r.id] > measured_concentration_coverage_criteria and \
                #dG_r_coverage[r.id]>measured_dG_f_coverage_criteria:
                dG_rv.lower_bound = dG_r[r.id]['dG_r_lb'];
                dG_rv.upper_bound = dG_r[r.id]['dG_r_ub'];
                #else: 
                #    dG_rv.lower_bound = self.dG_r_min;
                #    dG_rv.upper_bound = self.dG_r_max;
            else:
                dG_rv.lower_bound = self.dG_r_min;
                dG_rv.upper_bound = self.dG_r_max;
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
            cobra_model_irreversible.reactions.get_by_id('dG_rv_' + r.id).objective_coefficient = 0
            # check to see if the model broke
            cobra_model_irreversible.solver = 'glpk'
            cobra_model_irreversible.optimize()
            if not cobra_model_irreversible.objective.value\
                or cobra_model_irreversible.solver.status == 'infeasible':
                print(dG_rv.id + ' broke the model!');
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
            # record dG_rv variables
            dG_r_variables[r.id] = dG_rv;
        if return_dG_r_variables:
            return dG_r_variables;

    def tfba(self, cobra_model_irreversible, dG_r, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, 
    measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99, use_measured_dG_r=True, solver='glpk'):
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
        reaction_list = [r.id for r in cobra_model_irreversible.reactions]
        # add constraints
        dG_r_constraints = self._add_dG_r_constraints(cobra_model_irreversible,dG_r, dG_r_coverage, thermodynamic_consistency_check, use_measured_dG_r);
        # optimize
        cobra_model_irreversible.solver = solver
        solution = cobra_model_irreversible.optimize()
        self.tfba_data = {k:v for k,v in solution.fluxes.to_dict().items() if k in reaction_list}

    def tfba_conc_ln(self,cobra_model_irreversible, measured_concentration, estimated_concentration, 
        dG0_r, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check,
        measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,
        use_measured_concentrations=True,use_measured_dG0_r=True, solver='glpk'):
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
        # add constraints
        reaction_list = [r.id for r in cobra_model_irreversible.reactions]
        conc_ln_variables = self._add_conc_ln_constraints(cobra_model_irreversible,measured_concentration, estimated_concentration, dG0_r, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, 
            measured_concentration_coverage_criteria, measured_dG_f_coverage_criteria,
            use_measured_concentrations,use_measured_dG0_r);
        # optimize
        cobra_model_irreversible.solver = solver
        solution = cobra_model_irreversible.optimize()
        self.tfba_data = {k:v for k,v in solution.fluxes.to_dict().items() if k in reaction_list}

    def tfva(self, cobra_model_irreversible, dG_r, dG_r_coverage, thermodynamic_consistency_check, use_measured_dG_r=True,
             reaction_list=None,fraction_of_optimum=1.0, solver='glpk',
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
        reaction_list = [r for r in cobra_model_irreversible.reactions]
        # add dG_r constraints: # adding constraints here is slower!
        self._add_dG_r_constraints(cobra_model_irreversible,dG_r, dG_r_coverage, thermodynamic_consistency_check, use_measured_dG_r);

        from cobra.flux_analysis import flux_variability_analysis
        fva_data = flux_variability_analysis(cobra_model_irreversible, fraction_of_optimum=1.0,
                                        objective_sense='maximize',
                                        reaction_list=reaction_list,
                                        )
        simulatedData = cobra_simulatedData()
        self.tfva_data = simulatedData._convert_fluxBounds2var(dict(zip(list(fva_data.index),fva_data.to_dict('records'))))

    def tfva_dG_r(self, cobra_model_irreversible, dG_r,  
        dG_r_coverage, thermodynamic_consistency_check, 
        use_measured_dG_r=True,
        reaction_list=None,fraction_of_optimum=1.0, solver='glpk',
        objective_sense="maximize", **solver_args):
        """performs thermodynamic dG_r variability analysis to find max/min dG_r values

        cobra_model : :class:`~cobra.core.Model`:

        reaction_list : list of :class:`~cobra.core.Reaction`: or their id's
            The id's for which FVA should be run. If this is None, the bounds
            will be computed for all reactions in the model.

        fraction_of_optimum : fraction of optimum which must be maintained.
            The original objective reaction is constrained to be greater than
            maximal_value * fraction_of_optimum

        solver : string of solver name
            If None is given, the default solver will be used.

        """
        # add dG_r constraints: # adding constraints here is slower!
        dG_r_variables = self._add_dG_r_constraints(cobra_model_irreversible,dG_r, dG_r_coverage, thermodynamic_consistency_check, use_measured_dG_r,True);

        from cobra.flux_analysis import flux_variability_analysis
        fva_data = flux_variability_analysis(cobra_model_irreversible, fraction_of_optimum=0.9,
                                        objective_sense='maximize',
                                        reaction_list=list(dG_r_variables.values()),
                                        )
        simulatedData = cobra_simulatedData()
        self.tfva_dG_r_data = simulatedData._convert_fluxBounds2var(dict(zip(list(fva_data.index),fva_data.to_dict('records'))))

    def tfva_concentrations(self, cobra_model_irreversible, measured_concentration, estimated_concentration, 
        dG0_r, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check,
        measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,
        use_measured_concentrations=True,use_measured_dG0_r=True, reaction_list=None,fraction_of_optimum=1.0, solver='glpk',
        objective_sense="maximize", **solver_args):
        '''performs thermodynamic metabolite concentration variability analysis'''

        # add constraints
        conc_ln_variables = self._add_conc_ln_constraints(cobra_model_irreversible,measured_concentration, estimated_concentration, dG0_r, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check, measured_concentration_coverage_criteria, measured_dG_f_coverage_criteria,
                        use_measured_concentrations,use_measured_dG0_r,True,False);

        from cobra.flux_analysis import flux_variability_analysis
        fva_data = flux_variability_analysis(cobra_model_irreversible, fraction_of_optimum=0.9,
                                        objective_sense='maximize',
                                        reaction_list=list(conc_ln_variables.values()),
                                        )
        simulatedData = cobra_simulatedData()
        self.tfva_concentrations_data = simulatedData._convert_fluxBounds2var(dict(zip(list(fva_data.index),fva_data.to_dict('records'))))

    def analyze_tfva_results(self,threshold=1e-6):
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
        for k,v in self.tfva_data.items():
            self.tfva_analysis[k] = {
                'blocked':False,
                'essential':False,
                'substitutable':False,
                'constrained':False,
            }
            if v['flux_lb']<threshold and v['flux_ub']<threshold:
                self.tfva_analysis[k]['blocked'] = True;
                blocked_list.append(k);
                blocked_cnt+=1;
            if v['flux_lb']>threshold and v['flux_ub']>threshold:
                self.tfva_analysis[k]['essential'] = True;
                essential_list.append(k);
                essential_cnt+=1;
            if v['flux_lb']<threshold and v['flux_ub']>threshold:
                self.tfva_analysis[k]['substitutable'] = True;
                substitutable_list.append(k);
                substitutable_cnt+=1;
            if v['flux_lb']-v['flux_ub']<threshold:
                self.tfva_analysis[k]['constrained'] = True;
                constrained_list.append(k);
                constrained_cnt+=1;

        # print("blocked reactions (" + str(blocked_cnt) + "): " + str(blocked_list));
        # print("essential reactions (" + str(essential_cnt) + "): " + str(essential_list));
        # print("substitutable reactions (" + str(substitutable_cnt) + "): " + str(substitutable_list))
        # print("constrained reactions (" + str(constrained_cnt) + "): " + str(constrained_list));

    def tsampling_matlab_import(self,cobra_model_irreversible, dG_r, dG_r_coverage, thermodynamic_consistency_check, use_measured_dG_r=True, solver='glpk',
                  fraction_optimal = 1.0, matlab_data='data\\tsampling\\tsampler_out.mat',sampler_model_name='tsampler_out',plot_reactions=[]):
        '''performs sampling with bounds on dG_r values'''

        # copy the model:
        cobra_model_copy = cobra_model_irreversible.copy();
        # add constraints
        dG_r_variables = self._add_dG_r_constraints(cobra_model_copy,dG_r, dG_r_coverage, thermodynamic_consistency_check, use_measured_dG_r, True);
        # optimize
        cobra_model_copy.solver = 'glpk'
        cobra_model_copy.optimize()
        # confine the objective to a fraction of maximum optimal
        objective = [x.id for x in cobra_model_copy.reactions if x.objective_coefficient == 1]
        cobra_model_copy.reactions.get_by_id(objective[0]).upper_bound = fraction_optimal * cobra_model_copy.objective.value;
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
            for k,v in dG_r_variables.items():
                dG_r_constraints.append(v.id);
        tsampling.plot_points(dG_r_constraints);
        # extract out dG_r data
        for k,v in dG_r_variables.items():
            self.tsampling_dG_r_data[k] = {'dG_r_points':tsampling.points[v.id]['points'],
                                           'dG_r':tsampling.points[v.id]['mean'],
                                           'dG_r_var':pow(tsampling.points[v.id]['std'],2),
                                           'dG_r_lb':tsampling.points[v.id]['lb'],
                                           'dG_r_ub':tsampling.points[v.id]['ub']};

    def _add_conc_ln_constraints(self,cobra_model_irreversible, measured_concentration, estimated_concentration, 
        dG0_r, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check,
        measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,
        use_measured_concentrations=True,use_measured_dG0_r=True,
        return_concentration_variables=False,return_dG0_r_variables=False):
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
        original_metabolites = cobra_model_irreversible.metabolites[:];
        conc_lnv_dict = {}; # dictionary to record conc_ln variables
        dG0_r_dict = {};
        variables_break = [];
        for i,r in enumerate(cobra_model_irreversible.reactions[:]):
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
                        if p.id in list(measured_concentration.keys()):
                            # check that if the p_id has already been created:
                            if p.id in list(conc_lnv_dict.keys()):
                                # add constraints
                                conc_lnv_dict[p.id].add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)}); #reaction is also updated in model automatically
                                #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[p.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + p.id);
                                conc_lnv.lower_bound = log(measured_concentration[p.id]['concentration_lb']);
                                conc_lnv.upper_bound = log(measured_concentration[p.id]['concentration_ub']);
                                conc_lnv.variable_kind = 'continuous';
                                
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                                # record the new variable
                                conc_lnv_dict[p.id] = conc_lnv;
                                # add the new variable to the model
                                cobra_model_irreversible.add_reaction(conc_lnv);
                        elif p.id in list(estimated_concentration.keys()):
                            # check that if the p_id has already been created:
                            if p.id in list(conc_lnv_dict.keys()):
                                # add constraints
                                conc_lnv_dict[p.id].add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                                #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[p.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + r.id + '_' + p.id);
                                conc_lnv.lower_bound = log(estimated_concentration[p.id]['concentration_lb']);
                                conc_lnv.upper_bound = log(estimated_concentration[p.id]['concentration_ub']);
                                conc_lnv.variable_kind = 'continuous';
                                
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                                # record the new variable
                                conc_lnv_dict[p.id] = conc_lnv;
                                # add the new variable to the model
                                cobra_model_irreversible.add_reaction(conc_lnv);
                    else:
                        # check that if the p_id has already been created:
                        if p.id in list(conc_lnv_dict.keys()):
                            # add constraints
                            conc_lnv_dict[p.id].add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                            #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[p.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                        else:
                            conc_lnv = Reaction('conc_lnv_' + p.id);
                            conc_lnv.lower_bound = self.conc_min_ln;
                            conc_lnv.upper_bound = self.conc_max_ln;
                            conc_lnv.variable_kind = 'continuous';
                            
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
                        if react.id in list(measured_concentration.keys()):
                            # check that if the react_id has already been created:
                            if react.id in list(conc_lnv_dict.keys()):
                                # add constraints
                                conc_lnv_dict[react.id].add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                                #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[react.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + react.id);
                                conc_lnv.lower_bound = log(measured_concentration[react.id]['concentration_ub']);
                                conc_lnv.upper_bound = log(measured_concentration[react.id]['concentration_lb']);
                                conc_lnv.variable_kind = 'continuous';
                                
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                                # record the new variable
                                conc_lnv_dict[react.id] = conc_lnv;
                                # add the new variable to the model
                                cobra_model_irreversible.add_reaction(conc_lnv);
                        elif react.id in list(estimated_concentration.keys()):
                            # check that if the react_id has already been created:
                            if react.id in list(conc_lnv_dict.keys()):
                                # add constraints
                                conc_lnv_dict[react.id].add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                                #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[react.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + r.id + '_' + react.id);
                                conc_lnv.lower_bound = log(estimated_concentration[react.id]['concentration_ub']);
                                conc_lnv.upper_bound = log(estimated_concentration[react.id]['concentration_lb']);
                                conc_lnv.variable_kind = 'continuous';
                                
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                                # record the new variable
                                conc_lnv_dict[react.id] = conc_lnv;
                                # add the new variable to the model
                                cobra_model_irreversible.add_reaction(conc_lnv);
                    else:
                        # check that if the react_id has already been created:
                        if react.id in list(conc_lnv_dict.keys()):
                            # add constraints
                            conc_lnv_dict[react.id].add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                            #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[react.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                        else:
                            conc_lnv = Reaction('conc_lnv_' + react.id);
                            conc_lnv.lower_bound = self.conc_min_ln;
                            conc_lnv.upper_bound = self.conc_max_ln;
                            conc_lnv.variable_kind = 'continuous';
                            
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
            #                    
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
            #                    
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
            #                
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
            if use_measured_dG0_r and r.id in thermodynamic_consistency_check.keys() and thermodynamic_consistency_check[r.id] and r.id != 'NTD4': # ignore inconsistent reactions:
                if r.id in dG_r_coverage.keys() and dG_r_coverage[r.id]>measured_dG_f_coverage_criteria:
                    dG0_rv.lower_bound = dG0_r[r.id]['dG_r_lb'];
                    dG0_rv.upper_bound = dG0_r[r.id]['dG_r_ub'];
                    
                else:
                    dG0_rv.lower_bound = self.dG0_r_min;
                    dG0_rv.upper_bound = self.dG0_r_max;
                    
            else:
                dG0_rv.lower_bound = self.dG0_r_min;
                dG0_rv.upper_bound = self.dG0_r_max;
                
            dG0_rv.variable_kind = 'continuous';
            # add constraints
            dG0_rv.add_metabolites({conc_ln_constraint: 1.0})
            # record dG_rv variables
            dG0_r_dict[r.id] = dG0_rv
            # add indicator reactions to the model
            cobra_model_irreversible.add_reaction(dG0_rv);
            # check to see if the model broke
            cobra_model_irreversible.optimize();
            if not cobra_model_irreversible.objective.value or cobra_model_irreversible.solution.status == 'infeasible':
               print dG0_rv.id + ' broke the model!';
               variables_break.append(dG0_rv.id);
               #cobra_model_irreversible.remove_reactions(indicator)
               cobra_model_irreversible.remove_reactions(dG0_rv)
               #cobra_model_irreversible.reactions.get_by_id(r.id).subtract_metabolites({indicator_plus:1})
               dG0_rv.lower_bound = self.dG0_r_min;
               dG0_rv.upper_bound = self.dG0_r_max;
               dG0_rv.add_metabolites({conc_ln_constraint: 1.0})
               cobra_model_irreversible.add_reaction(dG0_rv);
        #for dG0_rv in dG0_r_dict.values():
        #    # add indicator reactions to the model
        #    cobra_model_irreversible.add_reaction(dG0_rv);
        #    # check to see if the model broke
        #    cobra_model_irreversible.optimize(solver='glpk');
        #    if not cobra_model_irreversible.objective.value or cobra_model_irreversible.solution.status == 'infeasible':
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

    def _add_conc_ln_constraints_transport(self,cobra_model_irreversible, measured_concentration, estimated_concentration,
        dG0_r, pH, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check,
        measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,
        use_measured_concentrations=True,use_measured_dG0_r=True,return_concentration_variables=False,return_dG0_r_variables=False,
        diagnose_I=False,diagnose_solver_I='glpk',diagnose_threshold_I=0.98,diagnose_break_I=0.1):
        '''
        Add thermodynamic constraints for concentrations and transport reactions

        Args
            cobra_model_irreversible = irreversible cobra model

            diagnose_I = check the growth rate after each constrain is added
            return_concentration_variables = return added concentration variables?
            return_dG0_r_variables = False add dG0_r variables?
            diagnose_solver_I = solver used in the diagnose FBA
            diagnose_threshold_I = % of orginal growth rate to flag a constrain
            diagnose_break_I = % of original growth rate to stop the diagnosis

        Returns
            cobra_model_irreversible = irreversible cobra model with dG0r and conc_ln constraints added
            conc_lnv_dict = dictionary of conc_ln variables
            dG0_r_dict = dictionary of dG0_r variables
            diagnosed_variables_O = dictionary of constraints that reduce the model solution by diagnose_threshold_I
        '''
        
        # pre-process the data
        dG0_r = self._scale_dG_r(dG0_r);
        #measured_concentration = self._scale_conc(measured_concentration);
        #estimated_concentration = self._scale_conc(estimated_concentration);
        if diagnose_I:
            diagnosed_variables_O = {};
            # original solution:
            cobra_model_irreversible.solver = diagnose_solver_I
            cobra_model_irreversible.optimize()
            sol_original = cobra_model_irreversible.objective.value
            diagnose_sol = cobra_model_irreversible.objective.value;
        # initialize hydrogens:
        hydrogens = [];
        compartments = list(set(cobra_model_irreversible.metabolites.list_attr('compartment')));
        for compart in compartments:
             hydrogens.append('h_' + compart);
        # find system boundaries and the objective reaction
        system_boundaries = [x.id for x in cobra_model_irreversible.reactions if x.boundary == 'system_boundary'];
        objectives = [x.id for x in cobra_model_irreversible.reactions if x.objective_coefficient == 1];
        transporters = find_transportRxns(cobra_model_irreversible);
        # adjustment for transport reactions (Henry et al, 2007, Biophysical Journal 92(5) 1792?1805)
        mets_trans = find_transportMetsAndRxns(cobra_model_irreversible);
        # bounds
        dG_r_indicator_constraint = 1-1e-6;
        # add variables and constraints to model for tfba
        reactions = [r for r in cobra_model_irreversible.reactions];
        conc_lnv_dict = {}; # dictionary to record conc_ln variables
        dG_r_variables = {};
        dG0_r_dict = {};
        dG_r_mem_dict = {};
        dG_r_pH_dict = {};
        variables_break = [];
        for i,r in enumerate(reactions[:]):
            if r.id in system_boundaries or r.id in objectives or r.id in transporters:
                continue;
            # create a constraint for vi-zi*vmax<=0
            indicator_plus = Metabolite(r.id + '_plus');
            indicator_plus._constraint_sense = 'L';
            indicator_plus._bound = 0;
            ## create a constraint for vi+zi*vmax>=0
            #indicator_minus = Metabolite(r.id + '_minus');
            #indicator_minus._constraint_sense = 'G';
            # create additional constraint for dG0_ri + RT*SUM[sij*ln(xj)] + K*zi<=K-1e-4
            conc_ln_constraint = Metabolite(r.id + '_conc');
            conc_ln_constraint._constraint_sense = 'L';
            conc_ln_constraint._bound = self.K - dG_r_indicator_constraint
            # make continuous variables for conc
            products = [p for p in r.products];
            # Method 2:
            for p in products:
                if not(p.id in hydrogens): # exclude hydrogen because it has already been accounted for when adjusting for the pH
                    if use_measured_concentrations:# and thermodynamic_consistency_check[r.id]: # ignore inconsistent reactions:
                        if p.id in list(measured_concentration.keys()):
                            # check that if the p_id has already been created:
                            if p.id in list(conc_lnv_dict.keys()):
                                # add constraints
                                conc_lnv_dict[p.id].add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)}); #reaction is also updated in model automatically
                                #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[p.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + p.id);
                                conc_lnv.lower_bound = log(measured_concentration[p.id]['concentration_lb']);
                                conc_lnv.upper_bound = log(measured_concentration[p.id]['concentration_ub']);
                                conc_lnv.variable_kind = 'continuous';
                                
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                                # record the new variable
                                conc_lnv_dict[p.id] = conc_lnv;
                                # add the new variable to the model
                                cobra_model_irreversible.add_reaction(conc_lnv);
                        elif p.id in list(estimated_concentration.keys()):
                            # check that if the p_id has already been created:
                            if p.id in list(conc_lnv_dict.keys()):
                                # add constraints
                                conc_lnv_dict[p.id].add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                                #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[p.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + p.id);
                                conc_lnv.lower_bound = log(estimated_concentration[p.id]['concentration_lb']);
                                conc_lnv.upper_bound = log(estimated_concentration[p.id]['concentration_ub']);
                                conc_lnv.variable_kind = 'continuous';
                                
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                                # record the new variable
                                conc_lnv_dict[p.id] = conc_lnv;
                                # add the new variable to the model
                                cobra_model_irreversible.add_reaction(conc_lnv);
                    else:
                        # check that if the p_id has already been created:
                        if p.id in list(conc_lnv_dict.keys()):
                            # add constraints
                            conc_lnv_dict[p.id].add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                            #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[p.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                        else:
                            conc_lnv = Reaction('conc_lnv_' + p.id);
                            conc_lnv.lower_bound = self.conc_min_ln;
                            conc_lnv.upper_bound = self.conc_max_ln;
                            conc_lnv.variable_kind = 'continuous';
                            # add constraints
                            conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)});
                            # record the new variable
                            conc_lnv_dict[p.id] = conc_lnv;
                            # add the new variable to the model
                            cobra_model_irreversible.add_reaction(conc_lnv);
                            cobra_model_irreversible.reactions.get_by_id('conc_lnv_' + p.id).objective_coefficient = 0
                    if r.id in mets_trans and p.name in mets_trans[r.id]: 
                        if p.id in list(dG_r_mem_dict.keys()):
                            # add constraints
                            dG_r_mem_dict[p.id].add_metabolites({conc_ln_constraint:1.0});
                        else:
                            dG_r_mem = Reaction('dG_r_mem_' + p.id)
                            dG_r_mem.lower_bound = fabs(r.get_coefficient(p.id))/2.0*p.charge/2.0*self.F*(33.3*r.get_coefficient(p.id)/fabs(r.get_coefficient(p.id))*pH[p.compartment]['pH']-143.33/2.0);
                            dG_r_mem.upper_bound = fabs(r.get_coefficient(p.id))/2.0*p.charge/2.0*self.F*(33.3*r.get_coefficient(p.id)/fabs(r.get_coefficient(p.id))*pH[p.compartment]['pH']-143.33/2.0);
                            dG_r_mem.variable_kind = 'continuous';
                            # add constraints
                            dG_r_mem.add_metabolites({conc_ln_constraint:1.0});
                            # record the new variable
                            dG_r_mem_dict[p.id] = dG_r_mem;
                            # add the new variable to the model
                            cobra_model_irreversible.add_reaction(dG_r_mem);
                            cobra_model_irreversible.reactions.get_by_id('dG_r_mem_' + p.id).objective_coefficient = 0
                else: 
                    if r.id in mets_trans and p.name in mets_trans[r.id]: 
                        if p.id in list(dG_r_mem_dict.keys()):
                            # add constraints
                            dG_r_mem_dict[p.id].add_metabolites({conc_ln_constraint:1.0});
                        else:
                            dG_r_mem = Reaction('dG_r_mem_' + p.id)
                            dG_r_mem.lower_bound = fabs(r.get_coefficient(p.id))/2.0*p.charge/2.0*self.F*(33.3*r.get_coefficient(p.id)/fabs(r.get_coefficient(p.id))*pH[p.compartment]['pH']-143.33/2.0);
                            dG_r_mem.upper_bound = fabs(r.get_coefficient(p.id))/2.0*p.charge/2.0*self.F*(33.3*r.get_coefficient(p.id)/fabs(r.get_coefficient(p.id))*pH[p.compartment]['pH']-143.33/2.0);
                            dG_r_mem.variable_kind = 'continuous';
                            # add constraints
                            dG_r_mem.add_metabolites({conc_ln_constraint:1.0});
                            # record the new variable
                            dG_r_mem_dict[p.id] = dG_r_mem;
                            # add the new variable to the model
                            cobra_model_irreversible.add_reaction(dG_r_mem);
                            cobra_model_irreversible.reactions.get_by_id('dG_r_mem_' + p.id).objective_coefficient = 0
                        if p.id in list(dG_r_pH_dict.keys()):
                            # add constraints
                            dG_r_pH_dict[p.id].add_metabolites({conc_ln_constraint:1.0});
                        else:
                            dG_r_pH = Reaction('dG_r_pH_' + p.id)
                            dG_r_pH.lower_bound = log(10)*self.R*temperature[p.compartment]['temperature']*pH[p.compartment]['pH']*r.get_coefficient(p.id)/2.0;
                            dG_r_pH.upper_bound = log(10)*self.R*temperature[p.compartment]['temperature']*pH[p.compartment]['pH']*r.get_coefficient(p.id)/2.0;
                            dG_r_pH.variable_kind = 'continuous';
                            # add constraints
                            dG_r_pH.add_metabolites({conc_ln_constraint:1.0});
                            # record the new variable
                            dG_r_pH_dict[p.id] = dG_r_pH;
                            # add the new variable to the model
                            cobra_model_irreversible.add_reaction(dG_r_pH);
                            cobra_model_irreversible.reactions.get_by_id('dG_r_pH_' + p.id).objective_coefficient = 0
            reactants = [react for react in r.reactants];
            # Method 2:
            for react in reactants:
                if not(react.id in hydrogens): # exclude hydrogen because it has already been accounted for when adjusting for the pH
                    if use_measured_concentrations:# and thermodynamic_consistency_check[r.id]: # ignore inconsistent reactions:
                        if react.id in list(measured_concentration.keys()):
                            # check that if the react_id has already been created:
                            if react.id in list(conc_lnv_dict.keys()):
                                # add constraints
                                conc_lnv_dict[react.id].add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                                #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[react.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + react.id);
                                #conc_lnv.lower_bound = log(measured_concentration[react.id]['concentration_ub']);
                                #conc_lnv.upper_bound = log(measured_concentration[react.id]['concentration_lb']);
                                conc_lnv.lower_bound = log(measured_concentration[react.id]['concentration_lb']);
                                conc_lnv.upper_bound = log(measured_concentration[react.id]['concentration_ub']);
                                conc_lnv.variable_kind = 'continuous';
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                                # record the new variable
                                conc_lnv_dict[react.id] = conc_lnv;
                                # add the new variable to the model
                                cobra_model_irreversible.add_reaction(conc_lnv);
                                cobra_model_irreversible.reactions.get_by_id('conc_lnv_' + react.id).objective_coefficient = 0
                        elif react.id in list(estimated_concentration.keys()):
                            # check that if the react_id has already been created:
                            if react.id in list(conc_lnv_dict.keys()):
                                # add constraints
                                conc_lnv_dict[react.id].add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                                #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[react.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + react.id);
                                #conc_lnv.lower_bound = log(estimated_concentration[react.id]['concentration_ub']);
                                #conc_lnv.upper_bound = log(estimated_concentration[react.id]['concentration_lb']);
                                conc_lnv.lower_bound = log(estimated_concentration[react.id]['concentration_lb']);
                                conc_lnv.upper_bound = log(estimated_concentration[react.id]['concentration_ub']);
                                conc_lnv.variable_kind = 'continuous';
                                #  #Errors with optlang
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                                # record the new variable
                                conc_lnv_dict[react.id] = conc_lnv;
                                # add the new variable to the model
                                cobra_model_irreversible.add_reaction(conc_lnv);
                                cobra_model_irreversible.reactions.get_by_id('conc_lnv_' + react.id).objective_coefficient = 0
                    else:
                        # check that if the react_id has already been created:
                        if react.id in list(conc_lnv_dict.keys()):
                            # add constraints
                            conc_lnv_dict[react.id].add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                            #cobra_model_irreversible.reactions.get_by_id(conc_lnv_dict[react.id].id).add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                        else:
                            conc_lnv = Reaction('conc_lnv_' + react.id);
                            conc_lnv.lower_bound = self.conc_min_ln;
                            conc_lnv.upper_bound = self.conc_max_ln;
                            conc_lnv.variable_kind = 'continuous';
                            #  #Errors with optlang
                            # add constraints
                            conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)});
                            # record the new variable
                            conc_lnv_dict[react.id] = conc_lnv;
                            # add the new variable to the model
                            cobra_model_irreversible.add_reaction(conc_lnv);
                            cobra_model_irreversible.reactions.get_by_id('conc_lnv_' + react.id).objective_coefficient = 0
                    if r.id in mets_trans and react.name in mets_trans[r.id]: 
                        if react.id in list(dG_r_mem_dict.keys()):
                            # add constraints
                            dG_r_mem_dict[react.id].add_metabolites({conc_ln_constraint:1.0});
                        else:
                            dG_r_mem = Reaction('dG_r_mem_' + react.id)
                            dG_r_mem.lower_bound = fabs(r.get_coefficient(react.id))/2.0*react.charge/2.0*self.F*(33.3*r.get_coefficient(react.id)/fabs(r.get_coefficient(react.id))*pH[react.compartment]['pH']-143.33/2.0);
                            dG_r_mem.upper_bound = fabs(r.get_coefficient(react.id))/2.0*react.charge/2.0*self.F*(33.3*r.get_coefficient(react.id)/fabs(r.get_coefficient(react.id))*pH[react.compartment]['pH']-143.33/2.0);
                            dG_r_mem.variable_kind = 'continuous';
                            # add constraints
                            dG_r_mem.add_metabolites({conc_ln_constraint:1.0});
                            # record the new variable
                            dG_r_mem_dict[react.id] = dG_r_mem;
                            # add the new variable to the model
                            cobra_model_irreversible.add_reaction(dG_r_mem);
                            cobra_model_irreversible.reactions.get_by_id('dG_r_mem_' + react.id).objective_coefficient = 0
                else: 
                    if r.id in mets_trans and react.name in mets_trans[r.id]: 
                        if react.id in list(dG_r_mem_dict.keys()):
                            # add constraints
                            dG_r_mem_dict[react.id].add_metabolites({conc_ln_constraint:1.0});
                        else:
                            dG_r_mem = Reaction('dG_r_mem_' + react.id)
                            dG_r_mem.lower_bound = fabs(r.get_coefficient(react.id))/2.0*react.charge/2.0*self.F*(33.3*r.get_coefficient(react.id)/fabs(r.get_coefficient(react.id))*pH[react.compartment]['pH']-143.33/2.0);
                            dG_r_mem.upper_bound = fabs(r.get_coefficient(react.id))/2.0*react.charge/2.0*self.F*(33.3*r.get_coefficient(react.id)/fabs(r.get_coefficient(react.id))*pH[react.compartment]['pH']-143.33/2.0);
                            dG_r_mem.variable_kind = 'continuous';
                            # add constraints
                            dG_r_mem.add_metabolites({conc_ln_constraint:1.0});
                            # record the new variable
                            dG_r_mem_dict[react.id] = dG_r_mem;
                            # add the new variable to the model
                            cobra_model_irreversible.add_reaction(dG_r_mem);
                            cobra_model_irreversible.reactions.get_by_id('dG_r_mem_' + react.id).objective_coefficient = 0
                        if react.id in list(dG_r_pH_dict.keys()):
                            # add constraints
                            dG_r_pH_dict[react.id].add_metabolites({conc_ln_constraint:1.0});
                        else:
                            dG_r_pH = Reaction('dG_r_pH_' + react.id)
                            dG_r_pH.lower_bound = log(10)*self.R*temperature[react.compartment]['temperature']*pH[react.compartment]['pH']*r.get_coefficient(react.id)/2.0;
                            dG_r_pH.upper_bound = log(10)*self.R*temperature[react.compartment]['temperature']*pH[react.compartment]['pH']*r.get_coefficient(react.id)/2.0;
                            dG_r_pH.variable_kind = 'continuous';
                            # add constraints
                            dG_r_pH.add_metabolites({conc_ln_constraint:1.0});
                            # record the new variable
                            dG_r_pH_dict[react.id] = dG_r_pH;
                            # add the new variable to the model
                            cobra_model_irreversible.add_reaction(dG_r_pH);
                            cobra_model_irreversible.reactions.get_by_id('dG_r_pH_' + react.id).objective_coefficient = 0
            # make a boolean indicator variable
            indicator = Reaction('indicator_' + r.id);
            indicator.lower_bound = 0;
            indicator.upper_bound = 1;
            indicator.variable_kind = 'integer';
            # make a continuous variable for dG0_r
            dG0_rv = Reaction('dG0_rv_' + r.id);
            dG0_rv.variable_kind = 'continuous';
            if use_measured_dG0_r and r.id in thermodynamic_consistency_check.keys() and thermodynamic_consistency_check[r.id]:# and r.id != 'NTD4': # ignore inconsistent reactions:
                if r.id in dG_r_coverage.keys() and dG_r_coverage[r.id]['measured_dG_f_coverage']>measured_dG_f_coverage_criteria:
                    dG0_rv.lower_bound = dG0_r[r.id]['dG_r_lb'];
                    dG0_rv.upper_bound = dG0_r[r.id]['dG_r_ub'];
                else:
                    dG0_rv.lower_bound = self.dG0_r_min;
                    dG0_rv.upper_bound = self.dG0_r_max;
            else:
                dG0_rv.lower_bound = self.dG0_r_min;
                dG0_rv.upper_bound = self.dG0_r_max;
            # add constraints to the variables
            indicator.add_metabolites({indicator_plus: -r.upper_bound});#,indicator_minus: -r.lower_bound});
            indicator.add_metabolites({conc_ln_constraint: self.K});#,indicator_minus: -r.lower_bound});
            dG0_rv.add_metabolites({conc_ln_constraint: 1.0})
            cobra_model_irreversible.reactions.get_by_id(r.id).add_metabolites({indicator_plus: 1.0});#,indicator_minus: 1});
            # add indicator reactions to the model
            cobra_model_irreversible.add_reaction(indicator);
            cobra_model_irreversible.add_reaction(dG0_rv);
            cobra_model_irreversible.reactions.get_by_id('dG0_rv_' + r.id).objective_coefficient = 0
            # record dG_rv variables
            dG0_r_dict[r.id] = dG0_rv
            if diagnose_I:
                # check to see if the model broke
                cobra_model_irreversible.solver = diagnose_solver_I
                cobra_model_irreversible.optimize();
                print(r.id + " solution: "+ str(cobra_model_irreversible.objective.value));
                if cobra_model_irreversible.objective.value<sol_original*diagnose_break_I:
                    diagnosed_variables_O[r.id]={'solution_before':diagnose_sol,
                                                 'solution_after':cobra_model_irreversible.objective.value};
                    return diagnosed_variables_O;
                elif cobra_model_irreversible.objective.value<diagnose_sol*diagnose_threshold_I:
                    diagnosed_variables_O[r.id]={'solution_before':diagnose_sol,
                                                 'solution_after':cobra_model_irreversible.objective.value};
                    diagnose_sol=cobra_model_irreversible.objective.value;

        #Returns
        if diagnose_I:
            return diagnosed_variables_O;
        if return_concentration_variables and not return_dG0_r_variables:
            return conc_lnv_dict;
        if return_dG0_r_variables and not return_concentration_variables:
            return dG0_r_dict;
        if return_concentration_variables and return_dG0_r_variables:
            return conc_lnv_dict,dG0_r_dict;

    def check_conc_ln_constraints_transport(self,cobra_model_irreversible, measured_concentration, estimated_concentration, dG0_r, pH, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check,
                              measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,
                              n_checks_I = 5,
                              diagnose_solver_I='glpk',diagnose_threshold_I=0.98,diagnose_break_I=0.1):
        '''Check conc_ln_constraints_trasport
        1. check without using measured concentrations or measured dG0_r
        2. check without using measured concentrations, but using measured dG0_r
        3. check using both measured concentrations and measured dG0_r
        reactions involved with variables found to break the model are
        changed from "feasible:True" to "feasible:False"

        Args
            n_checks_I = number of loops per check
            diagnose_solver_I = solver used in the diagnose FBA
            diagnose_threshold_I = % of orginal growth rate to flag a constrain
            diagnose_break_I = % of original growth rate to stop the diagnosis

        Returns
            thermodynamic_constraints_check = thermodynamic_consistency_check updated from the check
            inconsistent_tcc = list of feasible reactions that break the model when tfba constraints are added
            diagnose_variables_1 = results of check 1
            diagnose_variables_2 = results of check 2
            diagnose_variables_3 = results of check 3
        '''
        thermodynamic_constraints_check = thermodynamic_consistency_check;
        # check 1:
        diagnose_variables_1 = {}
        for i in range(n_checks_I):
            cobra_model_check = cobra_model_irreversible.copy();
            diagnose_variables_tmp = {}
            diagnose_variables_tmp = self._add_conc_ln_constraints_transport(cobra_model_check, measured_concentration, estimated_concentration, dG0_r, pH, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check,
                                  measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,
                                  use_measured_concentrations=False,use_measured_dG0_r=False,
                                  diagnose_I=True,diagnose_solver_I=diagnose_solver_I,
                                  diagnose_threshold_I=diagnose_threshold_I,diagnose_break_I=diagnose_break_I);
            if diagnose_variables_tmp:
                for k,v in diagnose_variables_tmp.items():
                    thermodynamic_constraints_check[k] = False;
                    diagnose_variables_1.update(diagnose_variables_tmp)
            else:
                break;
        # check 2:
        diagnose_variables_2 = {}
        for i in range(n_checks_I):
            cobra_model_check = cobra_model_irreversible.copy();
            diagnose_variables_tmp = {}
            diagnose_variables_tmp = self._add_conc_ln_constraints_transport(cobra_model_check, measured_concentration, estimated_concentration, dG0_r, pH, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check,
                                  measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,
                                  use_measured_concentrations=False,use_measured_dG0_r=True,
                                  diagnose_I=True,diagnose_solver_I=diagnose_solver_I,
                                  diagnose_threshold_I=diagnose_threshold_I,diagnose_break_I=diagnose_break_I);
            if diagnose_variables_tmp:
                for k,v in diagnose_variables_tmp.items():
                    thermodynamic_constraints_check[k] = False;
                    diagnose_variables_2.update(diagnose_variables_tmp)
            else:
                break;
        # check 3:
        diagnose_variables_3 = {}
        for i in range(n_checks_I):
            cobra_model_check = cobra_model_irreversible.copy();
            diagnose_variables_tmp = {}
            diagnose_variables_tmp = self._add_conc_ln_constraints_transport(cobra_model_check, measured_concentration, estimated_concentration, dG0_r, pH, temperature, metabolomics_coverage, dG_r_coverage, thermodynamic_consistency_check,
                                      measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,
                                      use_measured_concentrations=True,use_measured_dG0_r=True,
                                      diagnose_I=True,diagnose_solver_I=diagnose_solver_I,
                                      diagnose_threshold_I=diagnose_threshold_I,diagnose_break_I=diagnose_break_I);
            if diagnose_variables_tmp:
                for k,v in diagnose_variables_tmp.items():
                    thermodynamic_constraints_check[k] = False;
                    diagnose_variables_3.update(diagnose_variables_tmp)
            else:
                break;

        # list out all identified reactions
        # for k,v in thermodynamic_constraints_check.items(): 
        #     if not v: 
        #         print(k)
        inconsistent_tcc = [k for k,v in thermodynamic_constraints_check.items() if not v];
        return thermodynamic_constraints_check,inconsistent_tcc,diagnose_variables_1,diagnose_variables_2,diagnose_variables_3

    def get_variableTypeAndUnits(self,rxn_id):
        '''return the variable type and units based on the name of the rxn_id
        Args
        rxn_id = string
        OUTPUT
        type_O = string
        units_O = string
        '''
        type_O = None;
        units_O = None;
        if 'conc_lnv_' in rxn_id:
            type_O = 'metabolite_activity';
            units_O = 'ln(M)';
        elif 'dG0_rv_' in rxn_id:
            type_O = 'dG0_r';
            units_O = 'kJ*mol-1';
        elif 'dG_rv_' in rxn_id:
            type_O = 'dG_r';
            units_O = 'kJ*mol-1';
        elif 'indicator_' in rxn_id:
            type_O = 'indicator';
            units_O = '';
        elif 'dG_r_mem_' in rxn_id:
            type_O = 'dG_mem';
            units_O = 'kJ*mol-1';
        elif 'dG_r_pH_' in rxn_id:
            type_O = 'dG_pH';
            units_O = 'kJ*mol-1';
        else:
            type_O = 'flux';
            units_O = 'mmol*gDCW-1*hr-1';
        return type_O,units_O;

    def _add_osmoticPressure_constraints(self,cobra_model_irreversible):
        """In development..."""
        return
    def _add_ionicStrength_constraints(self,cobra_model_irreversible, measured_concentration):
        """In development..."""
        return
    def _add_solubility_constraints(self,cobra_model_irreversible, measured_concentration):
        """In development..."""
        return
    def _add_pH_constraints(self,cobra_model_irreversible, measured_concentration):
        """In development..."""
        return
    def _add_membranePotential_constraints(self,cobra_model_irreversible, measured_concentration):
        """In development..."""
        return
    def _add_electroneutrality_constraints(self,cobra_model_irreversible, measured_concentration, compartments):
        """In development..."""
        return