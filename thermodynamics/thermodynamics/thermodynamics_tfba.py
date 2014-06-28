from __future__ import with_statement
from math import floor,ceil,log,sqrt,pow,exp,fabs
from copy import deepcopy
from cobra.core.Metabolite import Metabolite
from cobra.core.Reaction import Reaction
from collections import Counter
from math import log, exp
from warnings import warn

from six import iteritems, string_types
from cobra.solvers import solver_dict, get_solver_name

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

    def _scale_dG_r(self,dG_r):
        '''scale dG_r lb/ub to be within pre-defined bounds'''
        # scale down the magnitude of dG_r
        for k,v in dG_r.iteritems():
            if v['dG_r_lb']<self.dG_r_min:
                dG_r[k]['dG_r_lb']=self.dG_r_min;
            if v['dG_r_ub']>self.dG_r_max:
                dG_r[k]['dG_r_ub']=self.dG_r_max;  
        return dG_r;
    def _scale_conc(self,concentrations):
        '''scale conc lb/ub to be within pre-defined bounds'''
        # scale down the magnitude of conc
        for k,v in concentrations.iteritems():
            if v['concentration_lb']<self.conc_min:
                concentrations[k]['concentration_lb']=self.conc_min;
            if v['concentration_ub']>self.conc_max:
                concentrations[k]['concentration_ub']=self.conc_max;  
        return concentrations;
    def _add_dG_r_constraints(self, cobra_model_irreversible, dG_r, use_measured_dG_r=True, return_dG_r_variables=False):
        '''add constraints for dG_r to the model'''
        # pre-process the data
        dG_r = self._scale_dG_r(dG_r);
        # bounds
        dG_r_indicator_constraint = 1-1e-6;
        # find system boundaries and the objective reaction
        system_boundaries = [x.id for x in cobra_model_irreversible.reactions if x.boundary == 'system_boundary'];
        objectives = [x.id for x in cobra_model_irreversible.reactions if x.objective_coefficient == 1];
        # add variables and constraints to model for tfba
        reactions = [r for r in cobra_model_irreversible.reactions];
        dG_r_variables = {};
        for i,r in enumerate(reactions):
            if r.id in system_boundaries or r.id in objectives:
                continue;
            # make a boolean indicator variable
            indicator = Reaction('indicator_' + r.id);
            indicator.lower_bound = 0;
            indicator.upper_bound = 1;
            indicator.variable_kind = 'integer';
            # make a continuous variable for dG_r
            dG_rv = Reaction('dG_rv_' + r.id);
            if use_measured_dG_r:
                dG_rv.lower_bound = dG_r[r.id]['dG_r_lb'];
                dG_rv.upper_bound = dG_r[r.id]['dG_r_ub'];
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
            dG_r_constraint._bound = dG_r_indicator_constraint
            # add constraints to the variables
            indicator.add_metabolites({indicator_plus: -r.upper_bound});#,indicator_minus: -r.lower_bound});
            indicator.add_metabolites({dG_r_constraint: 1});#,indicator_minus: -r.lower_bound});
            dG_rv.add_metabolites({dG_r_constraint: 1.0/self.K})
            cobra_model_irreversible.reactions.get_by_id(r.id).add_metabolites({indicator_plus: 1});#,indicator_minus: 1});
            # add indicator reactions to the model
            cobra_model_irreversible.add_reaction(indicator);
            cobra_model_irreversible.add_reaction(dG_rv);
            # record dG_rv variables
            dG_r_variables[r.id] = dG_rv;
        if return_dG_r_variables:
            return dG_r_variables;
    def _add_conc_ln_constraints(self,cobra_model_irreversible, measured_concentration, estimated_concentration, dG0_r, temperature,
                     use_measured_concentrations=True,use_measured_dG0_r=True,return_concentration_variables=False,return_dG0_r_variables=False):
        # pre-process the data
        dG0_r = self._scale_dG_r(dG0_r);
        measured_concentration = self._scale_conc(measured_concentration);
        estimated_concentration = self._scale_conc(estimated_concentration);
        # initialize hydrogens:
        hydrogens = [];
        compartments = list(set(cobra_model_irreversible.metabolites.list_attr('compartment')));
        for compart in compartments:
             hydrogens.append('h_' + compart);
        # find system boundaries and the objective reaction
        system_boundaries = [x.id for x in cobra_model_irreversible.reactions if x.boundary == 'system_boundary'];
        objectives = [x.id for x in cobra_model_irreversible.reactions if x.objective_coefficient == 1];
        # bounds
        dG_r_indicator_constraint = 1-1e-6;
        # add variables and constraints to model for tfba
        reactions = [r for r in cobra_model_irreversible.reactions];
        conc_lnv_dict = {}; # dictionary to record conc_ln variables
        dG0_r_dict = {};
        for i,r in enumerate(reactions):
            if r.id in system_boundaries or r.id in objectives:
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
            conc_ln_constraint._bound = dG_r_indicator_constraint
            # make continuous variables for conc
            metabolites = [p for p in r.products] + [react for react in r.reactants];
            for met in metabolites:
                if not(met.id in hydrogens): # exclude hydrogen because it has already been accounted for when adjusting for the pH
                    if use_measured_concentrations:
                        if met.id in measured_concentration.keys():
                            # check that if the met_id has already been created:
                            if met.id in conc_lnv_dict.keys():
                                # add constraints
                                conc_lnv_dict[met.id].add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)/self.K});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + met.id);
                                conc_lnv.lower_bound = log(measured_concentration[met.id]['concentration_lb']);
                                conc_lnv.upper_bound = log(measured_concentration[met.id]['concentration_ub']);
                                conc_lnv.variable_kind = 'continuous';
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)/self.K});
                                # record the new variable
                                conc_lnv_dict[met.id] = conc_lnv;
                        elif met.id in estimated_concentration.keys():
                            # check that if the met_id has already been created:
                            if met.id in conc_lnv_dict.keys():
                                # add constraints
                                conc_lnv_dict[met.id].add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)/self.K});
                            else:
                                conc_lnv = Reaction('conc_lnv_' + r.id + '_' + met.id);
                                conc_lnv.lower_bound = log(estimated_concentration[met.id]['concentration_lb']);
                                conc_lnv.upper_bound = log(estimated_concentration[met.id]['concentration_ub']);
                                conc_lnv.variable_kind = 'continuous';
                                # add constraints
                                conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)/self.K});
                                # record the new variable
                                conc_lnv_dict[met.id] = conc_lnv;
                    else:
                        # check that if the met_id has already been created:
                        if met.id in conc_lnv_dict.keys():
                            # add constraints
                            conc_lnv_dict[met.id].add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)/self.K});
                        else:
                            conc_lnv = Reaction('conc_lnv_' + met.id);
                            conc_lnv.lower_bound = self.conc_min_ln;
                            conc_lnv.upper_bound = self.conc_max_ln;
                            conc_lnv.variable_kind = 'continuous';
                            # add constraints
                            conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[met.compartment]['temperature']*r.get_coefficient(met.id)/self.K});
                            # record the new variable
                            conc_lnv_dict[met.id] = conc_lnv;
            # make a boolean indicator variable
            indicator = Reaction('indicator_' + r.id);
            indicator.lower_bound = 0;
            indicator.upper_bound = 1;
            indicator.variable_kind = 'integer';
            # add constraints
            indicator.add_metabolites({indicator_plus: -r.upper_bound});#,indicator_minus: -r.lower_bound});
            indicator.add_metabolites({conc_ln_constraint: 1});#,indicator_minus: -r.lower_bound});
            # add indicator reactions to the model
            cobra_model_irreversible.add_reaction(indicator);
            cobra_model_irreversible.reactions.get_by_id(r.id).add_metabolites({indicator_plus: 1});#,indicator_minus: 1});
            # add constraints to the boolean indicator variable
            ## make a continuous variable for dG_r
            #dG_rv = Reaction('dG_rv_' + r.id);
            #dG_rv.variable_kind = 'continuous';
            # make a continuous variable for dG0_r
            dG0_rv = Reaction('dG_rv_' + r.id);
            if use_measured_dG0_r:
                dG0_rv.lower_bound = dG0_r[r.id]['dG_r_lb'];
                dG0_rv.upper_bound = dG0_r[r.id]['dG_r_ub'];
            else:
                dG0_rv.lower_bound = self.dG0_r_min;
                dG0_rv.upper_bound = self.dG0_r_max;
            dG0_rv.variable_kind = 'continuous';
            # add constraints
            dG0_rv.add_metabolites({conc_ln_constraint: 1.0/self.K})
            # add indicator reactions to the model
            cobra_model_irreversible.add_reaction(dG0_rv);
            dG0_r_dict[r.id] = dG0_rv
        # add all indicator reactions for conc_ln to the model
        cobra_model_irreversible.add_reactions(conc_lnv_dict.values());
        if return_concentration_variables and not return_dG0_r_variables:
            return conc_lnv_dict;
        if return_dG0_r_variables and not return_concentration_variables:
            return dG0_r_dict;
        if return_concentration_variables and return_dG0_r_variables:
            return conc_lnv_dict,dG0_r_dict;

    def tfba(self, cobra_model_irreversible, dG_r, use_measured_dG_r=True, solver=None):
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
        self._add_dG_r_constraints(cobra_model_copy,dG_r,False);
        # optimize
        cobra_model_copy.optimize(solver='gurobi');
    def tfba_conc_ln(self,cobra_model_irreversible, measured_concentration, estimated_concentration, dG0_r, temperature,
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
        self._add_conc_ln_constraints(cobra_model_copy,measured_concentration, estimated_concentration, dG0_r, temperature,
                     use_measured_concentrations,use_measured_dG0_r);
        # optimize
        cobra_model_copy.optimize(solver='gurobi');
        print cobra_model_copy.solution.f;
    def tfva(self, cobra_model_irreversible, dG_r, use_measured_dG_r=True,
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
        self._add_dG_r_constraints(cobra_model_copy,dG_r,False);

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
        #self._add_dG_r_constraints(cobra_model_copy,dG_r,False);
        # perform fva
        for r in reaction_list:
            ## print reaction for debugging
            #print 'TFVA for rxn ' + r.id
            i = cobra_model_copy.reactions.index(r)
            self.tfva_data[r.id] = {}
            solver.change_variable_objective(lp, i, 1.)
            solver.solve_problem(lp, objective_sense="maximize", **solver_args)
            self.tfva_data[r.id]["flux_lb"] = solver.get_objective_value(lp)
            solver.solve_problem(lp, objective_sense="minimize", **solver_args)
            self.tfva_data[r.id]["flux_ub"] = solver.get_objective_value(lp)
            self.tfva_data[r.id]['flux_units']= 'mmol*gDW-1*hr-1';
            # revert the problem to how it was before
            solver.change_variable_objective(lp, i, 0.)

    def tfva_dG_r(self, cobra_model_irreversible, dG_r, use_measured_dG_r=True,
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
        dG_r_variables = self._add_dG_r_constraints(cobra_model_copy,dG_r,False,True);

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
        #dG_r_variables = self._add_dG_r_constraints(cobra_model_copy,dG_r,False,True);
        # perform tfva on dG_r
        for r in dG_r_variables.values():
            # print reaction for debugging
            print 'TFVA_dG_r for rxn ' + r.id
            i = cobra_model_copy.reactions.index(r)
            self.tfva_dG_r_data[r.id] = {}
            solver.change_variable_objective(lp, i, 1.)
            solver.solve_problem(lp, objective_sense="maximize", **solver_args)
            self.tfva_dG_r_data[r.id]["dG_r_lb"] = solver.get_objective_value(lp)
            solver.solve_problem(lp, objective_sense="minimize", **solver_args)
            self.tfva_dG_r_data[r.id]["dG_r_ub"] = solver.get_objective_value(lp)
            self.tfva_dG_r_data[r.id]['flux_units']= 'mmol*gDW-1*hr-1';
            # revert the problem to how it was before
            solver.change_variable_objective(lp, i, 0.)
    def tfva_concentrations(self, cobra_model_irreversible, measured_concentration, estimated_concentration, dG0_r, temperature,
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
        conc_ln_variables = self._add_conc_ln_constraints(cobra_model_copy,measured_concentration, estimated_concentration, dG0_r, temperature,
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
        #conc_ln_variables = self._add_conc_ln_constraints(cobra_model_copy,measured_concentration, estimated_concentration, dG0_r, temperature,
        #             use_measured_concentrations,use_measured_dG0_r, True, False);
        # perform tfva on dG_r
        for r in conc_ln_variables.values():
            # print reaction for debugging
            print 'TFVA_dG_r for met ' + r.id
            i = cobra_model_copy.reactions.index(r)
            self.tfva_concentration_data[r.id] = {}
            solver.change_variable_objective(lp, i, 1.)
            solver.solve_problem(lp, objective_sense="maximize", **solver_args)
            self.tfva_concentration_data[r.id]["concentration_lb"] = exp(solver.get_objective_value(lp)) #convert from ln
            solver.solve_problem(lp, objective_sense="minimize", **solver_args)
            self.tfva_concentration_data[r.id]["concentration_ub"] = exp(solver.get_objective_value(lp)) #convert from ln
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
        