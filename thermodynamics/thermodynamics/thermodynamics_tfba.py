from __future__ import with_statement
from math import floor,ceil,log,sqrt,pow,exp,fabs
from copy import deepcopy
from cobra.core.Metabolite import Metabolite
from cobra.core.Reaction import Reaction
from collections import Counter
from math import log

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
        # set min/max for metabolite activity
        y = 1000.0
        conc = 5.0e-5
        self.conc_lb = 1/sqrt(y)*conc;
        self.conc_ub = sqrt(y)*conc;
        self.conc_lb_ln = log(self.conc_lb);
        self.conc_ub_ln = log(self.conc_ub);
        # initialize constants
        self.R = 8.314e-3; # gas constant 8.3144621 [kJ/K/mol]
        self.K=10*self.dG_r_max; # an arbitrary constant larger than dG_r_max

    def _scale_dG_r(self,dG_r):
        '''scale dG_r lb/ub to be within pre-defined bounds'''
        # scale down the magnitude of dG_r
        for k,v in dG_r.iteritems():
            if v['dG_r_lb']<self.dG_r_min:
                dG_r[k]['dG_r_lb']=self.dG_r_min;
            if v['dG_r_ub']>self.dG_r_max:
                dG_r[k]['dG_r_ub']=self.dG_r_max;  
        return dG_r;
    def _scale_conc(self,measured_concentrations):
        '''scale conc lb/ub to be within pre-defined bounds'''
        # scale down the magnitude of conc
        for k,v in conc.iteritems():
            if v['concentration_lb']<self.conc_min:
                conc[k]['concentration_lb']=self.conc_min;
            if v['concentration_ub']>self.conc_max:
                conc[k]['concentration_ub']=self.conc_max;  
        return conc;
    def tfba(self, cobra_model_irreversible, dG_r):
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
        # pre-process the data
        dG_r = self._scale_dG_r(dG_r);
        # bounds
        dG_r_indicator_constraint = 1-1e-6;
        # add variables and constraints to model for tfba
        reactions = [r for r in cobra_model_irreversible.reactions];
        for i,r in enumerate(reactions):
            ## temperature
            #T = temperature[r.get_compartments()[0]]['temperature'];
            # make a boolean indicator variable
            indicator = Reaction('indicator_' + r.id);
            indicator.lower_bound = 0;
            indicator.upper_bound = 1;
            indicator.variable_kind = 'integer';
            # make a continuous variable for dG_r
            dG_rv = Reaction('dG_rv_' + r.id);
            #dG_rv.lower_bound = self.dG_r_min;
            #dG_rv.upper_bound = self.dG_r_max;
            dG_rv.lower_bound = dG_r[r.id]['dG_r_lb'];
            dG_rv.upper_bound = dG_r[r.id]['dG_r_ub'];
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
        # optimize
        cobra_model_irreversible.optimize(solver='gurobi');

        return;
    def tfba_conc_ln(self,cobra_model_irreversible, measured_concentration, estimated_concentration, dG0_r, pH, ionic_strength, temperature):
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
        # pre-process the data
        dG0_r = self._scale_dG_r(dG_r);
        measured_concentration = self._scale_conc(measured_concentration);
        estimated_concentration = self._scale_conc(estimated_concentration);
        # bounds
        dG_r_indicator_constraint = 1-1e-6;
        # add variables and constraints to model for tfba
        reactions = [r for r in cobra_model_irreversible.reactions];
        for i,r in enumerate(reactions):
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
            #dG_rv.lower_bound = self.dG_r_min;
            #dG_rv.upper_bound = self.dG_r_max;
            #dG_rv.variable_kind = 'continuous';
            # make a continuous variable for dG0_r
            dG0_rv = Reaction('dG_rv_' + r.id);
            dG0_rv.lower_bound = dG0_r[r.id]['dG_r_lb'];
            dG0_rv.upper_bound = dG0_r[r.id]['dG_r_ub'];
            dG0_rv.variable_kind = 'continuous';
            # add constraints
            dG0_rv.add_metabolites({conc_ln_constraint: 1.0/self.K})
            # add indicator reactions to the model
            cobra_model_irreversible.add_reaction(dG0_rv);
            # make continuous variables for conc
            for p in r.products:
                if not(p.id in hydrogens): # exclude hydrogen because it has already been accounted for when adjusting for the pH
                    #if p.id in measured_concentration.keys():
                    #    conc_lnv = Reaction('conc_lnv_' + r.id);
                    #    conc_lnv.lower_bound = log(measured_concentration[p.id]['concentration_lb']);
                    #    conc_lnv.upper_bound = log(measured_concentration[p.id]['concentration_ub']);
                    #    conc_lnv.variable_kind = 'continuous';
                    #    # add constraints
                    #    concv.add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)/self.K});
                    #    # add indicator reactions to the model
                    #    cobra_model_irreversible.add_reaction(conc_lnv);
                    #elif p.id in estimated_concentration.keys():
                    #    conc_lnv = Reaction('conc_lnv_' + r.id);
                    #    conc_lnv.lower_bound = log(estimated_concentration[p.id]['concentration_lb']);
                    #    conc_lnv.upper_bound = log(estimated_concentration[p.id]['concentration_ub']);
                    #    conc_lnv.variable_kind = 'continuous';
                    #    # add constraints
                    #    conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)/self.K});
                    #    # add indicator reactions to the model
                    #    cobra_model_irreversible.add_reaction(conc_lnv);
                    conc_lnv = Reaction('concv_' + r.id);
                    conc_lnv.lower_bound = self.conc_ln_min;
                    conc_lnv.upper_bound = self.conc_ln_max;
                    conc_lnv.variable_kind = 'continuous';
                    # add constraints
                    conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[p.compartment]['temperature']*r.get_coefficient(p.id)/self.K});
                    # add indicator reactions to the model
                    cobra_model_irreversible.add_reaction(conc_lnv);
            for react in r.products:
                if not(react.id in hydrogens): # exclude hydrogen because it has already been accounted for when adjusting for the pH
                    #if react.id in measured_concentration.keys():
                    #    conc_lnv = Reaction('conc_lnv_' + r.id);
                    #    conc_lnv.lower_bound = log(measured_concentration[react.id]['concentration_lb']);
                    #    conc_lnv.upper_bound = log(measured_concentration[react.id]['concentration_ub']);
                    #    conc_lnv.variable_kind = 'continuous';
                    #    # add constraints
                    #    concv.add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)/self.K});
                    #    # add indicator reactions to the model
                    #    cobra_model_irreversible.add_reaction(conc_lnv);
                    #elif react.id in estimated_concentration.keys():
                    #    conc_lnv = Reaction('conc_lnv_' + r.id);
                    #    conc_lnv.lower_bound = log(estimated_concentration[react.id]['concentration_lb']);
                    #    conc_lnv.upper_bound = log(estimated_concentration[react.id]['concentration_ub']);
                    #    conc_lnv.variable_kind = 'continuous';
                    #    # add constraints
                    #    conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)/self.K});
                    #    # add indicator reactions to the model
                    #    cobra_model_irreversible.add_reaction(conc_lnv);
                    conc_lnv = Reaction('concv_' + r.id);
                    conc_lnv.lower_bound = self.conc_ln_min;
                    conc_lnv.upper_bound = self.conc_ln_max;
                    conc_lnv.variable_kind = 'continuous';
                    # add constraints
                    conc_lnv.add_metabolites({conc_ln_constraint:self.R*temperature[react.compartment]['temperature']*r.get_coefficient(react.id)/self.K});
                    # add indicator reactions to the model
                    cobra_model_irreversible.add_reaction(conc_lnv);
        # optimize
        cobra_model_irreversible.optimize(solver='gurobi');

    def tfva(self, cobra_model_irreversible, measured_concentration, estimated_concentration, dG0_r, dG_r, temperature):
        '''performs thermodynamic flux variability analysis'''
        return;
    def tfva_dG_r(self, cobra_model_irreversible, measured_concentration, estimated_concentration, dG0_r, dG_r, temperature):
        '''performs thermodynamic dG_r variability analysis'''
        return;
    def tfva_concentrations(self, cobra_model_irreversible, measured_concentration, estimated_concentration, dG0_r, dG_r, temperature):
        '''performs thermodynamic metabolite concentration variability analysis'''
        return;