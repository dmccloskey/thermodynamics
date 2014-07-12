
import json, csv
from math import sqrt,exp,pow
from numpy import average, var, log
from copy import copy

from thermodynamics_io import thermodynamics_io

class thermodynamics_metabolomicsData(thermodynamics_io):
    """Class to handle metabolomics data"""

    def __init__(self):
        
        self.measured_concentrations = {}
        self.estimated_concentrations = {}
        self.measured_concentrations_intracellular = {}
        self.measured_concentrations_extracellular = {}

    def import_metabolomics_data(self, concentration_filename_I): #deprecated
        '''import measured values required for analysis'''

        self.measured_concentrations = self._checkInput_concentrations(self.import_values_json(concentration_filename_I));

    def import_metabolomics_data_intracellular(self, concentration_filename_I):
        '''import measured values required for analysis'''

        self.measured_concentrations_intracellular = self._checkInput_concentrations(self.import_values_json(concentration_filename_I));

    def import_metabolomics_data_extracellular(self, concentration_filename_I):
        '''import measured values required for analysis'''

        self.measured_concentrations_extracellular = self._checkInput_concentrations(self.import_values_json(concentration_filename_I));
        
    def format_metabolomics_data(self): #deprecated
        '''format data'''
        
        self.measured_concentrations = self._convert_metabolomics_names(self.measured_concentrations);
        self.measured_concentrations = self._compartementalize_concentrations(self.measured_concentrations);
        
    def format_metabolomics_data_intracellular(self):
        '''format data'''
        
        self.measured_concentrations_intracellular = self._convert_metabolomics_names(self.measured_concentrations_intracellular);
        self.measured_concentrations_intracellular = self._compartementalize_concentrations(self.measured_concentrations_intracellular);
        
    def format_metabolomics_data_extracellular(self):
        '''format data'''
        
        self.measured_concentrations_extracellular = self._convert_metabolomics_names(self.measured_concentrations_extracellular);
        self.measured_concentrations_extracellular = self._compartementalize_concentrations_extracellular(self.measured_concentrations_extracellular);
        
    def generate_estimated_metabolomics_data(self, cobra_model):
        '''generate estimated values'''

        self.estimated_concentrations = self._generalize_compartment2all_concentration(cobra_model);
        self.estimated_concentrations.update(self._generalize_externalCompartments2all_concentration(cobra_model));

    def check_metabolomics_data(self):
        '''check data integrity'''
        return

    def remove_measured_concentrations(self,mets_I):
        '''Remove measured metabolite concentrations'''

        for met in mets_I:
            v=self.measured_concentrations.pop(met);
        
    def _checkInput_concentrations(self, measured_values_I):
        """
        check concentration data input

        measured_values: measured values with variances
                                 {metabolite.id: {'concentration': float,
                                                 'concentration_var': float,
                                                 'concentration_units': 'M'}
        returns a dictionary: measured_values_O:
                                 {metabolite.id: {'concentration': float,
                                                 'concentration_var': float,
                                                 'concentration_units': 'M'}
        """
        # check units
        measured_values_O = {};
        for k,v in measured_values_I.iteritems():
            if v['concentration_units'] == 'M':
                measured_values_O[k] = v
            else:
                print (str(k) + ' has invalid units of ' + str(v) + ' and will be ignored')
        return measured_values_O;

    def _convert_cv2varAndmM2M_concentrations(self, measured_values):
        """
        convered measured concentration values from mM to M

        convert measured concentration values with a coefficient of varation
        (CV = SD/Ave*100; SD = CV/100*AVE) 

        measured_values: measured values with variances
                                 {metabolite.id: {'concentration': float,
                                                 'concentration_cv': float,
                                                 'concentration_units': 'mM'}
        returns a dictionary: measured_values_O:
                                 {metabolite.id: {'concentration': float,
                                                 'concentration_var': float,
                                                 'concentration_units': 'M'}
        """
        measured_values_O = {};
        for k,v in measured_values.iteritems():
             concM = v['concentration']*1e-3;
             concMvar = v['concentration_cv']/100*(v['concentration']*1e-3)*v['concentration_cv']/100*(v['concentration']*1e-3);
             measured_values_O[k + '_c'] = {'concentration': concM,
                                                 'concentration_var': concMvar,
                                                 'concentration_units': 'M'}
        return measured_values_O

    def _compartementalize_concentrations(self, measured_values):
        """
        add a compartment identifier to intracellular metabolites

        measured_values: measured values with variances
                                 {metabolite.id: {'concentration': float,
                                                 'concentration_var': float,
                                                 'concentration_units': 'M'}
        returns a dictionary: measured_values_O:
                                 {metabolite.id: {'concentration': float,
                                                 'concentration_var': float,
                                                 'concentration_units': 'M'}
        """
        measured_values_O = {};
        for k,v in measured_values.iteritems():
             measured_values_O[k + '_c'] = v

        return measured_values_O

    def _compartementalize_concentrations_extracellular(self, measured_values):
        """
        add a compartment identifier to extracellular metabolites

        measured_values: measured values with variances
                                 {metabolite.id: {'concentration': float,
                                                 'concentration_var': float,
                                                 'concentration_units': 'M'}
        returns a dictionary: measured_values_O:
                                 {metabolite.id: {'concentration': float,
                                                 'concentration_var': float,
                                                 'concentration_units': 'M'}
        """
        measured_values_O = {};
        for k,v in measured_values.iteritems():
             measured_values_O[k + '_e'] = v
             measured_values_O[k + '_p'] = v

        return measured_values_O

    def _convert_cv2lbubAndmM2M_concentrations(self, measured_values,min_value):
        """
        convered measured concentration values from mM to M

        convert measured concentration values with a coefficient of varation
        (CV = SD/Ave*100; SD = CV/100*AVE) to lb and ub
        currently use +/- SD
        lb = ave - sqrt(var) NOTE: if < 0 min_value is used instead
        ub = ave + sqrt(var)

        measured_values: measured values with variances
                                 {metabolite.id: {'concentration': float,
                                                 'concentration_cv': float,
                                                 'concentration_units': 'mM'}
        returns a dictionary: measured_values_O:
                                 {metabolite.id: {'concentration_lb': float,
                                                 'concentration_ub': float,
                                                 'concentration_units': 'M'}
        """
        measured_values_O = {};
        for k,v in measured_values.iteritems():
             concMlb = 0.0;
             concMub = 0.0;
             concMlb = v['concentration']*1e-3 - v['concentration_cv']/100*(v['concentration']*1e-3);
             concMub = v['concentration']*1e-3 + v['concentration_cv']/100*(v['concentration']*1e-3);
             if concMlb<0: concMlb = min_value;
             measured_values_O[k + '_c'] = {'concentration_lb': concMlb,
                                                 'concentration_ub': concMub,
                                                 'concentration_units': 'M'}
        return measured_values_O

    def _generalize_compartmentLBUB2all_concentration(self, cobra_model, lbub=None, exceptions=None):
        """
        takes a compartment and lb/ub for that compartment
        and updates each metabolite in that compartment

        allows for exceptions to the generalization as input

        cobra_model: a Model object

        lbub: 'compartment': {'concentration_lb': float,
                               'concentration_ub': float,
                               'concentration_units': 'M'}

        exceptions: metabolite.id: {'concentration_lb': float,
                               'concentration_ub': float,
                               'concentration_units': string}

        returns a dictionary: metabolite.id {'concentration_lb': float,
                               'concentration_ub': float,
                               'concentration_units': string}
        """
        if not(exceptions):
            exceptions = {};
            exceptions['pi_'] = {'concentration_lb': 1.0e-3,
                                   'concentration_ub': 1.0e-3,
                                   'concentration_units': 'M'};
            exceptions['h2o_'] = {'concentration_lb': 55.0,
                                   'concentration_ub': 55.0,
                                   'concentration_units': 'M'};
            exceptions['h2_'] = {'concentration_lb': 0.034e-3,
                                   'concentration_ub': 0.034e-3,
                                   'concentration_units': 'M'};
            exceptions['o2_'] = {'concentration_lb': 0.055e-3,
                                   'concentration_ub': 0.055e-3,
                                   'concentration_units': 'M'};
            exceptions['co2_'] = {'concentration_lb': 1.4e-3,
                                   'concentration_ub': 1.4e-3,
                                   'concentration_units': 'M'};

        if not(lbub):
            lbub = {};
            compartments = list(set(cobra_model.metabolites.list_attr('compartment')));
            for c in compartments:
                lbub[c] = {'concentration_lb':0.001e-6,
                                       'concentration_ub':0.01,
                                       'concentration_units':'M'};

        default_values = {};
        for m in cobra_model.metabolites:
            default_values[m.id] = {'concentration_lb': lbub[m.compartment]['concentration_lb'],
                                   'concentration_ub': lbub[m.compartment]['concentration_ub'],
                                   'concentration_units': lbub[m.compartment]['concentration_units']};
            if exceptions:        
                for k,v in exceptions.iteritems():
                    if k in m.id:
                        default_values[m.id] = {'concentration_lb':v['concentration_lb'],
                                       'concentration_ub':v['concentration_ub'],
                                       'concentration_units':v['concentration_units']};
    
        return default_values;

    def _generalize_compartment2all_concentration(self, cobra_model, concentration=None, exceptions=None):
        """
        takes a compartment and concentration for that compartment
        and updates each metabolite in that compartment

        allows for exceptions to the generalization as input

        cobra_model: a Model object

        concentration: 'compartment': {'concentration': float,
                               'concentration_var': float,
                               'concentration_lb': float,
                               'concentration_ub': float,
                               'concentration_units': string}}
            Note: the concentration is given in ln space

            Note that the concentration_var for the extimated concentration is the span (y)
            upper bound = sqrt(y)*concentration
            lower bound = 1/sqrt(y)*concentration
                where log(sqrt(y)*concentration)- log(1/sqrt(y)*concentration) = log(y)
            e.g. a span of y=1000 corresponds to approximately 3 orders of magnitude difference
                C = 50
                y = 1000
                C_ub = sqrt(y)*C = 1581.11388300841895
                C_lb = 1/sqrt(y)*C = 1.5811388300841894
            
                estimated C_var = (C-1/sqrt(y)*C)^2
            
            adapted from doi:  10.1093/bioinformatics/bts317

                               'concentration_units': 'M'}

        exceptions: metabolite.id: {'concentration': float,
                               'concentration_var': float,
                               'concentration_lb': float,
                               'concentration_ub': float,
                               'concentration_units': 'M'}

        returns a dictionary: metabolite.id {'concentration': float,
                               'concentration_var': float,
                               'concentration_lb': float,
                               'concentration_ub': float,
                               'concentration_units': 'M'}
        """
        if not(exceptions):
            exceptions = {};
            #exceptions['pi_'] = {'concentration': 1.0e-3,
            #                       'concentration_var': 0.0,
            #                       'concentration_lb':1.0e-3,
            #                       'concentration_ub':1.0e-3,
            #                       'concentration_units': 'M'};
            #exceptions['h2o_'] = {'concentration': 55.0,
            #                       'concentration_var': 0.0,
            #                       'concentration_lb':55.0,
            #                       'concentration_ub':55.0,
            #                       'concentration_units': 'M'};
            #exceptions['h2_'] = {'concentration': 0.034e-3,
            #                       'concentration_var': 0.0,
            #                       'concentration_lb':0.034e-3,
            #                       'concentration_ub':0.034e-3,
            #                       'concentration_units': 'M'};
            #exceptions['o2_'] = {'concentration': 0.055e-3,
            #                       'concentration_var': 0.0,
            #                       'concentration_lb':0.055e-3,
            #                       'concentration_ub':0.055e-3,
            #                       'concentration_units': 'M'};
            #exceptions['co2_'] = {'concentration': 1.4e-3,
            #                       'concentration_var': 0.0,
            #                       'concentration_lb':1.4e-3,
            #                       'concentration_ub':1.4e-3,
            #                       'concentration_units': 'M'};
            exceptions['pi_'] = {'concentration': 1.0e-3,
                                   'concentration_var': 1.0e-6,
                                   'concentration_lb':1.0e-3,
                                   'concentration_ub':1.0e-3,
                                   'concentration_units': 'M'};
            exceptions['h2o_'] = {'concentration': 55.0,
                                   'concentration_var': 1.0e-6,
                                   'concentration_lb':55.0,
                                   'concentration_ub':55.0,
                                   'concentration_units': 'M'};
            exceptions['h2_'] = {'concentration': 0.034e-3,
                                   'concentration_var': 1.0e-6,
                                   'concentration_lb':0.034e-3,
                                   'concentration_ub':0.034e-3,
                                   'concentration_units': 'M'};
            exceptions['o2_'] = {'concentration': 0.055e-3,
                                   'concentration_var': 1.0e-6,
                                   'concentration_lb':0.055e-3,
                                   'concentration_ub':0.055e-3,
                                   'concentration_units': 'M'};
            exceptions['co2_'] = {'concentration': 1.4e-3,
                                   'concentration_var': 1.0e-6,
                                   'concentration_lb':1.4e-3,
                                   'concentration_ub':1.4e-3,
                                   'concentration_units': 'M'};
        ## estimation implementation #1
        #if not(concentration):
        #    concentration = {};
        #    compartments = list(set(cobra_model.metabolites.list_attr('compartment')));
        #    for c in compartments:
        #        concentration[c] = {'concentration':5.0e-5,
        #                               'concentration_var':4.99e-5,
        #                               'concentration_lb':2e-8,
        #                               'concentration_ub':0.02,
        #                               'concentration_units':'M'};
        # estimation implementation #2
        # define the span
        y = 1000.0
        conc = 5.0e-5
        # create a dummy data set to estimate the geometric mean and variance
        geoconc = [sqrt(y)*conc,conc*exp(1),conc,conc/exp(1),1/sqrt(y)*conc]
        geomean = exp(average(log(geoconc)))
        geovar = exp(var(log(geoconc)))

        if not(concentration):
            concentration = {};
            compartments = list(set(cobra_model.metabolites.list_attr('compartment')));
            for c in compartments:
                #concentration[c] = {'concentration':conc,
                #                       'concentration_var':pow((conc-1/sqrt(y)*conc),2),
                #                       'concentration_lb':1/sqrt(y)*conc,
                #                       'concentration_ub':sqrt(y)*conc,
                #                       'concentration_units':'M'};
                concentration[c] = {'concentration':geomean,
                                       'concentration_var':geovar,
                                       'concentration_lb':1/sqrt(y)*conc,
                                       'concentration_ub':sqrt(y)*conc,
                                       'concentration_units':'M'};

        default_values = {};
        for m in cobra_model.metabolites:
            default_values[m.id] = {'concentration': concentration[m.compartment]['concentration'],
                                   'concentration_var': concentration[m.compartment]['concentration_var'],
                                   'concentration_lb': concentration[m.compartment]['concentration_lb'],
                                   'concentration_ub': concentration[m.compartment]['concentration_ub'],
                                   'concentration_units': concentration[m.compartment]['concentration_units']};
            if exceptions:        
                for k,v in exceptions.iteritems():
                    if k in m.id:
                        default_values[m.id] = {'concentration':v['concentration'],
                                       'concentration_var':v['concentration_var'],
                                       'concentration_lb':v['concentration_lb'],
                                       'concentration_ub':v['concentration_ub'],
                                       'concentration_units':v['concentration_units']};
    
        return default_values;

    def _generalize_externalCompartments2all_concentration(self, cobra_model, concentration=None, exceptions=None):
        """
        takes a compartment and concentration for that compartment
        and updates each metabolite in that compartment

        allows for exceptions to the generalization as input

        cobra_model: a Model object

        concentration: 'compartment': {'concentration': float,
                               'concentration_var': float,
                               'concentration_lb': float,
                               'concentration_ub': float,
                               'concentration_units': string}}
            Note: the concentration is given in ln space

            Note that the concentration_var for the extimated concentration is the span (y)
            upper bound = sqrt(y)*concentration
            lower bound = 1/sqrt(y)*concentration
                where log(sqrt(y)*concentration)- log(1/sqrt(y)*concentration) = log(y)
            e.g. a span of y=1000 corresponds to approximately 3 orders of magnitude difference
                C = 50
                y = 1000
                C_ub = sqrt(y)*C = 1581.11388300841895
                C_lb = 1/sqrt(y)*C = 1.5811388300841894
            
                estimated C_var = (C-1/sqrt(y)*C)^2
            
            adapted from doi:  10.1093/bioinformatics/bts317

                               'concentration_units': 'M'}

        exceptions: metabolite.id: {'concentration': float,
                               'concentration_var': float,
                               'concentration_lb': float,
                               'concentration_ub': float,
                               'concentration_units': 'M'}

        returns a dictionary: metabolite.id {'concentration': float,
                               'concentration_var': float,
                               'concentration_lb': float,
                               'concentration_ub': float,
                               'concentration_units': 'M'}
        """
        if not(exceptions):
            exceptions = {};
            #exceptions['pi_'] = {'concentration': 1.0e-3,
            #                       'concentration_var': 0.0,
            #                       'concentration_lb':1.0e-3,
            #                       'concentration_ub':1.0e-3,
            #                       'concentration_units': 'M'};
            #exceptions['h2o_'] = {'concentration': 55.0,
            #                       'concentration_var': 0.0,
            #                       'concentration_lb':55.0,
            #                       'concentration_ub':55.0,
            #                       'concentration_units': 'M'};
            #exceptions['h2_'] = {'concentration': 0.034e-3,
            #                       'concentration_var': 0.0,
            #                       'concentration_lb':0.034e-3,
            #                       'concentration_ub':0.034e-3,
            #                       'concentration_units': 'M'};
            #exceptions['o2_'] = {'concentration': 0.055e-3,
            #                       'concentration_var': 0.0,
            #                       'concentration_lb':0.055e-3,
            #                       'concentration_ub':0.055e-3,
            #                       'concentration_units': 'M'};
            #exceptions['co2_'] = {'concentration': 1.4e-3,
            #                       'concentration_var': 0.0,
            #                       'concentration_lb':1.4e-3,
            #                       'concentration_ub':1.4e-3,
            #                       'concentration_units': 'M'};
            exceptions['pi_'] = {'concentration': 1.0e-3,
                                   'concentration_var': 1.0e-6,
                                   'concentration_lb':1.0e-3,
                                   'concentration_ub':1.0e-3,
                                   'concentration_units': 'M'};
            exceptions['h2o_'] = {'concentration': 55.0,
                                   'concentration_var': 1.0e-6,
                                   'concentration_lb':55.0,
                                   'concentration_ub':55.0,
                                   'concentration_units': 'M'};
            exceptions['h2_'] = {'concentration': 0.034e-3,
                                   'concentration_var': 1.0e-6,
                                   'concentration_lb':0.034e-3,
                                   'concentration_ub':0.034e-3,
                                   'concentration_units': 'M'};
            exceptions['o2_'] = {'concentration': 0.055e-3,
                                   'concentration_var': 1.0e-6,
                                   'concentration_lb':0.055e-3,
                                   'concentration_ub':0.055e-3,
                                   'concentration_units': 'M'};
            exceptions['co2_'] = {'concentration': 1.4e-3,
                                   'concentration_var': 1.0e-6,
                                   'concentration_lb':1.4e-3,
                                   'concentration_ub':1.4e-3,
                                   'concentration_units': 'M'};
        ## estimation implementation #1
        #if not(concentration):
        #    concentration = {};
        #    compartments = list(set(cobra_model.metabolites.list_attr('compartment')));
        #    for c in compartments:
        #        concentration[c] = {'concentration':5.0e-5,
        #                               'concentration_var':4.99e-5,
        #                               'concentration_lb':2e-8,
        #                               'concentration_ub':0.02,
        #                               'concentration_units':'M'};
        # estimation implementation #2
        # define the span
        y = 1.0e9
        conc = 1.0e-5
        # create a dummy data set to estimate the geometric mean and variance
        geoconc = [1e-12,sqrt(y)*conc,conc*exp(1),conc,conc/exp(1),1/sqrt(y)*conc,1.0]
        geomean = exp(average(log(geoconc)))
        geovar = exp(var(log(geoconc)))

        compartments = ['e','p'];

        if not(concentration):
            concentration = {};
            #compartments = list(set(cobra_model.metabolites.list_attr('compartment')));
            for c in compartments:
                concentration[c] = {'concentration':geomean,
                                       'concentration_var':geovar,
                                       'concentration_lb':1e-12,
                                       'concentration_ub':1.0,
                                       'concentration_units':'M'};

        default_values = {};
        for m in cobra_model.metabolites:
            if m.compartment in compartments:
                default_values[m.id] = {'concentration': concentration[m.compartment]['concentration'],
                                       'concentration_var': concentration[m.compartment]['concentration_var'],
                                       'concentration_lb': concentration[m.compartment]['concentration_lb'],
                                       'concentration_ub': concentration[m.compartment]['concentration_ub'],
                                       'concentration_units': concentration[m.compartment]['concentration_units']};
                if exceptions:        
                    for k,v in exceptions.iteritems():
                        if k in m.id:
                            default_values[m.id] = {'concentration':v['concentration'],
                                           'concentration_var':v['concentration_var'],
                                           'concentration_lb':v['concentration_lb'],
                                           'concentration_ub':v['concentration_ub'],
                                           'concentration_units':v['concentration_units']};
    
        return default_values;

    def _convert_metabolomics_names(self, measured_values):
        '''Format measured metabolite names'''
        # convert names:
        measured_concentrations_converted = {};
        for met in measured_values.keys():
            met_conv = copy(met);
            met_conv = met_conv.replace('-','_DASH_');
            met_conv = met_conv.replace('(','_LPAREN_');
            met_conv = met_conv.replace(')','_RPAREN_');
            measured_concentrations_converted[met_conv] = measured_values[met];
        #self.measured_concentrations = measured_concentrations_converted;
        return measured_concentrations_converted;

    def combine_measured_concentrations(self):
        '''combine intracellular and extracellular concentration measurements'''
        self.measured_concentrations = self.measured_concentrations_intracellular;
        self.measured_concentrations.update(self.measured_concentrations_extracellular);