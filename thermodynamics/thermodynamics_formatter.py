# -*- coding: utf-8 -*-
"""helper functions for generating thermodynamic analysis input"""

import json, csv
from math import sqrt,exp,pow
from numpy import average, var, log

def convert_fluxBounds2var(flux_bounds):
    """
    convert flux bounds from FVA to median and variance
    
    variance = (max - median)^2

    flux_bounds: {reaction.id: {'maximum': float, 'minimum': float}}

    returns a dictionary: {reaction.id: {'flux': float, 'flux_var': float, 'flux_units': 'mmol*gDW-1*hr-1'}

    """

    flux_bounds_O = {};
    for k,v in flux_bounds.items():
        median = (v['maximum'] - v['minimum'])/2;
        variance = (v['maximum'] - median)*(v['maximum'] - median);
        flux_bounds_O[k] = {'flux': median, 'flux_var': variance, 'flux_units': 'mmol*gDW-1*hr-1',
                            'flux_lb': v['minimum'], 'flux_ub': v['maximum']};

    return flux_bounds_O

def convert_cv2varAndmM2M_concentrations(measured_values):
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
    for k,v in measured_values.items():
         concM = v['concentration']*1e-3;
         concMvar = v['concentration_cv']/100*(v['concentration']*1e-3)*v['concentration_cv']/100*(v['concentration']*1e-3);
         measured_values_O[k + '_c'] = {'concentration': concM,
                                             'concentration_var': concMvar,
                                             'concentration_units': 'M'}
    return measured_values_O

def compartementalize_concentrations(measured_values):
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
    for k,v in measured_values.items():
         measured_values_O[k + '_c'] = v

    return measured_values_O

def convert_cv2lbubAndmM2M_concentrations(measured_values,min_value):
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
    for k,v in measured_values.items():
         concMlb = 0.0;
         concMub = 0.0;
         concMlb = v['concentration']*1e-3 - v['concentration_cv']/100*(v['concentration']*1e-3);
         concMub = v['concentration']*1e-3 + v['concentration_cv']/100*(v['concentration']*1e-3);
         if concMlb<0: concMlb = min_value;
         measured_values_O[k + '_c'] = {'concentration_lb': concMlb,
                                             'concentration_ub': concMub,
                                             'concentration_units': 'M'}
    return measured_values_O

def convert_var2lbub_dG_f(measured_values):
    """
    convert measured dG_f values with a stadard deviation
    currently use +/- SD
    lb = ave - sqrt(var)
    ub = ave + sqrt(var)

    measured_values: measured values with variances
                             {metabolite.id: {'dG_f': float,
                                             'dG_f_var': float,
                                             'dG_f_units': 'kJ/mol'}
    returns a dictionary: measured_values_O:
                             {metabolite.id: {'dG_f': float,
                                             'dG_f_var': float,
                                             'dG_f_lb': float,
                                             'dG_f_ub': float,
                                             'dG_f_units': 'kJ/mol'}
    """
    measured_values_O = {};
    for k,v in measured_values.items():
         concMlb = 0.0;
         concMub = 0.0;
         if v['dG_f_var']:
             concMlb = v['dG_f'] - sqrt(v['dG_f_var']);
             concMub = v['dG_f'] + sqrt(v['dG_f_var']);
         else:
             concMlb = v['dG_f'];
             concMub = v['dG_f'];
         #if concMlb<0: concMlb = min_value;
         measured_values_O[k] = {'dG_f': v['dG_f'],
                                 'dG_f_var': v['dG_f_var'],
                                 'dG_f_lb': concMlb,
                                             'dG_f_ub': concMub,
                                             'dG_f_units': 'kJ/mol'}
    return measured_values_O

def convert_std2lbub_dG_f(measured_values,min_value):
    """
    convert measured values with a standard deviation
    (CV = SD/Ave*100; SD = CV/100*AVE) to lb and ub
    currently use +/- SD
    lb = ave - sqrt(var) NOTE: if < 0 min_value is used instead
    ub = ave + sqrt(var)

    measured_values: measured values with variances
                             {metabolite.id: {'dG_f': float,
                                             'dG_f_cv': float,
                                             'dG_f_units': 'mM'}
    returns a dictionary: measured_values_O:
                             {metabolite.id: {'dG_f_lb': float,
                                             'dG_f_ub': float,
                                             'dG_f_units': 'M'}
    """
    measured_values_O = {};
    for k,v in measured_values.items():
         concMlb = 0.0;
         concMub = 0.0;
         concMlb = v['dG_f'] - v['dG_f_var'];
         concMub = v['dG_f'] + v['dG_f_var'];
         if concMlb<0: concMlb = min_value;
         measured_values_O[k] = {'dG_f_lb': concMlb,
                                             'dG_f_ub': concMub,
                                             'dG_f_units': 'M'}
    return measured_values_O

def generalize_compartmentLBUB2all_concentration(cobra_model, lbub=None, exceptions=None):
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
            for k,v in exceptions.items():
                if k in m.id:
                    default_values[m.id] = {'concentration_lb':v['concentration_lb'],
                                   'concentration_ub':v['concentration_ub'],
                                   'concentration_units':v['concentration_units']};
    
    return default_values;

def generalize_compartment2all_concentration(cobra_model, concentration=None, exceptions=None):
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
            for k,v in exceptions.items():
                if k in m.id:
                    default_values[m.id] = {'concentration':v['concentration'],
                                   'concentration_var':v['concentration_var'],
                                   'concentration_lb':v['concentration_lb'],
                                   'concentration_ub':v['concentration_ub'],
                                   'concentration_units':v['concentration_units']};
    
    return default_values;

def generalize_compartmentLBUB2all_dG_f(cobra_model, lbub=None, exceptions=None):
    """
    takes a compartment and lb/ub for that compartment
    and updates each metabolite in that compartment

    allows for exceptions to the generalization as input

    cobra_model: a Model object

    lbub: metabolite.compartment: {'dG_f_lb': float,
                           'dG_f_ub': float,
                           'dG_f_units': 'kJ/mol'}

    returns a dictionary: metabolite.id {'dG_f_lb': float,
                           'dG_f_ub': float,
                           'dG_f_units': string}
    """
    if not(lbub):
        lbub = {};
        compartments = list(set(cobra_model.metabolites.list_attr('compartment')));
        for c in compartments:
            lbub[c] = {'dG_f_lb':0.0,
                                   'dG_f_ub':0.0,
                                   'dG_f_units':'kJ/mol'};
    default_values = {};
    for m in cobra_model.metabolites:
        default_values[m.id] = {'dG_f_lb': lbub[m.compartment]['dG_f_lb'],
                               'dG_f_ub': lbub[m.compartment]['dG_f_ub'],
                               'dG_f_units': lbub[m.compartment]['dG_f_units']};
        if exceptions:    
            for k,v in exceptions.items():
                if k in m.id:
                    default_values[m.id] = {'dG_f_lb':v['dG_f_lb'],
                                   'dG_f_ub':v['dG_f_ub'],
                                   'dG_f_units':v['dG_f_units']};
    
    return default_values;

def generalize_compartment2all_dG_f(cobra_model, dG_f=None, exceptions=None):
    """
    takes a compartment and lb/ub for that compartment
    and updates each metabolite in that compartment

    allows for exceptions to the generalization as input

    the lower bounds and upper bounds are estimated as
    lb = mean - sqrt(var)
    ub = mean + sqrt(var)

    cobra_model: a Model object

    lbub: metabolite.compartment: {'dG_f': float,
                           'dG_f_var': float,
                           'dG_f_units': 'kJ/mol'}

    returns a dictionary: metabolite.id {'dG_f': float,
                           'dG_f_var': float,
                           'dG_f_lb': float,
                           'dG_f_ub': float,
                           'dG_f_units': 'kJ/mol'}
    """
    if not(dG_f):
        dG_f = {};
        compartments = list(set(cobra_model.metabolites.list_attr('compartment')));
        for c in compartments:
            dG_f[c] = {'dG_f':0.0,
                                   'dG_f_var':1e12, # based on the convention described in
                                                   # doi:10.1371/journal.pcbi.1003098
                                   'dG_f_units':'kJ/mol'};
    default_values = {};
    for m in cobra_model.metabolites:
        default_values[m.id] = {'dG_f': dG_f[m.compartment]['dG_f'],
                               'dG_f_var': dG_f[m.compartment]['dG_f_var'],
                               'dG_f_lb': dG_f[m.compartment]['dG_f'] - sqrt(dG_f[m.compartment]['dG_f_var']),
                               'dG_f_ub': dG_f[m.compartment]['dG_f'] + sqrt(dG_f[m.compartment]['dG_f_var']),
                               'dG_f_units': dG_f[m.compartment]['dG_f_units']};
        if exceptions:    
            for k,v in exceptions.items():
                if k in m.id:
                    default_values[m.id] = {'dG_f':v['dG_f'],
                                   'dG_f_var':v['dG_f_var'],
                                   'dG_f_lb':v['dG_f_lb'],
                                   'dG_f_ub':v['dG_f_ub'],
                                   'dG_f_units':v['dG_f_units']};
    
    return default_values;