"""helper functions for generating thermodynamic analysis input"""

import json, csv
from math import sqrt,exp,pow
from numpy import average, var, log

def checkInput_concentrations(measured_values_I):
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

def checkInput_dG_f(measured_values_I):
    """
    check dG_O_f data input

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
    # check units
    
    measured_values_O = {};
    for k,v in measured_values_I.iteritems():
        if v['dG_f_units'] == 'kJ/mol':
            measured_values_O[k] = v
        else:
            print (str(k) + ' has invalid units of ' + str(v) + ' and will be ignored')
    return measured_values_O;

def import_values_json(filename):
    '''import values from a json file'''
    data = json.load(open(filename))
    return data;

