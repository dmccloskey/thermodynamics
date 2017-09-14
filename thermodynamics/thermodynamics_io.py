# -*- coding: utf-8 -*-
"""helper functions for generating thermodynamic analysis input"""

import json, csv
from math import sqrt,exp,pow

class thermodynamics_io:
    """base class for the input and output of data"""

    def import_values_json(self, filename):
        """import values from a json file"""
        data = json.load(open(filename))
        return data;

    def export_values_json(self, filename, data):
        """export values to a json file"""
        
        with open(filename, 'w') as outfile:
            json.dump(data, outfile, indent=4);

def checkInput_concentrations(measured_values_I):
    """check concentration data input

    Args:
        measured_values (dict): measured values with variances
            {metabolite.id: {'concentration': float,
            'concentration_var': float,
            'concentration_units': 'M'}
    Returns:
        dict: output: measured_values_O: {metabolite.id: {'concentration': float,
            'concentration_var': float,
            'concentration_units': 'M'}
    """
    # check units
    measured_values_O = {};
    for k,v in measured_values_I.items():
        if v['concentration_units'] == 'M':
            measured_values_O[k] = v
        else:
            print((str(k) + ' has invalid units of ' + str(v) + ' and will be ignored'))
    return measured_values_O;

def checkInput_dG_f(measured_values_I):
    """check dG_O_f data input

    Args:
        measured_values (dict): measured values with variances
            {metabolite.id: {'dG_f': float,
            'dG_f_var': float,
            'dG_f_units': 'kJ/mol'}

    Returns:
        dict: output: measured_values_O: {metabolite.id: {'dG_f': float,
            'dG_f_var': float,
            'dG_f_lb': float,
            'dG_f_ub': float,
            'dG_f_units': 'kJ/mol'}
    """
    # check units
    
    measured_values_O = {};
    for k,v in measured_values_I.items():
        if v['dG_f_units'] == 'kJ/mol':
            measured_values_O[k] = v
        else:
            print((str(k) + ' has invalid units of ' + str(v) + ' and will be ignored'))
    return measured_values_O;

def import_values_json(filename):
    """import values from a json file"""
    data = json.load(open(filename))
    return data;

