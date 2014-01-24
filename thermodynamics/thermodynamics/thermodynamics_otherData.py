
import json, csv
from math import sqrt,exp,pow
from numpy import average, var, log

from thermodynamics_io import thermodynamics_io

class thermodynamics_otherData(thermodynamics_io):
    """Class to handle other input data"""

    def __init__(self):
        
        self.pH = {}
        self.temperature = {}
        self.ionic_strength = {}

    def load_data(self,pH_I,ionic_strength_I,temperature_I):
        '''load pH, ionic_strength, and temperature'''
        self.pH = pH_I
        self.temperature = ionic_strength_I
        self.ionic_strength = temperature_I

    def check_data(self):
        '''check data integrity'''
        return

    def clear_data(self):
        '''clear data'''
        self.pH = {}
        self.temperature = {}
        self.ionic_strength = {}

    def load_defaultData(self):
        '''load default pH, ionic_strength, and temperature'''
        # define the pH, ionic strength, and temperature
        pH_I = {};
        pH_I['c'] = {'pH':7.5};
        pH_I['p'] = {'pH':7.0};
        pH_I['e'] = {'pH':7.0};

        ionic_strength_I = {};
        ionic_strength_I['c'] = {'ionic_strength': .2,'ionic_strength_units': 'M'}
        ionic_strength_I['p'] = {'ionic_strength': .2,'ionic_strength_units': 'M'}
        ionic_strength_I['e'] = {'ionic_strength': .2,'ionic_strength_units': 'M'}

        temperature_I = {};
        temperature_I['c'] = {'temperature': 310.15,'temperature_units': 'K'}
        temperature_I['p'] = {'temperature': 310.15,'temperature_units': 'K'}
        temperature_I['e'] = {'temperature': 310.15,'temperature_units': 'K'}

        self.pH = pH_I
        self.temperature = temperature_I
        self.ionic_strength = ionic_strength_I



