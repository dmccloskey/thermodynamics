import pytest
import numpy as np

# Dependencies from cobra
from cobra.io.json import load_json_model
from cobra.manipulation.modify import convert_to_irreversible

# Other dependencies
import csv,json,sys

# Dependencies from thermodynamics
# from thermodynamics.thermodynamics_simulatedData import thermodynamics_simulatedData
from cobra_utilities.cobra_simulatedData import cobra_simulatedData
from thermodynamics.thermodynamics_metabolomicsData import thermodynamics_metabolomicsData
from thermodynamics.thermodynamics_otherData import thermodynamics_otherData
from thermodynamics.thermodynamics_dG_f_data import thermodynamics_dG_f_data
from thermodynamics.thermodynamics_dG_r_data import thermodynamics_dG_r_data
from thermodynamics.thermodynamics_utility import load_thermoModel, simulate_thermoConstraints, add_pykA
from thermodynamics.thermodynamics_dG_p_data import thermodynamics_dG_p_data
from thermodynamics.thermodynamics_tfba import thermodynamics_tfba

from . import data_dir

class test_thermodynamics():
    
    def __init__(self):
        self.cobra_model = None
        self.simulated_data = None
        self.other_data = None
        self.metabolomics_data = None
        self.tcc = None
        self.dG_f_data = None
        self.tccp = None

    def test_cobra_model(self):
        cobra_model = load_json_model(data_dir + '/mini.json')
        convert_to_irreversible(cobra_model)
        solution = cobra_model.optimize()
        assert(solution.objective_value == 30.0)
        assert(solution.fluxes['ENO'] == 20.0)
        self.cobra_model = cobra_model

    def test_simulatedData(self):
        data_fva = data_dir + '/test_fva.json'
        data_srd = data_dir + '/test_srd.json'
        simulated_data = cobra_simulatedData()
        simulated_data.generate_sra_data(self.cobra_model) # perform single reaction deletion analysis
        simulated_data.generate_fva_data(self.cobra_model) # perform flux variability analysis
        assert(simulated_data.fva_data['ENO']['maximum'] == 1000.0)
        assert(simulated_data.fva_data['ENO']['minimum'] == 18.0)
        assert(np.isnan(simulated_data.sra_data['ENO']['gr']))
        assert(np.isnan(simulated_data.sra_data['ENO']['gr_ratio']))
        assert(simulated_data.sra_data['H2Ot']['gr'] == 30.0)
        assert(simulated_data.sra_data['H2Ot']['gr_ratio'] == 1.0)
        test_fva_data = simulated_data._convert_fluxBounds2var(simulated_data.fva_data)
        assert(test_fva_data['ENO']['flux'] == 491.0)
        assert(test_fva_data['ENO']['flux_var'] == 259081.0)
        assert(test_fva_data['ENO']['flux_units'] == 'mmol*gDW-1*hr-1')
        assert(test_fva_data['ENO']['flux_ub'] == 1000.0)
        assert(test_fva_data['ENO']['flux_lb'] == 18.0)
        simulated_data.export_sra_data(data_srd) # save results for later use
        simulated_data.export_fva_data(data_fva) # save results for later use
        simulated_data.import_sra_data(data_srd)
        simulated_data.import_fva_data(data_fva)
        assert(simulated_data.fva_data['ENO']['flux_ub'] == 1000.0)
        assert(simulated_data.fva_data['ENO']['flux_lb'] == 18.0)
        assert(np.isnan(simulated_data.sra_data['ENO']['gr']))
        assert(np.isnan(simulated_data.sra_data['ENO']['gr_ratio']))
        assert(simulated_data.sra_data['H2Ot']['gr'] == 30.0)
        assert(simulated_data.sra_data['H2Ot']['gr_ratio'] == 1.0)
        simulated_data.check_data()
        self.simulated_data = simulated_data

    def test_otherData(self):
        # load pH, ionic_strength, and temperature parameters
        other_data = thermodynamics_otherData()
        other_data.load_defaultData()
        other_data.check_data()
        assert(other_data.pH['c']['pH'] == 7.5)
        assert(other_data.pH['e']['pH'] == 7.0)
        assert(other_data.pH['p']['pH'] == 7.0)
        assert(other_data.ionic_strength['c']['ionic_strength'] == 0.2)
        assert(other_data.ionic_strength['c']['ionic_strength_units'] == 'M')
        assert(other_data.ionic_strength['e']['ionic_strength'] == 0.2)
        assert(other_data.ionic_strength['e']['ionic_strength_units'] == 'M')
        assert(other_data.ionic_strength['p']['ionic_strength'] == 0.2)
        assert(other_data.ionic_strength['p']['ionic_strength_units'] == 'M')        
        assert(other_data.temperature['c']['temperature'] == 310.15)
        assert(other_data.temperature['c']['temperature_units'] == 'K')
        assert(other_data.temperature['e']['temperature'] == 310.15)
        assert(other_data.temperature['e']['temperature_units'] == 'K')
        assert(other_data.temperature['p']['temperature'] == 310.15)
        assert(other_data.temperature['p']['temperature_units'] == 'K')
        self.other_data = None

    def test_dG_f_data(self):
        # generate dG_f data for all model compounds
        # calculate the dG_f for each compound in each compartment
        # export/load and check the dG_f data
        data_dG0_transformed = data_dir + '/test_dG_f01.json'
        dG_f_data = thermodynamics_dG_f_data(id2KEGGID_filename_I=data_dir + '/id2KEGGID.csv')
        dG_f_data.get_transformed_dG_f(data_dir + '/compounds_dG0_f.json',
            self.cobra_model,self.other_data.pH,
            self.other_data.temperature,self.other_data.ionic_strength
            ) # adjust the non-transformed dG0_f data to physiological pH, temperature, and ionic strength
        dG_f_data.export_dG_f(data_dG0_transformed) # save results for later use
        dG_f_data.import_dG_f(data_dG0_transformed)
        dG_f_data.format_dG_f()
        dG_f_data.generate_estimated_dG_f(self.cobra_model)
        dG_f_data.check_data()
        self.dG_f_data = dG_f_data

    def test_metabolomicsData(self):
        # load metabolomics data
        data_concentrations = data_dir + '/test_geo01.json'
        metabolomics_data = thermodynamics_metabolomicsData()
        metabolomics_data.import_metabolomics_data(data_concentrations)
        metabolomics_data.format_metabolomics_data() # add compartment identifiers to metabolite ids
        metabolomics_data.generate_estimated_metabolomics_data(self.cobra_model)
        self.metabolomics_data = metabolomics_data

    def test_dG_r_data(self):
        # calculate dG_r and perform a consistency check based on model simulations
        data_ta = data_dir + '/test_ta.csv'
        data_dG0 = data_dir + '/test_dG0.json'
        data_dG = data_dir + '/test_dG.json'
        data_tcc = data_dir + '/test_tcc.json'
        tcc = thermodynamics_dG_r_data()
        # tcc.import_dG0_r_json(data_dG0)
        # tcc.import_dG_r_json(data_dG)
        # tcc.import_tcc_json(data_tcc)
        tcc.calculate_dG0_r(self.cobra_model, dG_f_data.measured_dG_f, dG_f_data.estimated_dG_f, other_data.temperature) # calculate the change in free energy of reaction without accounting for metabolite concentrations
        tcc.calculate_dG_r(self.cobra_model,metabolomics_data.measured_concentrations, metabolomics_data.estimated_concentrations,
                        other_data.pH, other_data.ionic_strength, other_data.temperature) # adjust the change in free energy of reaction for intracellular metabolite concentrations
        tcc.check_thermodynamicConsistency(self.cobra_model,simulated_data.fva_data,
                        metabolomics_data.measured_concentrations,
                        metabolomics_data.estimated_concentrations,
                        other_data.pH,other_data.ionic_strength,other_data.temperature) # check the thermodynamic consistency of the data
        tcc.export_dG0_r_json(data_dG0) # save for later use
        tcc.export_dG_r_json(data_dG) # save for later use
        tcc.export_tcc_json(data_tcc) # save for later use
        tcc.export_summary(self.cobra_model,simulated_data.fva_data,data_ta) # write summary of the analysis to csv file
        self.tcc = tcc

    def test_dG_p_data(self):        
        # calculate the dG for biosynthetic pathways
        tccp = thermodynamics_dG_p_data()
        tccp.calculate_dG_p(self.cobra_model,self.tcc.dG0_r,self.tcc.dG_r)
        self.tccp = tccp 