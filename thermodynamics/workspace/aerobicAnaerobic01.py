# Dependencies from cobra
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.flux_analysis.objective import update_objective
from cobra.manipulation.modify import convert_to_irreversible
from cobra.flux_analysis.parsimonious import optimize_minimal_flux

# Other dependencies
import csv
import json

# Dependencies from thermodynamics
from thermodynamics.thermodynamics_simulatedData import thermodynamics_simulatedData
from thermodynamics.thermodynamics_metabolomicsData import thermodynamics_metabolomicsData
from thermodynamics.thermodynamics_otherData import thermodynamics_otherData
from thermodynamics.thermodynamics_dG_f_data import thermodynamics_dG_f_data
from thermodynamics.thermodynamics_dG_r_data import thermodynamics_dG_r_data
from thermodynamics.load_model import load_thermoModel
from thermodynamics.thermodynamics_dG_p_data import thermodynamics_dG_p_data

# Read in the sbml file and define the model conditions
# Anaerobic specific changes:
cobra_model_anoxic = load_thermoModel(anoxic = True);
convert_to_irreversible(cobra_model_anoxic);
# Aerobic specific changes
cobra_model_oxic = load_thermoModel(anoxic = False);
convert_to_irreversible(cobra_model_oxic);

data_fva_anoxic = 'data\\ijo1366_fva_anoxic.json';
data_srd_anoxic = 'data\\ijo1366_srd_anoxic.json';
# make/load simulated data for anaerobic conditions
simulated_data_anoxic = thermodynamics_simulatedData();
simulated_data_anoxic.generate_sra_data(cobra_model_anoxic);
simulated_data_anoxic.generate_fva_data(cobra_model_anoxic);
simulated_data_anoxic.export_sra_data(data_srd_anoxic);
simulated_data_anoxic.export_fva_data(data_fva_anoxic);
simulated_data_anoxic.import_sra_data(data_srd_anoxic);
simulated_data_anoxic.import_fva_data(data_fva_anoxic);
simulated_data_anoxic.check_data()

data_fva_oxic = 'data\\ijo1366_fva_oxic.json';
data_srd_oxic = 'data\\ijo1366_srd_oxic.json';
# make/load simulated data for aerobic conditions
simulated_data_oxic = thermodynamics_simulatedData();
simulated_data_oxic.generate_sra_data(cobra_model_oxic);
simulated_data_oxic.generate_fva_data(cobra_model_oxic);
simulated_data_oxic.export_sra_data(data_srd_oxic);
simulated_data_oxic.export_fva_data(data_fva_oxic);
simulated_data_oxic.import_sra_data(data_srd_oxic);
simulated_data_oxic.import_fva_data(data_fva_oxic);
simulated_data_oxic.check_data()

# load pH, ionic_strength, and temperature parameters
other_data = thermodynamics_otherData();
other_data.load_defaultData();
other_data.check_data();

# generate dG_f data for all model compounds
# calculate the dG_f for each compound in each compartment
# export/load and check the dG_f data
data_dG0_transformed = 'data\\ijo1366_dG_f.json';
dG_f_data = thermodynamics_dG_f_data(id2KEGGID_filename_I='data\\id2KEGGID.csv');
#dG_f_data.make_dG0_f_pH0(); # only if the data has not been generated previously!
dG_f_data.get_transformed_dG_f('data\\compounds_dG0_f.json',cobra_model,other_data.pH,other_data.temperature,other_data.ionic_strength);
dG_f_data.export_dG_f(data_dG0_transformed);
dG_f_data.import_dG_f(data_dG0_transformed);
dG_f_data.format_dG_f();
dG_f_data.generate_estimated_dG_f(cobra_model)
dG_f_data.check_data()

# load metabolomics data for anoxic conditions
data_concentrations_anoxic = 'data\\aerobicAnaerobic01_anoxic_geo.json';
metabolomics_data_anoxic = thermodynamics_metabolomicsData();
metabolomics_data_anoxic.import_metabolomics_data(data_concentrations_anoxic);
metabolomics_data_anoxic.format_metabolomics_data();
metabolomics_data_anoxic.generate_estimated_metabolomics_data(cobra_model_anoxic);
# load metabolomics data for oxic conditions
data_concentrations_oxic = 'data\\aerobicAnaerobic01_oxic_geo.json';
metabolomics_data_oxic = thermodynamics_metabolomicsData();
metabolomics_data_oxic.import_metabolomics_data(data_concentrations_oxic);
metabolomics_data_oxic.format_metabolomics_data();
metabolomics_data_oxic.generate_estimated_metabolomics_data(cobra_model_oxic);

data_ta_anoxic = 'data\\aerobicAnaerobic01_anoxic_ta.json';
data_dG0_anoxic = 'data\\aerobicAnaerobic01_anoxic_dG0.json';
data_dG_anoxic = 'data\\aerobicAnaerobic01_anoxic_dG.json';
# calculate dG_r and perform a consistency check based on model simulations for anoxic conditions
tcc_anoxic = thermodynamics_dG_r_data();
tcc_anoxic.calculate_dG0_r(cobra_model_anoxic, dG_f_data.measured_dG_f, dG_f_data.estimated_dG_f);
tcc_anoxic.calculate_dG_r(cobra_model_anoxic,metabolomics_data_anoxic.measured_concentrations, metabolomics_data_anoxic.estimated_concentrations,
                   other_data.pH, other_data.ionic_strength, other_data.temperature)
tcc_anoxic.check_thermodynamicConsistency(cobra_model_anoxic,simulated_data_anoxic.fva_data,
                   metabolomics_data_anoxic.measured_concentrations,
                   metabolomics_data_anoxic.estimated_concentrations,
                   other_data.pH,other_data.ionic_strength,other_data.temperature)

data_ta_oxic = 'data\\aerobicAnaerobic01_oxic_ta.json';
data_dG0_oxic = 'data\\aerobicAnaerobic01_oxic_dG0.json';
data_dG_oxic = 'data\\aerobicAnaerobic01_oxic_dG.json';
# calculate dG_r and perform a consistency check based on model simulations for oxic conditions
tcc_oxic = thermodynamics_dG_r_data();
tcc_oxic.calculate_dG0_r(cobra_model_oxic, dG_f_data.measured_dG_f, dG_f_data.estimated_dG_f);
tcc_oxic.calculate_dG_r(cobra_model_oxic,metabolomics_data_oxic.measured_concentrations, metabolomics_data_oxic.estimated_concentrations,
                   other_data.pH, other_data.ionic_strength, other_data.temperature)
tcc_oxic.check_thermodynamicConsistency(cobra_model_oxic,simulated_data_oxic.fva_data,
                   metabolomics_data_oxic.measured_concentrations,
                   metabolomics_data_oxic.estimated_concentrations,
                   other_data.pH,other_data.ionic_strength,other_data.temperature)

# constrain the model solution and simulate the reaction using pFBA

# calculate the dG_r for biosynthetic pathways for anoxic conditions
tccp_anoxic = thermodynamics_dG_p_data();
###TODO
#tccp_anoxic.calculate_dG_p(cobra_model_anoxic,tcc_anoxic.dG0_r,tcc_anoxic.dG_r);
# calculate the dG_r for biosynthetic pathways for oxic conditions
tccp_oxic = thermodynamics_dG_p_data();
###TODO
#tccp_oxic.calculate_dG_p(cobra_model_oxic,tcc_oxic.dG0_r,tcc_oxic.dG_r);
# visualize the results
###TODO