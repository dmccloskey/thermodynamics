# Dependencies from cobra
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.flux_analysis.objective import update_objective

# Other dependencies
import csv
import json

# Dependencies from thermodynamics
from thermodynamics.thermodynamics_simulatedData import thermodynamics_simulatedData
from thermodynamics.thermodynamics_metabolomicsData import thermodynamics_metabolomicsData
from thermodynamics.thermodynamics_otherData import thermodynamics_otherData
from thermodynamics.thermodynamics_dG_f_data import thermodynamics_dG_f_data
from thermodynamics.thermodynamics_dG_r_data import thermodynamics_dG_r_data

# Read in the sbml file and define the model conditions
ijo1366_sbml = "data\\iJO1366.xml"
cobra_model = create_cobra_model_from_sbml_file(ijo1366_sbml, print_time=True)
# Change the objective
update_objective(cobra_model,{'Ec_biomass_iJO1366_WT_53p95M':1.0})
# Change uptake reactions for growth on glycerol
cobra_model.reactions.get_by_id('EX_glc_LPAREN_e_RPAREN_').lower_bound = -10.0;

data_fva = 'data\\ijo1366_fva_glc.json';
data_srd = 'data\\ijo1366_srd_glc.json';
# make/load simulated data
#simulated_data = thermodynamics_simulatedData();
#simulated_data.generate_sra_data(cobra_model);
#simulated_data.generate_fva_data(cobra_model);
#simulated_data.export_sra_data(data_srd);
simulated_data.export_fva_data(data_fva);
simulated_data.import_sra_data(data_srd);
simulated_data.import_fva_data(data_fva);
simulated_data.check_data()

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

data_concentrations = 'data\\ALEWt01_OxicEvo04Glc_0_geo.json';
metabolomics_data = thermodynamics_metabolomicsData();
metabolomics_data.import_metabolomics_data(data_concentrations);
metabolomics_data.format_metabolomics_data();
metabolomics_data.generate_estimated_metabolomics_data(cobra_model);

data_ta = 'data\\ALEWt01_OxicEvo04Glc_ta.json';
data_dG0 = 'data\\ALEWt01_OxicEvo04Glc_dG0.json';
data_dG = 'data\\ALEWt01_OxicEvo04Glc_dG.json';
# calculate dG_r and perform a consistency check based on model simulations
tcc = thermodynamics_dG_r_data();
tcc.calculate_dG0_r(cobra_model, dG_f_data.measured_dG_f, dG_f_data.estimated_dG_f);
tcc.calculate_dG_r(cobra_model,metabolomics_data.measured_concentrations, metabolomics_data.estimated_concentrations,
                   other_data.pH, other_data.ionic_strength, other_data.temperature)
tcc.check_thermodynamicConsistency(cobra_model,simulated_data.fva_data,
                   metabolomics_data.measured_concentrations,
                   metabolomics_data.estimated_concentrations,
                   other_data.pH,other_data.ionic_strength,other_data.temperature)
print ('thermodynamics complete')

# data files:
#data_concentrations = 'data\\ALEWt01_OxicWtGlc_0_geo.json';
#data_ta = 'data\\ALEWt01_OxicWtGlc_ta.json';
#data_dG0 = 'data\\ALEWt01_OxicWtGlc_dG0.json';
#data_dG = 'data\\ALEWt01_OxicWtGlc_dG.json';
#data_concentrations = 'data\\ALEWt01_OxicEvo03Glc_0_geo.json';
#data_ta = 'data\\ALEWt01_OxicEvo03Glc_ta.json';
#data_dG0 = 'data\\ALEWt01_OxicEvo03Glc_dG0.json';
#data_dG = 'data\\ALEWt01_OxicEvo03Glc_dG.json';
#data_concentrations = 'data\\ALEWt01_OxicEvo0Glc_0_geo.json';
#data_ta = 'data\\ALEWt01_OxicEvo04Glc_ta.json';
#data_dG0 = 'data\\ALEWt01_OxicEvo04Glc_dG0.json';
#data_dG = 'data\\ALEWt01_OxicEvo04Glc_dG.json';
#data_concentrations = 'data\\ALEWt01_OxicEvo08Glc_0_geo.json';
#data_ta = 'data\\ALEWt01_OxicEvo08Glc_ta.json';
#data_dG0 = 'data\\ALEWt01_OxicEvo08Glc_dG0.json';
#data_dG = 'data\\ALEWt01_OxicEvo08Glc_dG.json';
#data_concentrations = 'data\\ALEWt01_OxicEvo09Glc_0_geo.json';
#data_ta = 'data\\ALEWt01_OxicEvo09Glc_ta.json';
#data_dG0 = 'data\\ALEWt01_OxicEvo09Glc_dG0.json';
#data_dG = 'data\\ALEWt01_OxicEvo09Glc_dG.json';