# Dependencies from cobra
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.flux_analysis.objective import update_objective
from cobra.manipulation.modify import convert_toersible
from cobra.flux_analysis.parsimonious import optimize_minimal_flux

# Other dependencies
import csv,json,sys

# Dependencies from thermodynamics
from thermodynamics.thermodynamics_simulatedData import thermodynamics_simulatedData
from thermodynamics.thermodynamics_metabolomicsData import thermodynamics_metabolomicsData
from thermodynamics.thermodynamics_otherData import thermodynamics_otherData
from thermodynamics.thermodynamics_dG_f_data import thermodynamics_dG_f_data
from thermodynamics.thermodynamics_dG_r_data import thermodynamics_dG_r_data
from thermodynamics.thermodynamics_utility import load_thermoModel, simulate_thermoConstraints, add_pykA
from thermodynamics.thermodynamics_dG_p_data import thermodynamics_dG_p_data
from thermodynamics.thermodynamics_tfba import thermodynamics_tfba

# Read in the model sbml file and define the model conditions
# wild-type specific changes:
cobra_model_wt01 = create_cobra_model_from_sbml_file('data\\WtEcoli01_geo\\iteration3_140601_iDMisotopomer_WTEColi_113C80_U13C20_01_OxicWtGlc.xml');
# evolved wild-type specific changes
cobra_model_evo04 = load_thermoModel(wt01 = False);

# make/load simulated data for wild-type conditions
data_fva_wt01 = 'data\\WtEcoli01_geo\\WtEcoli01_fva_wt01.json';
data_srd_wt01 = 'data\\WtEcoli01_geo\\WtEcoli01_srd_wt01.json';
simulated_data_wt01 = thermodynamics_simulatedData();
simulated_data_wt01.generate_sra_data(cobra_model_wt01); # perform single reaction deletion analysis
simulated_data_wt01.generate_fva_data(cobra_model_wt01); # perform flux variability analysis
simulated_data_wt01.export_sra_data(data_srd_wt01); # save results for later use
simulated_data_wt01.export_fva_data(data_fva_wt01); # save results for later use
simulated_data_wt01.import_sra_data(data_srd_wt01);
simulated_data_wt01.import_fva_data(data_fva_wt01);
simulated_data_wt01.check_data();

# make/load simulated data for evolved wild-type conditions
data_fva_evo04 = 'data\\WtEcoli01_geo\\WtEcoli01_fva_evo04.json';
data_srd_evo04 = 'data\\WtEcoli01_geo\\WtEcoli01_srd_evo04.json';
simulated_data_evo04 = thermodynamics_simulatedData();
simulated_data_evo04.generate_sra_data(cobra_model_evo04); # perform single reaction deletion analysis
simulated_data_evo04.generate_fva_data(cobra_model_evo04); # perform flux variability analysis
simulated_data_evo04.export_sra_data(data_srd_evo04); # save results for later use
simulated_data_evo04.export_fva_data(data_fva_evo04); # save results for later use
simulated_data_evo04.import_sra_data(data_srd_evo04);
simulated_data_evo04.import_fva_data(data_fva_evo04);
simulated_data_evo04.check_data();

# load pH, ionic_strength, and temperature parameters
other_data = thermodynamics_otherData();
other_data.load_defaultData();
other_data.check_data();

# generate dG_f data for all model compounds
# calculate the dG_f for each compound in each compartment
# export/load and check the dG_f data
#data_dG0_transformed = 'data\\ijo1366_dG_f.json';
data_dG0_transformed = 'data\\WtEcoli01_geo\\WtEcoli01_dG_f01.json';
dG_f_data = thermodynamics_dG_f_data(id2KEGGID_filename_I='data\\id2KEGGID.csv');
##dG_f_data.make_dG0_f_pH0(); # only if the data has not been generated previously!
#dG_f_data.get_transformed_dG_f('data\\compounds_dG0_f.json',cobra_model_evo04,other_data.pH,other_data.temperature,other_data.ionic_strength); # adjust the non-transformed dG0_f data to physiological pH, temperature, and ionic strength (this step has already been completed)
#dG_f_data.export_dG_f(data_dG0_transformed); # save results for later use
dG_f_data.import_dG_f(data_dG0_transformed);
dG_f_data.format_dG_f();
dG_f_data.generate_estimated_dG_f(cobra_model_evo04)
dG_f_data.check_data()

# load metabolomics data for wt01 conditions
data_concentrations_wt01 = 'data\\WtEcoli01_geo\\ALEWt01_OxicWtGlc_0_geo.json';
metabolomics_data_wt01 = thermodynamics_metabolomicsData();
metabolomics_data_wt01.import_metabolomics_data(data_concentrations_wt01);
metabolomics_data_wt01.format_metabolomics_data(); # add compartment identifiers to metabolite ids
metabolomics_data_wt01.generate_estimated_metabolomics_data(cobra_model_wt01);
# load metabolomics data for evo04 conditions
data_concentrations_evo04 = 'data\\WtEcoli01_geo\\ALEWt01_OxicEvo04Glc_0_geo.json';
metabolomics_data_evo04 = thermodynamics_metabolomicsData();
metabolomics_data_evo04.import_metabolomics_data(data_concentrations_evo04);
metabolomics_data_evo04.format_metabolomics_data(); # add compartment identifiers to metabolite ids
metabolomics_data_evo04.generate_estimated_metabolomics_data(cobra_model_evo04);

# calculate dG_r and perform a consistency check based on model simulations for wt01 conditions
data_ta_wt01 = 'data\\WtEcoli01_geo\\WtEcoli01_wt01_ta.csv';
data_dG0_wt01 = 'data\\WtEcoli01_geo\\WtEcoli01_wt01_dG0.json';
data_dG_wt01 = 'data\\WtEcoli01_geo\\WtEcoli01_wt01_dG.json';
data_tcc_wt01 = 'data\\WtEcoli01_geo\\WtEcoli01_wt01_tcc.json';
data_dG_escher_wt01 = 'data\\WtEcoli01_geo\\WtEcoli01_wt01_dG_escher.json';
data_displacement_escher_wt01 = 'data\\WtEcoli01_geo\\WtEcoli01_wt01_displacement_escher.json';
data_concentrations_escher_wt01 = 'data\\WtEcoli01_geo\\WtEcoli01_wt01_concentrations_escher.json';
tcc_wt01 = thermodynamics_dG_r_data();
tcc_wt01.calculate_dG0_r_v2(cobra_model_wt01, dG_f_data.measured_dG_f, dG_f_data.estimated_dG_f, other_data.temperature); # calculate the change in free energy of reaction without accounting for metabolite concentrations
tcc_wt01.calculate_dG_r_v2(cobra_model_wt01,metabolomics_data_wt01.measured_concentrations, metabolomics_data_wt01.estimated_concentrations,
                   other_data.pH, other_data.ionic_strength, other_data.temperature); # adjust the change in free energy of reaction for intracellular metabolite concentrations
tcc_wt01.check_thermodynamicConsistency(cobra_model_wt01,simulated_data_wt01.fva_data,
                   metabolomics_data_wt01.measured_concentrations,
                   metabolomics_data_wt01.estimated_concentrations,
                   other_data.pH,other_data.ionic_strength,other_data.temperature); # check the thermodynamic consistency of the data
tcc_wt01.calculate_displacement(cobra_model_wt01,metabolomics_data_wt01.measured_concentrations, metabolomics_data_wt01.estimated_concentrations); # calculate the displacements from equillibrium
tcc_wt01.export_dG0_r_json(data_dG0_wt01); # save for later use
tcc_wt01.export_dG_r_json(data_dG_wt01); # save for later use
tcc_wt01.export_tcc_json(data_ta_wt01); # save for later use
tcc_wt01.export_summary(cobra_model_wt01,simulated_data_wt01.fva_data,data_ta_wt01); # write summary of the analysis to csv file
tcc_wt01.export_dG_r_escher(data_dG_escher_wt01); # save for later use
tcc_wt01.export_displacement_escher(data_displacement_escher_wt01); # save for later use
tcc_wt01.export_concentrations_escher(data_concentrations_escher_wt01,metabolomics_data_wt01.measured_concentrations); # save for later use
#tcc_wt01.import_dG0_r_json(data_dG0_wt01); 
#tcc_wt01.import_dG_r_json(data_dG_wt01);
#tcc_wt01.import_tcc_json(data_ta_wt01);

# calculate dG_r and perform a consistency check based on model simulations for evo04 conditions
data_ta_evo04 = 'data\\WtEcoli01_geo\\WtEcoli01_evo04_ta.csv';
data_dG0_evo04 = 'data\\WtEcoli01_geo\\WtEcoli01_evo04_dG0.json';
data_dG_evo04 = 'data\\WtEcoli01_geo\\WtEcoli01_evo04_dG.json';
data_tcc_evo04 = 'data\\WtEcoli01_geo\\WtEcoli01_evo04_tcc.json';
data_dG_escher_evo04 = 'data\\WtEcoli01_geo\\WtEcoli01_evo04_dG_escher.json';
data_displacement_escher_evo04 = 'data\\WtEcoli01_geo\\WtEcoli01_evo04_displacement_escher.json';
data_concentrations_escher_evo04 = 'data\\WtEcoli01_geo\\WtEcoli01_evo04_concentrations_escher.json';
tcc_evo04 = thermodynamics_dG_r_data();
tcc_evo04.calculate_dG0_r_v2(cobra_model_evo04, dG_f_data.measured_dG_f, dG_f_data.estimated_dG_f, other_data.temperature); # calculate the change in free energy of reaction without accounting for metabolite concentrations
tcc_evo04.calculate_dG_r_v2(cobra_model_evo04,metabolomics_data_evo04.measured_concentrations, metabolomics_data_evo04.estimated_concentrations,
                   other_data.pH, other_data.ionic_strength, other_data.temperature); # adjust the change in free energy of reaction for intracellular metabolite concentrations
tcc_evo04.check_thermodynamicConsistency(cobra_model_evo04,simulated_data_evo04.fva_data,
                   metabolomics_data_evo04.measured_concentrations,
                   metabolomics_data_evo04.estimated_concentrations,
                   other_data.pH,other_data.ionic_strength,other_data.temperature); # check the thermodynamic consistency of the data
tcc_evo04.calculate_displacement(cobra_model_evo04,metabolomics_data_evo04.measured_concentrations, metabolomics_data_evo04.estimated_concentrations); # calculate the displacements from equillibrium
tcc_evo04.export_dG0_r_json(data_dG0_evo04); # save for later use
tcc_evo04.export_dG_r_json(data_dG_evo04); # save for later use
tcc_evo04.export_tcc_json(data_tcc_evo04); # save for later use
tcc_evo04.export_summary(cobra_model_evo04,simulated_data_evo04.fva_data,data_ta_evo04); # write summary of the analysis to csv file
tcc_evo04.export_dG_r_escher(data_dG_escher_evo04); # save for later use
tcc_evo04.export_displacement_escher(data_displacement_escher_evo04); # save for later use
tcc_evo04.export_concentrations_escher(data_concentrations_escher_evo04,metabolomics_data_evo04.measured_concentrations); # save for later use
#tcc_evo04.import_dG0_r_json(data_dG0_evo04);
#tcc_evo04.import_dG_r_json(data_dG_evo04);
#tcc_evo04.import_tcc_json(data_tcc_evo04);

# inspect the thermodynamic analysis results

## constrain the model solution and simulate optimal growth
#gr_analysis_wt01 = simulate_thermoConstraints(cobra_model_wt01,['PGCD','ACACT1r','NDPK2']);
#gr_analysis_evo04 = simulate_thermoConstraints(cobra_model_evo04,['PGCD','ACACT1r']);

## calculate the dG for biosynthetic pathways for wt01 conditions
#tccp_wt01 = thermodynamics_dG_p_data();
#tccp_wt01.calculate_dG_p(cobra_model_wt01,tcc_wt01.dG0_r,tcc_wt01.dG_r);
## calculate the dG for biosynthetic pathways for evo04 conditions
#tccp_evo04 = thermodynamics_dG_p_data();
#tccp_evo04.calculate_dG_p(cobra_model_evo04,tcc_evo04.dG0_r,tcc_evo04.dG_r);

## expand the reaction set of the wt01 model to reflect the enzyme permiscuity of pykA
#add_pykA(cobra_model_wt01);
#gr_analysis_wt01 = simulate_thermoConstraints(cobra_model_wt01,['PGCD','ACACT1r','NDPK2']);

# visualize the results
###TODO

#tfba test
tfba = thermodynamics_tfba()
#tfba.tfba_looplaw(cobra_model_wt01); # constraints need some work
tfba.tsampling(cobra_model_wt01,tcc_wt01.dG_r, tcc_wt01.metabolomics_coverage, tcc_wt01.dG_r_coverage, tcc_wt01.thermodynamic_consistency_check, 0.5, 0.99,use_measured_dG_r=True);
tfba.tsampling_conc_ln(cobra_model_wt01, metabolomics_data_wt01.measured_concentrations, metabolomics_data_wt01.estimated_concentrations, tcc_wt01.dG0_r,
                  other_data.temperature, tcc_wt01.dG_r_coverage, 0.99, use_measured_concentrations=False,use_measured_dG0_r=True, solver='gurobi');
#tfba.tfba(cobra_model_wt01,tcc_wt01.dG_r, tcc_wt01.metabolomics_coverage, tcc_wt01.dG_r_coverage, tcc_wt01.thermodynamic_consistency_check, 0.5, 0.99,use_measured_dG_r=True);
#tfba.tfba_conc_ln(cobra_model_wt01, metabolomics_data_wt01.measured_concentrations, metabolomics_data_wt01.estimated_concentrations, tcc_wt01.dG0_r,
#                  other_data.temperature, tcc_wt01.dG_r_coverage, 0.99, use_measured_concentrations=False,use_measured_dG0_r=True, solver='gurobi');
#tfba.tfva(cobra_model_wt01, tcc_wt01.dG_r, tcc_wt01.metabolomics_coverage, tcc_wt01.dG_r_coverage, tcc_wt01.thermodynamic_consistency_check, 0.5, 0.99, use_measured_dG_r=True, solver='gurobi');
#tfba.tfva_dG_r(cobra_model_wt01, tcc_wt01.dG_r, tcc_wt01.metabolomics_coverage, tcc_wt01.dG_r_coverage, tcc_wt01.thermodynamic_consistency_check, 0.5, 0.99, use_measured_dG_r=True, solver='gurobi');
tfba.tfva_concentrations(cobra_model_wt01, metabolomics_data_wt01.measured_concentrations, metabolomics_data_wt01.estimated_concentrations,
                  tcc_wt01.dG0_r, other_data.temperature, tcc_wt01.dG_r_coverage, 0.99, use_measured_concentrations=False,use_measured_dG0_r=True,solver='gurobi');

rxn_ids = tfba.tfva_dG_r_data.keys();
print ('rxn_id\tfva_min\tfva_max\ttfva_min\ttfvamax')
for rxn in rxn_ids:
    print ('%s\t%s\t%s\t%s\t%s' %(rxn,simulated_data_wt01.fva_data[rxn]['flux_lb'],simulated_data_wt01.fva_data[rxn]['flux_ub'],fva_results[rxn]['flux_lb'],fva_results[rxn]['flux_ub']))
rxn_ids = fva_results.keys();
print ('rxn_id\tdG_r_min\tdG_r_max')
for rxn in rxn_ids:
    print ('%s\t%s\t%s' %(rxn,tfba.tfva_dG_r_data[rxn]['dG_r_lb'],tfba.tfva_dG_r_data[rxn]['dG_r_ub']))
met_ids = tfba.tfva_concentration_data.keys();
print ('met_id\tconcentration_min\tconcentration_max')
for met in met_ids:
    print ('%s\t%s\t%s' %(met,tfba.tfva_concentration_data[met]['concentration_lb'],tfba.tfva_concentration_data[met]['concentration_ub']))