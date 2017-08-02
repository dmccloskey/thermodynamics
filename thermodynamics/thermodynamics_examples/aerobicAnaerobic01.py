# Dependencies from cobra
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
from cobra.flux_analysis.parsimonious import optimize_minimal_flux

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

def _main_():
    ##PART 1: Working
    #-------
    # Read in the model sbml file and define the model conditions
    # Anaerobic specific changes:
    cobra_model_anoxic = load_thermoModel(anoxic = True)
    #convert_to_irreversible(cobra_model_anoxic)
    # Aerobic specific changes
    cobra_model_oxic = load_thermoModel(anoxic = False)
    #convert_to_irreversible(cobra_model_oxic)

    ##PART 2: Working
    #-------
    # make/load simulated data for anaerobic conditions
    data_fva_anoxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_fva_anoxic.json'
    data_srd_anoxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_srd_anoxic.json'
    simulated_data_anoxic = cobra_simulatedData()
    # simulated_data_anoxic.generate_sra_data(cobra_model_anoxic) # perform single reaction deletion analysis
    # simulated_data_anoxic.generate_fva_data(cobra_model_anoxic) # perform flux variability analysis
    # simulated_data_anoxic.export_sra_data(data_srd_anoxic) # save results for later use
    # simulated_data_anoxic.export_fva_data(data_fva_anoxic) # save results for later use
    simulated_data_anoxic.import_sra_data(data_srd_anoxic)
    simulated_data_anoxic.import_fva_data(data_fva_anoxic)
    simulated_data_anoxic.check_data()

    # make/load simulated data for aerobic conditions
    data_fva_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_fva_oxic.json'
    data_srd_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_srd_oxic.json'
    simulated_data_oxic = cobra_simulatedData()
    # simulated_data_oxic.generate_sra_data(cobra_model_oxic) # perform single reaction deletion analysis
    # simulated_data_oxic.generate_fva_data(cobra_model_oxic) # perform flux variability analysis
    # simulated_data_oxic.export_sra_data(data_srd_oxic) # save results for later use
    # simulated_data_oxic.export_fva_data(data_fva_oxic) # save results for later use
    simulated_data_oxic.import_sra_data(data_srd_oxic)
    simulated_data_oxic.import_fva_data(data_fva_oxic)
    simulated_data_oxic.check_data()

    ##PART 3: Working
    #-------
    # load pH, ionic_strength, and temperature parameters
    other_data = thermodynamics_otherData()
    other_data.load_defaultData()
    other_data.check_data()

    # generate dG_f data for all model compounds
    # calculate the dG_f for each compound in each compartment
    # export/load and check the dG_f data
    #data_dG0_transformed = 'thermodynamics_data/ijo1366_dG_f.json'
    data_dG0_transformed = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_dG_f01.json'
    dG_f_data = thermodynamics_dG_f_data(id2KEGGID_filename_I='/home/user/code/thermodynamics/thermodynamics_data/id2KEGGID.csv')
    # # # dG_f_data.make_dG0_f_pH0() # only if the data has not been generated previously!
    # dG_f_data.get_transformed_dG_f('/home/user/code/thermodynamics/thermodynamics_data/compounds_dG0_f.json',cobra_model_oxic,other_data.pH,other_data.temperature,other_data.ionic_strength) # adjust the non-transformed dG0_f data to physiological pH, temperature, and ionic strength (this step has already been completed)
    #dG_f_data.export_dG_f(data_dG0_transformed) # save results for later use
    dG_f_data.import_dG_f(data_dG0_transformed)
    dG_f_data.format_dG_f()
    dG_f_data.generate_estimated_dG_f(cobra_model_oxic)
    dG_f_data.check_data()

    ##PART 4: Working
    #-------
    # load metabolomics data for anoxic conditions
    data_concentrations_anoxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_anoxic_geo01.json'
    metabolomics_data_anoxic = thermodynamics_metabolomicsData()
    metabolomics_data_anoxic.import_metabolomics_data(data_concentrations_anoxic)
    metabolomics_data_anoxic.format_metabolomics_data() # add compartment identifiers to metabolite ids
    metabolomics_data_anoxic.generate_estimated_metabolomics_data(cobra_model_anoxic)
    # load metabolomics data for oxic conditions
    data_concentrations_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_oxic_geo01.json'
    metabolomics_data_oxic = thermodynamics_metabolomicsData()
    metabolomics_data_oxic.import_metabolomics_data(data_concentrations_oxic)
    metabolomics_data_oxic.format_metabolomics_data() # add compartment identifiers to metabolite ids
    metabolomics_data_oxic.generate_estimated_metabolomics_data(cobra_model_oxic)

    #PART 5: Working
    -------
    # calculate dG_r and perform a consistency check based on model simulations for anoxic conditions
    data_ta_anoxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_anoxic_ta.csv'
    data_dG0_anoxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_anoxic_dG0.json'
    data_dG_anoxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_anoxic_dG.json'
    data_tcc_anoxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_anoxic_tcc.json'
    tcc_anoxic = thermodynamics_dG_r_data()
    tcc_anoxic.calculate_dG0_r(cobra_model_anoxic, dG_f_data.measured_dG_f, dG_f_data.estimated_dG_f, other_data.temperature) # calculate the change in free energy of reaction without accounting for metabolite concentrations
    tcc_anoxic.calculate_dG_r(cobra_model_anoxic,metabolomics_data_anoxic.measured_concentrations, metabolomics_data_anoxic.estimated_concentrations,
                       other_data.pH, other_data.ionic_strength, other_data.temperature) # adjust the change in free energy of reaction for intracellular metabolite concentrations
    tcc_anoxic.check_thermodynamicConsistency(cobra_model_anoxic,simulated_data_anoxic.fva_data,
                       metabolomics_data_anoxic.measured_concentrations,
                       metabolomics_data_anoxic.estimated_concentrations,
                       other_data.pH,other_data.ionic_strength,other_data.temperature) # check the thermodynamic consistency of the data
    tcc_anoxic.export_dG0_r_json(data_dG0_anoxic) # save for later use
    tcc_anoxic.export_dG_r_json(data_dG_anoxic) # save for later use
    tcc_anoxic.export_tcc_json(data_ta_anoxic) # save for later use
    tcc_anoxic.export_summary(cobra_model_anoxic,simulated_data_anoxic.fva_data,data_ta_anoxic) # write summary of the analysis to csv file

    # calculate dG_r and perform a consistency check based on model simulations for oxic conditions
    data_ta_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_oxic_ta.csv'
    data_dG0_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_oxic_dG0.json'
    data_dG_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_oxic_dG.json'
    data_tcc_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_oxic_tcc.json'
    tcc_oxic = thermodynamics_dG_r_data()
    # tcc_oxic.import_dG0_r_json(data_dG0_oxic)
    # tcc_oxic.import_dG_r_json(data_dG_oxic)
    # tcc_oxic.import_tcc_json(data_tcc_oxic)
    tcc_oxic.calculate_dG0_r(cobra_model_oxic, dG_f_data.measured_dG_f, dG_f_data.estimated_dG_f, other_data.temperature) # calculate the change in free energy of reaction without accounting for metabolite concentrations
    tcc_oxic.calculate_dG_r(cobra_model_oxic,metabolomics_data_oxic.measured_concentrations, metabolomics_data_oxic.estimated_concentrations,
                       other_data.pH, other_data.ionic_strength, other_data.temperature) # adjust the change in free energy of reaction for intracellular metabolite concentrations
    tcc_oxic.check_thermodynamicConsistency(cobra_model_oxic,simulated_data_oxic.fva_data,
                       metabolomics_data_oxic.measured_concentrations,
                       metabolomics_data_oxic.estimated_concentrations,
                       other_data.pH,other_data.ionic_strength,other_data.temperature) # check the thermodynamic consistency of the data
    tcc_oxic.export_dG0_r_json(data_dG0_oxic) # save for later use
    tcc_oxic.export_dG_r_json(data_dG_oxic) # save for later use
    tcc_oxic.export_tcc_json(data_tcc_oxic) # save for later use
    tcc_oxic.export_summary(cobra_model_oxic,simulated_data_oxic.fva_data,data_ta_oxic) # write summary of the analysis to csv file

    ##PART 6: Working
    #-------
    # inspect the thermodynamic analysis results

    # constrain the model solution and simulate optimal growth
    gr_analysis_anoxic = simulate_thermoConstraints(cobra_model_anoxic,['PGCD','ACACT1r','NDPK2'])
    gr_analysis_oxic = simulate_thermoConstraints(cobra_model_oxic,['PGCD','ACACT1r'])

    # expand the reaction set of the anoxic model to reflect the enzyme permiscuity of pykA
    add_pykA(cobra_model_anoxic)
    gr_analysis_anoxic = simulate_thermoConstraints(cobra_model_anoxic,['PGCD','ACACT1r','NDPK2'])

    ##PART 7: Working
    #-------
    # calculate the dG for biosynthetic pathways

    # calculate the dG for biosynthetic pathways for anoxic conditions
    tccp_anoxic = thermodynamics_dG_p_data()
    tccp_anoxic.calculate_dG_p(cobra_model_anoxic,tcc_anoxic.dG0_r,tcc_anoxic.dG_r)
    # calculate the dG for biosynthetic pathways for oxic conditions
    tccp_oxic = thermodynamics_dG_p_data()
    tccp_oxic.calculate_dG_p(cobra_model_oxic,tcc_oxic.dG0_r,tcc_oxic.dG_r)    

    # visualize the results
    ##TODO: methods are defined, but an example is not yet given