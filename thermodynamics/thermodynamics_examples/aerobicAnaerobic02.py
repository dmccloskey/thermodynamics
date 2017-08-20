# -*- coding: utf-8 -*-
# Dependencies from cobra
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.manipulation.modify import convert_to_irreversible
from cobra.flux_analysis.parsimonious import optimize_minimal_flux

# Other dependencies
import csv,json,sys

# Dependencies from thermodynamics
# from thermodynamics.thermodynamics_simulatedData import thermodynamics_simulatedData
from thermodynamics.thermodynamics_simulatedData import thermodynamics_simulatedData
from cobra_sampling.optGpSampler_sampling import  optGpSampler_sampling
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
    # Aerobic specific changes
    cobra_model_oxic = load_thermoModel(anoxic = False)
    convert_to_irreversible(cobra_model_oxic)

    ##PART 2: Working
    #-------
    # make/load simulated data for aerobic conditions
    data_fva_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_fva_irrev_oxic.json'
    data_srd_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_srd_irrev_oxic.json'
    simulated_data_oxic = thermodynamics_simulatedData()
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
    # load metabolomics data for oxic conditions
    data_concentrations_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_oxic_geo01.json'
    metabolomics_data_oxic = thermodynamics_metabolomicsData()
    metabolomics_data_oxic.import_metabolomics_data(data_concentrations_oxic)
    metabolomics_data_oxic.format_metabolomics_data() # add compartment identifiers to metabolite ids
    metabolomics_data_oxic.generate_estimated_metabolomics_data(cobra_model_oxic)

    ##PART 5: Working
    #-------
    # calculate dG_r and perform a consistency check based on model simulations for oxic conditions
    data_ta_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_oxic_irrev_ta.csv'
    data_dG0_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_oxic_orrev_dG0.json'
    data_dG_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_oxic_irrev_dG.json'
    data_tcc_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_oxic_irrev_tcc.json'
    tcc_oxic = thermodynamics_dG_r_data()
    tcc_oxic.import_dG0_r_json(data_dG0_oxic)
    tcc_oxic.import_dG_r_json(data_dG_oxic)
    tcc_oxic.import_tcc_json(data_tcc_oxic)
    # tcc_oxic.calculate_dG0_r(cobra_model_oxic, dG_f_data.measured_dG_f, dG_f_data.estimated_dG_f, other_data.temperature) # calculate the change in free energy of reaction without accounting for metabolite concentrations
    # tcc_oxic.calculate_dG_r(cobra_model_oxic,metabolomics_data_oxic.measured_concentrations, metabolomics_data_oxic.estimated_concentrations,
    #                    other_data.pH, other_data.ionic_strength, other_data.temperature) # adjust the change in free energy of reaction for intracellular metabolite concentrations
    # tcc_oxic.check_thermodynamicConsistency(cobra_model_oxic,simulated_data_oxic.fva_data,
    #                    metabolomics_data_oxic.measured_concentrations,
    #                    metabolomics_data_oxic.estimated_concentrations,
    #                    other_data.pH,other_data.ionic_strength,other_data.temperature) # check the thermodynamic consistency of the data
    # tcc_oxic.export_dG0_r_json(data_dG0_oxic) # save for later use
    # tcc_oxic.export_dG_r_json(data_dG_oxic) # save for later use
    # tcc_oxic.export_tcc_json(data_tcc_oxic) # save for later use
    # tcc_oxic.export_summary(cobra_model_oxic,simulated_data_oxic.fva_data,data_ta_oxic) # write summary of the analysis to csv file    
    
    ##PART 8:
    #-------
	# Diagnose model variables and constraints prior to FBA, FVA, and sampling (Oxic condition only)

    # identified inconsistent concentrations/dG_f/tcc values
    inconsistent_concentrations_I = []
    inconsistent_dG_f_I = []
    inconsistent_tcc_I = ['EX_glc_LPAREN_e_RPAREN__reverse']
    # remove an inconsistent dGf values
    dG_f_data.remove_measured_dG_f(inconsistent_dG_f_I)
    # remove an inconsistent concentration values
    metabolomics_data_oxic.remove_measured_concentrations(inconsistent_concentrations_I)
    # remove an inconcsistent tcc
    tcc_oxic.change_feasibleReactions(inconsistent_tcc_I)    
    # diagnose tfba constraints
    tfba = thermodynamics_tfba()
    # thermodynamic_constraints_check,\
    #     inconsistent_tcc,diagnose_variables_1,\
    #     diagnose_variables_2,\
    #     diagnose_variables_3 = tfba.check_conc_ln_constraints_transport(
    #         cobra_model_oxic,
    #         metabolomics_data_oxic.measured_concentrations, metabolomics_data_oxic.estimated_concentrations,
    #         tcc_oxic.dG0_r, other_data.pH,other_data.temperature,tcc_oxic.metabolomics_coverage,
    #         tcc_oxic.dG_r_coverage, tcc_oxic.thermodynamic_consistency_check,
    #         0.5, 0.99,n_checks_I = 0,
    #         diagnose_solver_I=None,diagnose_threshold_I=0.98,diagnose_break_I=0.1)        
    # tfba._add_conc_ln_constraints_transport(cobra_model_oxic,
    #         metabolomics_data_oxic.measured_concentrations, metabolomics_data_oxic.estimated_concentrations,
    #         tcc_oxic.dG0_r, other_data.pH,other_data.temperature,tcc_oxic.metabolomics_coverage,
    #         tcc_oxic.dG_r_coverage, tcc_oxic.thermodynamic_consistency_check,
    #         0.5, 0.99)
    # # cobra_model_oxic.optimize()
    # from cobra.flux_analysis import flux_variability_analysis
    # fva_data = flux_variability_analysis(cobra_model_oxic, fraction_of_optimum=0.9,
    #                                   objective_sense='maximize',
    #                                   reaction_list=[cobra_model_oxic.reactions.get_by_id('conc_lnv_fum_c')],
    #                                   )
    
    # #PART 9:
    # -------
	# perform thermodynamic FBA and FVA (Oxic condition only)

    # # run TFBA
    # cobra_model_copy = cobra_model_oxic.copy()
    # tfba.tfba(cobra_model_copy,
    #     tcc_oxic.dG0_r,other_data.temperature,
    #     tcc_oxic.dG_r_coverage, tcc_oxic.thermodynamic_consistency_check,
    #     use_measured_dG_r=True, solver='glpk',)

    # cobra_model_copy = cobra_model_oxic.copy()
    # tfba.tfba_conc_ln(cobra_model_copy, 
    #     metabolomics_data_oxic.measured_concentrations, metabolomics_data_oxic.estimated_concentrations,
    #     tcc_oxic.dG0_r,other_data.temperature,tcc_oxic.metabolomics_coverage,
    #     tcc_oxic.dG_r_coverage, tcc_oxic.thermodynamic_consistency_check,
    #     measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,
    #     use_measured_concentrations=True,use_measured_dG0_r=True, solver='glpk',)
    
    # run TFVA
    data_tfva_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_oxic_irrev_tfva.csv'
    data_tfva_analysis_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_oxic_irrev_tfva_analysis.csv'
    data_tfva_dG_r_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_oxic_orrev_tfva_dG_r.json'
    data_tfva_concentrations_oxic = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/aerobicAnaerobic01_oxic_irrev_tfva_concentrations.json'
    cobra_model_copy = cobra_model_oxic.copy()
    tfba.tfva(cobra_model_copy, 
        tcc_oxic.dG0_r,other_data.temperature,
        tcc_oxic.dG_r_coverage, tcc_oxic.thermodynamic_consistency_check,
        use_measured_dG0_r=True, reaction_list=None,fraction_of_optimum=1.0, solver='glpk',
        objective_sense="maximize")
    tfba.export_tfva_data(data_tfva_oxic)
    tfba.analyze_tfva_results(flux_threshold=1e-6)
    tfba.export_tfva_analysis(data_tfva_analysis_oxic)

    cobra_model_copy = cobra_model_oxic.copy()
    tfba.tfva_dG_r(cobra_model_copy, 
        tcc_oxic.dG0_r,other_data.temperature,
        tcc_oxic.dG_r_coverage, tcc_oxic.thermodynamic_consistency_check,
        use_measured_dG0_r=True, fraction_of_optimum=1.0, solver='glpk',
        objective_sense="maximize")
    tfba.export_tfva_dG_r_data(data_tfva_dG_r_oxic)

    cobra_model_copy = cobra_model_oxic.copy()
    tfba.tfva_concentrations(cobra_model_copy, 
        metabolomics_data_oxic.measured_concentrations, metabolomics_data_oxic.estimated_concentrations,
        tcc_oxic.dG0_r,other_data.temperature,tcc_oxic.metabolomics_coverage,
        tcc_oxic.dG_r_coverage, tcc_oxic.thermodynamic_consistency_check,
        measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,
        use_measured_concentrations=True,use_measured_dG0_r=True,fraction_of_optimum=1.0, solver='glpk',
        objective_sense="maximize")
    tfba.export_tfva_concentrations_data(data_tfva_concentrations_oxic)
    
    # ##PART 10:
    # #-------
	# # perform thermodynamic Tsampling (Oxic condition only)
    # # NOTE: requires optGpSampler
    
    # # run Tsampling
    # data_dir = '/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo'
    # sampling = optGpSampler_sampling(data_dir_I = data_dir);
    # simulation_id_I = 'aerobicAnaerobic01_oxic_tsampling'
    # filename_model = simulation_id_I + '.json';
    # filename_script = simulation_id_I + '.py';
    # filename_points = simulation_id_I + '_points' + '.json';
    # filename_warmup = simulation_id_I + '_warmup' + '.json';
    # sampling.export_sampling_optGpSampler(cobra_model=cobra_model_oxic,
    #     filename_model=filename_model,
    #     filename_script=filename_script,
    #     filename_points=filename_points,
    #     filename_warmup=filename_warmup,
    #     solver_id_I = 'optGpSampler',
    #     n_points_I = 2*len(cobra_model_oxic.reactions),
    #     n_steps_I = 5000,
    #     n_threads_I = 2)
    
    ##PART 11:
    #-------
	# Analyze thermodynamic sampling (Oxic condition only)

    sampling = optGpSampler_sampling(
        data_dir_I = data_dir,
        model_I=cobra_model_oxic);
    sampling.get_points_json(filename_points);
    sampling.get_warmup_json(filename_warmup);
    sampling.calculate_mixFraction();
    # check if the model contains loops
    #loops_bool = self.sampling.check_loops();
    sampling.simulate_loops(
        data_fva='/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/loops_fva_tmp.json',
        solver_I = 'optGpSampler');
    sampling.find_loops(data_fva='/home/user/code/thermodynamics/thermodynamics_data/aerobicAnaerobic01_geo/loops_fva_tmp.json');
    sampling.remove_loopsFromPoints();
    # calculate the flux descriptive statistics
    sampling.descriptive_statistics(points_I='flux');
    # calculate descriptive stats for metabolites
    sampling.convert_points2MetabolitePoints();
    sampling.descriptive_statistics(points_I='metabolite');
    # calculate descriptive stats for subsystems
    sampling.convert_points2SubsystemPoints();
    sampling.descriptive_statistics(points_I='subsystem');
    

    # visualize the results
    ##TODO: methods are defined, but an example is not yet given