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
        self.other_data = other_data

    def test_dG_f_data(self):
        # generate dG_f data for all model compounds
        # calculate the dG_f for each compound in each compartment
        # export/load and check the dG_f data
        data_dG0_transformed = data_dir + '/test_dG_f01.json'
        dG_f_data = thermodynamics_dG_f_data(id2KEGGID_filename_I=data_dir + '/id2KEGGID.csv')
        assert(dG_f_data.KEGGID2id['C00234'] == '10fthf')
        assert(dG_f_data.id2KEGGID['10fthf'] == 'C00234')
        dG_f_data.get_transformed_dG_f(data_dir + '/compounds_dG0_f.json',
            self.cobra_model,self.other_data.pH,
            self.other_data.temperature,self.other_data.ionic_strength
            ) # adjust the non-transformed dG0_f data to physiological pH, temperature, and ionic strength
        assert(dG_f_data.dG0_f['C00236']['dG_f'] == -2354.086788)
        assert(dG_f_data.dG0_f['C00236']['dG_f_units'] == 'kJ/mol')
        assert(dG_f_data.dG0_f['C00236']['dG_f_var'] == 17.8)
        assert(dG_f_data.dG_f['13dpg_c']['dG_f'] == -2094.9687955266982)
        assert(dG_f_data.dG_f['13dpg_c']['dG_f_units'] == 'kJ/mol')
        assert(dG_f_data.dG_f['13dpg_c']['dG_f_var'] == 17.8)
        dG_f_data.export_dG_f(data_dG0_transformed) # save results for later use
        dG_f_data.import_dG_f(data_dG0_transformed)
        dG_f_data.format_dG_f()
        assert(dG_f_data.measured_dG_f['13dpg_c']['dG_f'] == -2094.9687955266982)
        assert(dG_f_data.measured_dG_f['13dpg_c']['dG_f_units'] == 'kJ/mol')
        assert(dG_f_data.measured_dG_f['13dpg_c']['dG_f_var'] == 17.8)
        assert(dG_f_data.measured_dG_f['13dpg_c']['dG_f_lb'] == -2099.187800148644)
        assert(dG_f_data.measured_dG_f['13dpg_c']['dG_f_ub'] == -2090.7497909047524)        
        dG_f_data.generate_estimated_dG_f(self.cobra_model)
        assert(dG_f_data.estimated_dG_f['13dpg_c']['dG_f'] == 0)
        assert(dG_f_data.estimated_dG_f['13dpg_c']['dG_f_units'] == 'kJ/mol')
        assert(dG_f_data.estimated_dG_f['13dpg_c']['dG_f_var'] == 1000000000000.0)
        assert(dG_f_data.estimated_dG_f['13dpg_c']['dG_f_lb'] == -1000000.0)
        assert(dG_f_data.estimated_dG_f['13dpg_c']['dG_f_ub'] == 1000000.0)  
        dG_f_data.check_data()
        self.dG_f_data = dG_f_data

    def test_metabolomicsData(self):
        # load metabolomics data
        data_concentrations = data_dir + '/test_geo01.json'
        metabolomics_data = thermodynamics_metabolomicsData()
        metabolomics_data.import_metabolomics_data(data_concentrations)
        assert(metabolomics_data.measured_concentrations['pep']['concentration'] == 0.000122661)  
        assert(metabolomics_data.measured_concentrations['pep']['concentration_lb'] == 2.73e-05)  
        assert(metabolomics_data.measured_concentrations['pep']['concentration_ub'] == 0.000551246)  
        assert(metabolomics_data.measured_concentrations['pep']['concentration_var'] == 1.32798777)  
        assert(metabolomics_data.measured_concentrations['pep']['concentration_units'] == 'M')  
        metabolomics_data.format_metabolomics_data() # add compartment identifiers to metabolite ids
        assert(metabolomics_data.measured_concentrations['pep_c']['concentration'] == 0.000122661)  
        assert(metabolomics_data.measured_concentrations['pep_c']['concentration_lb'] == 2.73e-05)  
        assert(metabolomics_data.measured_concentrations['pep_c']['concentration_ub'] == 0.000551246)  
        assert(metabolomics_data.measured_concentrations['pep_c']['concentration_var'] == 1.32798777)  
        assert(metabolomics_data.measured_concentrations['pep_c']['concentration_units'] == 'M')  
        metabolomics_data.generate_estimated_metabolomics_data(self.cobra_model)
        assert(metabolomics_data.estimated_concentrations['pep_c']['concentration'] == 5.000000000000004e-05)  
        assert(metabolomics_data.estimated_concentrations['pep_c']['concentration_lb'] == 1.5811388300841896e-06)  
        assert(metabolomics_data.estimated_concentrations['pep_c']['concentration_ub'] == 0.0015811388300841897)  
        assert(metabolomics_data.estimated_concentrations['pep_c']['concentration_var'] == 176.21560953198042)  
        assert(metabolomics_data.estimated_concentrations['pep_c']['concentration_units'] == 'M')  
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
        tcc.calculate_dG0_r(self.cobra_model, 
            self.dG_f_data.measured_dG_f, self.dG_f_data.estimated_dG_f, 
            self.other_data.temperature) # calculate the change in free energy of reaction without accounting for metabolite concentrations
        assert(tcc.dG0_r['ENO']['dG_r'] == -3.8598069140609823)
        assert(tcc.dG0_r['ENO']['dG_r_var'] == 53.400000000000006)
        assert(tcc.dG0_r['ENO']['dG_r_lb'] == -8.078811536006697)
        assert(tcc.dG0_r['ENO']['dG_r_ub'] == 0.35919770788473215)
        assert(tcc.dG0_r['ENO']['dG_r_units'] == 'kJ/mol')
        assert(tcc.dG0_r['ENO']['Keq_lb'] == 22.943577827979983)
        assert(tcc.dG0_r['ENO']['Keq_ub'] == 0.8699668221299198)
        tcc.calculate_dG_r(self.cobra_model,self.metabolomics_data.measured_concentrations,
            self.metabolomics_data.estimated_concentrations,
            self.other_data.pH, self.other_data.ionic_strength, self.other_data.temperature) # adjust the change in free energy of reaction for intracellular metabolite concentrations
        assert(tcc.dG_r['ENO']['dG_r'] == 8.716270424972741)
        assert(tcc.dG_r['ENO']['dG_r_var'] == 53.400000000000006)
        assert(tcc.dG_r['ENO']['dG_r_lb'] == 4.793062495527913)
        assert(tcc.dG_r['ENO']['dG_r_ub'] == 12.628180518569243)
        assert(tcc.dG_r['ENO']['dG_r_units'] == 'kJ/mol')
        assert(tcc.dG_r['ENO']['Keq_lb'] == 0.15586046876150242)
        assert(tcc.dG_r['ENO']['Keq_ub'] == 0.007466525151002594)
        assert(tcc.dG_r_coverage['ENO'] == 1.0)
        tcc.check_thermodynamicConsistency(self.cobra_model,self.simulated_data.fva_data,
            self.metabolomics_data.measured_concentrations,
            self.metabolomics_data.estimated_concentrations,
            self.other_data.pH,self.other_data.ionic_strength,self.other_data.temperature) # check the thermodynamic consistency of the data  
        assert(not tcc.thermodynamic_consistency_check['ENO'])
        assert(tcc.thermodynamic_consistency_check['ATPM'])
        tcc.export_dG0_r_json(data_dG0) # save for later use
        tcc.export_dG_r_json(data_dG) # save for later use
        tcc.export_tcc_json(data_tcc) # save for later use
        tcc.export_summary(self.cobra_model,self.simulated_data.fva_data,data_ta) # write summary of the analysis to csv file
        tcc.import_dG0_r_json(data_dG0)
        tcc.import_dG_r_json(data_dG)
        tcc.import_tcc_json(data_tcc)
        assert(tcc.dG0_r['ENO']['dG_r'] == -3.8598069140609823)
        assert(not tcc.thermodynamic_consistency_check['ENO'])
        assert(tcc.dG_r['ENO']['dG_r'] == 8.716270424972741)
        assert(tcc.dG_r_coverage['ENO'] == 1.0)
        self.tcc = tcc

    def test_dG_p_data(self):        
        # calculate the dG for biosynthetic pathways
        tccp = thermodynamics_dG_p_data()
        assert(tccp.pathways['Glycolysis']['reactions'] == ['PGI', 'PFK', 'FBA', 'TPI', 'GAPD', 'PGK', 'PGM', 'ENO'])
        assert(tccp.pathways['Glycolysis']['stoichiometry'] == [1, 1, 1, 1, 1, -1, -1, 1])
        tccp.calculate_dG_p(self.cobra_model,self.tcc.dG0_r,self.tcc.dG_r)
        assert(tccp.dG0_p['Glycolysis']['dG0_p'] == 1.9688173094308468)
        assert(tccp.dG0_p['Glycolysis']['dG0_p_var'] == 267.0)
        assert(tccp.dG0_p['Glycolysis']['dG0_p_lb'] == -10.688196556406865)
        assert(tccp.dG0_p['Glycolysis']['dG0_p_ub'] == 14.625831175268445)
        assert(tccp.dG0_p['Glycolysis']['dG0_p_units'] == 'kJ/mol')
        assert(tccp.dG_p['Glycolysis']['dG_p'] == -2.93851827063736)
        assert(tccp.dG_p['Glycolysis']['dG_p_var'] == 267.0)
        assert(tccp.dG_p['Glycolysis']['dG_p_lb'] == -29.148908111427986)
        assert(tccp.dG_p['Glycolysis']['dG_p_ub'] == 23.27562458558087)
        assert(tccp.dG_p['Glycolysis']['dG_p_units'] == 'kJ/mol')
        self.tccp = tccp    
    
    def test_tfba_constraints(self):
        # Diagnose model variables and constraints prior to FBA, FVA

        # identified inconsistent concentrations/dG_f/tcc values
        inconsistent_concentrations_I = []
        inconsistent_dG_f_I = []
        inconsistent_tcc_I = []
        # remove an inconsistent dGf values
        self.dG_f_data.remove_measured_dG_f(inconsistent_dG_f_I)
        # remove an inconsistent concentration values
        self.metabolomics_data.remove_measured_concentrations(inconsistent_concentrations_I)
        # remove an inconcsistent tcc
        self.tcc.change_feasibleReactions(inconsistent_tcc_I)    
        # diagnose tfba constraints
        tfba = thermodynamics_tfba()
        # test concentration constraints
        thermodynamic_constraints_check,\
            inconsistent_tcc,diagnose_variables_1,\
            diagnose_variables_2,\
            diagnose_variables_3 = tfba.check_conc_ln_constraints_transport(
                self.cobra_model,
                self.metabolomics_data.measured_concentrations, self.metabolomics_data.estimated_concentrations,
                self.tcc.dG0_r, self.other_data.pH,self.other_data.temperature,
                self.tcc.metabolomics_coverage,
                self.tcc.dG_r_coverage, self.tcc.thermodynamic_consistency_check,
                n_checks_I = 2,
                diagnose_solver_I='glpk',diagnose_threshold_I=29,diagnose_break_I=29)
        assert(diagnose_variables_1['ENO']['solution_before'] == 30.0)
        assert(diagnose_variables_1['ENO']['solution_after'] == 30.0)
        assert(diagnose_variables_2['ENO']['solution_before'] == 30.0)
        assert(diagnose_variables_2['ENO']['solution_after'] == 30.0)
        assert(diagnose_variables_3['ENO']['solution_before'] == 30.0)
        assert(diagnose_variables_3['ENO']['solution_after'] == 30.0)
        assert(not thermodynamic_constraints_check['ENO'])
        assert(thermodynamic_constraints_check['ATPM'])
        assert('ENO' in inconsistent_tcc)  
    
    def test_tfba(self):
        #perform thermodynamic FBA
        tfba = thermodynamics_tfba()

        # # remove inconsistent reactions
        # inconsistent_tcc_I = ['ENO', 'EX_glc__D_e', 'EX_h_e', 'EX_lac__D_e', 'FBA', 'GAPD', 'GLCpts', 'LDH_D', 'PFK', 'PGI', 'PGK',
        # 'PGM', 'PYK', 'TPI', 'ENO_reverse', 'EX_glc__D_e_reverse', 'EX_h_e_reverse', 'FBA_reverse',
        # 'GAPD_reverse', 'LDH_D_reverse', 'PGI_reverse', 'PGK_reverse', 'PGM_reverse', 'TPI_reverse']
        # self.tcc.change_feasibleReactions(inconsistent_tcc_I)

        # run TFBA
        cobra_model_copy = self.cobra_model.copy()
        tfba.tfba(cobra_model_copy,
            self.tcc.dG_r,self.other_data.temperature,
            self.tcc.dG_r_coverage, self.tcc.thermodynamic_consistency_check,
            use_measured_dG_r=True, solver='glpk',)
        assert(cobra_model_copy.objective.value == 12.585)
        assert(tfba.tfba_data['ENO'] == 8.3900000000000006)

        cobra_model_copy = self.cobra_model.copy()
        tfba.tfba_conc_ln(cobra_model_copy, 
            self.metabolomics_data.measured_concentrations, self.metabolomics_data.estimated_concentrations,
            self.tcc.dG0_r,self.other_data.temperature,self.tcc.metabolomics_coverage,
            self.tcc.dG_r_coverage, self.tcc.thermodynamic_consistency_check,
            measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,
            use_measured_concentrations=True,use_measured_dG0_r=True, solver='glpk',)
        assert(cobra_model_copy.objective.value == 12.585)
        assert(tfba.tfba_data['ENO'] == 8.5962128757921494)
    
    def test_tfva(self):
        # run TFVA
        tfba = thermodynamics_tfba()

        data_tfva = data_dir + 'test_tfva.csv'
        data_tfva_analysis = data_dir + 'test_tfva_analysis.csv'
        data_tfva_dG_r = data_dir + 'test_tfva_dG_r.json'
        data_tfva_concentrations = data_dir + 'test_tfva_concentrations.json'
        cobra_model_copy = self.cobra_model.copy()
        tfba.tfva(cobra_model_copy, 
            self.tcc.dG0_r,self.other_data.temperature,
            self.tcc.dG_r_coverage, self.tcc.thermodynamic_consistency_check,
            use_measured_dG0_r=True, reaction_list=None,fraction_of_optimum=1.0, solver='glpk',
            objective_sense="maximize")
        assert(tfba.fva_data['ENO']['maximum'] == 1000.0)
        assert(tfba.fva_data['ENO']['minimum'] == 18.0)
        tfba.export_tfva_data(data_tfva)
        tfba.analyze_tfva_results(flux_threshold=1e-6)
        assert(tfba.tfva_analysis['ENO'] == 30.0)
        tfba.export_tfva_analysis(data_tfva_analysis)

        cobra_model_copy = self.cobra_model.copy()
        tfba.tfva_dG_r(cobra_model_copy, 
            self.tcc.dG0_r,other_data.temperature,
            self.tcc.dG_r_coverage, self.tcc.thermodynamic_consistency_check,
            use_measured_dG0_r=True, fraction_of_optimum=1.0, solver='glpk',
            objective_sense="maximize")
        assert(tfba.tfva_dG_r_data['dGr_rv_ENO']['maximum'] == 1000.0)
        assert(tfba.tfva_dG_r_data['dGr_rv_ENO']['minimum'] == 18.0)
        tfba.export_tfva_dG_r_data(data_tfva_dG_r)

        cobra_model_copy = self.cobra_model.copy()
        tfba.tfva_concentrations(cobra_model_copy, 
            self.metabolomics_data.measured_concentrations, self.metabolomics_data.estimated_concentrations,
            self.tcc.dG0_r,self.other_data.temperature,self.tcc.metabolomics_coverage,
            self.tcc.dG_r_coverage, self.tcc.thermodynamic_consistency_check,
            measured_concentration_coverage_criteria = 0.5, measured_dG_f_coverage_criteria = 0.99,
            use_measured_concentrations=True,use_measured_dG0_r=True,fraction_of_optimum=1.0, solver='glpk',
            objective_sense="maximize")
        assert(tfba.tfva_concentration_data['conc_lnv_pep_c']['maximum'] == 1000.0)
        assert(tfba.tfva_concentration_data['conc_lnv_pep_c']['minimum'] == 18.0)
        tfba.export_tfva_concentrations_data(data_tfva_concentrations)
    
    def test_tsampling(self):
        # perform thermodynamic Tsampling
        sampling = optGpSampler_sampling(data_dir_I = data_dir);
        simulation_id_I = 'test_tsampling'
        filename_model = simulation_id_I + '.json';
        filename_script = simulation_id_I + '.py';
        filename_points = simulation_id_I + '_points' + '.json';
        filename_warmup = simulation_id_I + '_warmup' + '.json';
        sampling.export_sampling_optGpSampler(cobra_model=self.cobra_model,
            filename_model=filename_model,
            filename_script=filename_script,
            filename_points=filename_points,
            filename_warmup=filename_warmup,
            solver_id_I = 'optGpSampler',
            n_points_I = 2*len(self.cobra_model.reactions),
            n_steps_I = 5000,
            n_threads_I = 2)
        assert(diagnose_variables_1['ENO']['solution_before'] == 30.0)
    
    def test_tsampling_analysis(self):
        # Analyze thermodynamic sampling
        simulation_id_I = 'test_tsampling'
        filename_points = simulation_id_I + '_points' + '.json';
        filename_warmup = simulation_id_I + '_warmup' + '.json';
        sampling = optGpSampler_sampling(
            data_dir_I = data_dir,
            model_I=self.cobra_model);
        sampling.get_points_json(filename_points);
        sampling.get_warmup_json(filename_warmup);
        sampling.calculate_mixFraction();
        assert(diagnose_variables_1['ENO']['solution_before'] == 30.0)
        # check if the model contains loops
        #loops_bool = self.sampling.check_loops();
        sampling.simulate_loops(
            data_fva = data_dir + 'test_loops_fva.json',
            solver_I = 'optGpSampler');
        sampling.find_loops(data_fva = data_dir + 'test_loops_fva.json');
        sampling.remove_loopsFromPoints();
        assert(diagnose_variables_1['ENO']['solution_before'] == 30.0)
        # calculate the flux descriptive statistics
        sampling.descriptive_statistics(points_I='flux');
        assert(diagnose_variables_1['ENO']['solution_before'] == 30.0)
        # calculate descriptive stats for metabolites
        sampling.convert_points2MetabolitePoints();
        sampling.descriptive_statistics(points_I='metabolite');
        assert(diagnose_variables_1['ENO']['solution_before'] == 30.0)
        # calculate descriptive stats for subsystems
        sampling.convert_points2SubsystemPoints();
        sampling.descriptive_statistics(points_I='subsystem');
        assert(diagnose_variables_1['ENO']['solution_before'] == 30.0)