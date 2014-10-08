#rpomut01_OxicWTGlyc_0

# Dependencies from cobra
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis import single_deletion
from cobra.flux_analysis.objective import update_objective
# Dependencies from thermodynamics
from thermodynamics.thermodynamics_dG0_f_data import thermodynamics_dG0_f_data
from thermodynamics.thermodynamics_analysis import thermodynamic_analysis
from thermodynamics.thermodynamics_analysis_helper import import_values_json,\
                                                                convert_fluxBounds2var,\
                                                                convert_cv2varAndmM2M_concentrations, \
                                                                compartementalize_concentrations, \
                                                                convert_var2lbub_dG_f, \
                                                                generalize_compartment2all_concentration,\
                                                                generalize_compartment2all_dG_f
# Other dependies
import csv
import json

def _main_():
    # Read in the sbml file and define the model conditions
    ijo1366_sbml = "data\\iJO1366.xml"
    cobra_model = create_cobra_model_from_sbml_file(ijo1366_sbml, print_time=True)
    # Change the objective
    update_objective(cobra_model,{'Ec_biomass_iJO1366_WT_53p95M':1.0})
    # Change uptake reactions for growth on glycerol
    cobra_model.reactions.get_by_id('EX_glc_LPAREN_e_RPAREN_').lower_bound = 0.0;
    cobra_model.reactions.get_by_id('EX_glyc_LPAREN_e_RPAREN_').lower_bound = -10.0;

    # define the pH, ionic strength, and temperature
    pH  = {};
    pH['c'] = {'pH':7.5};
    pH['p'] = {'pH':7.0};
    pH['e'] = {'pH':7.0};

    ionic_strength = {};
    ionic_strength['c'] = {'ionic_strength': .2,'ionic_strength_units': 'M'}
    ionic_strength['p'] = {'ionic_strength': .2,'ionic_strength_units': 'M'}
    ionic_strength['e'] = {'ionic_strength': .2,'ionic_strength_units': 'M'}

    temperature = {};
    temperature['c'] = {'temperature': 310.15,'temperature_units': 'K'}
    temperature['p'] = {'temperature': 310.15,'temperature_units': 'K'}
    temperature['e'] = {'temperature': 310.15,'temperature_units': 'K'}

    # calculate the dG_f for each compound in each compartment
    dG_f_data = thermodynamics_dG0_f_data(id2KEGGID_filename_I='data\\id2KEGGID.csv');
    #dG_f_data.make_dG0_f_pH0(); 
    dG_f_data.get_transformed_dG_f('data\\compounds_dG0_f.json',cobra_model,pH,temperature,ionic_strength);
    dG_f_data.export_dG_f('data\\ijo1366_dG_f.json');

    # calculate the reaction bounds using FVA
    reaction_bounds = flux_variability_analysis(cobra_model, fraction_of_optimum=0.9,
                                  objective_sense='maximize', the_reactions=None,
                                  allow_loops=True, solver='gurobi',
                                  the_problem='return', tolerance_optimality=1e-6,
                                  tolerance_feasibility=1e-6, tolerance_barrier=1e-8,
                                  lp_method=1, lp_parallel=0, new_objective=None,
                                  relax_b=None, error_reporting=None,
                                  number_of_processes=1, copy_model=True);
    # Update the data file
    with open('data\\ijo1366_fva_glyc.json', 'w') as outfile:
        json.dump(reaction_bounds, outfile, indent=4);

    # single reaction deletion
    single_reaction_deletions = single_deletion(cobra_model, element_list=None,
                        method='fba', the_problem='return',
                        element_type='reaction', solver='gurobi',
                        error_reporting=None);

    # FBA
    cobra_model.optimize(solver='gurobi');

    # Update the data file
    reaction_deletions = {};
    for k,v in single_reaction_deletions[0].iteritems():
	    reaction_deletions[k.id] = {'gr':None,'gr_ratio':None};
	    if v:
		    reaction_deletions[k.id] = {'gr':v,'gr_ratio':v/cobra_model.solution.f};

    with open('data\\ijo1366_srd_glyc.json', 'w') as outfile:
        json.dump(reaction_deletions, outfile, indent=4);

    # read in the data
    measured_concentrations_import = import_values_json('data\\rpomut01_OxicWTGlyc_0_geo.json');
    measured_concentrations = compartementalize_concentrations(measured_concentrations_import);

    measured_dG_fs_import = import_values_json('data\\ijo1366_dG_f.json');
    measured_dG_fs = convert_var2lbub_dG_f(measured_dG_fs_import);

    ijo1366_fva = json.load(open('data\\ijo1366_fva_glyc.json'));
    reaction_bounds = convert_fluxBounds2var(ijo1366_fva);

    reaction_deletions = json.load(open('data\\ijo1366_srd_glyc.json'));

    # generate estimates
    estimated_concentrations = generalize_compartment2all_concentration(cobra_model);
    estimated_dG_fs = generalize_compartment2all_dG_f(cobra_model);

    # perform thermodynamics_analysis
    thermodynamic_analysis(cobra_model, reaction_bounds,
                               measured_concentrations, measured_dG_fs,
                               estimated_concentrations, estimated_dG_fs,
                               pH, ionic_strength, temperature, 
                               measured_concentration_coverage_criteria = 0.5,
                               measured_dG_f_coverage_criteria = 1.0);


