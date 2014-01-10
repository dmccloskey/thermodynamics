# A sample script to illustrate the use of GMS 
#
# The tests were run with the 32-bit CPLEX 12.5.0 
# solver, with CoBRApy 0.2.x and Python 2.7 on MacOSX 10.8.4.
# There is a free academic CPLEX license available from IBM.
#
# If running cplex in Python/OSX, be sure to set 32 bit python at the shell
# command line since IBM hasn't released a 64-bit version for OSX.
# For macs:
# export VERSIONER_PYTHON_PREFER_32_BIT=yes
from gim3e.core import gim3e
from gim3e.sampling import gms
import pickle
from types import *
from cobra.io.sbml import create_cobra_model_from_sbml_file
import json

selected_tolerance = 1E-8
selected_growth = 0
selected_penalty = 0

gim3e_dir = gim3e.__file__
n_remove_chars = len('/core/gim3e.py')
gim3e_dir = gim3e_dir[:(-1 * (n_remove_chars))]
data_dir = gim3e_dir + "data/"

# develop the sampling algorithm with E coli core as a first approach
#sbml_file = 'E_coli_core_M9.xml'
#cobra_model = create_cobra_model_from_sbml_file(data_dir + 'E_coli_core_M9.xml', print_time=True)
ijo1366_sbml = "data\\iJO1366.xml"
cobra_model = create_cobra_model_from_sbml_file(ijo1366_sbml, print_time=True)

cobra_model.reactions.get_by_id('Ec_biomass_iJO1366_WT_53p95M').objective_coefficient = 1
cobra_model.optimize()
cobra_model.solution.f

# Make a sampling object
#sampling_object = gms.sample_container(gim3e_model)
sampling_object = gms.sample_container(cobra_model)
gms.sampling_optimize(sampling_object.cobra_model_full, objective_sense = 'maximize', 
                   the_problem = None, solver = 'gurobi',  
                   error_reporting = None,
                   tolerance_optimality = selected_tolerance, 
                   tolerance_feasibility = selected_tolerance,
                   tolerance_barrier = 0.0001 * selected_tolerance,
                   tolerance_integer = selected_tolerance)
# make warmup points
gms.create_warmup_points(sampling_object, solver = 'gurobi', solver_tolerance = 1E-8, force_vms = True, additional_forced_reactions = ['penalty'])
# get rid of redundant points
gms.reduce_warmup_points(sampling_object, solver_tolerance = selected_tolerance)
# sample
gms.achr_sampler(sampling_object, solver = "gurobi", solver_tolerance = 1E-8, max_time = 60 * 60 * 24, n_points = 1000)
gms.save_sampling_object(sampling_object, "sampling_trial_core")

the_converted_results, converted_reaction_list = gms.convert_sampling_results_to_reversible(sampling_object)
the_converted_result_dict = {}
n_reactions, n_points = the_converted_results.shape
for the_row_index, the_row_id in enumerate(converted_reaction_list):
    the_converted_result_dict[the_row_id] = {}
    for the_column_index in range(0, n_points):
        the_converted_result_dict[the_row_id][the_column_index] = the_converted_results[(the_row_index, the_column_index)]
mix_frac = gms.mix_fraction(sampling_object.sampled_points, sampling_object.initial_points, fixed = sampling_object.const_ind)
print(mix_frac)
with open('data\\test2_sampling.json', 'w') as outfile:
    json.dump(the_converted_result_dict, outfile, indent=4);