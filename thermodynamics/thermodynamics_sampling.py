# Dependencies
import operator, json, csv
# Dependencies from 3rd party
import scipy.io
from numpy import histogram, mean, std, loadtxt, savetxt
import matplotlib as mpl
import matplotlib.pyplot as plt
# Dependencies from cobra
from cobra.io.mat import load_matlab_model
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.flux_analysis import flux_variability_analysis, single_deletion
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
from cobra.flux_analysis.objective import update_objective
# Dependencies from thermodynamics
from .thermodynamics_io import thermodynamics_io
from .thermodynamics_utility import mean_confidence_interval

class thermodynamics_sampling(thermodynamics_io):

    def __init__(self):
        self.points = {};
        self.loops = {};

    def get_points_matlab(self,matlab_data,sampler_model_name):
        '''load sampling points from MATLAB'''

        # load model from MATLAB file
        model = load_matlab_model(matlab_data,sampler_model_name);

        # load sample points from MATLAB file into numpy array
        points = scipy.io.loadmat(matlab_data)[sampler_model_name]['points'][0][0];
        #mat = scipy.io.loadmat('data\\EvoWt.mat')
        #points = mat['model_WT_sampler_out']['points'][0][0]

        points_dict = {};
        for i,r in enumerate(model.reactions):
            # convert names:
            r_id_conv = r.id.replace('-','_DASH_');
            r_id_conv = r_id_conv.replace('(','_LPAREN_');
            r_id_conv = r_id_conv.replace(')','_RPAREN_');
            # extract points
            m,lb,ub = mean_confidence_interval(points[i,:],confidence = 0.95)
            points_dict[r_id_conv] = {'points':points[i,:],
                                 'mean':mean(points[i,:]),
                                 'std':std(points[i,:]),
                                 'lb':lb,
                                 'ub':ub}

        self.points = points_dict;

    def get_points_numpy(self,numpy_data,ijo1366_sbml):
        '''load sampling points from numpy file'''

        # load points from numpy file
        points = loadtxt(numpy_data);

        # Read in the sbml file and define the model conditions
        cobra_model = create_cobra_model_from_sbml_file(ijo1366_sbml, print_time=True)

        points_dict = {};
        for i,r in enumerate(cobra_model.reactions):
            # extract points
            m,lb,ub = mean_confidence_interval(points[i,:],confidence = 0.95)
            points_dict[r.id] = {'points':points[i,:],
                                 'mean':mean(points[i,:]),
                                 'std':std(points[i,:]),
                                 'lb':lb,
                                 'ub':ub}

        self.points = points_dict;

    def plot_points(self,reaction_lst=None):
        '''plot sampling points from MATLAB'''
        if not reaction_lst:
            reaction_lst = ['ENO','FBA','FBP','G6PP','GAPD','GLBRAN2',
                        'GLCP','GLCP2','GLDBRAN2','HEX1','PDH','PFK',
                        'PGI','PGK','PGM','PPS','PYK','TPI','ENO_reverse',
                        'FBA_reverse','GAPD_reverse','PGI_reverse',
                        'PGK_reverse','PGM_reverse','TPI_reverse']
        for r in reaction_lst:
            # loop through each reaction in the list
            plt.figure()
            n, bins, patches = plt.hist(self.points[r]['points'],50,label = [r])
            plt.legend()
            plt.show()

    def check_loops(self,cobra_model):
        '''Check if the model contains loops'''

        # Change all uptake reactions to 0
        system_boundaries = [x.id for x in cobra_model.reactions if x.boundary == 'system_boundary'];
        for rxn in cobra_model.reactions:
            if rxn.id in system_boundaries:
            #if 'EX_' in rxn.id and '_LPAREN_e_RPAREN_' in rxn.id:
                cobra_model.reactions.get_by_id(rxn.id).lower_bound = 0.0;
        # Set ATPM to 0
        cobra_model.reactions.get_by_id('ATPM').lower_bound = 0.0

        loops_bool = True;
        cobra_model.optimize(solver='gurobi');
        if not cobra_model.solution.f:
            loops_bool = False;

        return loops_bool;

    def simulate_loops(self,cobra_model,data_fva):
        '''Simulate FVA after closing exchange reactions and setting ATPM to 0
        reactions with flux will be involved in loops'''

        # Change all uptake reactions to 0
        system_boundaries = [x.id for x in cobra_model.reactions if x.boundary == 'system_boundary'];
        for rxn in cobra_model.reactions:
            if rxn.id in system_boundaries:
            #if 'EX_' in rxn.id and '_LPAREN_e_RPAREN_' in rxn.id:
                cobra_model.reactions.get_by_id(rxn.id).lower_bound = 0.0;
        # Set ATPM to 0
        cobra_model.reactions.get_by_id('ATPM').lower_bound = 0.0;

        # calculate the reaction bounds using FVA
        reaction_bounds = flux_variability_analysis(cobra_model, fraction_of_optimum=1.0,
                                          the_reactions=None, solver='gurobi');

        # Update the data file
        with open(data_fva, 'wb') as outfile:
            json.dump(reaction_bounds, outfile, indent=4);

    def simulate_loops_sbml(self,ijo1366_sbml,data_fva):
        '''Simulate FVA after closing exchange reactions and setting ATPM to 0
        reactions with flux will be involved in loops'''

        # Read in the sbml file and define the model conditions
        cobra_model = create_cobra_model_from_sbml_file(ijo1366_sbml, print_time=True)
        # Change all uptake reactions to 0
        for rxn in cobra_model.reactions:
            if 'EX_' in rxn.id and '_LPAREN_e_RPAREN_' in rxn.id:
                rxn.lower_bound = 0.0;
        # Set ATPM to 0
        cobra_model.reactions.get_by_id('ATPM').lower_bound = 0.0

        # calculate the reaction bounds using FVA
        reaction_bounds = flux_variability_analysis(cobra_model, fraction_of_optimum=0.9,
                                          objective_sense='maximize', the_reactions=None,
                                          allow_loops=True, solver='gurobi',
                                          the_problem='return', tolerance_optimality=1e-6,
                                          tolerance_feasibility=1e-6, tolerance_barrier=1e-8,
                                          lp_method=1, lp_parallel=0, new_objective=None,
                                          relax_b=None, error_reporting=None,
                                          number_of_processes=1, copy_model=False);

        # Update the data file
        with open(data_fva, 'wb') as outfile:
            json.dump(reaction_bounds, outfile, indent=4);

    def find_loops(self,data_fva):
        '''extract out loops from simulate_loops'''

        data_loops = json.load(open(data_fva))
        rxn_loops = [];
        for k,v in data_loops.items():
            if abs(v['minimum'])>1.0 or abs(v['maximum'])>1.0:
                rxn_loops.append(k);
        #return rxn_loops
        self.loops = rxn_loops;

    def remove_loopsFromPoints(self):
        '''remove reactions with loops from sampling points'''

        points_loopless = {};
        for k,v in self.points.items():
            if k in self.loops: continue
            else: 
                points_loopless[k] = {'points':v,
                                 'mean':v['mean'],
                                 'std':v['std']};

        #return points_loopless_mean;
        self.points = points_loopless;

    def export_points_numpy(self,filename):
        '''export sampling points'''

        savetxt(filename,self.points);

    def export_sampling_matlab(self,cobra_model,filename_model='C:\\Users\\Public\\Documents\\sample_model.mat',filename_script='C:\\Users\\Public\\Documents\\sample_script.m'):
        '''export model and script for sampling using matlab cobra_toolbox'''
        # copy the model:
        cobra_model_copy = cobra_model.copy();
        # optimize
        cobra_model_copy.optimize(solver='gurobi');
        # confine the objective to a fraction of maximum optimal
        objective = [x.id for x in cobra_model_copy.reactions if x.objective_coefficient == 1]
        cobra_model_copy.reactions.get_by_id(objective[0]).upper_bound = fraction_optimal * cobra_model_copy.solution.f;
        # write model to mat
        save_matlab_model(cobra_model_copy,filename_model);
        ## write model to xml
        #write_sbml_model(cobra_model_copy,'data\\tsampling\\tsampler_conc_ln.xml');
        # write the sampling script to file\
        mat_script = "% initialize with Tomlab_CPLEX\n"+\
                      "load('" + filename_model + "')\n"+\
                      "initCobraToolbox();\n"+\
                      "% sample\n"+\
                      "[tsampler_out, mixedFrac] = gpSampler(" + cobra_model_copy.description + ", [], [], [], [], [], true);\n"+\
                      "[tsampler_out, mixedFrac] = gpSampler(tsampler_out, [], [], [], 20000, [], true);";
        with open(filename_script,'w') as f:
            f.write(mat_script);