# Dependencies from cobra
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis import single_deletion
from cobra.flux_analysis.objective import update_objective
# Dependencies from thermodynamics
from thermodynamics.thermodynamics_dG0_f_data import thermodynamics_dG0_f_data
from thermodynamics.thermodynamics_analysis import thermodynamic_analysis
from thermodynamics.thermodynamics_io import import_values_json, checkInput_concentrations, \
                                                                checkInput_dG_f
from thermodynamics.thermodynamics_formatter import convert_fluxBounds2var,\
                                                                convert_cv2varAndmM2M_concentrations, \
                                                                compartementalize_concentrations, \
                                                                convert_var2lbub_dG_f, \
                                                                generalize_compartment2all_concentration,\
                                                                generalize_compartment2all_dG_f
# Other dependencies
import csv
import json

# data files:
data_dG0_transformed = 'data\\ijo1366_dG_f.json';

data_concentrations = 'data\\ALEWt01_OxicEvo04Glc_0_geo.json';

data_ta1 = json.load(open('data\\ALEWt01_OxicEvo04Glc_ta.json'));
data_dG01 = json.load(open('data\\ALEWt01_OxicEvo04Glc_dG0.json'));
data_dG1 = json.load(open('data\\ALEWt01_OxicEvo04Glc_dG.json'));

data_ta2 = json.load(open('data\\ALEWt01_OxicEvo04Glc_ta2.json'));
data_dG02 = json.load(open('data\\ALEWt01_OxicEvo04Glc_dG02.json'));
data_dG2 = json.load(open('data\\ALEWt01_OxicEvo04Glc_dG2.json'));

data_ta3 = json.load(open('data\\ALEWt01_OxicEvo04Glc_ta3.json'));
data_dG03 = json.load(open('data\\ALEWt01_OxicEvo04Glc_dG03.json'));
data_dG3 = json.load(open('data\\ALEWt01_OxicEvo04Glc_dG3.json'));

print 'thermodynamics_analysis completed'


