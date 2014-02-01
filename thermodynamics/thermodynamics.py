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

# pathways:
pathways = {};
pathways['de novo purine biosynthesis'] = {'reactions':['GLUPRT','PRAGSr','GARFT','PRFGS','PRAIS','AIRC2','AIRC3','PRASCSi','ADSL2r','AICART','IMPC'],
                                           'stoichiometry':[1,1,1,1,1,1,-1,1,1,1,-1]};
pathways['do novo pyrimidine biosynthesis (aerobic)'] = {'reactions':['ASPCT','DHORTS','DHORD2','ORPT','OMPDC','UMPK','NDPK2','CTPS2'],
                                               'stoichiometry':[1,-1,1,-1,1,1,1,1]}
pathways['do novo pyrimidine biosynthesis (anaerobic)'] = {'reactions':['ASPCT','DHORTS','DHORD5','ORPT','OMPDC','UMPK','NDPK2','CTPS2'],
                                               'stoichiometry':[1,-1,1,-1,1,1,1,1]}
pathways['fad biosynthesis'] = {'reactions':['GTPCII2','DHPPDA2','APRAUR','PMDPHT','DB4PS','RBFSa','RBFSb','RBFK','FMNAT'],
                                               'stoichiometry':[1,1,1,1,1,1,1,1,1]}
pathways['nad biosynthesis (aerobic)'] = {'reactions':['ASPO6','QULNS','NNDPR','NNATr','NADS1','NADK','NADPPPS'],
                                               'stoichiometry':[1,1,1,1,1,1,1]}
pathways['nad biosynthesis (anaerobic)'] = {'reactions':['ASPO5','QULNS','NNDPR','NNATr','NADS1','NADK','NADPPPS'],
                                               'stoichiometry':[1,1,1,1,1,1,1]}
pathways['gth biosynthesis'] = {'reactions':['GLUCYS','GTHS','GTHOr'],
                                               'stoichiometry':[1,1,1]}
pathways['nad salvage'] = {'reactions':['NADN','NNAM','NMNN','NAMNPP','NMNDA','NMNAT','NADDP'],
                           'stoichiometry':[1,1,1,1,1,1,1]}
pathways['e4p and prpp to skm'] = {'reactions':['DDPA', 'DHQS', 'SHK3Dr'],
                                    'stoichiometry':[1,1,1]}
pathways['asp-L to thr-L'] = {'reactions':['ASPK', 'ASAD', 'HSDy','HSK','THRS'],
                              'stoichiometry':[1,-1,-1,1,1]}
pathways['hom-L and cyst-L to met-L'] = {'reactions':['HSST', 'SHSL1', 'CYSTL', 'METS'],
                              'stoichiometry':[1,1,1,1]}
pathways['methylglyoxal bypass'] = {'reactions':['MGSA', 'LGTHL', 'GLYOX'],
                              'stoichiometry':[1,1,1]}
pathways['hom-L to thr-L'] = {'reactions':['HSK','THRS'],
                              'stoichiometry':[1,1]}
pathways['3pg to ser-L'] = {'reactions':['PGCD', 'PSERT', 'PSP_L'],
                              'stoichiometry':[1,1,1]}
pathways['prpp to his-L'] = {'reactions':['ATPPRT', 'PRATPP', 'PRAMPC', 'PRMICI', 'IG3PS', 'IGPDH', 'HSTPT', 'HISTP', 'HISTD'],
                              'stoichiometry':[1,1,1,1,1,1,1,1,1]}
pathways['dttp synthesis (anaerobic)'] = {'reactions':['RNTR4c', 'DUTPDP', 'TMDS', 'DTMPK', 'NDPK4'],
                              'stoichiometry':[1,1,1,1,1]}
pathways['dttp synthesis (aerobic)'] = {'reactions':['RNDR4', 'URIDK2r', 'TMDS', 'DTMPK', 'NDPK4'],
                              'stoichiometry':[1,-1,1,1,1]}
pathways['coa biosynthesis'] = {'reactions':['MOHMT', 'DPR', 'PANTS', 'PNTK', 'PPNCL2', 'PPCDC', 'PTPATi', 'DPCOAK'],
                              'stoichiometry':[1,1,1,1,1,1,1,1]}
pathways['thf biosynthesis'] = {'reactions':['GTPCI', 'DNTPPA', 'DNMPPA', 'DHNPA2', 'GCALDD', 'DHNPA2', 'HPPK2', 'ADCS', 'ADCL', 'DHPS2', 'DHFS', 'DHFR'],
                              'stoichiometry':[1,1,1,1,1,1,1,1,1,1,1,1]}
pathways['thr-L to ile-L'] = {'reactions':['THRD_L','ACHBS','KARA2','DHAD2','ILETA'],
                              'stoichiometry':[1,1,1,1,-1]}
pathways['akg to glu-L'] = {'reactions':['GLNS','GLUSy'],
                              'stoichiometry':[1,1]}
pathways['chor to trp-L'] = {'reactions':['ANPRT', 'ANS', 'IGPS', 'PRAIi', 'TRPS2', 'TRPS3'],
                              'stoichiometry':[1,1,1,1,1,1]}
pathways['chor to tyr-L'] = {'reactions':['CHORM', 'PPND', 'TYRTA'],
                              'stoichiometry':[1,1,-1]}
pathways['chor to phe-L'] = {'reactions':['CHORM', 'PPNDH', 'PHETA1'],
                              'stoichiometry':[1,1,-1]}
pathways['OPPP'] = {'reactions':['G6PDH2r', 'GND', 'PGL', 'RPI'],
                              'stoichiometry':[1,1,1,-1]}
pathways['Glycolysis'] = {'reactions':['PGI', 'PFK', 'FBA', 'TPI', 'GAPD', 'PGK', 'PGM', 'ENO'],
                              'stoichiometry':[1,1,1,1,1,-1,-1,1]}
pathways['glu-L to arg-L'] = {'reactions':['ACGS','ACGK','AGPR','ACOTA','ACODA','OCBT','ARGSS','ARGS','ARGSL'],
                              'stoichiometry':[1,1,-1,-1,1,1,1,1,-1]}

for r in pathways['nad salvage']['reactions']:
    print r, cobra_model.reactions.get_by_id(r).build_reaction_string()
for r in pathways['gth biosynthesis']['reactions']:
    print r, cobra_model.reactions.get_by_id(r).build_reaction_string()

data_fva = 'data\\ijo1366_fva_glc.json';
data_srd = 'data\\ijo1366_srd_glc.json';
# make/load simulated data
simulated_data = thermodynamics_simulatedData();
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