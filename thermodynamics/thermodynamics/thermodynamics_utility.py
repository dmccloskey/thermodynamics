# Dependencies
import operator, json, csv
from collections import Counter
# Dependencies from cobra
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.io.sbml import write_cobra_model_to_sbml_file
from cobra.io.mat import save_matlab_model
from cobra.manipulation.modify import convert_to_irreversible, revert_to_reversible
from cobra.flux_analysis.objective import update_objective
from cobra.core.Model import Model
from cobra.core.Metabolite import Metabolite
from cobra.core.Reaction import Reaction

def find_transportMets(cobra_model_I, reaction_id_I):
    # transport metabolite definition:
    #	1. different id (same base id but different compartment)
    #	2. same name
    #	3. different compartment

    ''' Not full proof
    import re
    compartments = list(set(cobra_model.metabolites.list_attr('compartment')));
    met_ID
    for compart in compartments:
	    comp_tmp = '_' + compart;
	    met_tmp = re.sub(comp_tmp,'',met_ID)
    met_ids.append(met_tmp)'''

    met_ids = cobra_model_I.reactions.get_by_id(reaction_id_I).get_products() + cobra_model_I.reactions.get_by_id(reaction_id_I).get_reactants()
    met_names = [];
    mets_O = [];
    for m in met_ids:
        met_names.append(m.name);
    met_O = [k for k,v in Counter(met_names).items() if v>1]
    return met_O;

def load_thermoModel(anoxic = False):
    ijo1366_sbml = "data\\iJO1366.xml"
    # Read in the sbml file and define the model conditions
    cobra_model = create_cobra_model_from_sbml_file(ijo1366_sbml, print_time=True)
    # Update AMPMS2
    coc = Metabolite('co_c','CO','carbon monoxide','c');
    coc.charge = 0;
    cop = Metabolite('co_p','CO','carbon monoxide','p');
    cop.charge = 0;
    coe = Metabolite('co_e','CO','carbon monoxide','e');
    coe.charge = 0;
    cobra_model.add_metabolites([coc,cop,coe])
    ampms2_mets = {};
    ampms2_mets[cobra_model.metabolites.get_by_id('air_c')] = -1;
    ampms2_mets[cobra_model.metabolites.get_by_id('amet_c')] = -1;
    ampms2_mets[cobra_model.metabolites.get_by_id('dad_DASH_5_c')] = 1;
    ampms2_mets[cobra_model.metabolites.get_by_id('met_DASH_L_c')] = 1;
    ampms2_mets[cobra_model.metabolites.get_by_id('4ampm_c')] = 1;
    ampms2_mets[cobra_model.metabolites.get_by_id('h_c')] = 3;
    ampms2_mets[cobra_model.metabolites.get_by_id('for_c')] = 1;
    ampms2_mets[cobra_model.metabolites.get_by_id('co_c')] = 1;
    ampms2 = Reaction('AMPMS3');
    ampms2.add_metabolites(ampms2_mets);
    copp_mets = {};
    copp_mets[cobra_model.metabolites.get_by_id('co_c')] = -1;
    copp_mets[cobra_model.metabolites.get_by_id('co_p')] = 1;
    copp = Reaction('COtpp');
    copp.add_metabolites(copp_mets);
    coex_mets = {};
    coex_mets[cobra_model.metabolites.get_by_id('co_p')] = -1;
    coex_mets[cobra_model.metabolites.get_by_id('co_e')] = 1;
    coex = Reaction('COtex');
    coex.add_metabolites(coex_mets);
    cotrans_mets = {};
    cotrans_mets[cobra_model.metabolites.get_by_id('co_e')] = -1;
    cotrans = Reaction('EX_co_LPAREN_e_RPAREN_');
    cotrans.add_metabolites(cotrans_mets);
    cobra_model.add_reactions([ampms2,copp,coex,cotrans]);
    cobra_model.remove_reactions(['AMPMS2']);
    # Define the model conditions:
    system_boundaries = [x.id for x in cobra_model.reactions if x.boundary == 'system_boundary'];
    for b in system_boundaries:
            cobra_model.reactions.get_by_id(b).lower_bound = 0.0;
            cobra_model.reactions.get_by_id(b).upper_bound = 0.0;
    # Reset demand reactions
    demand = ['DM_4CRSOL',
            'DM_5DRIB',
            'DM_AACALD',
            'DM_AMOB',
            'DM_MTHTHF',
            'DM_OXAM'];
    for d in demand:
            cobra_model.reactions.get_by_id(d).lower_bound = 0.0;
            cobra_model.reactions.get_by_id(d).upper_bound = 1000.0;
    # Change the objective
    update_objective(cobra_model,{'Ec_biomass_iJO1366_WT_53p95M':1.0})
    # Assign KOs

    # Specify media composition (M9 glucose):
    cobra_model.reactions.get_by_id('EX_glc_LPAREN_e_RPAREN_').lower_bound = -10.0;
    cobra_model.reactions.get_by_id('EX_o2_LPAREN_e_RPAREN_').lower_bound = -18.0;
    uptake = ['EX_cl_LPAREN_e_RPAREN_',
                'EX_so4_LPAREN_e_RPAREN_',
                'EX_ca2_LPAREN_e_RPAREN_',
                'EX_pi_LPAREN_e_RPAREN_',
                'EX_fe2_LPAREN_e_RPAREN_',
                'EX_cu2_LPAREN_e_RPAREN_',
                'EX_zn2_LPAREN_e_RPAREN_',
                'EX_cbl1_LPAREN_e_RPAREN_',
                'EX_mobd_LPAREN_e_RPAREN_',
                'EX_ni2_LPAREN_e_RPAREN_',
                'EX_mn2_LPAREN_e_RPAREN_',
                'EX_k_LPAREN_e_RPAREN_',
                'EX_nh4_LPAREN_e_RPAREN_',
                'EX_cobalt2_LPAREN_e_RPAREN_',
                'EX_mg2_LPAREN_e_RPAREN_'];
    for u in uptake:
        cobra_model.reactions.get_by_id(u).lower_bound = -1000.0;
    # Specify allowed secretion products
    secrete = ['EX_meoh_LPAREN_e_RPAREN_',
                'EX_5mtr_LPAREN_e_RPAREN_',
                'EX_h_LPAREN_e_RPAREN_',
                'EX_co2_LPAREN_e_RPAREN_',
                'EX_co_LPAREN_e_RPAREN_',
                'EX_h2o_LPAREN_e_RPAREN_',
                'EX_ac_LPAREN_e_RPAREN_',
                'EX_fum_LPAREN_e_RPAREN_',
                'EX_for_LPAREN_e_RPAREN_',
                'EX_etoh_LPAREN_e_RPAREN_',
                'EX_lac_DASH_L_LPAREN_e_RPAREN_',
                'EX_pyr_LPAREN_e_RPAREN_',
                'EX_succ_LPAREN_e_RPAREN_'];
    for s in secrete:
        cobra_model.reactions.get_by_id(s).upper_bound = 1000.0;
    # Constrain specific reactions
    noFlux = ['F6PA', 'DHAPT'];
    ammoniaExcess = ['GLUDy'];
    aerobic = ['RNDR1', 'RNDR2', 'RNDR3', 'RNDR4', 'DHORD2', 'ASPO6','LCARR'];
    anaerobic = ['RNTR1c2', 'RNTR2c2', 'RNTR3c2', 'RNTR4c2', 'DHORD5', 'ASPO5'];
    if anaerobic:
        rxnList = noFlux + ammoniaExcess + anaerobic;
        for rxn in rxnList:
            cobra_model.reactions.get_by_id(rxn).lower_bound = 0.0;
            cobra_model.reactions.get_by_id(rxn).upper_bound = 0.0;
    else:
        rxnList = noFlux + ammoniaExcess + aerobic;
        for rxn in rxnList:
            cobra_model.reactions.get_by_id(rxn).lower_bound = 0.0;
            cobra_model.reactions.get_by_id(rxn).upper_bound = 0.0;
    ## Set the direction for specific reactions
    #fattyAcidSynthesis = ['ACCOAC', 'ACOATA', 'HACD1', 'HACD2', 'HACD3', 'HACD4', 'HACD5', 'HACD6', 'HACD7', 'HACD8', 'KAS14', 'KAS15', 'MACPD', 'MCOATA', '3OAR100', '3OAR120', '3OAR121', '3OAR140', '3OAR141', '3OAR160', '3OAR161', '3OAR180', '3OAR181', '3OAR40', '3OAR60', '3OAR80']
    #fattyAcidOxidation = ['ACACT1r', 'ACACT2r', 'ACACT3r', 'ACACT4r', 'ACACT5r', 'ACACT6r', 'ACACT7r', 'ACACT8r', 'ACOAD1f', 'ACOAD2f', 'ACOAD3f', 'ACOAD4f', 'ACOAD5f', 'ACOAD6f', 'ACOAD7f', 'ACOAD8f', 'CTECOAI6', 'CTECOAI7', 'CTECOAI8', 'ECOAH1', 'ECOAH2', 'ECOAH3', 'ECOAH4', 'ECOAH5', 'ECOAH6', 'ECOAH7', 'ECOAH8']
    #ndpk = ['NDPK1','NDPK2','NDPK3','NDPK4','NDPK5','NDPK7','NDPK8'];
    #rxnList = fattyAcidSynthesis + fattyAcidOxidation;
    #for rxn in rxnList:
    #    cobra_model.reactions.get_by_id(rxn).lower_bound = 0.0;
    #    cobra_model.reactions.get_by_id(rxn).upper_bound = 1000.0;

    return cobra_model;

def simulate_thermoConstraints(cobra_model_I,reactions_id_I):
    '''simulate the effect of constraining a list of model reactions to
    thermodynamically determined directions'''
    # Input:
    #   cobra_model_I
    #   reactions_id_I = cobra model reaction ids
    # Output:
    #   gr_O = {'original':gr,
    #           'reaction_id 1':[gr, % change in growth],
    #           'reaction_id 2':[gr, % change in growth],...}

    gr_O = {};
    # determine the orginal growth rate of the model
    cobra_model_I.optimize();
    gr_original = cobra_model_I.solution.f;
    gr_O['original'] = gr_original;
    # iterate through each reaction
    for rxn in reactions_id_I:
        gr_O[rxn] = None;
        # constrain the upper reaction bounds of the model
        ub = cobra_model_I.reactions.get_by_id(rxn).upper_bound;
        cobra_model_I.reactions.get_by_id(rxn).upper_bound = 0.0;
        # simulate growth with the constraint
        cobra_model_I.optimize();
        gr = cobra_model_I.solution.f
        gr_O[rxn] = [gr, gr/gr_original*100];
        # reset reaction constraint
        cobra_model_I.reactions.get_by_id(rxn).upper_bound = ub;

    return gr_O;

def add_pykA(cobra_model_I):
    '''adds the reactions corresponding to the permiscuous activity of pykA to
    the model'''
    # Input:
    #   cobra_model = cobra model
    # Output
    #   cobra_model = cobra model (input model is manipulated in place)

    # add conversion of udp to utp
    pyka_mets = {};
    pyka_mets[cobra_model_I.metabolites.get_by_id('atp_c')] = -1;
    pyka_mets[cobra_model_I.metabolites.get_by_id('udp_c')] = -1;
    pyka_mets[cobra_model_I.metabolites.get_by_id('adp_c')] = 1;
    pyka_mets[cobra_model_I.metabolites.get_by_id('utp_c')] = 1;
    pyka = Reaction('pyk_utp');
    pyka.add_metabolites(pyka_mets);
    cobra_model_I.add_reactions([pyka]);
    # add conversion of cdp to ctp
    pyka_mets = {};
    pyka_mets[cobra_model_I.metabolites.get_by_id('atp_c')] = -1;
    pyka_mets[cobra_model_I.metabolites.get_by_id('cdp_c')] = -1;
    pyka_mets[cobra_model_I.metabolites.get_by_id('adp_c')] = 1;
    pyka_mets[cobra_model_I.metabolites.get_by_id('ctp_c')] = 1;
    pyka = Reaction('pyk_ctp');
    pyka.add_metabolites(pyka_mets);
    cobra_model_I.add_reactions([pyka]);
    # add conversion of gdp to gtp
    pyka_mets = {};
    pyka_mets[cobra_model_I.metabolites.get_by_id('atp_c')] = -1;
    pyka_mets[cobra_model_I.metabolites.get_by_id('gdp_c')] = -1;
    pyka_mets[cobra_model_I.metabolites.get_by_id('adp_c')] = 1;
    pyka_mets[cobra_model_I.metabolites.get_by_id('gtp_c')] = 1;
    pyka = Reaction('pyk_gtp');
    pyka.add_metabolites(pyka_mets);
    cobra_model_I.add_reactions([pyka]);
    # add conversion of dtdp to dttp
    pyka_mets = {};
    pyka_mets[cobra_model_I.metabolites.get_by_id('atp_c')] = -1;
    pyka_mets[cobra_model_I.metabolites.get_by_id('dtdp_c')] = -1;
    pyka_mets[cobra_model_I.metabolites.get_by_id('adp_c')] = 1;
    pyka_mets[cobra_model_I.metabolites.get_by_id('dttp_c')] = 1;
    pyka = Reaction('pyk_dttp');
    pyka.add_metabolites(pyka_mets);
    cobra_model_I.add_reactions([pyka]);