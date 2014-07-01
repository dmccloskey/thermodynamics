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

def find_transportRxns(cobra_model_I):
    # transport reaction definition:
    # 1. the metabolite satisfies the critera for a transport metabolite
    #	1. different id (same base id but different compartment)
    #	2. same name
    #	3. different compartment
    # 2. the reaction is not a system boundary reaction

    rxn_O = [];
    for rxn in cobra_model_I.reactions:
        mets = [];
        mets = find_transportMets(cobra_model_I, rxn.id);
        #products = rxn.get_products();
        #reactants = rxn.get_reactants();
        #product_names = [x.name for x in products];
        #reactant_names = [x.name for x in reactants];
        if mets:
            rxn_O.append(rxn.id);
    return rxn_O;

def load_thermoModel(anoxic = False):
    ijo1366_sbml = "data\\iJO1366.xml"
    # Read in the sbml file and define the model conditions
    cobra_model = create_cobra_model_from_sbml_file(ijo1366_sbml, print_time=True)
    ## Update AMPMS2
    #coc = Metabolite('co_c','CO','carbon monoxide','c');
    #coc.charge = 0;
    #cop = Metabolite('co_p','CO','carbon monoxide','p');
    #cop.charge = 0;
    #coe = Metabolite('co_e','CO','carbon monoxide','e');
    #coe.charge = 0;
    #cobra_model.add_metabolites([coc,cop,coe])
    #ampms2_mets = {};
    #ampms2_mets[cobra_model.metabolites.get_by_id('air_c')] = -1;
    #ampms2_mets[cobra_model.metabolites.get_by_id('amet_c')] = -1;
    #ampms2_mets[cobra_model.metabolites.get_by_id('dad_DASH_5_c')] = 1;
    #ampms2_mets[cobra_model.metabolites.get_by_id('met_DASH_L_c')] = 1;
    #ampms2_mets[cobra_model.metabolites.get_by_id('4ampm_c')] = 1;
    #ampms2_mets[cobra_model.metabolites.get_by_id('h_c')] = 3;
    #ampms2_mets[cobra_model.metabolites.get_by_id('for_c')] = 1;
    #ampms2_mets[cobra_model.metabolites.get_by_id('co_c')] = 1;
    #ampms2 = Reaction('AMPMS3');
    #ampms2.add_metabolites(ampms2_mets);
    #copp_mets = {};
    #copp_mets[cobra_model.metabolites.get_by_id('co_c')] = -1;
    #copp_mets[cobra_model.metabolites.get_by_id('co_p')] = 1;
    #copp = Reaction('COtpp');
    #copp.add_metabolites(copp_mets);
    #coex_mets = {};
    #coex_mets[cobra_model.metabolites.get_by_id('co_p')] = -1;
    #coex_mets[cobra_model.metabolites.get_by_id('co_e')] = 1;
    #coex = Reaction('COtex');
    #coex.add_metabolites(coex_mets);
    #cotrans_mets = {};
    #cotrans_mets[cobra_model.metabolites.get_by_id('co_e')] = -1;
    #cotrans = Reaction('EX_co_LPAREN_e_RPAREN_');
    #cotrans.add_metabolites(cotrans_mets);
    #cobra_model.add_reactions([ampms2,copp,coex,cotrans]);
    #cobra_model.remove_reactions(['AMPMS2']);
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
    cobra_model.remove_reactions(['Ec_biomass_iJO1366_core_53p95M'])
    # Assign KOs

    # Specify media composition (M9 glucose):
    if anoxic:
        cobra_model.reactions.get_by_id('EX_glc_LPAREN_e_RPAREN_').lower_bound = -18.0;
        cobra_model.reactions.get_by_id('EX_o2_LPAREN_e_RPAREN_').lower_bound = 0.0;
    else:
        cobra_model.reactions.get_by_id('EX_glc_LPAREN_e_RPAREN_').lower_bound = -10.0;
        cobra_model.reactions.get_by_id('EX_o2_LPAREN_e_RPAREN_').lower_bound = -18.0;
    #uptake = ['EX_cl_LPAREN_e_RPAREN_',
    #            'EX_h_LPAREN_e_RPAREN_',
    #            'EX_so4_LPAREN_e_RPAREN_',
    #            'EX_ca2_LPAREN_e_RPAREN_',
    #            'EX_pi_LPAREN_e_RPAREN_',
    #            'EX_fe2_LPAREN_e_RPAREN_',
    #            'EX_cu2_LPAREN_e_RPAREN_',
    #            'EX_zn2_LPAREN_e_RPAREN_',
    #            'EX_cbl1_LPAREN_e_RPAREN_',
    #            'EX_mobd_LPAREN_e_RPAREN_',
    #            'EX_ni2_LPAREN_e_RPAREN_',
    #            'EX_mn2_LPAREN_e_RPAREN_',
    #            'EX_k_LPAREN_e_RPAREN_',
    #            'EX_nh4_LPAREN_e_RPAREN_',
    #            'EX_cobalt2_LPAREN_e_RPAREN_',
    #            'EX_mg2_LPAREN_e_RPAREN_'];
    uptake = ['EX_ca2_LPAREN_e_RPAREN_',
                'EX_cbl1_LPAREN_e_RPAREN_',
                'EX_cl_LPAREN_e_RPAREN_',
                'EX_co2_LPAREN_e_RPAREN_',
                'EX_cobalt2_LPAREN_e_RPAREN_',
                'EX_cu2_LPAREN_e_RPAREN_',
                'EX_fe2_LPAREN_e_RPAREN_',
                'EX_fe3_LPAREN_e_RPAREN_',
                'EX_h_LPAREN_e_RPAREN_',
                'EX_h2o_LPAREN_e_RPAREN_',
                'EX_k_LPAREN_e_RPAREN_',
                'EX_mg2_LPAREN_e_RPAREN_',
                'EX_mn2_LPAREN_e_RPAREN_',
                'EX_mobd_LPAREN_e_RPAREN_',
                'EX_na1_LPAREN_e_RPAREN_',
                'EX_nh4_LPAREN_e_RPAREN_',
                'EX_ni2_LPAREN_e_RPAREN_',
                'EX_pi_LPAREN_e_RPAREN_',
                'EX_sel_LPAREN_e_RPAREN_',
                'EX_slnt_LPAREN_e_RPAREN_',
                'EX_so4_LPAREN_e_RPAREN_',
                'EX_tungs_LPAREN_e_RPAREN_',
                'EX_zn2_LPAREN_e_RPAREN_'];
    for u in uptake:
        cobra_model.reactions.get_by_id(u).lower_bound = -1000.0;
    # Specify allowed secretion products
    secrete = ['EX_meoh_LPAREN_e_RPAREN_',
                'EX_5mtr_LPAREN_e_RPAREN_',
                'EX_h_LPAREN_e_RPAREN_',
                'EX_co2_LPAREN_e_RPAREN_',
                #'EX_co_LPAREN_e_RPAREN_',
                'EX_h2o_LPAREN_e_RPAREN_',
                'EX_ac_LPAREN_e_RPAREN_',
                #'EX_glyclt_LPAREN_e_RPAREN_',
                'EX_fum_LPAREN_e_RPAREN_',
                'EX_for_LPAREN_e_RPAREN_',
                'EX_etoh_LPAREN_e_RPAREN_',
                'EX_lac_DASH_L_LPAREN_e_RPAREN_',
                'EX_pyr_LPAREN_e_RPAREN_',
                'EX_succ_LPAREN_e_RPAREN_'];
    #secrete = ['EX_12ppd_DASH_R_LPAREN_e_RPAREN_',
    #            'EX_12ppd_DASH_S_LPAREN_e_RPAREN_',
    #            'EX_14glucan_LPAREN_e_RPAREN_',
    #            'EX_15dap_LPAREN_e_RPAREN_',
    #            'EX_23camp_LPAREN_e_RPAREN_',
    #            'EX_23ccmp_LPAREN_e_RPAREN_',
    #            'EX_23cgmp_LPAREN_e_RPAREN_',
    #            'EX_23cump_LPAREN_e_RPAREN_',
    #            'EX_23dappa_LPAREN_e_RPAREN_',
    #            'EX_26dap_DASH_M_LPAREN_e_RPAREN_',
    #            'EX_2ddglcn_LPAREN_e_RPAREN_',
    #            'EX_34dhpac_LPAREN_e_RPAREN_',
    #            'EX_3amp_LPAREN_e_RPAREN_',
    #            'EX_3cmp_LPAREN_e_RPAREN_',
    #            'EX_3gmp_LPAREN_e_RPAREN_',
    #            'EX_3hcinnm_LPAREN_e_RPAREN_',
    #            'EX_3hpp_LPAREN_e_RPAREN_',
    #            'EX_3hpppn_LPAREN_e_RPAREN_',
    #            'EX_3ump_LPAREN_e_RPAREN_',
    #            'EX_4abut_LPAREN_e_RPAREN_',
    #            'EX_4hoxpacd_LPAREN_e_RPAREN_',
    #            'EX_5dglcn_LPAREN_e_RPAREN_',
    #            'EX_5mtr_LPAREN_e_RPAREN_',
    #            'EX_LalaDglu_LPAREN_e_RPAREN_',
    #            'EX_LalaDgluMdap_LPAREN_e_RPAREN_',
    #            'EX_LalaDgluMdapDala_LPAREN_e_RPAREN_',
    #            'EX_LalaLglu_LPAREN_e_RPAREN_',
    #            'EX_ac_LPAREN_e_RPAREN_',
    #            'EX_acac_LPAREN_e_RPAREN_',
    #            'EX_acald_LPAREN_e_RPAREN_',
    #            'EX_acgal_LPAREN_e_RPAREN_',
    #            'EX_acgal1p_LPAREN_e_RPAREN_',
    #            'EX_acgam_LPAREN_e_RPAREN_',
    #            'EX_acgam1p_LPAREN_e_RPAREN_',
    #            'EX_acmana_LPAREN_e_RPAREN_',
    #            'EX_acmum_LPAREN_e_RPAREN_',
    #            'EX_acnam_LPAREN_e_RPAREN_',
    #            'EX_acolipa_LPAREN_e_RPAREN_',
    #            'EX_acser_LPAREN_e_RPAREN_',
    #            'EX_ade_LPAREN_e_RPAREN_',
    #            'EX_adn_LPAREN_e_RPAREN_',
    #            'EX_adocbl_LPAREN_e_RPAREN_',
    #            'EX_ag_LPAREN_e_RPAREN_',
    #            'EX_agm_LPAREN_e_RPAREN_',
    #            'EX_akg_LPAREN_e_RPAREN_',
    #            'EX_ala_DASH_B_LPAREN_e_RPAREN_',
    #            'EX_ala_DASH_D_LPAREN_e_RPAREN_',
    #            'EX_ala_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_alaala_LPAREN_e_RPAREN_',
    #            'EX_all_DASH_D_LPAREN_e_RPAREN_',
    #            'EX_alltn_LPAREN_e_RPAREN_',
    #            'EX_amp_LPAREN_e_RPAREN_',
    #            'EX_anhgm_LPAREN_e_RPAREN_',
    #            'EX_arab_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_arbt_LPAREN_e_RPAREN_',
    #            'EX_arbtn_LPAREN_e_RPAREN_',
    #            'EX_arbtn_DASH_fe3_LPAREN_e_RPAREN_',
    #            'EX_arg_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_ascb_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_asn_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_aso3_LPAREN_e_RPAREN_',
    #            'EX_asp_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_btn_LPAREN_e_RPAREN_',
    #            'EX_but_LPAREN_e_RPAREN_',
    #            'EX_butso3_LPAREN_e_RPAREN_',
    #            'EX_ca2_LPAREN_e_RPAREN_',
    #            'EX_cbi_LPAREN_e_RPAREN_',
    #            'EX_cbl1_LPAREN_e_RPAREN_',
    #            'EX_cd2_LPAREN_e_RPAREN_',
    #            'EX_cgly_LPAREN_e_RPAREN_',
    #            'EX_chol_LPAREN_e_RPAREN_',
    #            'EX_chtbs_LPAREN_e_RPAREN_',
    #            'EX_cit_LPAREN_e_RPAREN_',
    #            'EX_cl_LPAREN_e_RPAREN_',
    #            'EX_cm_LPAREN_e_RPAREN_',
    #            'EX_cmp_LPAREN_e_RPAREN_',
    #            'EX_co2_LPAREN_e_RPAREN_',
    #            'EX_cobalt2_LPAREN_e_RPAREN_',
    #            'EX_colipa_LPAREN_e_RPAREN_',
    #            'EX_colipap_LPAREN_e_RPAREN_',
    #            'EX_cpgn_LPAREN_e_RPAREN_',
    #            'EX_cpgn_DASH_un_LPAREN_e_RPAREN_',
    #            'EX_crn_LPAREN_e_RPAREN_',
    #            'EX_crn_DASH_D_LPAREN_e_RPAREN_',
    #            'EX_csn_LPAREN_e_RPAREN_',
    #            'EX_cu_LPAREN_e_RPAREN_',
    #            'EX_cu2_LPAREN_e_RPAREN_',
    #            'EX_cyan_LPAREN_e_RPAREN_',
    #            'EX_cynt_LPAREN_e_RPAREN_',
    #            'EX_cys_DASH_D_LPAREN_e_RPAREN_',
    #            'EX_cys_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_cytd_LPAREN_e_RPAREN_',
    #            'EX_dad_DASH_2_LPAREN_e_RPAREN_',
    #            'EX_damp_LPAREN_e_RPAREN_',
    #            'EX_dca_LPAREN_e_RPAREN_',
    #            'EX_dcmp_LPAREN_e_RPAREN_',
    #            'EX_dcyt_LPAREN_e_RPAREN_',
    #            'EX_ddca_LPAREN_e_RPAREN_',
    #            'EX_dgmp_LPAREN_e_RPAREN_',
    #            'EX_dgsn_LPAREN_e_RPAREN_',
    #            'EX_dha_LPAREN_e_RPAREN_',
    #            'EX_dimp_LPAREN_e_RPAREN_',
    #            'EX_din_LPAREN_e_RPAREN_',
    #            'EX_dms_LPAREN_e_RPAREN_',
    #            'EX_dmso_LPAREN_e_RPAREN_',
    #            'EX_dopa_LPAREN_e_RPAREN_',
    #            'EX_doxrbcn_LPAREN_e_RPAREN_',
    #            'EX_dtmp_LPAREN_e_RPAREN_',
    #            'EX_dump_LPAREN_e_RPAREN_',
    #            'EX_duri_LPAREN_e_RPAREN_',
    #            'EX_eca4colipa_LPAREN_e_RPAREN_',
    #            'EX_enlipa_LPAREN_e_RPAREN_',
    #            'EX_enter_LPAREN_e_RPAREN_',
    #            'EX_etha_LPAREN_e_RPAREN_',
    #            'EX_ethso3_LPAREN_e_RPAREN_',
    #            'EX_etoh_LPAREN_e_RPAREN_',
    #            'EX_f6p_LPAREN_e_RPAREN_',
    #            'EX_fald_LPAREN_e_RPAREN_',
    #            'EX_fe2_LPAREN_e_RPAREN_',
    #            'EX_fe3_LPAREN_e_RPAREN_',
    #            'EX_fe3dcit_LPAREN_e_RPAREN_',
    #            'EX_fe3dhbzs_LPAREN_e_RPAREN_',
    #            'EX_fe3hox_LPAREN_e_RPAREN_',
    #            'EX_fe3hox_DASH_un_LPAREN_e_RPAREN_',
    #            'EX_fecrm_LPAREN_e_RPAREN_',
    #            'EX_fecrm_DASH_un_LPAREN_e_RPAREN_',
    #            'EX_feenter_LPAREN_e_RPAREN_',
    #            'EX_feoxam_LPAREN_e_RPAREN_',
    #            'EX_feoxam_DASH_un_LPAREN_e_RPAREN_',
    #            'EX_for_LPAREN_e_RPAREN_',
    #            'EX_fru_LPAREN_e_RPAREN_',
    #            'EX_frulys_LPAREN_e_RPAREN_',
    #            'EX_fruur_LPAREN_e_RPAREN_',
    #            'EX_fuc_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_fum_LPAREN_e_RPAREN_',
    #            'EX_fusa_LPAREN_e_RPAREN_',
    #            'EX_g1p_LPAREN_e_RPAREN_',
    #            'EX_g3pc_LPAREN_e_RPAREN_',
    #            'EX_g3pe_LPAREN_e_RPAREN_',
    #            'EX_g3pg_LPAREN_e_RPAREN_',
    #            'EX_g3pi_LPAREN_e_RPAREN_',
    #            'EX_g3ps_LPAREN_e_RPAREN_',
    #            'EX_g6p_LPAREN_e_RPAREN_',
    #            'EX_gal_LPAREN_e_RPAREN_',
    #            'EX_gal_DASH_bD_LPAREN_e_RPAREN_',
    #            'EX_gal1p_LPAREN_e_RPAREN_',
    #            'EX_galct_DASH_D_LPAREN_e_RPAREN_',
    #            'EX_galctn_DASH_D_LPAREN_e_RPAREN_',
    #            'EX_galctn_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_galt_LPAREN_e_RPAREN_',
    #            'EX_galur_LPAREN_e_RPAREN_',
    #            'EX_gam_LPAREN_e_RPAREN_',
    #            'EX_gam6p_LPAREN_e_RPAREN_',
    #            'EX_gbbtn_LPAREN_e_RPAREN_',
    #            'EX_gdp_LPAREN_e_RPAREN_',
    #            'EX_glc_LPAREN_e_RPAREN_',
    #            'EX_glcn_LPAREN_e_RPAREN_',
    #            'EX_glcr_LPAREN_e_RPAREN_',
    #            'EX_glcur_LPAREN_e_RPAREN_',
    #            'EX_glcur1p_LPAREN_e_RPAREN_',
    #            'EX_gln_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_glu_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_gly_LPAREN_e_RPAREN_',
    #            'EX_glyald_LPAREN_e_RPAREN_',
    #            'EX_glyb_LPAREN_e_RPAREN_',
    #            'EX_glyc_LPAREN_e_RPAREN_',
    #            'EX_glyc_DASH_R_LPAREN_e_RPAREN_',
    #            'EX_glyc2p_LPAREN_e_RPAREN_',
    #            'EX_glyc3p_LPAREN_e_RPAREN_',
    #            'EX_glyclt_LPAREN_e_RPAREN_',
    #            'EX_gmp_LPAREN_e_RPAREN_',
    #            'EX_gsn_LPAREN_e_RPAREN_',
    #            'EX_gthox_LPAREN_e_RPAREN_',
    #            'EX_gthrd_LPAREN_e_RPAREN_',
    #            'EX_gtp_LPAREN_e_RPAREN_',
    #            'EX_gua_LPAREN_e_RPAREN_',
    #            'EX_h_LPAREN_e_RPAREN_',
    #            'EX_h2_LPAREN_e_RPAREN_',
    #            'EX_h2o_LPAREN_e_RPAREN_',
    #            'EX_h2o2_LPAREN_e_RPAREN_',
    #            'EX_h2s_LPAREN_e_RPAREN_',
    #            'EX_hacolipa_LPAREN_e_RPAREN_',
    #            'EX_halipa_LPAREN_e_RPAREN_',
    #            'EX_hdca_LPAREN_e_RPAREN_',
    #            'EX_hdcea_LPAREN_e_RPAREN_',
    #            'EX_hg2_LPAREN_e_RPAREN_',
    #            'EX_his_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_hom_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_hxa_LPAREN_e_RPAREN_',
    #            'EX_hxan_LPAREN_e_RPAREN_',
    #            'EX_idon_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_ile_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_imp_LPAREN_e_RPAREN_',
    #            'EX_indole_LPAREN_e_RPAREN_',
    #            'EX_inost_LPAREN_e_RPAREN_',
    #            'EX_ins_LPAREN_e_RPAREN_',
    #            'EX_isetac_LPAREN_e_RPAREN_',
    #            'EX_k_LPAREN_e_RPAREN_',
    #            'EX_kdo2lipid4_LPAREN_e_RPAREN_',
    #            'EX_lac_DASH_D_LPAREN_e_RPAREN_',
    #            'EX_lac_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_lcts_LPAREN_e_RPAREN_',
    #            'EX_leu_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_lipa_LPAREN_e_RPAREN_',
    #            'EX_lipa_cold_LPAREN_e_RPAREN_',
    #            'EX_lipoate_LPAREN_e_RPAREN_',
    #            'EX_lys_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_lyx_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_mal_DASH_D_LPAREN_e_RPAREN_',
    #            'EX_mal_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_malt_LPAREN_e_RPAREN_',
    #            'EX_malthx_LPAREN_e_RPAREN_',
    #            'EX_maltpt_LPAREN_e_RPAREN_',
    #            'EX_malttr_LPAREN_e_RPAREN_',
    #            'EX_maltttr_LPAREN_e_RPAREN_',
    #            'EX_man_LPAREN_e_RPAREN_',
    #            'EX_man6p_LPAREN_e_RPAREN_',
    #            'EX_manglyc_LPAREN_e_RPAREN_',
    #            'EX_melib_LPAREN_e_RPAREN_',
    #            'EX_meoh_LPAREN_e_RPAREN_',
    #            'EX_met_DASH_D_LPAREN_e_RPAREN_',
    #            'EX_met_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_metsox_DASH_R_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_metsox_DASH_S_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_mg2_LPAREN_e_RPAREN_',
    #            'EX_mincyc_LPAREN_e_RPAREN_',
    #            'EX_minohp_LPAREN_e_RPAREN_',
    #            'EX_mmet_LPAREN_e_RPAREN_',
    #            'EX_mn2_LPAREN_e_RPAREN_',
    #            'EX_mnl_LPAREN_e_RPAREN_',
    #            'EX_mobd_LPAREN_e_RPAREN_',
    #            'EX_mso3_LPAREN_e_RPAREN_',
    #            'EX_n2o_LPAREN_e_RPAREN_',
    #            'EX_na1_LPAREN_e_RPAREN_',
    #            'EX_nac_LPAREN_e_RPAREN_',
    #            'EX_nh4_LPAREN_e_RPAREN_',
    #            'EX_ni2_LPAREN_e_RPAREN_',
    #            'EX_nmn_LPAREN_e_RPAREN_',
    #            'EX_no_LPAREN_e_RPAREN_',
    #            'EX_no2_LPAREN_e_RPAREN_',
    #            'EX_no3_LPAREN_e_RPAREN_',
    #            'EX_novbcn_LPAREN_e_RPAREN_',
    #            'EX_o16a4colipa_LPAREN_e_RPAREN_',
    #            'EX_o2_LPAREN_e_RPAREN_',
    #            'EX_o2s_LPAREN_e_RPAREN_',
    #            'EX_ocdca_LPAREN_e_RPAREN_',
    #            'EX_ocdcea_LPAREN_e_RPAREN_',
    #            'EX_octa_LPAREN_e_RPAREN_',
    #            'EX_orn_LPAREN_e_RPAREN_',
    #            'EX_orot_LPAREN_e_RPAREN_',
    #            'EX_pacald_LPAREN_e_RPAREN_',
    #            'EX_peamn_LPAREN_e_RPAREN_',
    #            'EX_phe_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_pheme_LPAREN_e_RPAREN_',
    #            'EX_pi_LPAREN_e_RPAREN_',
    #            'EX_pnto_DASH_R_LPAREN_e_RPAREN_',
    #            'EX_ppa_LPAREN_e_RPAREN_',
    #            'EX_ppal_LPAREN_e_RPAREN_',
    #            'EX_pppn_LPAREN_e_RPAREN_',
    #            'EX_ppt_LPAREN_e_RPAREN_',
    #            'EX_pro_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_progly_LPAREN_e_RPAREN_',
    #            'EX_psclys_LPAREN_e_RPAREN_',
    #            'EX_pser_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_ptrc_LPAREN_e_RPAREN_',
    #            'EX_pydam_LPAREN_e_RPAREN_',
    #            'EX_pydx_LPAREN_e_RPAREN_',
    #            'EX_pydxn_LPAREN_e_RPAREN_',
    #            'EX_pyr_LPAREN_e_RPAREN_',
    #            'EX_quin_LPAREN_e_RPAREN_',
    #            'EX_r5p_LPAREN_e_RPAREN_',
    #            'EX_rfamp_LPAREN_e_RPAREN_',
    #            'EX_rib_DASH_D_LPAREN_e_RPAREN_',
    #            'EX_rmn_LPAREN_e_RPAREN_',
    #            'EX_sbt_DASH_D_LPAREN_e_RPAREN_',
    #            'EX_sel_LPAREN_e_RPAREN_',
    #            'EX_ser_DASH_D_LPAREN_e_RPAREN_',
    #            'EX_ser_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_skm_LPAREN_e_RPAREN_',
    #            'EX_slnt_LPAREN_e_RPAREN_',
    #            'EX_so2_LPAREN_e_RPAREN_',
    #            'EX_so3_LPAREN_e_RPAREN_',
    #            'EX_so4_LPAREN_e_RPAREN_',
    #            'EX_spmd_LPAREN_e_RPAREN_',
    #            'EX_succ_LPAREN_e_RPAREN_',
    #            'EX_sucr_LPAREN_e_RPAREN_',
    #            'EX_sulfac_LPAREN_e_RPAREN_',
    #            'EX_tartr_DASH_D_LPAREN_e_RPAREN_',
    #            'EX_tartr_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_taur_LPAREN_e_RPAREN_',
    #            'EX_tcynt_LPAREN_e_RPAREN_',
    #            'EX_thm_LPAREN_e_RPAREN_',
    #            'EX_thr_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_thrp_LPAREN_e_RPAREN_',
    #            'EX_thym_LPAREN_e_RPAREN_',
    #            'EX_thymd_LPAREN_e_RPAREN_',
    #            'EX_tma_LPAREN_e_RPAREN_',
    #            'EX_tmao_LPAREN_e_RPAREN_',
    #            'EX_tre_LPAREN_e_RPAREN_',
    #            'EX_trp_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_tsul_LPAREN_e_RPAREN_',
    #            'EX_ttdca_LPAREN_e_RPAREN_',
    #            'EX_ttdcea_LPAREN_e_RPAREN_',
    #            'EX_ttrcyc_LPAREN_e_RPAREN_',
    #            'EX_tungs_LPAREN_e_RPAREN_',
    #            'EX_tym_LPAREN_e_RPAREN_',
    #            'EX_tyr_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_tyrp_LPAREN_e_RPAREN_',
    #            'EX_uacgam_LPAREN_e_RPAREN_',
    #            'EX_udpacgal_LPAREN_e_RPAREN_',
    #            'EX_udpg_LPAREN_e_RPAREN_',
    #            'EX_udpgal_LPAREN_e_RPAREN_',
    #            'EX_udpglcur_LPAREN_e_RPAREN_',
    #            'EX_ump_LPAREN_e_RPAREN_',
    #            'EX_ura_LPAREN_e_RPAREN_',
    #            'EX_urea_LPAREN_e_RPAREN_',
    #            'EX_uri_LPAREN_e_RPAREN_',
    #            'EX_val_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_xan_LPAREN_e_RPAREN_',
    #            'EX_xmp_LPAREN_e_RPAREN_',
    #            'EX_xtsn_LPAREN_e_RPAREN_',
    #            'EX_xyl_DASH_D_LPAREN_e_RPAREN_',
    #            'EX_xylu_DASH_L_LPAREN_e_RPAREN_',
    #            'EX_zn2_LPAREN_e_RPAREN_'];
    for s in secrete:
        cobra_model.reactions.get_by_id(s).upper_bound = 1000.0;
    # Constrain specific reactions
    noFlux = ['F6PA', 'DHAPT'];
    ammoniaExcess = ['GLUDy'];
    aerobic = ['RNDR1', 'RNDR2', 'RNDR3', 'RNDR4','RNDR1b', 'RNDR2b', 'RNDR3b', 'RNDR4b', 'DHORD2', 'ASPO6','LCARR'];
    anaerobic = ['RNTR1c2', 'RNTR2c2', 'RNTR3c2', 'RNTR4c2','RNDR1b', 'RNDR2b', 'RNDR3b', 'RNDR4b', 'DHORD5', 'ASPO5'];
    if anoxic:
        rxnList = noFlux + ammoniaExcess + aerobic;
        for rxn in rxnList:
            cobra_model.reactions.get_by_id(rxn).lower_bound = 0.0;
            cobra_model.reactions.get_by_id(rxn).upper_bound = 0.0;
    else:
        rxnList = noFlux + ammoniaExcess + anaerobic;
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