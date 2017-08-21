
# -*- coding: utf-8 -*-
from math import floor,ceil,log,sqrt,pow,exp,fabs
from copy import deepcopy
from collections import Counter

class thermodynamics_dG_p_data():
    """Runs pathway thermodynamic analysis analysis on a cobra.Model object

    #1 Calculate pathway thermodynamics

    cobra_model: a Model object (irreversible)
    
    dG0_r: calculated standard Gibbs energies of reaction:
                    reaction.id: {'dG_r': float,
                          'dG_r_var': float,
                          'dG_r_lb': float,
                          'dG_r_ub': float,
                          'Keq': float,
                          'dG_r_units': string}}

    dG_r: calculated in vivo Gibbs energies of reaction:
                    reaction.id: {'dG_r': float,
                          'dG_r_var': float,
                          'dG_r_lb': float,
                          'dG_r_ub': float,
                          'Keq': float,
                          'dG_r_units': string}}

    dG0_p: {pathway.id: {'dG_p': float,
                          'dG_p_var': float,
                          'dG_p_lb': float,
                          'dG_p_ub': float,
                          'Keq': float,
                          'dG_p_units': string,
                          'reactions':[string],
                          'stoichiometry':[int]}}

    dG_p: {pathway.id: {'dG_p': float,
                          'dG_p_var': float,
                          'dG_p_lb': float,
                          'dG_p_ub': float,
                          'Keq': float,
                          'dG_p_units': string},
                          'reactions':[string],
                          'stoichiometry':[int]}

    thermodynamic_consistency_check: {pathway.id: {'feasible': boolean, NOTE: or None if the below criterion
                                                                    were not met
                                         'measured_concentration_coverage': float,
                                         'measured_dG_p_coverage': float}
    """

    def __init__(self,pathways_I = {}, dG0_p_I = {}, dG_p_I = {}):
        # initialize pathways
        if pathways_I:
            self.pathways = pathways_I;
        else:
            self.pathways = self.init_pathways();

        # output
        if dG0_p_I:
            self.dG0_p = dG0_p_I;
        else:
            self.dG0_p = {};
        if dG_p_I:
            self.dG_p = dG_p_I;
        else:
            self.dG_p = {};
        self.dG_p_coverage = {}
        self.metabolomics_coverage = {}
        self.thermodynamic_consistency_check = {}

    def export_dG0_p_json(self, filename_I):
        # save the results to json file
        with open(filename_I, 'w') as outfile:
            json.dump(self.dG0_p, outfile, indent=4);

    def export_dG_p_json(self, filename_I):
        # save the results to json file
        with open(filename_I, 'w') as outfile:
            json.dump(self.dG_p, outfile, indent=4);

    def export_tcc_json(self, filename_I):
        # save the results to json file
        with open(filename_I, 'w') as outfile:
            json.dump(self.thermodynamic_consistency_check, outfile, indent=4);

    def calculate_dG_p(self, cobra_model, dG0_r, dG_r):
        """calculate the Gibbs free energy of an entire pathway

        Args:
            dG0_r
            dG_r

        Returns:
            dG0_p
            dG_p
            thermodynamic_consistency_check

        """

        dG0_p = {};
        dG_p = {};
        # iterate through each pathway
        for path,v in self.pathways.items():
            dG0_p[path] = None;
            dG_p[path] = None;
            dG0_p_dict = {};
            dG_p_dict = {};
            dG0_p_tmp = 0.0;
            dG0_p_var_tmp = 0.0;
            dG0_p_lb_tmp = 0.0;
            dG0_p_ub_tmp = 0.0;
            dG_p_tmp = 0.0;
            dG_p_var_tmp = 0.0;
            dG_p_lb_tmp = 0.0;
            dG_p_ub_tmp = 0.0;
            # iterate through each reaction
            for i,rxn in enumerate(v['reactions']):
                if rxn in dG0_r:
                    # dG0 calculations
                    dG0_p_tmp += dG0_r[rxn]['dG_r']*v['stoichiometry'][i];
                    dG0_p_var_tmp += dG0_r[rxn]['dG_r_var']*v['stoichiometry'][i];
                    # ensure bounds are correct
                    if dG0_r[rxn]['dG_r_lb']>dG0_r[rxn]['dG_r_ub']:
                        dG0_p_lb_tmp += dG0_r[rxn]['dG_r_ub']*v['stoichiometry'][i];
                        dG0_p_ub_tmp += dG0_r[rxn]['dG_r_lb']*v['stoichiometry'][i];
                    else:
                        dG0_p_lb_tmp += dG0_r[rxn]['dG_r_lb']*v['stoichiometry'][i];
                        dG0_p_ub_tmp += dG0_r[rxn]['dG_r_ub']*v['stoichiometry'][i];
                else: print(('dG0_r not calculated for reaction ' + rxn + '!'));
                if rxn in dG_r:
                    # dG calculations
                    dG_p_tmp += dG_r[rxn]['dG_r']*v['stoichiometry'][i];
                    dG_p_var_tmp += dG_r[rxn]['dG_r_var']*v['stoichiometry'][i];
                    # ensure bounds are correct
                    if dG_r[rxn]['dG_r_lb']>dG_r[rxn]['dG_r_ub']:
                        dG_p_lb_tmp += dG_r[rxn]['dG_r_ub']*v['stoichiometry'][i];
                        dG_p_ub_tmp += dG_r[rxn]['dG_r_lb']*v['stoichiometry'][i];
                    else:
                        dG_p_lb_tmp += dG_r[rxn]['dG_r_lb']*v['stoichiometry'][i];
                        dG_p_ub_tmp += dG_r[rxn]['dG_r_ub']*v['stoichiometry'][i];
                else: print(('dG0_r not calculated for reaction ' + rxn + '!'));
            # copy information into dG0
            dG0_p_dict['reactions'] = v['reactions'];
            dG0_p_dict['stoichiometry'] = v['stoichiometry'];
            dG0_p_dict['dG0_p'] = dG0_p_tmp;
            dG0_p_dict['dG0_p_var'] = dG0_p_var_tmp;
            dG0_p_dict['dG0_p_lb'] = dG0_p_lb_tmp;
            dG0_p_dict['dG0_p_ub'] = dG0_p_ub_tmp;
            dG0_p_dict['dG0_p_units'] = 'kJ/mol';
            dG0_p[path] = dG0_p_dict;
            # copy information into dG
            dG_p_dict['reactions'] = v['reactions'];
            dG_p_dict['stoichiometry'] = v['stoichiometry'];
            dG_p_dict['dG_p'] = dG_p_tmp;
            dG_p_dict['dG_p_var'] = dG_p_var_tmp;
            dG_p_dict['dG_p_lb'] = dG_p_lb_tmp;
            dG_p_dict['dG_p_ub'] = dG_p_ub_tmp;
            dG_p_dict['dG_p_units'] = 'kJ/mol';
            dG_p[path] = dG_p_dict;
        # return dG0 and dG for all pathways
        self.dG0_p = dG0_p;
        self.dG_p = dG_p;
    
    def init_pathways(self):
        # pathways:
        pathways_irreversible = {
        'ptrc_to_4abut_1':{'reactions':['PTRCTA','ABUTD'],
                           'stoichiometry':[1,1]},
        'ptrc_to_4abut_2':{'reactions':['GGPTRCS','GGPTRCO','GGGABADr','GGGABAH'],
                           'stoichiometry':[1,1,1,1]},
        'glu_DASH_L_to_acg5p':{'reactions':['ACGS','ACGK'],
                           'stoichiometry':[1,1]},
        '2obut_and_pyr_to_3mop':{'reactions':['ACHBS','KARA2','DHAD2'],
                           'stoichiometry':[1,1,1]},
        'pyr_to_23dhmb':{'reactions':['ACLS','KARA1_reverse'],
                           'stoichiometry':[1,1]},
        #'met_DASH_L_and_ptrc_to_spmd_and_5mta':{'reactions':['METAT','ADMDC','SPMS'],
        #                   'stoichiometry':[1,1,1]}, #cannot be lumped
        'chor_and_prpp_to_3ig3p':{'reactions':['ANS','ANPRT','PRAIi','IGPS'],
                           'stoichiometry':[1,1,1,1]},
        'hom_DASH_L_and_cyst_DASH_L_to_pyr_hcys_DASH_L':{'reactions':['HSST','SHSL1','CYSTL'],
                           'stoichiometry':[1,1,1]},
        'e4p_and_pep_to_3dhq':{'reactions':['DDPA','DHQS'],
                           'stoichiometry':[1,1]},
        'aspsa_to_sl2a6o':{'reactions':['DHDPS','DHDPRy','THDPS'],
                           'stoichiometry':[1,1,1]},
        'glu_DASH_L_to_glu5sa':{'reactions':['GLU5K','G5SD'],
                           'stoichiometry':[1,1]},
        'g1p_to_glycogen':{'reactions':['GLGC','GLCS1'],
                           'stoichiometry':[1,1]},
        'thr_DASH_L_to_gly':{'reactions':['THRD','GLYAT_reverse'],
                           'stoichiometry':[1,1]},
        'dhap_to_lac_DASH_D':{'reactions':['MGSA','LGTHL','GLYOX'],
                           'stoichiometry':[1,1,1]},
        'hom_DASH_L_to_thr_DASH_L':{'reactions':['HSK','THRS'],
                           'stoichiometry':[1,1]},
        '3pg_to_ser_DASH_L':{'reactions':['PGCD','PSERT','PSP_L'],
                           'stoichiometry':[1,1,1]},
        'prpp_to_his_DASH_L':{'reactions':['ATPPRT','PRATPP','PRAMPC','PRMICI','IG3PS','IGPDH','HSTPT','HISTP','HISTD'],
                           'stoichiometry':[1,1,1,1,1,1,1,1,1]},
        'UMPSYN_aerobic':{'reactions':['ASPCT','DHORTS_reverse','DHORD2','ORPT_reverse','OMPDC'],
                           'stoichiometry':[1,1,1,1,1]},
        'UMPSYN_anaerobic':{'reactions':['ASPCT','DHORTS_reverse','DHORD5','ORPT_reverse','OMPDC'],
                           'stoichiometry':[1,1,1,1,1]},
        'IMPSYN_1':{'reactions':['GLUPRT','PRAGSr','PRFGS','PRAIS'],
                           'stoichiometry':[1,1,1,1]},
        'IMPSYN_2':{'reactions':['AIRC2','AIRC3_reverse','PRASCSi','ADSL2r'],
                           'stoichiometry':[1,1,1,1]},
        'IMPSYN_3':{'reactions':['AICART','IMPC_reverse'],
                           'stoichiometry':[1,1]},
        'imp_to_gmp':{'reactions':['IMPD','GMPS2'],
                           'stoichiometry':[1,1]},
        'imp_to_amp':{'reactions':['ADSS','ADSL1r'],
                           'stoichiometry':[1,1]},
        'utp_to_dump_anaerobic':{'reactions':['RNTR4c2','DUTPDP'],
                           'stoichiometry':[1,1]},
        'udp_to_dump_aerobic':{'reactions':['RNDR4','NDPK6','DUTPDP'],
                           'stoichiometry':[1,1,1]},
        'dtmp_to_dttp':{'reactions':['DTMPK','NDPK4'],
                           'stoichiometry':[1,1]}, #cannot be lumped
        'COASYN':{'reactions':['ASP1DC','MOHMT','DPR','PANTS','PNTK','PPNCL2','PPCDC','PTPATi','DPCOAK'],
                           'stoichiometry':[1,1,1,1,1,1,1,1,1]},
        'FADSYN_1':{'reactions':['GTPCII2','DHPPDA2','APRAUR','PMDPHT','RBFSb'],
                           'stoichiometry':[1,1,1,1,1]},
        'FADSYN_2':{'reactions':['RBFSa','DB4PS'],
                           'stoichiometry':[1,1]},
        'FADSYN_3':{'reactions':['RBFK','FMNAT'],
                           'stoichiometry':[1,1]},
        'NADSYN_aerobic':{'reactions':['ASPO6','QULNS','NNDPR','NNATr','NADS1','NADK'],
                           'stoichiometry':[1,1,1,1,1,1]},
        'NADSYN_anaerobic':{'reactions':['ASPO5','QULNS','NNDPR','NNATr','NADS1','NADK'],
                           'stoichiometry':[1,1,1,1,1,1]},
        #'NADSALVAGE':{'reactions':['NADPPPS','NADN','NNAM','NAMNPP','NMNN','NMNDA','NMNAT','NADDP','ADPRDP'],
        #                   'stoichiometry':[1,1,1,1,1,1,1,1,1]}, #cannot be lumped
        'THFSYN':{'reactions':['GTPCI','DNTPPA','DNMPPA','DHNPA2r','HPPK2','ADCS','ADCL','DHPS2','DHFS'],
                           'stoichiometry':[1,1,1,1,1,1,1,1,1]},
        'GTHSYN':{'reactions':['GLUCYS','GTHS'],
                           'stoichiometry':[1,1]},
        'GLYCPHOSPHOLIPID_1':{'reactions':['DASYN181','AGPAT181','G3PAT181'],'stoichiometry':[1,1,1]},
        'GLYCPHOSPHOLIPID_2':{'reactions':['PSSA181','PSD181'],'stoichiometry':[1,1]},
        'GLYCPHOSPHOLIPID_3':{'reactions':['PGSA160','PGPP160'],'stoichiometry':[1,1]},
        'GLYCPHOSPHOLIPID_4':{'reactions':['DASYN161','AGPAT161','G3PAT161'],'stoichiometry':[1,1,1]},
        'GLYCPHOSPHOLIPID_5':{'reactions':['PGSA181','PGPP181'],'stoichiometry':[1,1]},
        'GLYCPHOSPHOLIPID_6':{'reactions':['PSD161','PSSA161'],'stoichiometry':[1,1]},
        'GLYCPHOSPHOLIPID_7':{'reactions':['PSSA160','PSD160'],'stoichiometry':[1,1]},
        'GLYCPHOSPHOLIPID_8':{'reactions':['DASYN160','AGPAT160','G3PAT160'],'stoichiometry':[1,1,1]},
        'GLYCPHOSPHOLIPID_9':{'reactions':['PGSA161','PGPP161'],'stoichiometry':[1,1]},
        'MOLYBDOPTERIN_1':{'reactions':['MPTAT','MPTS','CPMPS'],'stoichiometry':[1,1,1]},
        'MOLYBDOPTERIN_2':{'reactions':['MOCDS','MOGDS'],'stoichiometry':[1,1]},
        'MOLYBDOPTERIN_3':{'reactions':['MOADSUx','MPTSS'],'stoichiometry':[1,1]},
        'COFACTOR_1':{'reactions':['GLUTRR','G1SAT','GLUTRS'],'stoichiometry':[1,1,1]},
        'COFACTOR_2':{'reactions':['DHNAOT4','UPPDC1','DHNCOAT','DHNCOAS','SEPHCHCS','SUCBZS','SUCBZL','PPPGO3','FCLT','CPPPGO','SHCHCS3'],'stoichiometry':[1,1,1,1,1,1,1,1,1,1,1]},
        'COFACTOR_3':{'reactions':['TYRL','AMMQLT8','HEMEOS','UPP3MT','SHCHD2','SHCHF','ENTCS','CBLAT'],'stoichiometry':[1,1,1,1,1,1,1,1]},
        'VITB6':{'reactions':['E4PD','PERD','OHPBAT','PDX5PS','PDX5PO2'],'stoichiometry':[1,1,1,1,1]},
        #'THIAMIN':{'reactions':['AMPMS2','PMPK','THZPSN3','TMPPP','TMPK'],'stoichiometry':[1,1,1,1,1]}, # original pathway without correction
        'THIAMIN':{'reactions':['AMPMS3','PMPK','THZPSN3','TMPPP','TMPK'],'stoichiometry':[1,1,1,1,1]},
        'COFACTOR_4':{'reactions':['I4FE4ST','I4FE4SR','I2FE2SS2'],'stoichiometry':[1,1,1]},
        'COFACTOR_5':{'reactions':['BMOGDS1','BMOGDS2','BMOCOS'],'stoichiometry':[1,1,1]},
        'COFACTOR_6':{'reactions':['DMPPS','GRTT','DMATT'],'stoichiometry':[1,1,1]},
        'COFACTOR_7':{'reactions':['MECDPS','DXPRIi','MEPCT','CDPMEK','MECDPDH5'],'stoichiometry':[1,1,1,1,1]},
        'COFACTOR_8':{'reactions':['LIPOS','LIPOCT'],'stoichiometry':[1,1]},
        'COFACTOR_9':{'reactions':['OMMBLHX','OMPHHX','OPHHX','HBZOPT','DMQMT','CHRPL','OMBZLM','OPHBDC','OHPHM'],'stoichiometry':[1,1,1,1,1,1,1,1,1]},
        'COFACTOR_10':{'reactions':['SERASr','DHBD','UPP3S','HMBS','ICHORT','DHBS'],'stoichiometry':[1,1,1,1,1,1]},
        'COFACTOR_11':{'reactions':['PMEACPE','EGMEACPR','DBTS','AOXSr2','I2FE2SR','OPMEACPD','MALCOAMT','AMAOTr','OPMEACPS','OPMEACPR','OGMEACPD','OGMEACPR','OGMEACPS','EPMEACPR','BTS5'],'stoichiometry':[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]},
        'CELLENV_1':{'reactions':['UAMAGS','UAPGR','UAGPT3','PAPPT3','GLUR_reverse','UAGCVT','UAMAS','UDCPDP','UGMDDS','UAAGDS'],'stoichiometry':[1,1,1,1,1,1,1,1,1,1]},
        'CELLENV_2':{'reactions':['3HAD181','3OAR181','3OAS181','EAR181x'],'stoichiometry':[1,1,1,1]},
        'CELLENV_3':{'reactions':['3HAD160','3OAR160','EAR160x','3OAS160'],'stoichiometry':[1,1,1,1]},
        'CELLENV_4':{'reactions':['EAR120x','3OAR120','3HAD120','3OAS120','EAR100x'],'stoichiometry':[1,1,1,1,1]},
        'CELLENV_5':{'reactions':['G1PACT','UAGDP','PGAMT_reverse','GF6PTA'],'stoichiometry':[1,1,1,1]},
        'CELLENV_6':{'reactions':['3OAR40','EAR40x','3OAS60','3OAR60','3HAD80','3OAS80','3OAR80','EAR60x','3HAD60','EAR80x','3HAD40'],'stoichiometry':[1,1,1,1,1,1,1,1,1,1,1]},
        'CELLENV_7':{'reactions':['3HAD161','EAR161x','3OAS161','3OAR161','3OAS141','3HAD141','3OAR121','EAR121x','3HAD121','EAR141x','T2DECAI','3OAR141','3OAS121'],'stoichiometry':[1,1,1,1,1,1,1,1,1,1,1,1,1]},
        'CELLENV_8':{'reactions':['TDPGDH','TDPDRR','TDPDRE','G1PTT'],'stoichiometry':[1,1,1,1]},
        'CELLENV_9':{'reactions':['3OAS140','3OAR140'],'stoichiometry':[1,1]},
        'CELLENV_10':{'reactions':['3HAD140','EAR140x'],'stoichiometry':[1,1]},
        'CELLENV_11':{'reactions':['3OAR100','3HAD100','3OAS100'],'stoichiometry':[1,1,1]},
        'LIPOPOLYSACCHARIDE_1':{'reactions':['COLIPAabcpp','COLIPAabctex','EDTXS1','EDTXS2','GALT1','GLCTR1','GLCTR2','GLCTR3','HEPK1','HEPK2','HEPT1','HEPT2','HEPT3','HEPT4','LPADSS','MOAT','MOAT2','MOAT3C','RHAT1','TDSK','USHD'],'stoichiometry':[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]},
        'LIPOPOLYSACCHARIDE_2':{'reactions':['AGMHE','GMHEPAT','GMHEPK','GMHEPPA','S7PI'],'stoichiometry':[1,1,1,1,1]},
        'LIPOPOLYSACCHARIDE_3':{'reactions':['U23GAAT','UHGADA','UAGAAT'],'stoichiometry':[1,1,1]},
        'LIPOPOLYSACCHARIDE_4':{'reactions':['KDOPP','KDOCT2','KDOPS'],'stoichiometry':[1,1,1]},
        'ASTPathway':{'reactions':['AST','SADH','SGDS','SGSAD','SOTA'],'stoichiometry':[1,1,1,1,1]},
        'OPPP':{'reactions':['G6PDH2r', 'GND', 'PGL', 'RPI_reverse'],'stoichiometry':[1,1,1,1]},
        'Glycolysis':{'reactions':['PGI', 'PFK', 'FBA', 'TPI', 'GAPD', 'PGK_reverse', 'PGM_reverse', 'ENO'],
                              'stoichiometry':[1,1,1,1,1,1,1,1]},
        'ED_pathway':{'reactions':['G6PDH2r', 'PGL', 'EDD', 'EDA'],
                                      'stoichiometry':[1,1,1,1]},
        'ED_branch':{'reactions':['EDD', 'EDA'],
                                      'stoichiometry':[1,1]}
        };
        pathways = {
        'ptrc_to_4abut_1':{'reactions':['PTRCTA','ABUTD'],
                           'stoichiometry':[1,1]},
        'ptrc_to_4abut_2':{'reactions':['GGPTRCS','GGPTRCO','GGGABADr','GGGABAH'],
                           'stoichiometry':[1,1,1,1]},
        'glu_DASH_L_to_acg5p':{'reactions':['ACGS','ACGK'],
                           'stoichiometry':[1,1]},
        '2obut_and_pyr_to_3mop':{'reactions':['ACHBS','KARA2','DHAD2'],
                           'stoichiometry':[1,1,1]},
        'pyr_to_23dhmb':{'reactions':['ACLS','KARA1'],
                           'stoichiometry':[1,-1]},
        #'met_DASH_L_and_ptrc_to_spmd_and_5mta':{'reactions':['METAT','ADMDC','SPMS'],
        #                   'stoichiometry':[1,1,1]}, #cannot be lumped
        'chor_and_prpp_to_3ig3p':{'reactions':['ANS','ANPRT','PRAIi','IGPS'],
                           'stoichiometry':[1,1,1,1]},
        'hom_DASH_L_and_cyst_DASH_L_to_pyr_hcys_DASH_L':{'reactions':['HSST','SHSL1','CYSTL'],
                           'stoichiometry':[1,1,1]},
        'e4p_and_pep_to_3dhq':{'reactions':['DDPA','DHQS'],
                           'stoichiometry':[1,1]},
        'aspsa_to_sl2a6o':{'reactions':['DHDPS','DHDPRy','THDPS'],
                           'stoichiometry':[1,1,1]},
        'glu_DASH_L_to_glu5sa':{'reactions':['GLU5K','G5SD'],
                           'stoichiometry':[1,1]},
        'g1p_to_glycogen':{'reactions':['GLGC','GLCS1'],
                           'stoichiometry':[1,1]},
        'thr_DASH_L_to_gly':{'reactions':['THRD','GLYAT'],
                           'stoichiometry':[1,-1]},
        'dhap_to_lac_DASH_D':{'reactions':['MGSA','LGTHL','GLYOX'],
                           'stoichiometry':[1,1,1]},
        'hom_DASH_L_to_thr_DASH_L':{'reactions':['HSK','THRS'],
                           'stoichiometry':[1,1]},
        '3pg_to_ser_DASH_L':{'reactions':['PGCD','PSERT','PSP_L'],
                           'stoichiometry':[1,1,1]},
        'prpp_to_his_DASH_L':{'reactions':['ATPPRT','PRATPP','PRAMPC','PRMICI','IG3PS','IGPDH','HSTPT','HISTP','HISTD'],
                           'stoichiometry':[1,1,1,1,1,1,1,1,1]},
        'UMPSYN_aerobic':{'reactions':['ASPCT','DHORTS','DHORD2','ORPT','OMPDC'],
                           'stoichiometry':[1,-1,1,-1,1]},
        'UMPSYN_anaerobic':{'reactions':['ASPCT','DHORTS','DHORD5','ORPT','OMPDC'],
                           'stoichiometry':[1,-1,1,-1,1]},
        'IMPSYN_1':{'reactions':['GLUPRT','PRAGSr','PRFGS','PRAIS'],
                           'stoichiometry':[1,1,1,1]},
        'IMPSYN_2':{'reactions':['AIRC2','AIRC3','PRASCSi','ADSL2r'],
                           'stoichiometry':[1,-1,1,1]},
        'IMPSYN_3':{'reactions':['AICART','IMPC'],
                           'stoichiometry':[1,-1]},
        'imp_to_gmp':{'reactions':['IMPD','GMPS2'],
                           'stoichiometry':[1,1]},
        'imp_to_amp':{'reactions':['ADSS','ADSL1r'],
                           'stoichiometry':[1,1]},
        'utp_to_dump_anaerobic':{'reactions':['RNTR4c2','DUTPDP'],
                           'stoichiometry':[1,1]},
        'udp_to_dump_aerobic':{'reactions':['RNDR4','NDPK6','DUTPDP'],
                           'stoichiometry':[1,1,1]},
        'dtmp_to_dttp':{'reactions':['DTMPK','NDPK4'],
                           'stoichiometry':[1,1]}, #cannot be lumped
        'COASYN':{'reactions':['ASP1DC','MOHMT','DPR','PANTS','PNTK','PPNCL2','PPCDC','PTPATi','DPCOAK'],
                           'stoichiometry':[1,1,1,1,1,1,1,1,1]},
        'FADSYN_1':{'reactions':['GTPCII2','DHPPDA2','APRAUR','PMDPHT','RBFSb'],
                           'stoichiometry':[1,1,1,1,1]},
        'FADSYN_2':{'reactions':['RBFSa','DB4PS'],
                           'stoichiometry':[1,1]},
        'FADSYN_3':{'reactions':['RBFK','FMNAT'],
                           'stoichiometry':[1,1]},
        'NADSYN_aerobic':{'reactions':['ASPO6','QULNS','NNDPR','NNATr','NADS1','NADK'],
                           'stoichiometry':[1,1,1,1,1,1]},
        'NADSYN_anaerobic':{'reactions':['ASPO5','QULNS','NNDPR','NNATr','NADS1','NADK'],
                           'stoichiometry':[1,1,1,1,1,1]},
        #'NADSALVAGE':{'reactions':['NADPPPS','NADN','NNAM','NAMNPP','NMNN','NMNDA','NMNAT','NADDP','ADPRDP'],
        #                   'stoichiometry':[1,1,1,1,1,1,1,1,1]}, #cannot be lumped
        'THFSYN':{'reactions':['GTPCI','DNTPPA','DNMPPA','DHNPA2r','HPPK2','ADCS','ADCL','DHPS2','DHFS'],
                           'stoichiometry':[1,1,1,1,1,1,1,1,1]},
        'GTHSYN':{'reactions':['GLUCYS','GTHS'],
                           'stoichiometry':[1,1]},
        'GLYCPHOSPHOLIPID_1':{'reactions':['DASYN181','AGPAT181','G3PAT181'],'stoichiometry':[1,1,1]},
        'GLYCPHOSPHOLIPID_2':{'reactions':['PSSA181','PSD181'],'stoichiometry':[1,1]},
        'GLYCPHOSPHOLIPID_3':{'reactions':['PGSA160','PGPP160'],'stoichiometry':[1,1]},
        'GLYCPHOSPHOLIPID_4':{'reactions':['DASYN161','AGPAT161','G3PAT161'],'stoichiometry':[1,1,1]},
        'GLYCPHOSPHOLIPID_5':{'reactions':['PGSA181','PGPP181'],'stoichiometry':[1,1]},
        'GLYCPHOSPHOLIPID_6':{'reactions':['PSD161','PSSA161'],'stoichiometry':[1,1]},
        'GLYCPHOSPHOLIPID_7':{'reactions':['PSSA160','PSD160'],'stoichiometry':[1,1]},
        'GLYCPHOSPHOLIPID_8':{'reactions':['DASYN160','AGPAT160','G3PAT160'],'stoichiometry':[1,1,1]},
        'GLYCPHOSPHOLIPID_9':{'reactions':['PGSA161','PGPP161'],'stoichiometry':[1,1]},
        'MOLYBDOPTERIN_1':{'reactions':['MPTAT','MPTS','CPMPS'],'stoichiometry':[1,1,1]},
        'MOLYBDOPTERIN_2':{'reactions':['MOCDS','MOGDS'],'stoichiometry':[1,1]},
        'MOLYBDOPTERIN_3':{'reactions':['MOADSUx','MPTSS'],'stoichiometry':[1,1]},
        'COFACTOR_1':{'reactions':['GLUTRR','G1SAT','GLUTRS'],'stoichiometry':[1,1,1]},
        'COFACTOR_2':{'reactions':['DHNAOT4','UPPDC1','DHNCOAT','DHNCOAS','SEPHCHCS','SUCBZS','SUCBZL','PPPGO3','FCLT','CPPPGO','SHCHCS3'],'stoichiometry':[1,1,1,1,1,1,1,1,1,1,1]},
        'COFACTOR_3':{'reactions':['TYRL','AMMQLT8','HEMEOS','UPP3MT','SHCHD2','SHCHF','ENTCS','CBLAT'],'stoichiometry':[1,1,1,1,1,1,1,1]},
        'VITB6':{'reactions':['E4PD','PERD','OHPBAT','PDX5PS','PDX5PO2'],'stoichiometry':[1,1,1,1,1]},
        #'THIAMIN':{'reactions':['AMPMS2','PMPK','THZPSN3','TMPPP','TMPK'],'stoichiometry':[1,1,1,1,1]}, # original pathway without correction
        'THIAMIN':{'reactions':['AMPMS3','PMPK','THZPSN3','TMPPP','TMPK'],'stoichiometry':[1,1,1,1,1]},
        'COFACTOR_4':{'reactions':['I4FE4ST','I4FE4SR','I2FE2SS2'],'stoichiometry':[1,1,1]},
        'COFACTOR_5':{'reactions':['BMOGDS1','BMOGDS2','BMOCOS'],'stoichiometry':[1,1,1]},
        'COFACTOR_6':{'reactions':['DMPPS','GRTT','DMATT'],'stoichiometry':[1,1,1]},
        'COFACTOR_7':{'reactions':['MECDPS','DXPRIi','MEPCT','CDPMEK','MECDPDH5'],'stoichiometry':[1,1,1,1,1]},
        'COFACTOR_8':{'reactions':['LIPOS','LIPOCT'],'stoichiometry':[1,1]},
        'COFACTOR_9':{'reactions':['OMMBLHX','OMPHHX','OPHHX','HBZOPT','DMQMT','CHRPL','OMBZLM','OPHBDC','OHPHM'],'stoichiometry':[1,1,1,1,1,1,1,1,1]},
        'COFACTOR_10':{'reactions':['SERASr','DHBD','UPP3S','HMBS','ICHORT','DHBS'],'stoichiometry':[1,1,1,1,1,1]},
        'COFACTOR_11':{'reactions':['PMEACPE','EGMEACPR','DBTS','AOXSr2','I2FE2SR','OPMEACPD','MALCOAMT','AMAOTr','OPMEACPS','OPMEACPR','OGMEACPD','OGMEACPR','OGMEACPS','EPMEACPR','BTS5'],'stoichiometry':[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]},
        'CELLENV_1':{'reactions':['UAMAGS','UAPGR','UAGPT3','PAPPT3','GLUR','UAGCVT','UAMAS','UDCPDP','UGMDDS','UAAGDS'],'stoichiometry':[1,1,1,1,-1,1,1,1,1,1]},
        'CELLENV_2':{'reactions':['3HAD181','3OAR181','3OAS181','EAR181x'],'stoichiometry':[1,1,1,1]},
        'CELLENV_3':{'reactions':['3HAD160','3OAR160','EAR160x','3OAS160'],'stoichiometry':[1,1,1,1]},
        'CELLENV_4':{'reactions':['EAR120x','3OAR120','3HAD120','3OAS120','EAR100x'],'stoichiometry':[1,1,1,1,1]},
        'CELLENV_5':{'reactions':['G1PACT','UAGDP','PGAMT','GF6PTA'],'stoichiometry':[1,1,-1,1]},
        'CELLENV_6':{'reactions':['3OAR40','EAR40x','3OAS60','3OAR60','3HAD80','3OAS80','3OAR80','EAR60x','3HAD60','EAR80x','3HAD40'],'stoichiometry':[1,1,1,1,1,1,1,1,1,1,1]},
        'CELLENV_7':{'reactions':['3HAD161','EAR161x','3OAS161','3OAR161','3OAS141','3HAD141','3OAR121','EAR121x','3HAD121','EAR141x','T2DECAI','3OAR141','3OAS121'],'stoichiometry':[1,1,1,1,1,1,1,1,1,1,1,1,1]},
        'CELLENV_8':{'reactions':['TDPGDH','TDPDRR','TDPDRE','G1PTT'],'stoichiometry':[1,1,1,1]},
        'CELLENV_9':{'reactions':['3OAS140','3OAR140'],'stoichiometry':[1,1]},
        'CELLENV_10':{'reactions':['3HAD140','EAR140x'],'stoichiometry':[1,1]},
        'CELLENV_11':{'reactions':['3OAR100','3HAD100','3OAS100'],'stoichiometry':[1,1,1]},
        'LIPOPOLYSACCHARIDE_1':{'reactions':['COLIPAabcpp','COLIPAabctex','EDTXS1','EDTXS2','GALT1','GLCTR1','GLCTR2','GLCTR3','HEPK1','HEPK2','HEPT1','HEPT2','HEPT3','HEPT4','LPADSS','MOAT','MOAT2','MOAT3C','RHAT1','TDSK','USHD'],'stoichiometry':[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]},
        'LIPOPOLYSACCHARIDE_2':{'reactions':['AGMHE','GMHEPAT','GMHEPK','GMHEPPA','S7PI'],'stoichiometry':[1,1,1,1,1]},
        'LIPOPOLYSACCHARIDE_3':{'reactions':['U23GAAT','UHGADA','UAGAAT'],'stoichiometry':[1,1,1]},
        'LIPOPOLYSACCHARIDE_4':{'reactions':['KDOPP','KDOCT2','KDOPS'],'stoichiometry':[1,1,1]},
        'ASTPathway':{'reactions':['AST','SADH','SGDS','SGSAD','SOTA'],'stoichiometry':[1,1,1,1,1]},
        'OPPP':{'reactions':['G6PDH2r', 'GND', 'PGL', 'RPI'],'stoichiometry':[1,1,1,-1]},
        'Glycolysis':{'reactions':['PGI', 'PFK', 'FBA', 'TPI', 'GAPD', 'PGK', 'PGM', 'ENO'],
                              'stoichiometry':[1,1,1,1,1,-1,-1,1]},
        'ED_pathway':{'reactions':['G6PDH2r', 'PGL', 'EDD', 'EDA'],
                                      'stoichiometry':[1,1,1,1]},
        'ED_branch':{'reactions':['EDD', 'EDA'],
                                      'stoichiometry':[1,1]}
        };
        return pathways;

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
pathways['dttp synthesis (aerobic)'] = {'reactions':['RNDR4','NDPK6','DUTPDP', 'TMDS', 'DTMPK', 'NDPK4'],
                              'stoichiometry':[1,-1,1,1,1,1]}
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

#for r in pathways['nad salvage']['reactions']:
#    print r, cobra_model.reactions.get_by_id(r).build_reaction_string()
#for r in pathways['gth biosynthesis']['reactions']:
#    print r, cobra_model.reactions.get_by_id(r).build_reaction_string()