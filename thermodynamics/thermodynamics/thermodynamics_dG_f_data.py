from cobra.core.Reaction import Reaction
from cobra.core.Metabolite import Metabolite
from cobra.core.Model import Model
# dependencies from component-contribution
from component_contribution.python.training_data import TrainingData
from component_contribution.python.component_contribution import ComponentContribution
from component_contribution.python.kegg_model import KeggModel
from component_contribution.python.compound_cacher import CompoundCacher
from component_contribution.python.compound_model import compound_model

from .thermodynamics_io import thermodynamics_io

import csv
import re
from math import sqrt
import json

class thermodynamics_dG_f_data(thermodynamics_io):
    """Class for handling dG_f data
    #1: make the reactant contribution (RC) data from the component_contribution method
    #2: combine the RC data with the bibliomic and GC data taken from the equilibrator database
    #3: transform the thermodynamic data to the desired pH, ionic strength and temperature
    #4: load, format, and check the data"""

    def __init__(self,id2KEGGID_filename_I=None,id2KEGGID_I={},dG0_f_I={},dG_f_I={},measured_dG_f_I={},estimated_dG_f_I={}):
        '''
        cobra_model_I = cobra model object
                                                    
        returns a dictionary: measured_dG_fs: transformed compound Gibbs energies of formation
                                    metabolite.id: {'dG_f_lb': float,
                                     'dG_f_ub': float,
                                     'dG_f_units': string}
                                     
                                     '''

        if id2KEGGID_filename_I: 
            self.id2KEGGID = self._get_id2KEGGID_csv(id2KEGGID_filename_I);
            self.KEGGID2id = self._make_invDict(self.id2KEGGID);
        elif id2KEGGID_I: 
            self.id2KEGGID = id2KEGGID_I;
            self.KEGGID2id = self._make_invDict(self.id2KEGGID);
        else: print('no id2KEGGID mapping provided')

        if dG0_f_I:
            self.dG0_f = dG0_f_I
        else:
            self.dG0_f = {}; # units of kJ/mol
        if dG_f_I:
            self.dG_f = self._checkInput_dG_f(dG_f_I)
        else:
            self.dG_f = {}; # units of kJ/mol
        if measured_dG_f_I:
            self.measured_dG_f = measured_dG_f_I
        else:
            self.measured_dG_f = {}; # units of kJ/mol
        if estimated_dG_f_I:
            self.estimated_dG_f = estimated_dG_f_I
        else:
            self.estimated_dG_f = {}; # units of kJ/mol

    def _get_id2KEGGID_csv(self, id2KEGGID_filename_I):
        #Read in the id2KEGGID mapping
        
        id2KEGGID = {};
        with open(id2KEGGID_filename_I,mode='r') as infile:
            reader = csv.reader(infile);
            for i,r in enumerate(reader):
                if not(i==0):
                    id = r[0].replace('-','_DASH_')
                    id2KEGGID[id] = r[1];
        return id2KEGGID

    def _make_invDict(self,dict_I):
        #Make an inverse dictionary from the id2KEGGID mapping

        #check that the values are unique
        if len(dict_I) == len(list(set(dict_I.values()))):
            invDict = {}
            #create the inverse mapping
            invDict = {v: k for k, v in list(dict_I.items())}
            return invDict;
        else:
            raise Exception('Non-unique values present');
  
    ### 1 Start ###
    # implementation of the component contribution method
    # there appears to be numerical instability in calculation of the GC values
    # leading to erroneous dG_r values
    # consequently, RC values were manually taken during the calculation      
    def get_transformed_dG0_f(self,cobra_model_I, pH_I,temperature_I,ionic_strength_I):
        '''get the transformed Gibbs free energies of formation

        relies on component_contribution: doi:10.1371/journal.pcbi.1003098

        cobra_model_I: cobra_model

        pH: {metabolite.compartment {'pH': float}}

        temperature: {metabolite.compartment {'temperature': float,
                'temperature_units': K}}

        ionic strength: {metabolite.compartment {'ionic_strength': float,
                                                    'ionic_strength_units': units}}'''

        self._add_KEGGID(cobra_model_I);
        reaction_list = self._make_KEGG_reaction(cobra_model_I);

        compartments = list(set(cobra_model_I.metabolites.list_attr('compartment')));
        for c in compartments:
            dG_f_intkeys = {};
            dG_var_f_intkeys = {};
            dG_f_strkeys = {};
            dG_var_f_strkeys = {};
            dG_f_intkeys, dG_var_f_intkeys = self._component_contribution_wrapper(list(reaction_list.values()),
                                                                  pH_I[c]['pH'],
                                                                  temperature_I[c]['temperature'],
                                                                  ionic_strength_I[c]['ionic_strength']);
            # convert integer KEGGID to string KEGGID
            for key,value in dG_f_intkeys.items():
                dG_f_strkeys['C%05d' % key] = value;

            for key,value in dG_var_f_intkeys.items():
                dG_var_f_strkeys['C%05d' % key] = value;     
                
            # convert KEGGID to metabolite.id
            for key,value in dG_f_strkeys.items():
                id = self.KEGGID2id[key];
                if id:
                    id_compartment = id + '_' + c;
                    self.dG_f[id_compartment] = {'dG_f':value,
                                 'dG_f_var':dG_var_f_strkeys[key],
                                 'dG_f_units':'kJ/mol'};

            # works but far too slow, inverse dictionary used instead
            #for key,value in dG_f_strkeys.iteritems():
            #    id = self._get_union(cobra_model_I,key,c)
            #    if id: # some reactions may not be in the specified compartment b/c we are 
            #           # feeding in all reactions to component_contribution in order
            #           #  to predict all compounds in a given compartment
            #        if len(id)>2: continue
            #        self.dG_f[id[0].id] = {'dG_f':value,
            #                     'dG_f_var':dG_var_f_strkeys[key],
            #                     'dG_f_units':'kJ/mol'};

    def export_dG0_f(self,filename_I):
        # write json to file
        with open(filename_I, 'w') as outfile:
            json.dump(self.dG_f, outfile, indent=4);

    def _add_KEGGID(self, cobra_model_I):
        #Add KEGGID to notes section of model_metabolites
        for met in cobra_model_I.metabolites:
             id = met.id.replace('_DASH_','-')
             id = id.split('_')[0]
             if id in list(self.id2KEGGID.keys()):
                 cobra_model_I.metabolites.get_by_id(met.id).notes['KEGGID'] = self.id2KEGGID[id];
             else: cobra_model_I.metabolites.get_by_id(met.id).notes['KEGGID'] = None;

    def _make_KEGG_reaction(self, cobra_model_I):
        '''covert a cobra_model reaction list with BIGG ids
        to a reaction list with KEGG ids

        cobra_model_I: cobra model

        returns a dictionary: reaction.id: formula'''

        reaction_list_tmp = {};
        reaction_list = {}
        for rxn in cobra_model_I.reactions:
            r_tmp = '';
            add2rxnlist = True;
            # check for exchange reactions
            if len(rxn.get_reactants()) < 1: continue
            if len(rxn.get_products()) < 1: continue
            # check for diffusion reactions
            if len(rxn.get_reactants()) == len(rxn.get_products()) == 1:
                r = rxn.get_reactants();
                p = rxn.get_products();
                if r[0].id.split('_')[0] == p[0].id.split('_')[0]: continue
            for r in rxn.get_reactants():
                if r.notes['KEGGID']: 
                    r_tmp = r_tmp + str(int(abs(rxn.get_coefficient(r.id)))) + ' ' + r.notes['KEGGID'] + ' + ';
                else: add2rxnlist = False;
            r_tmp = r_tmp[:-2];
            r_tmp = r_tmp + ' <=>  ';
            for p in rxn.get_products():
                if p.notes['KEGGID']: 
                    r_tmp = r_tmp + str(int(abs(rxn.get_coefficient(p.id)))) + ' ' + p.notes['KEGGID'] + ' + ';
                else: add2rxnlist = False;
            r_tmp = r_tmp[:-2];
            if add2rxnlist: reaction_list_tmp[rxn.id] = r_tmp;

        # check for empty reaction
        for key,value in reaction_list_tmp.items():
            check = value.replace(" ", "");
            if not(check == '<=>') and len(check.split('<=>')[0])>0 and len(check.split('<=>')[1])>0:
                reaction_list[key] = value

        return reaction_list;

    def _component_contribution_wrapper(self,reaction_list_I,pH_I,temperature_I,ionic_strength_I):
        '''wrapper for component contribution'''

        # Orginal implementation involves an input of a reaction list
        # the dG_prime and dG_std for each reaction are then calculated.
        # Because we are interested in accounting for transportation
        # and metabolite concentration,
        # we need to extraction out the dG_prime and dG_std
        # of formation and calculated the dG_prime and
        # dG_std for each reaction after accounting for
        # transportation and metabolite concentration.
        # TODO: change input from reaction list to metabolite list
        # however this may be problematic due to the implementation
        # of component_contribution (i.e., the group contribution
        # and reactant contribution are calculated along orthoganol planes
        # in order to gain greater coverage and accuracy of predicted
        # dG_f values)

        # Debugging:
        #wolf = ['C00049 + C00002 <=> C03082 + C00008',
        #        'C00026 + C00049 <=> C00025 + C00036',
        #        'C00158 <=> C00311',
        #        'C00631 <=> C00001 + C00074',
        #        'C00354 <=> C00111 + C00118',
        #        'C00020 + C00002 <=> 2 C00008',
        #        'C00002 + C00001 <=> C00008 + C00009',
        #        'C00015 + C00002 <=> C00075 + C00008',
        #        'C00033 + C00002 + C00010 <=> C00024 + C00020 + C00013']
        #model = KeggModel.from_formulas(wolf) 

        model = KeggModel.from_formulas(reaction_list_I)    

        td = TrainingData()
        cc = ComponentContribution(td)
        model.add_thermo(cc)
    
        dG0_prime_f, dG0_var_f = model.get_transformed_dG0(pH=pH_I, I=ionic_strength_I, T=temperature_I)

        return dG0_prime_f, dG0_var_f

    def _convert_KEGGID2id(self, KEGGID_I):
        '''convert KEGGID to model id

        KEGGID_I = list of kegg ids

        returns a list of metabolite.ids
        '''
        return

    def _get_keggFromNotes(self, cobra_model_I,KEGGID_I):
        return cobra_model_I.metabolites.query(lambda x:x['KEGGID']==KEGGID_I,'notes')

    def _get_compartments(self, cobra_model_I,compartment_I):
        return cobra_model_I.metabolites.query(lambda x:x==compartment_I,'compartment')

    def _get_union(self, cobra_model_I,KEGGID_I,compartment_I):
        return cobra_model_I.metabolites.query(lambda x:x in self._get_keggFromNotes(cobra_model_I,KEGGID_I) and x in self._get_compartments(cobra_model_I,compartment_I))

    ### 1 End ###

    ### 2 Start ###
    # equates to a set of hacks to extraction out the RC data using the component contribution
    # method and combine with the bibliomic and GC data previously described using the
    # psuedoisomer contribution method that is stored on equilibrator
    def make_dG0_f_pH0(self):
        '''
        Set of functions to generate a combined RC, bibliomic, and GC set of dG0_f values
        priority 0: RC (from component contribution; see _make_ccache())
        priority 1: bibliomic (from pseudoisomer contriubtion; from equilibrator website)
        priority 2: GC (from pseudoisomer contribution; from equilibrator website)
        '''
        # get dG0_f (pH = 0, ionic strength = 0, temperature = 298.15 K) for priority 0, 1, and 2
        dG0_f_0 = {};
        dG0_f_1 = {};
        dG0_f_2 = {};
        dG0_f_0 = self._get_pseudoisomer_priority0_pH0();
        dG0_f_1 = self._get_pseudoisomer_priority1_pH0();
        dG0_f_2 = self._get_pseudoisomer_priority2_pH0();

        self._combine_dG0_f_pH0(dG0_f_0, dG0_f_1, dG0_f_2, 'data\\compounds_dG0_f.json')

        #>>> len(keys)
        #12841
        #>>> len(dG0_f_0)
        #671
        #>>> len(dG0_f_1)
        #191
        #>>> len(dG0_f_2)
        #12775
        #>>> 

    def _make_ccache(self):
        pseudo = json.load(open('data\\kegg_pseudoisomers.json'));
        rxns = [];
        cnt = 0;
        for i, item in enumerate(pseudo):
            rxns.append('1 ' + item['CID'] + ' <=> 1 ' + item['CID'])
            cnt = cnt + 1;
            if cnt>2000:
                model = KeggModel.from_formulas(rxns) 
                model.ccache.dump()
                rxns = [];
                cnt = 0;
        model = KeggModel.from_formulas(rxns) 
        model.ccache.dump()

    def _get_pseudoisomer_priority2_pH0(self):
       # compounds = json.load(open('cobra\\thermodynamics\\component_contribution\\cache\\compounds.json'))
        compounds = json.load(open('data\\compounds.json'))
        pseudo = json.load(open('data\\kegg_pseudoisomers.json'))
        # units: dG0_f: kJ/mol
        #        variance = (kJ/mol)^2
        # variance: we are estimating that the variance for the GC method is 62.0 (kJ/mol)^2,
        #           based on an estimate of the uncertainties in the groups in the GC method
        #           doi:10.1371/journal.pcbi.1003098 page 10

        zs_neg = [];
        zs_min = {};
        for item in compounds:
             pKas = item['pKas']
             if pKas:
                 min_pKas = min(pKas);
                 if min_pKas < 0:
                     pKas_neg.append(item['id']);
                     min_pKas = min([x for x in pKas if x>0]);
                 min_pKas_index = [i for i,v in enumerate(pKas) if v==min_pKas][0];
                 min_pKas_index = min_pKas_index + 1; # number of zs = number of pkas + 1
                 zs_min[item['id']] = {'pka_min':min_pKas,\
                                       'z_min':item['zs'][min_pKas_index]};
             else:
                zs_min[item['id']] = {'pka_min':None,\
                                       'z_min':item['zs'][0]};

        dG0_f = {};
        for item in pseudo:
            if 'pmaps' in list(item.keys()):
                pmaps = item['pmaps'];
                for p in pmaps:
                    if 'priority' in list(p.keys()):
                        if p['priority'] == 2:
                            species = p['species'];
                            for s in species:
                                if s['z']==zs_min[item['CID']]['z_min']:
                                    value = s['dG0_f'];
                                    if not(value==type(0.0)): value = float(value);                        
                                    if 'source' in list(p.keys()) and p['priority'] == 2:
                                        dG0_f[item['CID']] = {'priority':2,'source':p['source'], 'dG0_f':value, 'dG0_f_units': 'kJ/mol', 'dG0_f_var':62.0};
                                    else:
                                        dG0_f[item['CID']] = {'priority':2,'source':'', 'dG0_f':value, 'dG0_f_units': 'kJ/mol', 'dG0_f_var':62.0};
        return dG0_f;

    def _get_pseudoisomer_priority1_pH0(self):
        pseudo = json.load(open('data\\kegg_pseudoisomers.json'))
        # units: dG0_f: kJ/mol
        #        variance = (kJ/mol)^2
        # variance: we are estimating that the experimental variance for bibliomic data is 31.0 (kJ/mol)^2,
        #           which is half the variance reported for the group contribution method

        dG0_f = {};
        for item in pseudo:
            if 'pmaps' in list(item.keys()):
                pmaps = item['pmaps'];
                for p in pmaps:
                    if 'priority' in list(p.keys()):
                        if p['priority'] == 1:
                            species = p['species'];
                            zs_min = [];
                            for s in species:
                                zs_min.append(s['z']);
                            z_min = min(zs_min);
                            for s in species:
                                if s['z']==z_min:
                                    value = s['dG0_f'];
                                    if not(value==type(0.0)): value = float(value);
                                    if 'source' in list(p.keys()) and p['priority'] == 1:
                                        dG0_f[item['CID']] = {'priority':1,'source':p['source'], 'dG0_f':value, 'dG0_f_units': 'kJ/mol', 'dG0_f_var':31.0};
                                    else:
                                        dG0_f[item['CID']] = {'priority':1,'source':'', 'dG0_f':value, 'dG0_f_units': 'kJ/mol', 'dG0_f_var':31.0};
        return dG0_f;

    def _get_pseudoisomer_priority0_pH0(self):
        # units: dG0_f: kJ/mol
        #        variance = (kJ/mol)^2
        # variance: we are estimating that the variance for the RC method is 17.8 (kJ/mol)^2,
        #           doi:10.1371/journal.pcbi.1003098 page 10

        dG0_f = {};
        with open('data\\dG0_f_rc_v2.csv','r') as infile:
            reader = csv.reader(infile)
            headers = next(reader);
            for r in reader:
                value = r[1];
                if not(value==type(0.0)): value = float(value);
                dG0_f[r[0]] = {'priority':0,'source':'Reactant Contribution', 'dG0_f':value, 'dG0_f_units': 'kJ/mol', 'dG0_f_var':17.8};

        return dG0_f

    def _combine_dG0_f_pH0(self, dG0_f_0, dG0_f_1, dG0_f_2, filename):
        '''
        combine priority 0, 1, and 2 dG0_f values into a single json datafile
        '''
        compounds_dG0_f = {};
        keys = [];
        if dG0_f_0:keys.extend(list(dG0_f_0.keys()))
        if dG0_f_1:keys.extend(list(dG0_f_1.keys()))
        if dG0_f_2: keys.extend(list(dG0_f_2.keys()))
        keys = sorted(list(set(keys)));
        # add the data into a single json file
        for k in keys:
            compounds_dG0_f[k] = [];
            if k in list(dG0_f_0.keys()):
                compounds_dG0_f[k].append(dG0_f_0[k])
            if k in list(dG0_f_1.keys()):
                compounds_dG0_f[k].append(dG0_f_1[k])
            if k in list(dG0_f_2.keys()):
                compounds_dG0_f[k].append(dG0_f_2[k])
        with open(filename,'w') as outfile:
            json.dump(compounds_dG0_f, outfile, indent=4);

    ### 2 End ###

    ### 3 Start ###
    # upload the combined dG0_f data file and tranfsorm to the desired ph, temp, and ionic strength         
    def get_transformed_dG_f(self,dG0_f_I, cobra_model_I, pH_I,temperature_I,ionic_strength_I):
        '''get the transformed Gibbs free energies of formation

        relies on RC data taken from the component_contribution: doi:10.1371/journal.pcbi.1003098
        and bibliomic/GC data taken from the psuedoisomer contribution method: doi:10.1093/bioinformatics/bts317

        cobra_model_I: cobra_model

        pH: {metabolite.compartment {'pH': float}}

        temperature: {metabolite.compartment {'temperature': float,
                'temperature_units': K}}

        ionic strength: {metabolite.compartment {'ionic_strength': float,
                                                    'ionic_strength_units': units}}'''

        # upload dG0_f values
        dG0_f_KEGG_all = {};
        if type(dG0_f_I) == type({}): dG0_f_KEGG_all=dG0_f_I
        elif type(dG0_f_I) == type(''):
            try:
                dG0_f_KEGG_all = json.load(open(filename_dG0_f_I));
            except IOError as e:
                print(e);
                return;

        # extract out dG0_f values in the model
        for m in cobra_model_I.metabolites:
            m_id_model = m.id[:-2] # assuming that the compartment is a 1 letter abbreviation!
            m_id_kegg = '';
            if m_id_model in list(self.id2KEGGID.keys()):
                m_id_kegg = self.id2KEGGID[m_id_model];
            else: 
                continue;
            if m_id_kegg and m_id_kegg in list(dG0_f_KEGG_all.keys()) and not(m_id_model in list(self.dG0_f.keys())):
                # get minimum priority value
                priority = [];
                min_priority = 0;
                for p in dG0_f_KEGG_all[m_id_kegg]:
                    priority.append(p['priority'])
                min_priority = min(priority);
                # copy the minimum priority dG0_f data
                for p in dG0_f_KEGG_all[m_id_kegg]:
                    if p['priority'] == min_priority:
                        self.dG0_f[m_id_kegg] = {'dG_f':p['dG0_f'],
                                 'dG_f_var':p['dG0_f_var'],
                                 'dG_f_units':p['dG0_f_units']}; 

        # transform the data
        compartments = list(set(cobra_model_I.metabolites.list_attr('compartment')));
        for c in compartments:
            # get metabolite ids and dG0_f data of the particular compartment
            dG0_f_ids = [];
            dG0_f_data = [];
            dG0_f_var = [];
            dG0_f_units = [];

            for key,value in self.dG0_f.items():
                #id_compartment = key + '_' + c;
                #dG0_f_ids.append(id_compartment);
                dG0_f_ids.append(key);
                dG0_f_data.append(value['dG_f']);
                dG0_f_var.append(value['dG_f_var']);
                dG0_f_units.append(value['dG_f_units']);

            dG_f_data = [];
            dG_f_data = self._get_transformed_f(dG0_f_ids, dG0_f_data, 
                                                pH_I[c]['pH'],
                                                temperature_I[c]['temperature'],
                                                ionic_strength_I[c]['ionic_strength']);

            # convert KEGGID to metabolite.id
            for i,kegg_id in enumerate(dG0_f_ids):
                id = self.KEGGID2id[kegg_id];
                if id:
                    id_compartment = id + '_' + c;
                    self.dG_f[id_compartment] = {'dG_f':dG_f_data[i],
                                 'dG_f_var':dG0_f_var[i],
                                 'dG_f_units':dG0_f_units[i]};

    def export_dG_f(self,filename_I):
        # write json to file
        with open(filename_I, 'w') as outfile:
            json.dump(self.dG_f, outfile, indent=4);

    def _get_transformed_f(self, dG0_f_id_I, dG0_f_I, pH_I, temperature_I, ionic_strength_I):

        c_model = compound_model(dG0_f_id_I, dG0_f_I);
        dG_f = c_model.get_transformed_dG0_f(pH=pH_I, I=ionic_strength_I, T=temperature_I)

        return dG_f;
    ### 3 END ###

    ### 4 START ###
    # load, format, and check data

    def import_dG_f(self, dG_f_filename_I):
        '''import measured values required for analysis'''

        self.dG_f = {};
        self.dG_f = self._checkInput_dG_f(self.import_values_json(dG_f_filename_I));

    def format_dG_f(self):
        '''format data'''

        self.measured_dG_f = self._convert_var2lbub_dG_f(self.dG_f);

    def generate_estimated_dG_f(self, cobra_model):
        '''generate estimated values'''

        self.estimated_dG_f = self._generalize_compartment2all_dG_f(cobra_model);

    def check_data(self):
        '''check data integrity'''
        return;

    def _checkInput_dG_f(self, measured_values_I):
        """
        check dG_O_f data input

        measured_values: measured values with variances
                                 {metabolite.id: {'dG_f': float,
                                                 'dG_f_var': float,
                                                 'dG_f_lb': float,
                                                 'dG_f_ub': float,
                                                 'dG_f_units': 'kJ/mol'}
        returns a dictionary: measured_values_O:
                                 {metabolite.id: {'dG_f': float,
                                                 'dG_f_var': float,
                                                 'dG_f_lb': float,
                                                 'dG_f_ub': float,
                                                 'dG_f_units': 'kJ/mol'}
        """
        # check units
    
        measured_values_O = {};
        for k,v in measured_values_I.items():
            if v['dG_f_units'] == 'kJ/mol':
                measured_values_O[k] = v
            else:
                print((str(k) + ' has invalid units of ' + str(v) + ' and will be ignored'))
        return measured_values_O;

    def _convert_var2lbub_dG_f(self, measured_values):
        """
        convert measured dG_f values with a stadard deviation
        currently use +/- SD
        lb = ave - sqrt(var)
        ub = ave + sqrt(var)

        measured_values: measured values with variances
                                 {metabolite.id: {'dG_f': float,
                                                 'dG_f_var': float,
                                                 'dG_f_units': 'kJ/mol'}
        returns a dictionary: measured_values_O:
                                 {metabolite.id: {'dG_f': float,
                                                 'dG_f_var': float,
                                                 'dG_f_lb': float,
                                                 'dG_f_ub': float,
                                                 'dG_f_units': 'kJ/mol'}
        """
        measured_values_O = {};
        for k,v in measured_values.items():
             concMlb = 0.0;
             concMub = 0.0;
             if v['dG_f_var']:
                 concMlb = v['dG_f'] - sqrt(v['dG_f_var']);
                 concMub = v['dG_f'] + sqrt(v['dG_f_var']);
             else:
                 concMlb = v['dG_f'];
                 concMub = v['dG_f'];
             #if concMlb<0: concMlb = min_value;
             measured_values_O[k] = {'dG_f': v['dG_f'],
                                     'dG_f_var': v['dG_f_var'],
                                     'dG_f_lb': concMlb,
                                                 'dG_f_ub': concMub,
                                                 'dG_f_units': 'kJ/mol'}
        return measured_values_O

    def _convert_std2lbub_dG_f(self, measured_values,min_value):
        """
        convert measured values with a standard deviation
        (CV = SD/Ave*100; SD = CV/100*AVE) to lb and ub
        currently use +/- SD
        lb = ave - sqrt(var) NOTE: if < 0 min_value is used instead
        ub = ave + sqrt(var)

        measured_values: measured values with variances
                                 {metabolite.id: {'dG_f': float,
                                                 'dG_f_cv': float,
                                                 'dG_f_units': 'mM'}
        returns a dictionary: measured_values_O:
                                 {metabolite.id: {'dG_f_lb': float,
                                                 'dG_f_ub': float,
                                                 'dG_f_units': 'M'}
        """
        measured_values_O = {};
        for k,v in measured_values.items():
             concMlb = 0.0;
             concMub = 0.0;
             concMlb = v['dG_f'] - v['dG_f_var'];
             concMub = v['dG_f'] + v['dG_f_var'];
             if concMlb<0: concMlb = min_value;
             measured_values_O[k] = {'dG_f_lb': concMlb,
                                                 'dG_f_ub': concMub,
                                                 'dG_f_units': 'M'}
        return measured_values_O

    def _generalize_compartmentLBUB2all_dG_f(self, cobra_model, lbub=None, exceptions=None):
        """
        takes a compartment and lb/ub for that compartment
        and updates each metabolite in that compartment

        allows for exceptions to the generalization as input

        cobra_model: a Model object

        lbub: metabolite.compartment: {'dG_f_lb': float,
                               'dG_f_ub': float,
                               'dG_f_units': 'kJ/mol'}

        returns a dictionary: metabolite.id {'dG_f_lb': float,
                               'dG_f_ub': float,
                               'dG_f_units': string}
        """
        if not(lbub):
            lbub = {};
            compartments = list(set(cobra_model.metabolites.list_attr('compartment')));
            for c in compartments:
                lbub[c] = {'dG_f_lb':-1e6,
                                       'dG_f_ub':1e6,
                                       'dG_f_units':'kJ/mol'};
        default_values = {};
        for m in cobra_model.metabolites:
            default_values[m.id] = {'dG_f_lb': lbub[m.compartment]['dG_f_lb'],
                                   'dG_f_ub': lbub[m.compartment]['dG_f_ub'],
                                   'dG_f_units': lbub[m.compartment]['dG_f_units']};
            if exceptions:    
                for k,v in exceptions.items():
                    if k in m.id:
                        default_values[m.id] = {'dG_f_lb':v['dG_f_lb'],
                                       'dG_f_ub':v['dG_f_ub'],
                                       'dG_f_units':v['dG_f_units']};
    
        return default_values;

    def _generalize_compartment2all_dG_f(self, cobra_model, dG_f=None, exceptions=None):
        """
        takes a compartment and lb/ub for that compartment
        and updates each metabolite in that compartment

        allows for exceptions to the generalization as input

        the lower bounds and upper bounds are estimated as
        lb = mean - sqrt(var)
        ub = mean + sqrt(var)

        cobra_model: a Model object

        lbub: metabolite.compartment: {'dG_f': float,
                               'dG_f_var': float,
                               'dG_f_units': 'kJ/mol'}

        returns a dictionary: metabolite.id {'dG_f': float,
                               'dG_f_var': float,
                               'dG_f_lb': float,
                               'dG_f_ub': float,
                               'dG_f_units': 'kJ/mol'}
        """
        if not(dG_f):
            dG_f = {};
            compartments = list(set(cobra_model.metabolites.list_attr('compartment')));
            for c in compartments:
                dG_f[c] = {'dG_f':0.0,
                                       'dG_f_var':1e12, # based on the convention described in
                                                       # doi:10.1371/journal.pcbi.1003098
                                       'dG_f_units':'kJ/mol'};
        default_values = {};
        for m in cobra_model.metabolites:
            default_values[m.id] = {'dG_f': dG_f[m.compartment]['dG_f'],
                                   'dG_f_var': dG_f[m.compartment]['dG_f_var'],
                                   'dG_f_lb': dG_f[m.compartment]['dG_f'] - sqrt(dG_f[m.compartment]['dG_f_var']),
                                   'dG_f_ub': dG_f[m.compartment]['dG_f'] + sqrt(dG_f[m.compartment]['dG_f_var']),
                                   'dG_f_units': dG_f[m.compartment]['dG_f_units']};
            if exceptions:    
                for k,v in exceptions.items():
                    if k in m.id:
                        default_values[m.id] = {'dG_f':v['dG_f'],
                                       'dG_f_var':v['dG_f_var'],
                                       'dG_f_lb':v['dG_f_lb'],
                                       'dG_f_ub':v['dG_f_ub'],
                                       'dG_f_units':v['dG_f_units']};
    
        return default_values;    

    def remove_measured_dG_f(self,mets_I):
        '''Remove measured metabolite dG_f'''

        for met in mets_I:
            v=self.measured_dG_f.pop(met);

