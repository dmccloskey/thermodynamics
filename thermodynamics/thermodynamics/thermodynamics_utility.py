from cobra.core.Model import Model
from cobra.core.Metabolite import Metabolite
from cobra.core.Reaction import Reaction
from collections import Counter

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