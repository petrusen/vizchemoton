import sys
import numpy as np
import networkx as nx
from math import exp, log10, log
import scine_utilities as utils
import scine_database as db
from scine_chemoton.gears.pathfinder import Pathfinder as pf
from scine_database.energy_query_functions import (get_energy_change,
    get_barriers_for_elementary_step_by_type,
    rate_constant_from_barrier, get_energy_for_structure
)
from itertools import combinations
from kinetic_utility_library import * 
import json

manager = db.Manager()
db_name = "ozone_ethene_orca"
ip = 'localhost'
port = '27018' #'8889'
credentials = db.Credentials(ip, int(port), db_name)
manager.set_credentials(credentials)
manager.connect()

if True:

    model1 = db.Model("dft", "pbe-d3bj", "def2-tzvp")
    #model = db.Model("cc", "dlpno-ccsd(t)", "cc-pvtz")
    #model1.spin_mode = "any"
    model1.program = "orca"
    #model1.solvation = "none"


    #e_threshold = np.inf
    e_threshold = 150
    #energy_type = 'electronic_energy'
    fingerprint = 'test'
    path_edges = 'edges_' + fingerprint
    path_nodes = 'nodes_' + fingerprint
    path_labels = 'labels_' + fingerprint
    #path_edges = '/home/petrusen/PycharmProjects

#########################################################################

calculations = manager.get_collection("calculations")
structures = manager.get_collection("structures")
reactions = manager.get_collection("reactions")
flasks = manager.get_collection("flasks")
compounds = manager.get_collection("compounds")
properties = manager.get_collection('properties')
elementary_steps = manager.get_collection('elementary_steps')

# # # Load Pathfinder and assign NetworkX Digraph
pathfinder = pf(manager)
pathfinder.load_graph("orca_test.json")
if False:
  pathfinder.options.model = model1
  pathfinder.options.graph_handler = "barrier"
  pathfinder.options.use_structure_model = True
  pathfinder.options.structure_model = model1  # primitive_model # db.Model("gfn2", "gfn2", "")
  pathfinder.build_graph()
  pathfinder.export_graph("orca_test.json")


print(len([node for node in pathfinder.graph_handler.graph.nodes if ';' in node]) / 2)

# # # List of compounds and reactions
cmp_list = [node for node in pathfinder.graph_handler.graph.nodes if ";" not in node]
lhs_rxn_list = [node for node in pathfinder.graph_handler.graph.nodes if ";0;" in node]

cmp_dict = {}
stoich_id_dict = {}
excluded_rxn, activation_barriers = [], []
tetramolecular, l_new_reac, l_reac_id = 0, list(), list()
comp_idx = len(cmp_list)
trimolecular_constants = list()
my_list = []
cmp_idx = 1
html_reactions, html_compounds = list(), dict()

for rxn_ind, rxn_id in enumerate(lhs_rxn_list):
    rxn = db.Reaction(db.ID(rxn_id[:-3]), reactions)
    reactants = rxn.get_reactants(db.Side.BOTH)
    lhs, rhs = reactants
    if False:
        _idx_lhs = [(p,cmp_dict[p.string()]) for p in _lhs]
        _idx_rhs = [(p,cmp_dict[p.string()]) for p in _rhs]
        redundant_id = list()
        for o in _idx_lhs:
            redun = None
            if o in _idx_rhs:
                redundant += 1; redun = o[-1]; redundant_id.append(o)

        lhs = tuple(a for a,b in _idx_lhs if (a,b) not in redundant_id)
        rhs = tuple(a for a,b in _idx_rhs if (a,b) not in redundant_id)

    s_lhs, s_rhs = len(lhs), len(rhs)
    reactants = (lhs, rhs)

    if (s_lhs < 3 and s_rhs < 3): # and not any(item in rxn_id for item in my_list):  # Unimolecular and Bimolecular
        # # # Get reagent indices
        r_index = [0.0, 0.0]
        count = 0
        if len(reactants[0]) == 1:
           node_x = reactants[0][0].string()
           if node_x in cmp_dict.keys():
               r_index[count * -1 + 1] = float(cmp_dict[node_x])
           else:
               cmp_dict[node_x] = cmp_idx
               r_index[count * -1 + 1] = float(cmp_dict[node_x])
               cmp_idx = cmp_idx + 1
        
        elif len(reactants[0]) == 2:
            node_x = "//".join(sorted([o.string() for o in reactants[0]]))
            if node_x in cmp_dict.keys():
                r_index[count * -1 + 1] = float(cmp_dict[node_x])
            else:
                cmp_dict[node_x] = cmp_idx
                r_index[count * -1 + 1] = float(cmp_dict[node_x])
                cmp_idx = cmp_idx + 1
        
        
        p_index = [0.0, 0.0]
        if len(reactants[1]) == 1:
            node_y = reactants[1][0].string()
            if node_y in cmp_dict.keys():
                p_index[count * -1 + 1] = float(cmp_dict[node_y])
            else:
                cmp_dict[node_y] = cmp_idx
                p_index[count * -1 + 1] = float(cmp_dict[node_y])
                cmp_idx = cmp_idx + 1
        elif len(reactants[1]) > 1:
            node_y = "//".join(sorted([o.string() for o in reactants[1]]))
            if node_y in cmp_dict.keys():
                p_index[count * -1 + 1] = float(cmp_dict[node_y])
            else:
                cmp_dict[node_y] = cmp_idx
                p_index[count * -1 + 1] = float(cmp_dict[node_y])
                cmp_idx = cmp_idx + 1

        # # # Get selected Elementary Step if present in graph
        if "elementary_step_id" in pathfinder.graph_handler.graph.nodes(data=True)[rxn_id]:
            es_id = db.ID(
                pathfinder.graph_handler.graph.nodes(data=True)[rxn_id]["elementary_step_id"])
            es_from_graph = db.ElementaryStep(es_id, elementary_steps)
            _energy, barriers, not_None = get_energy_and_barriers('electronic_energy', es_id, elementary_steps, model1, structures, properties, es_from_graph)

            if es_from_graph.get_type() == db.ElementaryStepType.BARRIERLESS and not_None: #_energy != None:  # check as jacs
                html_reactions.append((cmp_dict[node_x], cmp_dict[node_y], None))
            elif not_None:
                b1, b2 = barriers
                node_ts = es_from_graph.get_transition_state().string() + ";"
                if node_ts in cmp_dict.keys():
                  pass
                else:
                  cmp_dict[node_ts] = cmp_idx
                  cmp_idx = cmp_idx + 1
                html_reactions.append((cmp_dict[node_x], cmp_dict[node_y], cmp_dict[node_ts]))


for compound_id in cmp_dict:
    html_compounds[cmp_dict[compound_id]] = {}
    html_compounds[cmp_dict[compound_id]]['id'] = compound_id
    html_compounds[cmp_dict[compound_id]]['xyz'] = []
    html_compounds[cmp_dict[compound_id]]['charge'] = []
    html_compounds[cmp_dict[compound_id]]['multiplicity'] = []
    html_compounds[cmp_dict[compound_id]]['energy'] = []

    if "//" in compound_id:
        ids = compound_id.split("//")
        for _ids in ids:
          type_object = pathfinder.graph_handler.graph.nodes(data=True)[_ids]["type"]
          if type_object == db.CompoundOrFlask.COMPOUND.name:
              compound = db.Compound(db.ID(_ids), compounds)
          else:
              compound = db.Flask(db.ID(_ids), flasks)
          structure = compound.get_centroid()
          structure_obj = db.Structure(structure, structures)
          xyz = [(str(o.element), tuple(o.position)) for o in structure_obj.get_atoms()]
          #print(structure_obj.get_charge(), structure_obj.multiplicity, dir(structure_obj))
          z, s = structure_obj.get_charge(), structure_obj.multiplicity
          e = get_energy_for_structure(structure_obj, 'electronic_energy', model1, structures, properties)
          e_kj = e * utils.KJPERMOL_PER_HARTREE
          html_compounds[cmp_dict[compound_id]]['xyz'].append(xyz)
          html_compounds[cmp_dict[compound_id]]['charge'].append(z)
          html_compounds[cmp_dict[compound_id]]['multiplicity'].append(s)         
          html_compounds[cmp_dict[compound_id]]['energy'].append(e_kj)     


    elif ";" in compound_id:  # its a ts
        structure = compound_id[0:-1]
        structure_obj = db.Structure(db.ID(structure), structures)
        xyz = [(str(o.element), tuple(o.position)) for o in structure_obj.get_atoms()]
        z, s = structure_obj.get_charge(), structure_obj.multiplicity
        e = get_energy_for_structure(structure_obj, 'electronic_energy', model1, structures, properties)
        e_kj = e * utils.KJPERMOL_PER_HARTREE
        #print(structure_obj.get_charge(), structure_obj.multiplicity, dir(structure_obj))
        html_compounds[cmp_dict[compound_id]]['xyz'].append(xyz)  
        html_compounds[cmp_dict[compound_id]]['charge'].append(z)
        html_compounds[cmp_dict[compound_id]]['multiplicity'].append(s)      
        html_compounds[cmp_dict[compound_id]]['energy'].append(e_kj)
    else:
        type_object = pathfinder.graph_handler.graph.nodes(data=True)[compound_id]["type"]
        if type_object == db.CompoundOrFlask.COMPOUND.name:
            compound = db.Compound(db.ID(compound_id), compounds)
        else:
            compound = db.Flask(db.ID(compound_id), flasks)
        structure = compound.get_centroid()
        structure_obj = db.Structure(structure, structures)
        xyz = [(str(o.element), tuple(o.position)) for o in structure_obj.get_atoms()]
        z, s = structure_obj.get_charge(), structure_obj.multiplicity
        e = get_energy_for_structure(structure_obj, 'electronic_energy', model1, structures, properties)
        e_kj = e * utils.KJPERMOL_PER_HARTREE
        html_compounds[cmp_dict[compound_id]]['xyz'].append(xyz)
        html_compounds[cmp_dict[compound_id]]['charge'].append(z)
        html_compounds[cmp_dict[compound_id]]['multiplicity'].append(s)
        html_compounds[cmp_dict[compound_id]]['energy'].append(e_kj) 

# Open a file in write mode
with open('reactions_epetrus_241010.csv', 'w') as f:
    # Loop through the list and write each tuple to the file
    for item in html_reactions:
        r, p, ts = item
        f.write("{r},{p},{ts}\n".format(r=r, p=p, ts=ts))

with open('compounds_epetrus_241010.txt', 'w') as file:
     file.write(json.dumps(html_compounds)) # use `json.loads` to do the revers






# print("Included Reactions:", len(lhs_rxn_list) - len(excluded_rxn))


