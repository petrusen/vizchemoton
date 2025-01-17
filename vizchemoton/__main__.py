'''
Enric Petrus, December 2024. Added SCINE helper functiosn to link with the amk-tools generation of html files.
Diego Garay-Ruiz, November 2023. Collection of helper functions to link amk-tools and grrm-tools, generating interactive
HTML dashboards to visualize GRRM-generated reaction networks.
'''

# Standard Library Imports
import sys
from collections import Counter
import argparse
import json

#Third-Party Library Imports
import numpy as np
import networkx as nx
import bokeh.plotting
import bokeh.models as bkm
import RXVisualizer as arxviz
import networkx as nx

# Project-Specific SCINE imports
import scine_utilities as utils
import scine_database as db
from scine_chemoton.gears.pathfinder import Pathfinder as pf
from scine_database.energy_query_functions import (get_energy_change,
    get_barriers_for_elementary_step_by_type,
    rate_constant_from_barrier, get_energy_for_structure
)

def get_energy_and_barriers(energy_type, es_id, elementary_steps, model1, structures, properties, es_from_graph):
    """
    Wrapper function Gets the elementary step ID with the lowest energy of the corresponding transition state of a
    reaction.

    Input:
      - energy_type (str): name of the energy property such as 'electronic_energy' or 'gibbs_free_energy'
      - es_id (str): id of the elementary_step
      - elementary_steps (db.Collection): the elementary step collection
      - model1 (dict): dictionary with the method_family, method, basis_set and program keys.
      - structures (db.Collection): the structures collection
      - properties (db.Collection): the properties step collection
      - es_from_graph (db.ElementaryStep): db object for the given es_id 

    Returns:
      - energy (float): energy state in the reaction.
      - barriers (tuple): forward and backward energy barriers 
      - not_None (bool): returns True if no None was found in the barriers tuple 
    """
    energy = get_energy_change(db.ElementaryStep(es_id, elementary_steps), energy_type, model1, structures, 
                               properties)
    barriers = get_barriers_for_elementary_step_by_type(es_from_graph, energy_type, model1, structures, properties)
    
    if None in barriers: 
        not_None = False
    else:
        not_None = True

    return energy, barriers, not_None


def get_reactions_and_compounds(db_name, ip, port, dict_method, read_pathfinder=False, write_pathfinder=False,
                                verbose=True):
    """
    Extract the chemical reactions, compounds and transition states from the Mongo-DB where the exploration
    with Chemoton was run.

    Input:
      - db_name (str): name of the database
      - ip (str): internet protocol to the database
      - port (int): port number to the database
      - read_pathfinder (bool, str, optional): either False if no file to read, or string with the path to the file
      - write_pathfinder (bool, str, optional): either False if no file to write, or string with the path to the new file
      - verbose (bool, optional): if True, enable verbose output. Default is True.

    Returns:
      - html_reactions (list): a list of tuples with the indexes of the reactant, product, and TS.
      - html_compounds (dict): a dictionary for each compound containing relevant information (charge, spin, xyz ...)
    """

    manager = db.Manager()
    credentials = db.Credentials(ip, int(port), db_name)
    manager.set_credentials(credentials)
    manager.connect()
    model1 = db.Model(dict_method["method_family"], dict_method["method"], dict_method["basis_set"])
    model1.program = dict_method["program"]

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

    if isinstance(read_pathfinder, str):
        if verbose: print("## reading pathfinder object with name "+read_pathfinder)
        pathfinder.load_graph(read_pathfinder)
    elif isinstance(write_pathfinder, str):
        if verbose: print("## writing new pathfinder object with name "+write_pathfinder)
        pathfinder.options.model = model1
        pathfinder.options.graph_handler = "barrier"
        pathfinder.options.use_structure_model = True
        pathfinder.options.structure_model = model1
        pathfinder.build_graph()
        pathfinder.export_graph(write_pathfinder)


    # print(len([node for node in pathfinder.graph_handler.graph.nodes if ';' in node]) / 2)

    # # # List of compounds and reactions
    lhs_rxn_list = [node for node in pathfinder.graph_handler.graph.nodes if ";0;" in node]
    cmp_idx = 1
    cmp_dict, html_reactions, html_compounds = dict(), list(), dict()

    if verbose: print("## iterating through reactions")
    for rxn_ind, rxn_id in enumerate(lhs_rxn_list):
        # Iterate through the reations of the network
        rxn = db.Reaction(db.ID(rxn_id[:-3]), reactions)
        reactants = rxn.get_reactants(db.Side.BOTH)
        lhs, rhs = reactants
        s_lhs, s_rhs = len(lhs), len(rhs)
        reactants = (lhs, rhs)

        if s_lhs < 3 and s_rhs < 3:

            # Get reactant indexes
            cmp_dict_keys = cmp_dict.keys()
            if len(reactants[0]) == 1:
                node_x = reactants[0][0].string()
                if node_x not in cmp_dict_keys:
                    cmp_dict[node_x] = cmp_idx
                    cmp_idx = cmp_idx + 1
            elif len(reactants[0]) == 2:
                for node_i in [o.string() for o in reactants[0]]:
                    if node_i not in cmp_dict_keys:
                        cmp_dict[node_i] = cmp_idx
                        cmp_idx = cmp_idx + 1
                # flasks (i.e., adducts) are depicted with //
                node_x = "//".join(sorted([o.string() for o in reactants[0]]))
                if node_x not in cmp_dict_keys:
                    cmp_dict[node_x] = cmp_idx
                    cmp_idx = cmp_idx + 1

            # Get product indexes
            if len(reactants[1]) == 1:
                node_y = reactants[1][0].string()
                if node_y not in cmp_dict_keys:
                    cmp_dict[node_y] = cmp_idx
                    cmp_idx = cmp_idx + 1
            elif len(reactants[1]) == 2:
                for node_i in [o.string() for o in reactants[1]]:
                    if node_i not in cmp_dict_keys:
                        cmp_dict[node_i] = cmp_idx
                        cmp_idx = cmp_idx + 1
                # flasks (i.e., adducts) are depicted with //
                node_y = "//".join(sorted([o.string() for o in reactants[1]]))
                if node_y not in cmp_dict_keys:
                    cmp_dict[node_y] = cmp_idx
                    cmp_idx = cmp_idx + 1

            # Get elementary steps and energies
            if "elementary_step_id" in pathfinder.graph_handler.graph.nodes(data=True)[rxn_id]:
                es_id = db.ID(
                    pathfinder.graph_handler.graph.nodes(data=True)[rxn_id]["elementary_step_id"])
                es_from_graph = db.ElementaryStep(es_id, elementary_steps)
                _energy, barriers, not_None = get_energy_and_barriers('electronic_energy', es_id, elementary_steps,
                                                                      model1, structures, properties, es_from_graph)

                if es_from_graph.get_type() == db.ElementaryStepType.BARRIERLESS and not_None and _energy != None:  # check as jacs
                    html_reactions.append([cmp_dict[node_x], cmp_dict[node_y], None])
                elif not_None:
                    b1, b2 = barriers
                    node_ts = es_from_graph.get_transition_state().string() + ";"
                    if node_ts not in cmp_dict.keys():
                        cmp_dict[node_ts] = cmp_idx
                        cmp_idx = cmp_idx + 1
                    html_reactions.append([cmp_dict[node_x], cmp_dict[node_y], cmp_dict[node_ts]])

    if verbose: print("## preparing compounds json object")
    for compound_id in cmp_dict:
        html_compounds[cmp_dict[compound_id]] = {}
        if "//" in compound_id:  # checking the flasks
            # if the user is interested in uploading the data in ioChem-BD, this conditional
            # block should be disregarded by deactivating the following line:
            # continue
            ids = compound_id.split("//")
            html_compounds[cmp_dict[compound_id]]['crn_id'] = list()
            html_compounds[cmp_dict[compound_id]]['_mongodb_id'] = list()
            html_compounds[cmp_dict[compound_id]]['mongodb_id'] = list()
            html_compounds[cmp_dict[compound_id]]['xyz'] = list()
            html_compounds[cmp_dict[compound_id]]['charge'] = list()
            html_compounds[cmp_dict[compound_id]]['multiplicity'] = list()
            html_compounds[cmp_dict[compound_id]]['energy'] = list()
            html_compounds[cmp_dict[compound_id]]['method'] = list()
            html_compounds[cmp_dict[compound_id]]['basis_set'] = list()
            html_compounds[cmp_dict[compound_id]]['program'] = list()
            html_compounds[cmp_dict[compound_id]]['solvent'] = list()
            html_compounds[cmp_dict[compound_id]]['solvation'] = list()
            for _ids in ids:
                type_object = pathfinder.graph_handler.graph.nodes(data=True)[_ids]["type"]
                if type_object == db.CompoundOrFlask.COMPOUND.name:
                    compound = db.Compound(db.ID(_ids), compounds)
                    html_compounds[cmp_dict[compound_id]]['crn_id'].append("c" + str(cmp_dict[_ids]))
                else:
                    compound = db.Flask(db.ID(_ids), flasks)
                    html_compounds[cmp_dict[compound_id]]['crn_id'].append("f" + str(cmp_dict[_ids]))
                structure = compound.get_centroid()
                structure_obj = db.Structure(structure, structures)
                xyz = [(str(o.element), tuple(o.position)) for o in structure_obj.get_atoms()]
                #print(structure_obj.get_charge(), structure_obj.multiplicity, dir(structure_obj))
                z, s = structure_obj.get_charge(), structure_obj.multiplicity
                e = get_energy_for_structure(structure_obj, 'electronic_energy', model1, structures, properties)
                e_kj = e * utils.KJPERMOL_PER_HARTREE
                html_compounds[cmp_dict[compound_id]]['_mongodb_id'].append(_ids)
                html_compounds[cmp_dict[compound_id]]['xyz'].append(xyz)
                html_compounds[cmp_dict[compound_id]]['charge'].append(z)
                html_compounds[cmp_dict[compound_id]]['multiplicity'].append(s)
                html_compounds[cmp_dict[compound_id]]['energy'].append(e_kj)
                html_compounds[cmp_dict[compound_id]]['method'].append(model1.method)
                html_compounds[cmp_dict[compound_id]]['basis_set'].append(model1.basis_set)
                html_compounds[cmp_dict[compound_id]]['program'].append(model1.program + " " + model1.version)
                html_compounds[cmp_dict[compound_id]]['solvent'].append(model1.solvent)
                html_compounds[cmp_dict[compound_id]]['solvation'].append(model1.solvation)
            html_compounds[cmp_dict[compound_id]]['mongodb_id'] = "//".join(
                html_compounds[cmp_dict[compound_id]]['_mongodb_id'])

        elif ";" in compound_id:  # checking the transitions states
            structure = compound_id[0:-1]
            structure_obj = db.Structure(db.ID(structure), structures)
            xyz = [(str(o.element), tuple(o.position)) for o in structure_obj.get_atoms()]
            z, s = structure_obj.get_charge(), structure_obj.multiplicity
            e = get_energy_for_structure(structure_obj, 'electronic_energy', model1, structures, properties)
            e_kj = e * utils.KJPERMOL_PER_HARTREE
            model_obj = structure_obj.get_model()
            # print(structure_obj.get_charge(), structure_obj.multiplicity, dir(structure_obj))
            html_compounds[cmp_dict[compound_id]]['crn_id'] = "ts" + str(cmp_dict[compound_id])
            html_compounds[cmp_dict[compound_id]]['mongodb_id'] = compound_id
            html_compounds[cmp_dict[compound_id]]['xyz'] = xyz
            html_compounds[cmp_dict[compound_id]]['charge'] = z
            html_compounds[cmp_dict[compound_id]]['multiplicity'] = s
            html_compounds[cmp_dict[compound_id]]['energy'] = e_kj
            html_compounds[cmp_dict[compound_id]]['method'] = model1.method
            html_compounds[cmp_dict[compound_id]]['basis_set'] = model1.basis_set
            html_compounds[cmp_dict[compound_id]]['program'] = model1.program + " 5.0.3"
            html_compounds[cmp_dict[compound_id]]['solvent'] = model1.solvent
            html_compounds[cmp_dict[compound_id]]['solvation'] = model1.solvation

        else:  # checking compounds
            type_object = pathfinder.graph_handler.graph.nodes(data=True)[compound_id]["type"]
            if type_object == db.CompoundOrFlask.COMPOUND.name:
                compound = db.Compound(db.ID(compound_id), compounds)
                html_compounds[cmp_dict[compound_id]]['crn_id'] = "c" + str(cmp_dict[compound_id])
            else:
                html_compounds[cmp_dict[compound_id]]['crn_id'] = "f" + str(cmp_dict[compound_id])
                compound = db.Flask(db.ID(compound_id), flasks)
            structure = compound.get_centroid()
            structure_obj = db.Structure(structure, structures)
            xyz = [(str(o.element), tuple(o.position)) for o in structure_obj.get_atoms()]
            z, s = structure_obj.get_charge(), structure_obj.multiplicity
            e = get_energy_for_structure(structure_obj, 'electronic_energy', model1, structures, properties)
            e_kj = e * utils.KJPERMOL_PER_HARTREE
            model_obj = structure_obj.get_model()
            html_compounds[cmp_dict[compound_id]]['mongodb_id'] = compound_id
            html_compounds[cmp_dict[compound_id]]['xyz'] = xyz
            html_compounds[cmp_dict[compound_id]]['charge'] = z
            html_compounds[cmp_dict[compound_id]]['multiplicity'] = s
            html_compounds[cmp_dict[compound_id]]['energy'] = e_kj
            html_compounds[cmp_dict[compound_id]]['method'] = model1.method  # model_obj.method
            html_compounds[cmp_dict[compound_id]]['basis_set'] = model1.basis_set  # model_obj.basis_set
            html_compounds[cmp_dict[compound_id]][
                'program'] = model1.program + " 5.0.3"  # model_obj.program+" "+model_obj.version
            html_compounds[cmp_dict[compound_id]]['solvent'] = model1.solvent  # model_obj.solvent
            html_compounds[cmp_dict[compound_id]]['solvation'] = model1.solvation  # model_obj.solvation

    return html_reactions, html_compounds


def write_compound_reactions_files(html_reactions, html_compounds, reaction_file, compound_file,
                                   verbose=True):
    """
    Helper function to write reaction and compound files parsed from Chemoton.

    Input:
    - html_reactions (list): list of tuples of integers of the form [n1,n2,ts] specifying the indices of nodes and
      transition states from the set of compounds to define all elementary reactions in the network.
    - html_compounds (dict): dictionary mapping node/ts indices to the different computed fields that are available

    Output:
    - reaction_file (str): path to the reactions file
    - compounds_file (str): path to the compounds file
    """

    if verbose: print("## writing reactions and compound files")
    # Open a file in write mode
    with open(reaction_file, 'w') as f:
        # Loop through the list and write each tuple to the file
        for item in html_reactions:
            r, p, ts = item
            f.write("{r},{p},{ts}\n".format(r=r, p=p, ts=ts))

    with open(compound_file, 'w') as file:
        file.write(json.dumps(html_compounds))  # use `json.loads` to do the revers

    return None

def read_compound_reactions_files(reaction_file,compounds_file):
    """
    Helper function to read reaction and compound files parsed from Chemoton.

    Input:
    - reaction_file (str): path to the reactions file
    - compounds_file (str): path to the compounds file

    Output:
    - reaction_tuples (list): list of tuples of integers of the form [n1,n2,ts] specifying the indices of nodes and
      transition states from the set of compounds to define all elementary reactions in the network.
    - compounds (dict): dictionary mapping node/ts indices to the different computed fields that are available.
    """

    with open(reaction_file,"r") as freac:
        reaction_tuples = [line.strip().split(",") for line in freac.readlines()]
    with open(compounds_file,"r") as fcomp:
        compounds = json.load(fcomp)

    return reaction_tuples,compounds


def build_dashboard(G,title,outfile,size=(1400,800), layout_function=nx.kamada_kawai_layout,  map_field="energy"):
    """
    Wrapper function to generate HTML visualizations for a given network.

    Input:
    - G (nx.Graph): object as generated from RXReader. For profile support, it should contain a graph["pathList"] property.
    - title (str): title for the visualization.
    - outfile (str): name of the output HTML file.
    - size (tuple): tuple of integers, size of the final visualization in pixels.
    - layout_function (nx.object, optional): Function to generate graph layout.
    - map_field (str): name of the field used for node coloring.

    Output:
    - lay (bokey.obj): Bokeh layout as generated by full_view_layout()
    """

    ### Define sizing
    w1 = int(size[0]*4/7)
    w2 = int(size[0]*3/7)
    wu = int(size[0]/7)
    h = int(size[1]*6/8)

    sizing_dict = {'w1':w1,'w2':w2,'wu':wu,'h':h}

    ### Define custom classes

    style_template = """
    {% block postamble %}
	<script type="text/javascript" src="https://cdn.jsdelivr.net/gh/dgarayr/jsmol_to_bokeh/jsmol_to_bokeh.min.js"></script>
    <style>
    .bk-root .bk-btn-default {
        font-size: 1.2vh;
    }
    .bk-root .bk-input {
        font-size: 1.2vh;
        padding-bottom: 5px;
        padding-top: 5px;
    }
    .bk-root .bk {
        font-size: 1.2vh;
    }
    .bk-root .bk-clearfix{
        padding-bottom: 0.8vh;
    }
    </style>
    {% endblock %}
    """
    posx = layout_function(G)
    # Add model field to all nodes and edges & also vibrations
    arxviz.add_models(G)

    # Bokeh-powered visualization via RXVisualizer
    bk_fig,bk_graph = arxviz.bokeh_network_view(G,positions=posx,graph_title=title,width=w1,height=h,
                                                map_field=map_field,hide_energy=True)

    # bk_graph.selection_policy = bkm.NodesAndLinkedEdges()
    bk_graph.selection_policy = bkm.EdgesAndLinkedNodes()

    ### Modify the hovering tools here to add additional fields, removing the previous ones first
    valid_tools = [tool for tool in bk_fig.tools if tool.description]
    old_hovers = [tool for tool in valid_tools if "hover" in tool.description]
    for tool in old_hovers:
        bk_fig.tools.remove(tool)


    # custom edge hovering to reduce noise
    #
    hover_edgeJS = '''
    var erend = graph.edge_renderer.data_source
    var label1 = String.fromCharCode(916).concat("E1")
    var label2 = String.fromCharCode(916).concat("E2")
    if (cb_data.index.indices.length > 0) {
        var ndx = cb_data.index.indices[0]
        var tsname = erend.data["name"][ndx]
        if (tsname.includes('TSb')){
            hover.tooltips = [["tag","@name"]]
        } else {
            hover.tooltips = [["tag","@name"],["charge","@charge"],
                                ["multiplicity","@multiplicity"],["formula","@formula"],
                                [label1,"@deltaE1"],[label2,"@deltaE2"]]
        }
    }
    '''

    hover_node = bkm.HoverTool(description="Node hover",renderers=[bk_graph.node_renderer],
                               tooltips=[("tag","@name"),("charge","@charge"),("multiplicity","@multiplicity"),
                                         ("formula","@formula")],
                               formatters={"@energy":"printf"})
    bk_fig.add_tools(hover_node)
    hover_edge = bkm.HoverTool(description="Edge hover",renderers=[bk_graph.edge_renderer],
                               formatters={"@energy":"printf"},line_policy="interp")
    hover_edge.callback = bkm.CustomJS(args={"hover":hover_edge,"graph":bk_graph},code=hover_edgeJS)
    bk_fig.add_tools(hover_edge)

    highl_callback = bkm.CustomJS(args={"graph":bk_graph}, code=arxviz.js_callback_dict["highlightNeighbors"])


    lay = arxviz.full_view_layout(bk_fig,bk_graph,sizing_dict=sizing_dict)

    # add a button to the layout
    b_highlight = bkm.Button(label="Highlight neighbors",max_width=int(w1/4),align="center")
    b_highlight.js_on_click(highl_callback)


    sel_row = lay.children[0][0].children[2]
    sel_row.children = sel_row.children[0:2] + [b_highlight] + [sel_row.children[-1]]
    bokeh.plotting.output_file(outfile,title=title,mode="cdn")
    bokeh.plotting.save(lay,template=style_template)

    return lay,bk_fig,bk_graph

def scale_xyz_list(xyz,displ_vector=np.zeros(3)):
    """
    Bohr-to-angstrom scaling of a list of XYZ coordinates of the form [atom, [x, y, z]].

    Input:
    - xyz (list): XYZ coordinates, containing a list [atom, [x,y,z]] with atom being a string and x,y,z floats.
    - displ_vector (np.Array, optional): for translating the geometry.

    Output:
    - xyz_nw (list): scaled XYZ coordinates in the same format as the input.
    """

    bohr_to_ang = 0.529
    xyz_arr = np.array([item[1] for item in xyz]) * bohr_to_ang + displ_vector
    xyz_nw = [[item[0],list(xyz_arr[ii])] for ii,item in enumerate(xyz)]
    return xyz_nw

def xyz_list_to_xyz_block(xyz):
    """
    Transform a list of xyz coordinates [atom, [x, y, z]] into a string block.

    Input:
    - xyz (list): XYZ coordinates, containing a list [atom, [x,y,z]] with atom being a string and x,y,z floats.

    Output:
    - xyz_block (str): newline-joined block of the form a1,x1,y1,z2\na2,x2,y2,z2...
    """

    xyz_block = "\n".join(["%s %.6f %.6f %.6f" % (item[0],*item[1]) for item in xyz])
    return xyz_block

def formula_from_xyz_block(xyz):
    """
    Generates the molecular formula for a given XYZ geometry.

    Input:
    - xyz (list): XYZ coordinates, containing a list [atom, [x,y,z]] with atom being a string and x,y,z floats.

    Output:
    - formula (str): molecular formula from the input geometry.
    """
    labels = [item[0] for item in xyz]
    counter_list = sorted(Counter(labels).items())
    formula = ""
    for atom,ct in counter_list:
        if ct == 1:
            formula += atom
        else:
            formula += "%s%d" % (atom,ct)
    return formula

def sort_edge_names(edge_tuple):
    """
    Helper function to sort edge tuples lexicographically.

    Input:
    - edge_tuple (tuple): edge specification as a pair of node names.

    Output:
    - lexico_tuple (tuple): lexicographically sorted tuple.
    """

    n1,n2 = [int(nd) for nd in edge_tuple]
    srt_pair = sorted([n1,n2])
    lexico_tuple = tuple([str(nd) for nd in srt_pair])
    return lexico_tuple

def process_graph(reaction_list,compounds,dist_adduct=3.0):
    """
    Wrapper function to generate a nx.Graph from a list of reactions and a dictionary of compounds,
    including XYZ-formatted geometries where individual geometries of the species forming adducts are joined.

    Input:
    - reaction_list (list): list of tuples of integers of the form [n1,n2,ts] specifying the indices of nodes and
    transition states from the set of compounds to define all elementary reactions in the network.
    - compounds (dict): dictionary mapping node/ts indices to the different computed fields that are available
    - dist_adduct (float, optional): float, distance in angstrom between the centers of mass of adduct fragments for the
    joined 3D geometry.

    Output:
    - G (nx.Graph): containing network structure and the information required by RXVisualizer module to build the final
     dashboard.
    """
    bohr_to_ang = 0.529177

    G = nx.Graph()
    edge_list = [(item[0],item[1],{"tsidx":item[2]}) for item in reaction_list]
    G.add_edges_from(edge_list)
    node_renaming = {}

    ### Preprocessing compounds: for consistency, convert single elements to 1-element lists
    tgt_vars = ["energy","charge","multiplicity"]
    for comp in compounds.values():
        for vv in tgt_vars:
            if not isinstance(comp[vv],list):
                comp[vv] = [comp[vv]]

    # add node information
    for nd in G.nodes(data=True):
        comp = compounds[nd[0]]
        # nodes will be renamed to allow compounds
        # check for adducts, where both molecules must be brought together -> list of IDs
        if isinstance(comp["crn_id"],list):
            node_name = "+".join(comp["crn_id"])

            xyz_list = comp["xyz"]
            xyz0_arr = np.array([item[1] for item in xyz_list[0]]) * bohr_to_ang
            cntr = xyz0_arr.mean(axis=0)
            xyz0 = [[item[0],list(xyz0_arr[ii])] for ii,item in enumerate(xyz_list[0])]
            xyz_full = xyz0

            for ii,xyz in enumerate(xyz_list[1:]):
                displ_vec = cntr + (ii+1)*dist_adduct
                xyz_arr = np.array([item[1] for item in xyz]) * bohr_to_ang + displ_vec
                xyz_nw = [[item[0],list(xyz_arr[ii])] for ii,item in enumerate(xyz)]
                xyz_full += xyz_nw

        else:
            node_name = comp["crn_id"]
            xyz_list = [comp["xyz"]]
            # scale to angstrom
            xyz_arr = np.array([item[1] for item in xyz_list[0]]) * bohr_to_ang
            xyz_full = [[item[0],list(xyz_arr[ii])] for ii,item in enumerate(xyz_list[0])]

        node_renaming[nd[0]] = node_name
        # add this to the graph, with xyz-block format
        xyz_block = "\n".join(["%s %.6f %.6f %.6f" % (item[0],*item[1]) for item in xyz_full])
        nd[1]["geometry"] = xyz_block

        nd[1]["energy"] = sum(comp["energy"])
        nd[1]["ZPVE"] = 0.0
        nd[1]["name"] = node_name
        nd[1]["degree"] = G.degree(nd[0])
        # handle charge and multiplicity as strings to properly treat fragments
        nd[1]["charge"] = ";".join([str(item) for item in comp["charge"]])
        nd[1]["multiplicity"] = ";".join([str(item) for item in comp["multiplicity"]])
        nd[1]["formula"] = ";".join([formula_from_xyz_block(xyz) for xyz in xyz_list])
        nd[1]["neighbors"] = list(G.neighbors(nd[0]))

    for ii,ed in enumerate(G.edges(data=True)):
        e1,e2 = [sum(compounds[nd]["energy"]) for nd in ed[0:2]]
        if ed[2]["tsidx"] == "None":
            e_ts = max(e1,e2)
            ed[2]["name"] = "TSb_%04d" % ii
            ed[2]["geometry"] = None
            ed[2]["energy"] = 0.0
            ed[2]["ZPVE"] = 0.0
            delta_e1 = (e_ts - e1,ed[0])
            delta_e2 = (e_ts - e2,ed[1])
            ed[2]["deltaE1"] = "%.2f (%s)" % delta_e1
            ed[2]["deltaE2"] = "%.2f (%s)" % delta_e2
            continue
        ts_compound = compounds[ed[2]["tsidx"]]
        xyz_list = [ts_compound["xyz"]]
        geom = scale_xyz_list(xyz_list[0])
        ed[2]["geometry"] = xyz_list_to_xyz_block(geom)
        #ed[2]["name"] = "TS_%04d" % int(ed[2]["tsidx"])
        ed[2]["name"] = ts_compound["crn_id"]

        ### compute activation energy
        e_ts = sum(ts_compound["energy"])
        delta_e1 = (e_ts - e1,ed[0])
        delta_e2 = (e_ts - e2,ed[1])
        ### save string representations
        ed[2]["deltaE1"] = "%.2f (%s)" % delta_e1
        ed[2]["deltaE2"] = "%.2f (%s)" % delta_e2
        ed[2]["energy"] = sum(ts_compound["energy"])
        ed[2]["ZPVE"] = 0.0
        # handle charge and multiplicity as strings to properly treat fragments
        ed[2]["charge"] = ";".join([str(item) for item in ts_compound["charge"]])
        ed[2]["multiplicity"] = ";".join([str(item) for item in ts_compound["multiplicity"]])

        ed[2]["formula"] = ";".join([formula_from_xyz_block(xyz) for xyz in xyz_list])

    ## Apply renaming
    nx.relabel_nodes(G,node_renaming,copy=False)
    # and add neighbors now to ensure right naming
    for nd in G.nodes(data=True):
        nd[1]["neighbors"] = list(G.neighbors(nd[0]))
    return G


def main():

    ####################################################################################################
    # Parameters of the Chemoton exploration
    db_name = "ozone_tme_lcpbe"     # name of the database
    ip = 'localhost'                # ip where the database is hosted
    port = '8889'                   # port to the database
    dict_method = {"method_family": "dft", "method": 'lc-pbe',
                   "basis_set": "def2-svp", "program": "orca", } # method of the exploration

    # Files that only need to be generated once per exploration
    pathfinder_file = "crn_pathfinder.json"             # pathfinder object of the CRN
    reactions_file = "/home/petrusen/reactions3.csv"    # list of reactions
    compounds_file = "/home/petrusen/compounds3.json"   # dictionary of the compounds

    # Parameters for the generation of the html CRN file
    output_file = "network2.html"
    title_html = "Chemoton graph"
    dist_adduct = 3.0
    ####################################################################################################

    # only needs to be run once, then data is stored in reaction and compound files (can be commented)
    reactions,compounds = get_reactions_and_compounds(db_name, ip, port, dict_method,
                        write_pathfinder=False, read_pathfinder=pathfinder_file)
    write_compound_reactions_files(reactions,compounds,reactions_file, compounds_file)

    # generation of the html from the previously generated file
    reactions, compounds = read_compound_reactions_files(reactions_file, compounds_file)
    G = process_graph(reactions,compounds,dist_adduct)
    build_dashboard(G,title_html,output_file,
                               size=(1400,800),
                               layout_function=nx.kamada_kawai_layout,
                               map_field="degree")

if __name__ == '__main__':
    main()
