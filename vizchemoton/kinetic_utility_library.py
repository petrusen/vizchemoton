import sys
import numpy as np
import networkx as nx
from math import exp, log10, log
import scine_utilities as utils
import scine_database as db

from scine_database.energy_query_functions import (get_energy_change,
    get_barriers_for_elementary_step_by_type,
    rate_constant_from_barrier
)

from itertools import combinations

def stoich_to_chemical_formula(stoich):
    """
    Converts a tuple with the atom stoichiometry of a Compound/Flask to a Chemical Formula
    """
    nO, nC, nH = stoich
    formula = list()
    for str_z, num_z in zip(['O', 'C', 'H'], stoich):
        if num_z == 0:
            form_z = ""
        elif num_z == 1:
            form_z = str_z
        elif num_z >= 2:
            form_z = "${z}_{a}$".format(z=str_z, a=num_z)
        formula.append(form_z)

    return "".join(formula)


def Erying_K_to_AG(K):
    T, R, = 298.15, 0.001987096
    kb, h = 3.29983 * 10 ** (-27), 1.58367 * 10 ** (-37)
    AG = -T * R * log(K * h / (kb * T))
    return AG


def Erying_AG_to_K(AG):
    T, R, = 298.15, 0.001987096
    kb, h = 3.29983 * 10 ** (-27), 1.58367 * 10 ** (-37)
    K = kb * T / h * exp(- AG / (T * R))
    return K

def get_xyz_from_structure_id(structure, structures, verbose=True):
    s = db.Structure(structure, structures)
    xyz = s.get_atoms()
    cartesian_data = list()
    charge, multiplicity = s.charge, s.get_multiplicity()
    cartesian_data.append([charge, multiplicity])
    for elem, xyzi in zip(xyz.elements, xyz.positions):
        x, y, z = xyzi
        cartesian_data.append([str(elem), x, y, z])


    if verbose: print(s.get_compound())
    return cartesian_data

def get_energy_and_barriers_gibbs(es_id, elementary_steps, model1, structures, properties, model2, es_from_graph, verbose=False):
    if verbose: print("=================================================")
    if verbose: print("Composite Gibbs energy for es_id: ", es_id)
    _energy_k, correc_energy, barriers_k, correc_barriers, not_None = None, None, None, None, None
    _energy_i = get_energy_change(db.ElementaryStep(es_id, elementary_steps), 'electronic_energy',
                                           model1, structures, properties)
    barriers_i = get_barriers_for_elementary_step_by_type(es_from_graph, 'electronic_energy',
                                                                    model1,
                                                                    structures,
                                                                    properties)

    _energy_j = get_energy_change(db.ElementaryStep(es_id, elementary_steps), 'gibbs_free_energy',
                                           model1, structures, properties)
    barriers_j = get_barriers_for_elementary_step_by_type(es_from_graph, 'gibbs_free_energy',
                                                                    model1,
                                                                    structures,
                                                                    properties)
    #if verbose: print("preliminar energies", _energy_k, _energy_i, _energy_j)
    if None in barriers_i or None in barriers_j: # or None == _energy_i or None == _energy_j:
        not_None, _energy, barriers = False, None, (None,None)
    else:
        correc_energy = _energy_j - _energy_i
        correc_barriers = tuple(a-b for a,b in zip(barriers_j, barriers_i))
        if verbose: print("Gibbs Barrier Correction (kj/mol)", correc_barriers)

        _energy_k = get_energy_change(db.ElementaryStep(es_id, elementary_steps), 'electronic_energy',
                                           model2, structures, properties)
        barriers_k = get_barriers_for_elementary_step_by_type(es_from_graph, 'electronic_energy',
                                                                    model2,
                                                                    structures,
                                                                    properties)
        if None in barriers_k:
                       not_None, _energy, barriers = False, None, (None,None)
        else:
                       not_None = True
                       _energy = _energy_k + correc_energy
                       barriers = tuple(a+b for a,b in zip(barriers_k, correc_barriers))
        if verbose: print("Elec Barrier - Model 1 (kj/mol)", barriers_k)
        if verbose: print("Gibbs Corrected - Model 1+2 (kj/mol)", barriers)
    return _energy, barriers, not_None



def get_energy_and_barriers(energy_type, es_id, elementary_steps, model1, structures, properties, es_from_graph):

    _energy = get_energy_change(db.ElementaryStep(es_id, elementary_steps), energy_type,
                                           model1, structures, properties)
    barriers = get_barriers_for_elementary_step_by_type(es_from_graph, energy_type,
                                                                    model1,
                                                                    structures,
                                                                    properties)

    if None in barriers: # or None == _energy_i or None == _energy_j:
        not_None = False
    else:
        not_None = True 

    return _energy, barriers, not_None


