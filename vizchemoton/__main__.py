'''
Enric Petrus, December 2024. Added SCINE helper functions to link with the amk-tools generation of html files.
Diego Garay-Ruiz, November 2023. Collection of helper functions to link amk-tools and grrm-tools, generating interactive
HTML dashboards to visualize GRRM-generated reaction networks.
'''

import networkx as nx
from .vizchemoton_module import get_reactions_and_compounds, write_compound_reactions_files, read_compound_reactions_files, process_graph, build_dashboard


def main():

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
