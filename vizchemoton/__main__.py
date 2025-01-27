'''
Enric Petrus, December 2024. Added SCINE helper functions to link with the amk-tools generation of html files.
Diego Garay-Ruiz, November 2023. Collection of helper functions to link amk-tools and grrm-tools, generating interactive
HTML dashboards to visualize GRRM-generated reaction networks.
'''

import networkx as nx
from .vizchemoton_module import get_reactions_and_compounds, write_compound_reactions_files, read_compound_reactions_files, process_graph, build_dashboard, load_config

def main():
    # Load configuration
    config = load_config()

    # Parameters from config
    db_active = config["db"]["active"]
    db_pathfinder = config["db"]["pathfinder"]
    db_name = config["db"]["name"]
    ip = config["db"]["ip"]
    port = config["db"]["port"]
    dict_method = config["method"]

    pathfinder_file = config["files"]["pathfinder"]["path"]
    pathfinder_mode = config["files"]["pathfinder"]["mode"]

    reactions_file = config["files"]["reactions"]["path"]
    reactions_mode = config["files"]["reactions"]["mode"]

    compounds_file = config["files"]["compounds"]["path"]
    compounds_mode = config["files"]["compounds"]["mode"]

    output_file = config["output"]["file"]
    title_html = config["output"]["title"]
    verbose = config["output"]["verbose"]

    dist_adduct = config["graph"]["dist_adduct"]
    size = tuple(config["graph"]["size"])
    layout_function = getattr(nx, f"{config['graph']['layout']}_layout")
    map_field = config["graph"]["map_field"]

    if db_active: # the Mongo-DB is reachable

        if pathfinder_mode == 'read': # read the pathfinder object (to speed-up the process)
             if verbose: print("reading the {f} file".format(f=pathfinder_file))
             reactions, compounds = get_reactions_and_compounds(db_name, ip, port, dict_method,
             write_pathfinder=False, read_pathfinder=pathfinder_file)

        elif pathfinder_mode == 'write':  # write the pathfinder object
             if verbose: print("writing the {f} file".format(f=pathfinder_file))
             reactions, compounds = get_reactions_and_compounds(db_name, ip, port, dict_method,
             write_pathfinder=pathfinder_file, read_pathfinder=False)

        # write the reactions and compounds
        if reactions_mode == 'write' and compounds_mode == 'write':
            if verbose: print("writing the {f1} and {f2} files".format(f1=reactions_file, f2=compounds_file))
            write_compound_reactions_files(reactions, compounds, reactions_file, compounds_file)

    else: # the Mongo-DB is not reachable, or not necessary as reactions and compounds are stored in separate files
        if reactions_mode == 'read' and compounds_mode == 'read':
            if verbose: print("reading the {f1} and {f2} files".format(f1=reactions_file, f2=compounds_file))
            reactions, compounds = read_compound_reactions_files(reactions_file, compounds_file)
            G = process_graph(reactions, compounds, dist_adduct)
            build_dashboard(G, title_html, output_file, size=size, layout_function=layout_function, map_field=map_field)


if __name__ == '__main__':
    main()
