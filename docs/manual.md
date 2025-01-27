## Manual for Configuration File

This manual explains the configuration parameters of the config.yaml file for an optimal use 
of VizChemoton.

### 1. Database (`db`)

- **active** (`bool`): Enables (`True`) or disables (`False`) the use of the Mongo-DB. If it is
disabled, then this section is omitted. This feature is meant for scenarios where either the 
Mongo-DB is not reachable (but one has already the reaction.csv and compound.json file), or the 
Mongo-DB is very large so that the querying process comes a serious bottleneck of the workflow. 
- **pathfinder** (`bool`): Enables (`True`) or disables (`False`) the use of json file generated
by PATHFinder for the CRN of choice. This feature aims at speeding up the generation of the html
file, even though it still requires to query the Mongo-DB. 
- **name** (`str`): Name of the Mongo-DB.
- **ip** (`str`): IP address of the Mongo-DB server.
- **port** (`str`): Port number for Mongo-DB communication.

### 2. Computational Method (`method`)

- **method_family** (`str`): Specifies the family of the computational method (e.g., `dft`).
- **method** (`str`): Name of the computational method (e.g., `lc-pbe`).
- **basis_set** (`str`): Basis set used in the calculation (e.g., `def2-svp`).
- **program** (`str`): Quantum chemistry program used (e.g., `orca`).

### 3. File Paths (`files`)

#### pathfinder
- **path** (`str`): Path to the json file containing CRN data.
- **mode** (`str`): Either read a preexisting file (`read`) or write a new file (`write`)

#### reactions
- **path** (`str`): Path to the csv file containing reaction data.
- **mode** (`str`): Either read a preexisting file (`read`) or write a new file (`write`)

#### compounds
- **path** (`str`): Path to the json file containing compound data.
- **mode** (`str`): Either read a preexisting file (`read`) or write a new file (`write`)

### 4. Graph Settings (`graph`)

- **dist_adduct** (`float`): Distance threshold for adduct detection.
- **size** (`list[int, int]`): Graph size in pixels (`[width, height]`).
- **layout** (`str`): Graph layout algorithm (e.g., `kamada_kawai`).
- **map_field** (`str`): Property used for node mapping (e.g., `degree`).


### 5. Output (`output`)

- **file** (`str`): Output file name for the generated network visualization.
- **title** (`str`): Title of the network visualization.


---


