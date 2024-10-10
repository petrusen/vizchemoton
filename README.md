# Vizchemoton

This code reads a chemical reaction network generated with Chemoton[ref], and generates and html file using
the amktools[ref] in order to depict the reaction network.

# Dependencies

* python                3.6.x
* amktools
* nodejs
* scine-database        1.4.0
* scine-heron           1.0.0
* scine-molassembler    2.0.1
* scine-sparrow         5.1.0
* scine-utilities       9.0.0

# Instalation

```bash

# Create a conda enviroment
conda create --name my_env python=3.6
conda install -c conda-forge nodejs=14

# Clone this repository
git clone <this-repository>

# Install the visualization repository
git clone https://gitlab.com/dgarayr/amk_tools.git
pip install -e ./amk_tools/ 

# Install the reaction exploration repository
git clone https://github.com/qcscine/heron.git
pip install ./heron

# Run the code of this repository to obtain the html
python3 viz_chemoton.py ./reactions_epetrus_241010.csv ./compounds_epetrus_241010.txt

```

# References

[1]- chemoton
[2]- amktools

