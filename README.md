# VizChemoton - Visualization of Reaction Networks Generated with Chemoton

This package allows the **visualization of chemical reaction networks (CRNs)** generated with [Chemoton](https://github.com/qcscine/chemoton)[1] 
through the generation of standalone HTML files with the [amk-tools](https://github.com/dgarayr/amk_tools),[2]ultimately allowing the user to interact 
with the network **just via browser**.


![Example Image](./docs/example_crn_html.png)

*The image above shows which type of visualization (network on the left, structure on 
the right) the HTML file permits*


## Introduction to Chemoton

The research group of Markus Reiher at ETH Zürich has developed a the software enviroment [SCINE](https://github.com/qcscine) ("Software for Chemical Interaction Networks")[3] which pursues the performance of quantum chemical calculations with special focus on algorithmic stability, automation, interactivity, efficiency and error control. In this context, Chemoton is the SCINE module in charge of exploring CRNs in a fully automated fashion based on first-principles. 

The reaction network results from Chemoton can be visualized with the graphical user interface developed in the same SCINE framework, named [Heron](https://github.com/qcscine/heron).[4] This GUI has many functionalities beyond visualization, such as monitoring and steering the exploration, as well as creating json files to make the whole process reproducible. 

> Notice that here we offer a complementary visualization tool aimed at providing an interactive visualization of the CRN without any previous dependencies. We highlight, however, that *vizchemoton* is just meant for creating a light-weight and convenient visualization format of the CRN, thus making it easily accessible to the scientific community.     


## Introduction to amk-tools

The research group of Carles Bo at ICIQ (Tarragona, Spain) developed a library, named *amk-tools* allowing the reading, processing and visualizing the reaction networks generated by [AutoMeKin](https://github.com/emartineznunez/AutoMeKin),[5] an automated reaction discovery program developed by Emilio Martinez Nuñez (Galicia, Spain).[5] 


> Here we tailor the application of *amk-tools* to the visualization of CRNs generated with *Chemoton*. 


## Code overview

In the following Figure we describe the different modules -described above- that VizChemoton
uses to generate the html files.

![Example Image](./docs/vizchemoton_architecture.png)

Notice that the chemical data present in the network.html can also be uploaded in the FAIR repository
[ioChem-BD](http://dx.doi.org/10.19061/iochem-bd-6-430) via using the wrapper [json2orca](https://github.com/gruberlopez/json2orca). 

## Dependencies

Dependencies are detailed in the requierements.txt file and in the table below.

| Package              | Version |
|----------------------|---------|
| python               | 3.6.x   |
| amktools             | x.x.x   |
| nodejs               | 14.0.0  |
| scine-database       | 1.4.0   |
| scine-heron          | 1.0.0   |
| scine-molassembler   | 2.0.1   |
| scine-sparrow        | 5.1.0   |
| scine-utilities      | 9.0.0   |


## Instalation 

```bash

# Create a conda enviroment
conda create --name my_env python=3.6
conda install -c conda-forge nodejs=14

# Install Chemoton
git clone https://github.com/qcscine/chemoton.git
cd chemoton
git checkout 3.1.0
pip install -r requirements.txt
pip install .
cd ..

# Install amk-tools
git clone https://gitlab.com/dgarayr/amk_tools.git
pip install -e ./amk_tools/

# Install VizChemoton
git clone https://github.com/petrusen/vizchemoton.git
pip install -r ./vizchemoton/requirements.txt
pip install ./vizchemoton

```

## Usage of the HTML files

```bash

# Edit the main file
vi ./vizchemoton/__main__.py

# Run the code of this repository to obtain the html
python3 -m viz_chemoton.py 

```

## References

1. a) Gregor N. Simm, Markus Reiher. *J. Chem. Theory Comput.* **2017**, 13, 12, 6108-6119.  b) Jan P. Unsleber, Stephanie A. Grimmel, Markus Reiher. *J. Chem. Theory Comput.* **2022**, 18, 9, 5393-5409
2. Diego Garay-Ruiz, Moises Alvarez-Moreno, Carles Bo, Emilio Martinez-Nunez. *ACS Phys. Chem Au* **2022**, 2, 3, 225-236.
3. T. Weymuth, J. P. Unsleber, P. L. Türtscher, M. Steiner, J.-G. Sobez, C. H. Müller, M. Mörchen,
V. Klasovita, S. A. Grimmel, M. Eckhoff, K.-S. Csizi, F. Bosia, M. Bensberg, M. Reiher, *J. Chem. Phys.*, **2024**, *160*, 222501.
4. C. H. Müller, M. Steiner, J. P. Unsleber, T. Weymuth, M. Bensberg, K.-S. Csizi, M. Mörchen, P. L. Türtscher, M. Reiher, *J. Phys. Chem. A*, **2024**, 128, 9028−9044.
(DOI: 10.48550/arXiv.2406.09541)
5. E. Martínez-Núñez, G. L. Barnes, D. R. Glowacki, S. Kopec, D. Peláez, A. Rodríguez, R. Rodríguez-Fernández, R. J. Shannon, J. J. P. Stewart, P. G. Tahoces, S. A. Vazquez, *J. Comput. Chem.* **2021**, 42(28), 2036.
