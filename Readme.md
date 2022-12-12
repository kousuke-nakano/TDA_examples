# Scripts and data for the paper "Topological data analysis for revealing structural origin of density anomalies in silica glass"

## How to prepare an environment for the script
```
conda create --name TDA python=3.8.11
conda activate TDA
conda install pip
pip install -r requirements.txt
```

## PD_compute.py
`PD_compute.py` is a python script to compute PDs using atomic configurations obtained from MD simulations (`MDs/quench_atoms_list.pickle`). Here, `quench_atoms_list.pickle` is a pickled file containing the atomic configurations as ASE `atoms` instances.

> Warning: PD calculations for the -O-O- and -Si-O- networks need a large amount of memory. One should try the -Si-Si- network first.

## TDA_analysys.ipynb
`TDA_analysys.ipynb` is a jupyter notebook implementing the TDA workflow. The steps of the workflow are consistent with `Fig. 2` in the paper.
