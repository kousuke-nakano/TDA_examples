import pickle
from time import time
import os
import warnings
warnings.filterwarnings('ignore')

import numpy as np

from gph import ripser_parallel
from gtda.homology import VietorisRipsPersistence
from gtda.homology._utils import _postprocess_diagrams

from ase.io import write

root_dir=os.getcwd()
md_dir=os.path.join(root_dir, "MDs")
pd_dir=os.path.join(root_dir, "PDs"); makedirs(pd_dir, exist_ok=True)

# read temperatures (averaged)
with open(os.path.join(md_dir, f"quench_temp_list.pickle"), mode='rb') as f:
    quench_temp_list=pickle.load(f)

# read atoms (ASE-atoms instances)
with open(os.path.join(md_dir, f"quench_atoms_list.pickle"), mode='rb') as f:
    quench_atoms_list=pickle.load(f)

_ = [atoms.wrap() for atoms in quench_atoms_list] # warp!!!!

# extracting only Si and O atoms
confs = {temp:atoms for temp, atoms in zip(quench_temp_list, quench_atoms_list)}

chem_spec_list = [quench_atoms.get_chemical_symbols() for quench_atoms in quench_atoms_list]

indices_Si_list = [[i for i in range(len(chem_spec)) if chem_spec[i]=='Si'] for chem_spec in chem_spec_list]
indices_O_list = [[i for i in range(len(chem_spec)) if chem_spec[i]=='O'] for chem_spec in chem_spec_list]

Si_atoms = {k:v[indices_Si_list[i]] for i, (k, v) in enumerate(confs.items())}
O_atoms = {k:v[indices_O_list[i]] for i, (k, v) in enumerate(confs.items())}
Si_O_atoms = {k:v for k, v in confs.items()}

Si_atom_pos = {k:v.get_positions() for k, v in Si_atoms.items()}
O_atom_pos = {k:v.get_positions() for k, v in O_atoms.items()}
Si_O_atom_pos = {k:v.get_positions() for k, v in Si_O_atoms.items()}

#"""
t1 = time()
print("dgmSi")
VR = VietorisRipsPersistence(homology_dimensions=[0, 1, 2], n_jobs=-1, collapse_edges=True)
Si_ordered_keys = sorted([float(x) for x in Si_atom_pos.keys()])
diagramSi = VR.fit_transform([Si_atom_pos[k] for k in Si_ordered_keys])
dgmSi = {k: v for k, v in zip(Si_ordered_keys, diagramSi)}
with open(os.path.join(pd_dir, f"diagrams_Sigtda_{prefix}.pickle"), mode='wb') as f:
    pickle.dump(dgmSi, f)
print(f'Time to compute dgmSi: {time()-t1}')
#"""


t2 = time()
print("dgmO")
VR = VietorisRipsPersistence(homology_dimensions=[0, 1, 2], n_jobs=-1, collapse_edges=True)
O_ordered_keys = sorted([float(x) for x in O_atom_pos.keys()])
diagramO = VR.fit_transform([O_atom_pos[k] for k in O_ordered_keys])
dgmO = {k: v for k, v in zip(O_ordered_keys, diagramO)}
with open(os.path.join(pd_dir, f"diagrams_Ogtda_{prefix}.pickle"), mode='wb') as f:
    pickle.dump(dgmO, f)
print(f'Time to compute dgmO: {time()-t2}')
"""

"""
t3 = time()
print("dgmSiO")
VR = VietorisRipsPersistence(homology_dimensions=[0, 1, 2], n_jobs=-1, collapse_edges=True)
SiO_ordered_keys = sorted([float(x) for x in Si_O_atom_pos.keys()])
diagramSiO = VR.fit_transform([Si_O_atom_pos[k] for k in SiO_ordered_keys])
dgmSiO = {k: v for k, v in zip(SiO_ordered_keys, diagramSiO)}
with open(os.path.join(pd_dir, f"diagrams_SiOgtda_{prefix}.pickle"), mode='wb') as f:
    pickle.dump(dgmSiO, f)
print(f'Time to compute dgmSiO: {time()-t3}')


print(f'END')
