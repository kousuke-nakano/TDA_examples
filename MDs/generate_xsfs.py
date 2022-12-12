import pickle
from ase.io import read, write

# read dumped files
with open(f"quench_temp_list.pickle", mode='rb') as f:
    quench_temp_list = pickle.load(f)
with open(f"quench_density_list.pickle", mode='rb') as f:
    quench_density_list = pickle.load(f)
with open(f"quench_atoms_list.pickle", mode='rb') as f:
    quench_atoms_list = pickle.load(f)

# glass structures at each temp
for temp, density, atoms in zip(quench_temp_list, quench_density_list, quench_atoms_list):
    print(f"Ave. Temp = {temp:.0f} (K)")
    print(f"Ave. Density = {density} (g/cm^3)")
    print("===structure xsf===")
    print(f"num. atoms = {len(atoms)}")
    print(f"structure is saved as sio2_at_temp_{temp:.0f}_K.xsf")
    write(f"sio2_at_temp_{temp:.0f}_K.xsf", atoms)