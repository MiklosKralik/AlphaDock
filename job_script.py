# general troubleshooting
from AlphaDock.DockJob import DJ
from AlphaDock import get_full_path, clear_cache_folder
import random
import os

# reproducible
seed = 10550962
protein_paths = ['MSM7_MCE.pdb']

ligand_paths = ['cholesterol.sdf']
j = 1

for protein in protein_paths:
    protein = 'structures/' + protein

    for ligand in ligand_paths:
        ligand = 'structures/' + ligand
        clear_cache_folder()
        Job = DJ(protein, ligand, str(j), vina='vina')
        Job.to_pdbqt()
        Job.strip_protein()
        Job.box(padding=1.0)
        Job.dock(num_modes=9, exhaustiveness=8, seed=seed, energy_range=3)
        Job.surroundings()
        Job.plot()
        Job.extract_cache('./outputs')
        j += 1


#clear_cache_folder()
#
# JobTest = DJ('structures/MTB4_MCE.pdb', 'structures/cholesterol.sdf', 'Final',
#                vina='AlphaDock/vina')
#
# JobTest.to_pdbqt()
# JobTest.strip_protein()
# JobTest.box(padding=1.0)
# JobTest.dock(num_modes=9, exhaustiveness=8, seed=seed, energy_range=3)
# JobTest.surroundings()
# JobTest.plot(show=True)
# JobTest.extract_cache('./outputs')



## TODO:
# parametrize docking & everything else including outputs
# create plotting functions
# incorperate motif search
