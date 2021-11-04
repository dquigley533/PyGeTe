import numpy as np
import ase.spacegroup as spg
import PyGeTe as en

############################################
# Create the alpha GeTe structure with ASE #
############################################
A = 4.1719  # Lattice constant
C = 10.710

atomsGeTe = spg.crystal(["Ge", "Te"], [(0.0, 0.0, 0.763), (0.0, 0.0, 0.237)],
                        spacegroup=160, cellpar=[A, A, C, 90, 90, 120],
                        size=(5, 5, 2))

###########################################
# Calculate energy per atom with PyGeTe   #
###########################################

# Set species
speciesGeTe = [en.symdict[s] for s in atomsGeTe.get_chemical_symbols()]

# Make sure initial neighbour list is up to date
en.compute_neighbour_list(atomsGeTe.positions, 
                          atomsGeTe.get_cell(), speciesGeTe)
                
# Compute initial energy                   
total_energy = en.compute_model_energy(atomsGeTe.positions, 
                          atomsGeTe.get_cell(), speciesGeTe)
                

print("Energy per atom of initial structure = %.3f eV" %(total_energy/atomsGeTe.get_global_number_of_atoms()))


