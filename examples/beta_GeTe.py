import numpy as np
import ase.spacegroup as spg
import PyGeTe.energy as en
from PyGeTe import symdict

###########################################
# Create the beta GeTe structure with ASE #
###########################################
LAT = 6.009  # Lattice constant
n = 4        # Repeats

betaGeTe = spg.crystal(["Ge", "Te"], [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)],
                       spacegroup=225, cellpar=[LAT, LAT, LAT, 90, 90, 90],
                       size=(n, n, n))

###########################################
# Calculate energy per atom with PyGeTe   #
###########################################

# Set species
betaGeTeSpecies = [symdict[s] for s in betaGeTe.get_chemical_symbols()]

# Make sure initial neighbour list is up to date
en.compute_neighbour_list(betaGeTe.positions, 
                          betaGeTe.get_cell(), betaGeTeSpecies)
                
# Compute initial energy                   
total_energy = en.compute_model_energy(betaGeTe.positions,
                                       betaGeTe.get_cell(), betaGeTeSpecies)


print("Energy per atom = %.3f eV" %(total_energy/betaGeTe.get_global_number_of_atoms()))
