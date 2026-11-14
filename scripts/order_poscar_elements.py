#!/home/vempatti/.conda/envs/ase/bin/python

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar

# Load structure
structure = Structure.from_file("POSCAR")

# Define your desired element order here
custom_element_order = ["Au", "Cu", "Rh", "Pt", "O", "H", "C"]

# Sort the structure based on the custom order
structure = structure.get_sorted_structure(
    key=lambda site: custom_element_order.index(site.specie.symbol)
    if site.specie.symbol in custom_element_order
    else len(custom_element_order)
)

# Overwrite the original POSCAR file
Poscar(structure).write_file("POSCAR")

