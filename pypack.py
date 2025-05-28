import numpy as np
import math
import os
import argparse
import json

from openbabel import pybel
from ase import Atoms
from ase.io import read as ase_read, write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from scipy.spatial import cKDTree

def generate_obabel_xyz(smiles, name):
    mol = pybel.readstring("smi", smiles)
    mol.addh()
    mol.make3D()
    mol.localopt(forcefield="mmff94")
    xyz_file = f"{name}.xyz"
    mol.write("xyz", xyz_file, overwrite=True)
    return xyz_file

def obabel_to_ase(xyz_file):
    return ase_read(xyz_file)

def get_molecule_mass(ase_atoms):
    return ase_atoms.get_masses().sum()

def estimate_box_size_from_density(mol_masses, total_counts, density_gcc):
    total_mass_amu = sum(m * c for m, c in zip(mol_masses, total_counts))
    total_mass_g = total_mass_amu * 1.66054e-24
    volume_cm3 = total_mass_g / density_gcc
    volume_A3 = volume_cm3 * 1e24
    return round(volume_A3 ** (1/3), 2)

def random_rotate(mol):
    mol.rotate(np.random.rand() * 360, 'x', rotate_cell=False)
    mol.rotate(np.random.rand() * 360, 'y', rotate_cell=False)
    mol.rotate(np.random.rand() * 360, 'z', rotate_cell=False)
    return mol

def place_molecules_kdtree(molecules, box_size, min_dist, max_attempts):
    placed_atoms = Atoms(cell=[box_size]*3, pbc=[True]*3)

    for mol_index, mol in enumerate(molecules):
        for attempt in range(max_attempts):
            mol_copy = mol.copy()
            random_rotate(mol_copy)
            displacement = np.random.rand(3) * box_size
            mol_copy.translate(displacement)
            positions = mol_copy.get_positions() % box_size
            mol_copy.set_positions(positions)

            if len(placed_atoms) == 0:
                placed_atoms += mol_copy
                break

            all_positions = placed_atoms.get_positions()
            new_positions = mol_copy.get_positions()
            tree = cKDTree(all_positions, boxsize=box_size)
            distances, _ = tree.query(new_positions, k=1, distance_upper_bound=min_dist)

            if np.all(distances > min_dist):
                placed_atoms += mol_copy
                break
        else:
            print(f" Failed to place molecule {mol_index+1} after {max_attempts} attempts.")
            return None

    return placed_atoms

def build_system(molecule_defs, density_gcc, temperature, min_dist=1.2, max_attempts=2000):
    ase_molecules = []
    mol_masses = []
    total_counts = []
    flat_molecule_list = []

    for mol_def in molecule_defs:
        name = mol_def["name"]
        smiles = mol_def.get("smiles", "")
        file_path = mol_def.get("file", "")
        count = mol_def["count"]

        if file_path:
            ase_mol = ase_read(file_path)
        elif smiles:
            xyz_file = generate_obabel_xyz(smiles, name)
            ase_mol = obabel_to_ase(xyz_file)
        else:
            raise ValueError(f"Neither SMILES nor file provided for molecule '{name}'")

        mass = get_molecule_mass(ase_mol)

        ase_molecules.append(ase_mol)
        mol_masses.append(mass)
        total_counts.append(count)
        flat_molecule_list += [ase_mol.copy() for _ in range(count)]

    box_size = estimate_box_size_from_density(mol_masses, total_counts, density_gcc)

    for attempt in range(5):
        print(f"Packing attempt {attempt+1}: box = {box_size:.2f} Å")
        system = place_molecules_kdtree(flat_molecule_list, box_size, min_dist, max_attempts)
        if system:
            system.set_initial_charges([0.0] * len(system))
            MaxwellBoltzmannDistribution(system, temperature_K=temperature)
            return system
        else:
            print("Expanding box and retrying...")
            box_size *= 1.2

    raise RuntimeError("Failed to pack molecules after 5 attempts.")

# ---------------------------
# CLI Support
# ---------------------------
def cli():
    parser = argparse.ArgumentParser(description="Run molecular packing via CLI.")
    parser.add_argument("--json", help="Path to JSON file with molecule definitions")
    parser.add_argument("--format", default="lammps-data", choices=["lammps-data", "xyz"], help="Output format")
    parser.add_argument("--output", default="reax_input.data", help="Output file name")
    parser.add_argument("--atom_style", default="charge", choices=["atomic", "charge", "full"], help="Atom style (LAMMPS only)")
    parser.add_argument("--density", type=float, default=0.1, help="Target density (g/cm³)")
    parser.add_argument("--temperature", type=float, default=300, help="Temperature (K)")
    parser.add_argument("--min_dist", type=float, default=1.2, help="Minimum interatomic distance (Å)")
    parser.add_argument("--max_attempts", type=int, default=2000, help="Max placement attempts")

    args = parser.parse_args()

    if args.json:
        with open(args.json) as f:
            molecules = json.load(f)
    else:
        molecules = [
            {"name": "methane", "smiles": "C", "count": 30},
            {"name": "oxygen", "smiles": "O=O", "count": 12}
        ]

    system = build_system(
        molecules,
        density_gcc=args.density,
        temperature=args.temperature,
        min_dist=args.min_dist,
        max_attempts=args.max_attempts
    )

    write_args = dict(format=args.format)
    if args.format == "lammps-data":
        write_args.update(atom_style=args.atom_style, masses=True)

    print(f"Writing to {args.output} in {args.format} format...")
    write(args.output, system, **write_args)
    print("Done.")

if __name__ == "__main__":
    cli()
