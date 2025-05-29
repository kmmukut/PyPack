import io
import os
import argparse
import json
import random
import numpy as np
from scipy.spatial import cKDTree

from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem
from ase import Atoms
from ase.io import read as ase_read, write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution


def generate_obabel_xyz(smiles, name):
    """
    Generate a 3D conformation using OpenBabel from a SMILES string,
    write it to `<name>.xyz`, and return an ASE Atoms object.
    """
    mol = pybel.readstring("smi", smiles)
    mol.addh()
    mol.make3D()
    mol.localopt(forcefield="mmff94")
    xyz_file = f"{name}.xyz"
    # write to file
    mol.write("xyz", xyz_file, overwrite=True)
    # read back into ASE
    return ase_read(xyz_file)


def generate_single_conformer(smiles):
    """
    Use RDKit to embed a single conformer, UFF-optimize it, then refine with OpenBabel.
    Returns an ASE Atoms object (does NOT write to disk).
    """
    rdmol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    cid = AllChem.EmbedMolecule(rdmol)
    if cid < 0:
        raise RuntimeError("RDKit embedding failed")
    AllChem.UFFOptimizeMolecule(rdmol, confId=cid)
    xyz_block = Chem.MolToXYZBlock(rdmol, confId=cid)
    ob_mol = pybel.readstring("xyz", xyz_block)
    ob_mol.addh()
    ob_mol.localopt(forcefield="mmff94")
    xyz_data = ob_mol.write("xyz")
    xyz_str = xyz_data.decode('utf-8') if isinstance(xyz_data, bytes) else xyz_data
    return ase_read(io.StringIO(xyz_str), format='xyz')


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


def place_molecules_kdtree(mol_defs, box_size, min_dist, max_attempts):
    placed = Atoms(cell=[box_size]*3, pbc=[True]*3)
    for idx, mol_def in enumerate(mol_defs):
        for attempt in range(max_attempts):
            if mol_def.get('smiles'):
                mol = generate_single_conformer(mol_def['smiles'])
            else:
                mol = ase_read(mol_def['file'])
            mol = random_rotate(mol)
            disp = np.random.rand(3) * box_size
            mol.translate(disp)
            coords = mol.get_positions() % box_size
            mol.set_positions(coords)
            if len(placed) == 0:
                placed += mol
                break
            tree = cKDTree(placed.get_positions(), boxsize=box_size)
            dists, _ = tree.query(coords, k=1, distance_upper_bound=min_dist)
            if np.all(dists > min_dist):
                placed += mol
                break
        else:
            print(f"Failed to place molecule {idx+1} after {max_attempts} attempts.")
            return None
    return placed


def build_system(molecule_defs, density_gcc, temperature, min_dist=1.2, max_attempts=2000):
    # Generate and save one template XYZ for each SMILES species
    mol_masses = []
    total_counts = []
    for mol_def in molecule_defs:
        if mol_def.get('smiles'):
            # writes <name>.xyz to disk
            sample = generate_obabel_xyz(mol_def['smiles'], mol_def['name'])
        else:
            sample = ase_read(mol_def['file'])
        mol_masses.append(get_molecule_mass(sample))
        total_counts.append(mol_def['count'])
    box_size = estimate_box_size_from_density(mol_masses, total_counts, density_gcc)
    # Flatten definitions list
    flat_defs = []
    for mol_def in molecule_defs:
        flat_defs += [mol_def] * mol_def['count']
    # Attempt packing
    for attempt in range(5):
        print(f"Packing attempt {attempt+1}: box = {box_size:.2f} Å")
        system = place_molecules_kdtree(flat_defs, box_size, min_dist, max_attempts)
        if system is not None:
            system.set_initial_charges([0.0] * len(system))
            MaxwellBoltzmannDistribution(system, temperature_K=temperature)
            return system
        print("Expanding box and retrying...")
        box_size *= 1.2
    raise RuntimeError("Failed to pack molecules after 5 attempts.")


def cli():
    parser = argparse.ArgumentParser(description="Pack molecules with random on-the-fly conformers.")
    parser.add_argument("--json", help="JSON with molecule definitions")
    parser.add_argument("--format", default="lammps-data", choices=["lammps-data", "xyz"])
    parser.add_argument("--output", default="reax_input.data")
    parser.add_argument("--atom_style", default="charge", choices=["atomic", "charge", "full"])
    parser.add_argument("--density", type=float, default=0.1, help="Target density (g/cm³)")
    parser.add_argument("--temperature", type=float, default=300, help="Temperature (K)")
    parser.add_argument("--min_dist", type=float, default=1.2, help="Min interatomic distance (Å)")
    parser.add_argument("--max_attempts", type=int, default=2000, help="Max attempts per molecule")
    args = parser.parse_args()
    if args.json:
        molecule_defs = json.load(open(args.json))
    else:
        molecule_defs = [
            {"name": "methane", "smiles": "C", "count": 30},
            {"name": "oxygen", "smiles": "O=O", "count": 12}
        ]
    system = build_system(
        molecule_defs,
        density_gcc=args.density,
        temperature=args.temperature,
        min_dist=args.min_dist,
        max_attempts=args.max_attempts
    )
    write_kwargs = {"format": args.format}
    if args.format == "lammps-data":
        write_kwargs.update(atom_style=args.atom_style, masses=True)
    print(f"Writing to {args.output} in {args.format} format...")
    write(args.output, system, **write_kwargs)
    print("Done.")

if __name__ == "__main__":
    cli()
