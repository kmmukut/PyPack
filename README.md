# 🧬 PyPack

A Python-based toolkit for constructing custom molecular systems from SMILES strings and packing them into simulation-ready 3D boxes for LAMMPS simulations. Comes with a Tkinter GUI for convenient setup and execution.

---

## 🚀 Features

* Convert SMILES strings to 3D molecules with Open Babel
* Estimate box size based on user-defined target density
* Randomly orient and place molecules in a non-overlapping periodic box
* Export `.data` files compatible with LAMMPS ReaxFF simulations
* GUI to simplify molecule setup and parameter input

---

## 🧰 Requirements

### 🐍 Python Modules

Install all required dependencies using pip:

```bash
pip install numpy scipy openbabel ase
```

> **Note:** All required tools are installable via Python; no system-level package managers (like `apt` or `brew`) are necessary.

---

## 🗂️ File Structure

```
.
├── pypack.py         # Core packing logic
├── gui.py            # Tkinter-based GUI
└── reax_input.data   # Example output file (created after run)
```

---

## 🧪 Workflow

### 1. Using the GUI (Recommended)

Launch the graphical interface:

```bash
python gui.py
```

#### GUI Features:

* Add molecules with name, SMILES, and count
* Set parameters: density, temperature, min distance, max placement attempts
* Choose LAMMPS atom style (`charge`, `atomic`, `full`)
* Specify output filename and run

### 2. Command Line Usage

Edit the `molecules` list and parameters inside `pypack.py`, then run:

```bash
python pypack.py
```

The script will:

* Generate and optimize molecule geometries
* Estimate box size from total molecular mass and desired density
* Attempt to place molecules without overlap using `cKDTree`
* Write the final system as a LAMMPS-compatible `.data` file

---

## 📅 Input Format

Each molecule entry requires:

* `name`: Label used in filenames
* `smiles`: SMILES string for molecular structure
* `count`: Number of molecules

Example:

```python
molecules = [
    {"name": "methane", "smiles": "C", "count": 30},
    {"name": "oxygen", "smiles": "O=O", "count": 12}
]
```

---

## 📤 Output

The output is a LAMMPS `.data` file compatible with ReaxFF and similar force fields. Atom styles supported:

* `charge` (default)
* `atomic`
* `full`

---

## 🧠 Core Algorithm & Collision-Aware Packing

PyPack uses a collision-aware placement strategy based on a **k-dimensional tree (cKDTree)** spatial index from SciPy. This efficiently checks for overlaps during molecule placement within a cubic periodic boundary box.

### Core Highlights:

* **Molecule Rotation**: Each molecule is randomly rotated along all three axes before placement.
* **Placement Attempts**: Each molecule gets up to `max_attempts` (default: 2000) tries to be placed without overlapping any other molecule.
* **Minimum Distance**: A user-defined `min_dist` (e.g., 1.2 Å) is enforced between any two atoms using distance thresholding.
* **Box Size Adjustment**: If placement fails after 5 packing attempts, the box size is increased by **20% incrementally** and retried.

This iterative strategy ensures valid, non-overlapping configurations for complex molecular systems.

---

## ⚠️ Notes

* High-density configurations may fail placement—box size auto-expands over 5 attempts.
* SMILES-based generation relies on Open Babel for geometry optimization.
* Ensure all input molecules have valid SMILES strings.

---

## 📜 License

MIT License

---

## 🙌 Acknowledgments

* [Open Babel](https://openbabel.org/)
* [ASE - Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/)
* [LAMMPS Molecular Dynamics Simulator](https://lammps.org)

---

Built with 💻 for molecular dynamics enthusiasts!
