# 🧬 PyPack

**PyPack** is a Python toolkit for building molecular simulation systems from SMILES strings or structure files. It enables collision-free packing into 3D periodic boxes with support for LAMMPS and XYZ outputs. Includes a modern Tkinter GUI with embedded console.

---

## ✨ Features

* 🧪 Add molecules via SMILES strings or structure files (`.xyz`, `.mol`, `.sdf`)
* 📦 3D structure generation using Open Babel
* 📐 Density-based box estimation
* 🧱 Collision-aware packing using `cKDTree`
* 💾 Export to LAMMPS `.data` or `.xyz`
* 🖥️ GUI with live log output and auto file extension handling
* 🧾 JSON-based CLI input support

---

## 🧰 Requirements

Install dependencies:

```bash
pip install numpy scipy openbabel ase
```

---

## 📁 File Structure

```
.
├── pypack.py               # Core packing + CLI support
├── gui.py                  # Tkinter-based GUI
├── example_molecules.json  # Example JSON molecule input
└── reax_input.data         # Sample LAMMPS output
└── oxygen.xyz              # Sample O2 xyz input
```

---

## 🚀 Workflow

### 1. Graphical Interface (Recommended)

```bash
python gui.py
```

**Features:**

* Add molecules with name, input mode (SMILES or file), and count
* Set parameters: density, temperature, minimum distance, and max attempts
* Choose output format and atom style (if applicable)
* Real-time console feedback
* Auto-adjusts file extension

---

### 2. Command Line Interface

Define molecules in a JSON file and run:

```bash
python pypack.py --json example_molecules.json --format xyz --output packed_system.xyz
```

Example JSON (`example_molecules.json`):

```json
[
  {"name": "methane", "smiles": "C", "count": 30},
  {"name": "oxygen", "file": "oxygen.xyz", "count": 12}
]
```

Or use the built-in list:

> If `--json` is not provided, PyPack falls back to a hardcoded default list inside `pypack.py`, useful for quick tests:
>
> ```python
> molecules = [
>   {"name": "methane", "smiles": "C", "count": 30},
>   {"name": "oxygen", "smiles": "O=O", "count": 12}
> ]
> ```

```bash
python pypack.py --format lammps-data --atom_style charge --output reax_input.data
```

**CLI Options:**

* `--json`: Path to JSON file (overrides default list)
* `--format`: Output format (`lammps-data`, `xyz`)
* `--atom_style`: For `.data` format: `charge`, `atomic`, `full`
* `--output`: Output file name
* `--density`: Target density (g/cm³)
* `--temperature`: Temperature (K)
* `--min_dist`: Minimum distance (Å)
* `--max_attempts`: Max placement tries per molecule

---

## 📄 Output

Exported system:

* `.data`: LAMMPS-compatible, includes atom style and mass
* `.xyz`: For visualization, no atom styles

Only `.data` supports atom_style and mass info.

---

## 🧠 Algorithm Highlights

* 🔄 Random 3D rotation of molecules
* 📏 Enforces user-defined minimum interatomic distance
* 🔁 Each molecule gets multiple placement attempts
* 📦 Auto-box resizing: expands 20% if packing fails (max 5 rounds)

---

## ⚠️ Notes

* SMILES input requires Open Babel
* File input must be readable by ASE
* High-density systems may fail after retries

---

## 📜 License

MIT License

---

## 🙌 Acknowledgments

* [Open Babel](https://openbabel.org/)
* [ASE](https://wiki.fysik.dtu.dk/ase/)
* [LAMMPS](https://lammps.org/)

---

Built with 🧬 for molecular simulation enthusiasts.
