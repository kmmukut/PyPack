import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import pypack
import os

class MoleculeRow:
    def __init__(self, parent, remove_callback):
        self.frame = tk.Frame(parent)
        self.name_var = tk.StringVar()
        self.smi_var = tk.StringVar()
        self.count_var = tk.IntVar(value=1)

        tk.Entry(self.frame, textvariable=self.name_var, width=15).pack(side=tk.LEFT, padx=2)
        tk.Entry(self.frame, textvariable=self.smi_var, width=20).pack(side=tk.LEFT, padx=2)
        tk.Spinbox(self.frame, from_=1, to=1000, textvariable=self.count_var, width=5).pack(side=tk.LEFT, padx=2)
        tk.Button(self.frame, text="Remove", command=self._remove).pack(side=tk.LEFT, padx=2)
        self.remove_callback = remove_callback

    def _remove(self):
        self.frame.destroy()
        self.remove_callback(self)

    def pack(self, **kwargs):
        self.frame.pack(**kwargs)

    def get(self):
        return {"name": self.name_var.get(), "smiles": self.smi_var.get(), "count": self.count_var.get()}

class PackGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Molecular Packer GUI")

        # Molecule list
        mol_frame = tk.LabelFrame(self, text="Molecules")
        mol_frame.pack(fill=tk.X, padx=10, pady=5)
        header = tk.Frame(mol_frame)
        tk.Label(header, text="Name", width=15).pack(side=tk.LEFT)
        tk.Label(header, text="SMILES", width=20).pack(side=tk.LEFT)
        tk.Label(header, text="Count", width=5).pack(side=tk.LEFT)
        header.pack()

        self.rows = []
        self.rows_container = tk.Frame(mol_frame)
        self.rows_container.pack()
        tk.Button(mol_frame, text="Add molecule", command=self.add_row).pack(pady=5)

        # Parameters
        param_frame = tk.LabelFrame(self, text="Parameters")
        param_frame.pack(fill=tk.X, padx=10, pady=5)
        self.density_var = tk.DoubleVar(value=0.1)
        self.temp_var = tk.DoubleVar(value=300)
        self.min_dist_var = tk.DoubleVar(value=1.2)
        self.max_att_var = tk.IntVar(value=2000)
        self.atom_style_var = tk.StringVar(value="charge")

        params = [
            ("Density (g/cm³)", self.density_var),
            ("Temp (K)", self.temp_var),
            ("Min dist (Å)", self.min_dist_var),
            ("Max attempts", self.max_att_var)
        ]
        for label, var in params:
            row = tk.Frame(param_frame)
            tk.Label(row, text=label, width=15).pack(side=tk.LEFT)
            tk.Entry(row, textvariable=var, width=10).pack(side=tk.LEFT)
            row.pack(pady=2)

        # Atom style dropdown
        style_row = tk.Frame(param_frame)
        tk.Label(style_row, text="Atom style", width=15).pack(side=tk.LEFT)
        ttk.Combobox(style_row, textvariable=self.atom_style_var,
                     values=["atomic", "charge", "molecular"], width=10).pack(side=tk.LEFT)
        style_row.pack(pady=2)

        # Output file selection
        out_frame = tk.Frame(self)
        out_frame.pack(fill=tk.X, padx=10, pady=5)
        tk.Label(out_frame, text="Output file:", width=15).pack(side=tk.LEFT)
        self.out_var = tk.StringVar(value="reax_input.data")
        tk.Entry(out_frame, textvariable=self.out_var, width=25).pack(side=tk.LEFT)
        tk.Button(out_frame, text="Browse", command=self.browse_file).pack(side=tk.LEFT, padx=2)

        # Run button
        tk.Button(self, text="Run", command=self.run).pack(pady=10)

        # Initialize with one row
        self.add_row()

    def add_row(self):
        row = MoleculeRow(self.rows_container, self.remove_row)
        row.pack(fill=tk.X, pady=2)
        self.rows.append(row)

    def remove_row(self, row):
        self.rows.remove(row)

    def browse_file(self):
        f = filedialog.asksaveasfilename(defaultextension=".data", filetypes=[("LAMMPS data", "*.data"), ("All files", "*")])
        if f:
            self.out_var.set(f)

    def run(self):
        try:
            mol_defs = [r.get() for r in self.rows if r.get()["name"] and r.get()["smiles"]]
            if not mol_defs:
                raise ValueError("Add at least one molecule definition.")

            # Call build_system with positional args matching pypack signature
            system = pypack.build_system(
                mol_defs,
                self.density_var.get(),
                self.temp_var.get(),
                self.min_dist_var.get(),
                self.max_att_var.get()
            )

            # Write output with selected atom_style
            pypack.write(self.out_var.get(), system, format="lammps-data", atom_style=self.atom_style_var.get(), masses=True)
            messagebox.showinfo("Success", f"LAMMPS input written to:\n{self.out_var.get()}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    app = PackGUI()
    app.mainloop()
