# --- Enhanced gui.py with ttk styling, grid layout, and cleaner spacing ---
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import pypack as pypack
import os
import sys
import threading

class ConsoleRedirect:
    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, message):
        self.text_widget.after(0, self.text_widget.insert, tk.END, message)
        self.text_widget.after(0, self.text_widget.see, tk.END)

    def flush(self):
        pass

class MoleculeRow:
    def __init__(self, parent, remove_callback):
        self.frame = ttk.Frame(parent)
        self.name_var = tk.StringVar()
        self.input_mode = tk.StringVar(value="SMILES")
        self.smi_var = tk.StringVar()
        self.file_var = tk.StringVar()
        self.count_var = tk.IntVar(value=1)

        ttk.Entry(self.frame, textvariable=self.name_var, width=15).grid(row=0, column=0, padx=5, pady=2)
        mode_menu = ttk.Combobox(self.frame, textvariable=self.input_mode, values=["SMILES", "File"], width=8, state="readonly")
        mode_menu.grid(row=0, column=1, padx=5)
        mode_menu.bind("<<ComboboxSelected>>", self._toggle_input_mode)

        self.smi_entry = ttk.Entry(self.frame, textvariable=self.smi_var, width=22)
        self.file_entry = ttk.Entry(self.frame, textvariable=self.file_var, width=22)
        self.browse_btn = ttk.Button(self.frame, text="Browse", command=self._browse_file)

        self.smi_entry.grid(row=0, column=2, padx=5)

        ttk.Spinbox(self.frame, from_=1, to=1000, textvariable=self.count_var, width=5).grid(row=0, column=3, padx=5)
        ttk.Button(self.frame, text="Remove", command=self._remove).grid(row=0, column=4, padx=5)
        self.remove_callback = remove_callback

    def _toggle_input_mode(self, event=None):
        self.smi_entry.grid_forget()
        self.file_entry.grid_forget()
        self.browse_btn.grid_forget()
        if self.input_mode.get() == "SMILES":
            self.smi_entry.grid(row=0, column=2, padx=5)
        else:
            self.file_entry.grid(row=0, column=2, padx=5)
            self.browse_btn.grid(row=0, column=5, padx=5)

    def _browse_file(self):
        filepath = filedialog.askopenfilename(filetypes=[("Molecule files", "*.xyz *.mol *.sdf *.pdb"), ("All files", "*.*")])
        if filepath:
            self.file_var.set(filepath)

    def _remove(self):
        self.frame.destroy()
        self.remove_callback(self)

    def pack(self, **kwargs):
        self.frame.pack(**kwargs)

    def get(self):
        smiles = self.smi_var.get().strip() if self.input_mode.get() == "SMILES" else ""
        file_path = self.file_var.get().strip() if self.input_mode.get() == "File" else ""
        return {
            "name": self.name_var.get().strip(),
            "smiles": smiles,
            "file": file_path,
            "count": self.count_var.get()
        }

class PackGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Molecular Packer GUI")
        style = ttk.Style()
        style.theme_use('clam')

        self.format_extensions = {
            "lammps-data": ".data",
            "xyz": ".xyz"
        }

        mol_frame = ttk.LabelFrame(self, text="Molecules")
        mol_frame.pack(fill=tk.X, padx=10, pady=5)
        header = ttk.Frame(mol_frame)
        ttk.Label(header, text="Name", width=15).grid(row=0, column=0, padx=5)
        ttk.Label(header, text="Mode", width=8).grid(row=0, column=1, padx=5)
        ttk.Label(header, text="Input", width=22).grid(row=0, column=2, padx=5)
        ttk.Label(header, text="Count", width=5).grid(row=0, column=3, padx=5)
        header.pack()

        self.rows = []
        self.rows_container = ttk.Frame(mol_frame)
        self.rows_container.pack()
        ttk.Button(mol_frame, text="Add molecule", command=self.add_row).pack(pady=5)

        param_frame = ttk.LabelFrame(self, text="Parameters")
        param_frame.pack(fill=tk.X, padx=10, pady=5)
        self.density_var = tk.DoubleVar(value=0.1)
        self.temp_var = tk.DoubleVar(value=300)
        self.min_dist_var = tk.DoubleVar(value=1.2)
        self.max_att_var = tk.IntVar(value=2000)
        self.atom_style_var = tk.StringVar(value="charge")
        self.output_format_var = tk.StringVar(value="lammps-data")
        self.output_format_var.trace("w", self.sync_extension)

        for label, var in [
            ("Density (g/cm³)", self.density_var),
            ("Temp (K)", self.temp_var),
            ("Min dist (Å)", self.min_dist_var),
            ("Max attempts", self.max_att_var)
        ]:
            row = ttk.Frame(param_frame)
            ttk.Label(row, text=label, width=15).pack(side=tk.LEFT, padx=5)
            ttk.Entry(row, textvariable=var, width=10).pack(side=tk.LEFT, padx=5)
            row.pack(pady=2)

        style_row = ttk.Frame(param_frame)
        ttk.Label(style_row, text="Atom style", width=15).pack(side=tk.LEFT, padx=5)
        ttk.Combobox(style_row, textvariable=self.atom_style_var,
                     values=["atomic", "charge", "full"], width=10).pack(side=tk.LEFT)
        style_row.pack(pady=2)

        format_row = ttk.Frame(param_frame)
        ttk.Label(format_row, text="Output format", width=15).pack(side=tk.LEFT, padx=5)
        ttk.Combobox(format_row, textvariable=self.output_format_var,
                     values=list(self.format_extensions.keys()), width=15).pack(side=tk.LEFT)
        format_row.pack(pady=2)

        out_frame = ttk.LabelFrame(self, text="Output")
        out_frame.pack(fill=tk.X, padx=10, pady=5)
        ttk.Label(out_frame, text="Output file:", width=15).pack(side=tk.LEFT, padx=5)
        self.out_var = tk.StringVar(value="reax_input.data")
        ttk.Entry(out_frame, textvariable=self.out_var, width=30).pack(side=tk.LEFT, padx=5)
        ttk.Button(out_frame, text="Browse", command=self.browse_file).pack(side=tk.LEFT, padx=5)

        ttk.Button(self, text="Run", command=self.run).pack(pady=10)

        console_frame = ttk.LabelFrame(self, text="Console Output")
        console_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        self.console = tk.Text(console_frame, height=10, bg="#f9f9f9", relief=tk.SUNKEN)
        self.console.pack(fill=tk.BOTH, padx=5, pady=5, expand=True)

        sys.stdout = sys.stderr = ConsoleRedirect(self.console)
        self.add_row()

    def add_row(self):
        row = MoleculeRow(self.rows_container, self.remove_row)
        row.pack(fill=tk.X, pady=2)
        self.rows.append(row)

    def remove_row(self, row):
        self.rows.remove(row)

    def browse_file(self):
        ext = self.format_extensions.get(self.output_format_var.get(), ".data")
        f = filedialog.asksaveasfilename(defaultextension=ext, filetypes=[("All files", "*.*")])
        if f:
            self.out_var.set(f)

    def sync_extension(self, *args):
        base, _ = os.path.splitext(self.out_var.get())
        ext = self.format_extensions.get(self.output_format_var.get(), ".data")
        self.out_var.set(base + ext)

    def run(self):
        def worker():
            try:
                mol_defs = [r.get() for r in self.rows if r.get()["name"] and (r.get()["smiles"] or r.get()["file"])]
                if not mol_defs:
                    raise ValueError("Add at least one molecule definition with SMILES or file.")

                system = pypack.build_system(
                    mol_defs,
                    self.density_var.get(),
                    self.temp_var.get(),
                    self.min_dist_var.get(),
                    self.max_att_var.get()
                )

                fmt = self.output_format_var.get()
                kwargs = dict(format=fmt)
                if fmt == "lammps-data":
                    kwargs.update(atom_style=self.atom_style_var.get(), masses=True)

                pypack.write(self.out_var.get(), system, **kwargs)
                print(f"\nOutput written to: {self.out_var.get()}\n")
                messagebox.showinfo("Success", f"Output written to:\n{self.out_var.get()}")
            except Exception as e:
                messagebox.showerror("Error", str(e))

        threading.Thread(target=worker).start()

if __name__ == "__main__":
    app = PackGUI()
    app.mainloop()
