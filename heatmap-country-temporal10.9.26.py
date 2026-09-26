"""
================================================================================
 ROTAVIRUS GENE GENOTYPE HEATMAP GENERATOR
================================================================================
Builds a country x year heatmap showing the FULL genotype composition
(every genotype present, not just the dominant one) for a single rotavirus
gene dataset. Blank cells (no data for that country/year) are left empty --
no gray filler box is drawn.

HOW TO USE FOR EACH OF THE 11 GENES:
  1. Just press Run (F5). A file-browser window will pop up -- select that
     gene's Excel file. (If no popup appears, e.g. on some Linux/Spyder
     setups, it will instead ask you to paste the file path directly in
     the console -- just type or paste it there and hit Enter.)
  2. If needed, change TYPE_COL below to the column you want to color the
     heatmap by (see note under TYPE_COL).
  3. Change OUTPUT_NAME so each gene's figure gets its own file.
  4. The heatmap pops up in the Plots pane, and PNG + PDF files are saved
     next to this script.

Gene <-> genotype letter reference:
  VP7=G   VP4=P   VP6=I   VP1=R   VP2=C   VP3=M
  NSP1=A  NSP2=N  NSP3=T  NSP4=E  NSP5/6=H
================================================================================
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ==============================================================================
# =====================        1. INPUT SETTINGS        ======================
# ==============================================================================

# Leave this as None to get a file-browser popup (or a console prompt) each
# time you run the script. Or, if you'd rather hardcode a path instead of
# browsing, paste it here between the quotes, e.g. r"C:\data\VP6.xlsx"
EXCEL_FILE = None

SHEET_NAME = 0                     # sheet name (string) or index (0 = first sheet)

ACCESSION_COL = "Accession_ID"     # column with the accession/sequence ID
COUNTRY_COL   = "Country"          # column with the country
YEAR_COL      = "Year"             # column with the collection year

# Which column to color the heatmap by. All 11 gene files currently only
# carry the parental strain's G/P metadata (no per-gene genotype yet), so
# leave this as "Combined_Type" until each gene has been through the
# MAFFT + reference-strain classification step and has its own genotype
# column (e.g. "I", "R", "C", "M", "N", "T", "E", "H") -- switch to that
# column name once it exists and has real diversity to show.
TYPE_COL = "Combined_Type"

YEAR_BIN_SIZE = 3                  # group years into windows of this many years (1 = no binning)
OUTPUT_NAME   = "gene_heatmap"     # output file name, no extension (change per gene, e.g. "VP6_heatmap")
OUTPUT_DIR    = r"."                # folder to save the output into

# ==============================================================================
# =====================     2. NOTHING BELOW NEEDS EDITING    ================
# ==============================================================================


def get_excel_file_path(preset_path):
    """Return a valid Excel file path -- via popup browser, console prompt, or the preset value."""
    if preset_path:
        if not os.path.isfile(preset_path):
            raise FileNotFoundError(f"EXCEL_FILE is set but not found on disk: {preset_path}")
        return preset_path

    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        path = filedialog.askopenfilename(
            title="Select the gene's Excel file",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")]
        )
        root.destroy()
        if not path:
            raise ValueError("No file selected (dialog was closed/cancelled).")
        return path
    except ImportError:
        print("(No popup file browser available in this environment -- paste the file path instead.)")
        path = input("Full path to the Excel file: ").strip().strip('"')
        if not os.path.isfile(path):
            raise FileNotFoundError(f"File not found: {path}")
        return path


def build_heatmap(excel_file, sheet_name, country_col, year_col, type_col,
                   year_bin_size, output_name, output_dir):

    # ---- Load & clean ----
    df = pd.read_excel(excel_file, sheet_name=sheet_name)
    sub = df.dropna(subset=[year_col, country_col, type_col]).copy()
    sub[year_col] = sub[year_col].astype(int)

    # ---- Bin years ----
    sub["YearBin"] = (sub[year_col] // year_bin_size * year_bin_size).astype(int)
    sub["YearBin_label"] = sub["YearBin"].astype(str) + "-" + (sub["YearBin"] + year_bin_size - 1).astype(str)
    cols_sorted = sorted(sub["YearBin_label"].unique(), key=lambda x: int(x.split("-")[0]))

    # ---- All countries, no grouping, sorted by total sequence count ----
    rows_order = sub[country_col].value_counts().index.tolist()

    # ---- Full per-cell breakdown: every type present, not just the dominant one ----
    full = {}
    for country in rows_order:
        full[country] = {}
        for yb in cols_sorted:
            cell = sub[(sub[country_col] == country) & (sub["YearBin_label"] == yb)]
            if len(cell) == 0:
                full[country][yb] = None
            else:
                vc = cell[type_col].value_counts().to_dict()
                full[country][yb] = {k: int(v) for k, v in vc.items()}

    # ---- Color palette: one color per distinct type, most frequent first ----
    all_types = sorted(set(t for r in full.values() for c in r.values() if c for t in c.keys()))
    freq = {}
    for r in full.values():
        for c in r.values():
            if c:
                for t, n in c.items():
                    freq[t] = freq.get(t, 0) + n
    types_sorted = sorted(all_types, key=lambda t: -freq[t])

    base = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors[:8]) + list(plt.get_cmap("tab20c").colors[:8])
    if len(types_sorted) > len(base):
        print(f"NOTE: {len(types_sorted)} distinct '{type_col}' values found but only {len(base)} "
              f"distinct palette colors -- some colors will repeat.")
    palette = {t: base[i % len(base)] for i, t in enumerate(types_sorted)}

    # ---- Plot ----
    fig, ax = plt.subplots(figsize=(13, 0.42 * len(rows_order) + 2), dpi=200)
    pad = 0.06
    for yi, country in enumerate(rows_order):
        y = len(rows_order) - 1 - yi
        for xi, yb in enumerate(cols_sorted):
            cell = full[country][yb]
            if not cell:
                continue  # blank cell: draw nothing, no filler box
            total = sum(cell.values())
            items = sorted(cell.items(), key=lambda kv: -kv[1])
            bar_x0 = xi + pad
            bar_w = 1.0 - 2 * pad
            bar_y0 = y + 0.20
            bar_h = 0.42
            xcur = bar_x0
            for t, n in items:
                seg_w = bar_w * n / total
                ax.add_patch(mpatches.Rectangle((xcur, bar_y0), seg_w, bar_h,
                                                 facecolor=palette[t], edgecolor="white", linewidth=0.3))
                xcur += seg_w
            ax.text(xi + 0.5, y + 0.78, f"{total}", ha="center", va="center", fontsize=6.2, color="#333")

    ax.set_xlim(0, len(cols_sorted))
    ax.set_ylim(0, len(rows_order))
    ax.set_xticks(np.arange(len(cols_sorted)) + 0.5)
    ax.set_xticklabels(cols_sorted, fontsize=8.5)
    ax.set_yticks(np.arange(len(rows_order)) + 0.5)
    ax.set_yticklabels(rows_order[::-1], fontsize=7.8)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for xi in range(len(cols_sorted) + 1):
        ax.axvline(xi, color="#eee", linewidth=0.7, zorder=0)
    for yi in range(len(rows_order) + 1):
        ax.axhline(yi, color="#eee", linewidth=0.7, zorder=0)

    ax.set_title(
        f"{type_col} composition by country and year  (n={len(sub)} sequences, "
        f"{len(rows_order)} countries)\n"
        f"Each bar segment = one {type_col} value, sized by its share of sequences in that cell; "
        f"number = total sequences; blank = no data",
        fontsize=10.5, pad=14
    )

    legend_handles = [mpatches.Patch(facecolor=palette[t], label=t) for t in types_sorted]
    ax.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.025),
              ncol=9, fontsize=7.2, frameon=False, handlelength=1.1, handleheight=1.1, columnspacing=1.0)

    plt.tight_layout()

    png_path = os.path.join(output_dir, f"{output_name}.png")
    pdf_path = os.path.join(output_dir, f"{output_name}.pdf")
    plt.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.savefig(pdf_path, bbox_inches="tight", facecolor="white")

    print(f"Saved:\n  {png_path}\n  {pdf_path}")
    print(f"Countries: {len(rows_order)} | Year bins: {len(cols_sorted)} | "
          f"Distinct '{type_col}' values: {len(types_sorted)}")

    plt.show()  # opens in Spyder's Plots pane


if __name__ == "__main__":
    excel_path = get_excel_file_path(EXCEL_FILE)
    print(f"Using file: {excel_path}")
    build_heatmap(excel_path, SHEET_NAME, COUNTRY_COL, YEAR_COL, TYPE_COL,
                  YEAR_BIN_SIZE, OUTPUT_NAME, OUTPUT_DIR)