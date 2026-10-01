#!/usr/bin/env python3
"""
PCA of amino-acid (physicochemical class) substitution profiles
Rotavirus A structural (VP) and non-structural (NSP) proteins.

INPUT:
    Select a folder containing two Excel workbooks.

    The script automatically identifies:
        1. Structural workbook:
           filename containing "structural"

        2. Non-structural workbook:
           filename containing "nonstructural"
           or "non-structural"

    Each sheet/workbook holds one 4x4 block per protein:

             <PROTEIN> | From_Group | Hydrophobic | Polar | Basic | Acidic
                       | Hydrophobic| n  n  n  n
                       | Polar       | n  n  n  n
                       | Basic       | n  n  n  n
                       | Acidic      | n  n  n  n

    Rows = From class
    Columns = To class
    Values = substitution counts

METHOD:
    1. Each protein -> vector of 16 counts
       (all From->To cells, diagonal included).

    2. Counts divided by the protein's own total -> relative frequencies.
       The 16 values sum to 1.

    3. Each of the 16 columns z-scored across proteins.

    4. PCA using scikit-learn.

    PCA is run separately for structural and non-structural proteins.

OUTPUT:
    Figure6A_Structural_PCA.png
    Figure6B_NonStructural_PCA.png
    PCA_results.xlsx

    All output files are saved in the selected input folder.

REQUIRES:
    pandas
    numpy
    scikit-learn
    matplotlib
    openpyxl
    tkinter
"""

import os
import sys
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

import tkinter as tk
from tkinter import filedialog, messagebox


# ============================================================
# SETTINGS
# ============================================================

CLASSES = [
    "Hydrophobic",
    "Polar",
    "Basic",
    "Acidic"
]

COLS = [
    f"{a[:3]}->{b[:3]}"
    for a in CLASSES
    for b in CLASSES
]


# ============================================================
# SELECT INPUT FOLDER
# ============================================================

def select_folder():

    root = tk.Tk()
    root.withdraw()

    root.attributes("-topmost", True)

    folder = filedialog.askdirectory(
        title="Select folder containing Structural and Non-structural Excel files"
    )

    root.destroy()

    if not folder:
        print("\nNo folder selected. Program stopped.")
        sys.exit()

    return folder


# ============================================================
# FIND EXCEL FILES
# ============================================================

def find_excel_files(folder):

    excel_files = []

    for file in os.listdir(folder):

        if file.lower().endswith((".xlsx", ".xls")):

            # Ignore temporary Excel files
            if file.startswith("~$"):
                continue

            excel_files.append(file)

    if len(excel_files) == 0:

        messagebox.showerror(
            "No Excel files found",
            "No Excel (.xlsx or .xls) files were found in the selected folder."
        )

        sys.exit()

    structural_file = None
    nonstructural_file = None

    for file in excel_files:

        name = file.lower().replace("_", "").replace(" ", "")

        # ----------------------------------------------------
        # Structural file
        # ----------------------------------------------------

        if (
            "structural" in name
            and "nonstructural" not in name
        ):

            structural_file = file

        # ----------------------------------------------------
        # Non-structural file
        # ----------------------------------------------------

        if (
            "nonstructural" in name
            or "non-structural" in file.lower()
            or "non_structural" in file.lower()
        ):

            nonstructural_file = file

    # --------------------------------------------------------
    # Check that both files were found
    # --------------------------------------------------------

    if structural_file is None:

        messagebox.showerror(
            "Structural file not found",
            "Could not find an Excel file containing 'structural' "
            "in its filename.\n\n"
            "Example:\n"
            "Structural.xlsx"
        )

        sys.exit()

    if nonstructural_file is None:

        messagebox.showerror(
            "Non-structural file not found",
            "Could not find an Excel file containing "
            "'nonstructural' or 'non-structural' in its filename.\n\n"
            "Example:\n"
            "NonStructural.xlsx"
        )

        sys.exit()

    return (
        os.path.join(folder, structural_file),
        os.path.join(folder, nonstructural_file)
    )


# ============================================================
# READ 4x4 BLOCKS
# ============================================================

def read_blocks(xlsx_path):

    """
    Return {protein: 4x4 numpy array}
    from a workbook with stacked 4x4 blocks.
    """

    raw = pd.read_excel(
        xlsx_path,
        header=None
    )

    blocks = {}

    for i in range(len(raw)):

        label = raw.iloc[i, 0]

        if (
            isinstance(label, str)
            and label.strip()
            and label != "RV genes"
        ):

            block = raw.iloc[
                i + 1:i + 5,
                2:6
            ].astype(float).values

            assert (
                block.shape == (4, 4)
                and not np.isnan(block).any()
            ), f"Bad block for {label}"

            blocks[label.strip()] = block

    return blocks


# ============================================================
# CONVERT TO 16-VARIABLE MATRIX
# ============================================================

def to_matrix(blocks):

    """
    Proteins x 16 matrix of raw counts.
    """

    return pd.DataFrame(
        {
            k: v.flatten()
            for k, v in blocks.items()
        },
        index=COLS
    ).T


# ============================================================
# RUN PCA
# ============================================================

def run_pca(counts):

    # --------------------------------------------------------
    # Step 2:
    # Convert counts to relative frequencies
    # --------------------------------------------------------

    freq = counts.div(
        counts.sum(axis=1),
        axis=0
    )

    # --------------------------------------------------------
    # Step 3:
    # Z-score the 16 variables across proteins
    # --------------------------------------------------------

    z = StandardScaler().fit_transform(
        freq.values
    )

    # --------------------------------------------------------
    # Step 4:
    # PCA
    # --------------------------------------------------------

    pca = PCA().fit(z)

    # --------------------------------------------------------
    # PCA scores
    # --------------------------------------------------------

    scores = pd.DataFrame(
        pca.transform(z),
        index=counts.index,
        columns=[
            f"PC{i+1}"
            for i in range(pca.n_components_)
        ]
    )

    # --------------------------------------------------------
    # Explained variance
    # --------------------------------------------------------

    ev = pd.DataFrame(
        {
            "PC": scores.columns,

            "Variance_%":
                pca.explained_variance_ratio_ * 100,

            "Cumulative_%":
                np.cumsum(
                    pca.explained_variance_ratio_
                ) * 100
        }
    )

    # --------------------------------------------------------
    # PCA loadings
    # --------------------------------------------------------

    loadings = pd.DataFrame(
        pca.components_.T
        * np.sqrt(pca.explained_variance_),

        index=COLS,

        columns=scores.columns
    )

    return (
        freq,
        scores,
        ev,
        loadings
    )


# ============================================================
# PLOT PCA
# ============================================================

def plot(
    scores,
    ev,
    title,
    label,
    out_png
):

    fig, ax = plt.subplots(
        figsize=(6.4, 5.2)
    )

    cmap = plt.get_cmap("tab10")

    for i, name in enumerate(scores.index):

        x = scores.loc[
            name,
            "PC1"
        ]

        y = scores.loc[
            name,
            "PC2"
        ]

        ax.scatter(
            x,
            y,
            s=90,
            color=cmap(i),
            edgecolor="k",
            lw=.6,
            zorder=3
        )

        ax.annotate(
            name,
            (x, y),
            xytext=(7, 6),
            textcoords="offset points",
            fontweight="bold"
        )

    # --------------------------------------------------------
    # Zero lines
    # --------------------------------------------------------

    ax.axhline(
        0,
        color="#ccc",
        lw=.7
    )

    ax.axvline(
        0,
        color="#ccc",
        lw=.7
    )

    # --------------------------------------------------------
    # Axis labels
    # --------------------------------------------------------

    ax.set_xlabel(
        f"PC1 ({ev.loc[0, 'Variance_%']:.1f}% of variance)"
    )

    ax.set_ylabel(
        f"PC2 ({ev.loc[1, 'Variance_%']:.1f}% of variance)"
    )

    # --------------------------------------------------------
    # Title
    # --------------------------------------------------------

    ax.set_title(
        f"{label}. PCA of amino acid substitution profiles\n"
        f"{title}",
        fontsize=10.5
    )

    ax.margins(.18)

    ax.spines[
        ["top", "right"]
    ].set_visible(False)

    plt.tight_layout()

    plt.savefig(
        out_png,
        dpi=220
    )

    plt.close()


# ============================================================
# MAIN ANALYSIS
# ============================================================

def main():

    print("\n==============================================")
    print(" Rotavirus A Substitution Profile PCA")
    print("==============================================\n")

    # --------------------------------------------------------
    # Select folder
    # --------------------------------------------------------

    print("Please select the folder containing the Excel files...")

    folder = select_folder()

    print("\nSelected folder:")
    print(folder)

    # --------------------------------------------------------
    # Find the two Excel files
    # --------------------------------------------------------

    f_struct, f_nonstruct = find_excel_files(
        folder
    )

    print("\nStructural file:")
    print(os.path.basename(f_struct))

    print("\nNon-structural file:")
    print(os.path.basename(f_nonstruct))

    # --------------------------------------------------------
    # Define groups
    # --------------------------------------------------------

    groups = {

        "Structural": (
            read_blocks(f_struct),
            "A",
            "Figure6A_Structural_PCA.png"
        ),

        "Non-structural": (
            read_blocks(f_nonstruct),
            "B",
            "Figure6B_NonStructural_PCA.png"
        )
    }

    # --------------------------------------------------------
    # Output Excel file
    # --------------------------------------------------------

    output_excel = os.path.join(
        folder,
        "PCA_results.xlsx"
    )

    # --------------------------------------------------------
    # Write PCA results
    # --------------------------------------------------------

    with pd.ExcelWriter(
        output_excel,
        engine="openpyxl"
    ) as xl:

        for (
            gname,
            (
                blocks,
                label,
                png
            )
        ) in groups.items():

            print("\n----------------------------------------------")
            print(f"Processing {gname}")
            print("----------------------------------------------")

            # ------------------------------------------------
            # Convert blocks to matrix
            # ------------------------------------------------

            counts = to_matrix(
                blocks
            )

            print(
                f"Number of proteins: {len(counts)}"
            )

            # ------------------------------------------------
            # Run PCA
            # ------------------------------------------------

            (
                freq,
                scores,
                ev,
                loadings
            ) = run_pca(
                counts
            )

            # ------------------------------------------------
            # Plot
            # ------------------------------------------------

            output_png = os.path.join(
                folder,
                png
            )

            plot(
                scores,
                ev,
                f"{gname} proteins",
                label,
                output_png
            )

            # ------------------------------------------------
            # Excel sheet prefix
            # ------------------------------------------------

            tag = (
                "S"
                if gname == "Structural"
                else "NS"
            )

            # ------------------------------------------------
            # Input counts
            # ------------------------------------------------

            counts.assign(
                Total=counts.sum(axis=1)
            ).to_excel(
                xl,
                sheet_name=f"{tag}_input_counts"
            )

            # ------------------------------------------------
            # Relative frequencies
            # ------------------------------------------------

            freq.to_excel(
                xl,
                sheet_name=f"{tag}_frequencies"
            )

            # ------------------------------------------------
            # Explained variance
            # ------------------------------------------------

            ev.to_excel(
                xl,
                sheet_name=f"{tag}_explained_var",
                index=False
            )

            # ------------------------------------------------
            # Loadings
            # ------------------------------------------------

            loadings.to_excel(
                xl,
                sheet_name=f"{tag}_loadings"
            )

            # ------------------------------------------------
            # PCA scores
            # ------------------------------------------------

            scores.to_excel(
                xl,
                sheet_name=f"{tag}_scores"
            )

            # ------------------------------------------------
            # Print summary
            # ------------------------------------------------

            print(
                f"{gname}: "
                f"n={len(counts)} proteins | "
                f"PC1 "
                f"{ev.loc[0, 'Variance_%']:.1f}% | "
                f"PC2 "
                f"{ev.loc[1, 'Variance_%']:.1f}%"
            )

    # ========================================================
    # FINISHED
    # ========================================================

    print("\n==============================================")
    print(" PCA ANALYSIS COMPLETED")
    print("==============================================")

    print("\nOutput folder:")
    print(folder)

    print("\nGenerated files:")

    print(
        "1. Figure6A_Structural_PCA.png"
    )

    print(
        "2. Figure6B_NonStructural_PCA.png"
    )

    print(
        "3. PCA_results.xlsx"
    )

    print("\n==============================================\n")

    # --------------------------------------------------------
    # Show completion message
    # --------------------------------------------------------

    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)

    messagebox.showinfo(
        "PCA Analysis Completed",
        "PCA analysis completed successfully!\n\n"
        "The following files were saved in the selected folder:\n\n"
        "• Figure6A_Structural_PCA.png\n"
        "• Figure6B_NonStructural_PCA.png\n"
        "• PCA_results.xlsx"
    )

    root.destroy()


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    main()