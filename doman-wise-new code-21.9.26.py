#!/usr/bin/env python
# =========================================================
# Domain Variability Analysis Pipeline
# Input columns:
# Gene
# Residue range
# Domain/Motif
# Highly variable residues
# =========================================================

import os
import sys
import re
import pandas as pd
import matplotlib.pyplot as plt

try:
    import seaborn as sns
    HAS_SNS = True
except ImportError:
    HAS_SNS = False
    print("Seaborn not installed. Heatmaps will use matplotlib only.")


# -------------------------
# Select file
# -------------------------
def get_file_path():
    """
    Get the Excel file path across environments
    (CLI, Tkinter dialog, or console input).
    """

    # 1. Command-line argument
    if len(sys.argv) > 1:
        candidate = sys.argv[1]

        if os.path.isfile(candidate):
            return candidate
        else:
            print(f"Warning: Command line path does not exist: {candidate}")

    # 2. Try Tkinter dialog
    try:
        from tkinter import Tk, filedialog

        root = Tk()
        root.withdraw()
        root.attributes("-topmost", True)
        root.update()

        file = filedialog.askopenfilename(
            title="Select Excel File",
            filetypes=[("Excel files", "*.xlsx *.xls")],
            parent=root
        )

        root.destroy()

        if file:
            return file
        else:
            print("No file selected via Tkinter dialog.")

    except Exception as e:
        print(f"Tkinter file dialog unavailable or failed ({e}).")

    # 3. Fallback to default filename
    default_filename = "domain wise- 26.7.26.xlsx"

    if os.path.isfile(default_filename):
        print(f"Using default file found in directory: {default_filename}")
        return default_filename

    # 4. Manual entry
    print("Falling back to manual file path entry.")

    file = input(
        "Enter the full path to your Excel file: "
    ).strip().strip('"')

    if not file or not os.path.isfile(file):
        raise SystemExit("No valid file provided. Exiting.")

    return file


# -------------------------
# Robust Range Parsing
# -------------------------
def parse_domain_length(range_str):
    """
    Calculates total length from range strings like:
    '1–332'
    '333–488 & 524–594'
    '141–150, 208–221'
    """

    if pd.isna(range_str):
        return 0

    clean_str = str(range_str).replace('–', '-').strip()

    total_length = 0

    # Split across commas or ampersands
    parts = re.split(r'[,&]', clean_str)

    for part in parts:

        part = part.strip()

        # Look for start-end patterns
        match = re.search(r'(\d+)\s*-\s*(\d+)', part)

        if match:
            start, end = int(match.group(1)), int(match.group(2))
            total_length += max(0, end - start + 1)

        else:
            # Look for single integer residues
            match_single = re.search(r'^\d+$', part)

            if match_single:
                total_length += 1

    return total_length


# -------------------------
# Load and Preprocess Data
# -------------------------
file = get_file_path()

print(f"\nProcessing file: {file}")

# Load Excel
df = pd.read_excel(file)


# -------------------------
# Handle Gene column
# -------------------------
if "Unnamed: 0" in df.columns:

    if (
        df["Unnamed: 0"].isnull().all()
        or df["Unnamed: 0"].dtype == 'object'
    ):

        df.rename(
            columns={"Unnamed: 0": "Gene"},
            inplace=True
        )


# Remove accidental spaces from column names
df.columns = df.columns.str.strip()


# Forward fill Gene names
if "Gene" in df.columns:
    df["Gene"] = df["Gene"].ffill()


# -------------------------
# Required columns
# -------------------------
required = [
    "Gene",
    "Residue range",
    "Domain/Motif",
    "Highly variable residues"
]


missing = [
    c for c in required
    if c not in df.columns
]

if missing:

    print("Detected columns:", list(df.columns))

    raise ValueError(
        f"Missing required column(s): {missing}"
    )


# -------------------------
# Handle swapped columns
# -------------------------
for idx, row in df.iterrows():

    val_res = str(row["Residue range"])
    val_dom = str(row["Domain/Motif"])

    # If Domain/Motif contains digits with dashes/commas
    # and Residue range contains letters
    if (
        re.search(r'\d', val_dom)
        and re.search(r'[a-zA-Z]', val_res)
    ):

        df.loc[idx, "Residue range"] = val_dom
        df.loc[idx, "Domain/Motif"] = val_res


# -------------------------
# Calculate Domain Length
# -------------------------
df["Domain Length"] = (
    df["Residue range"]
    .apply(parse_domain_length)
)


# -------------------------
# Calculate Highly Variable %
# -------------------------
df["Highly Variable (%)"] = df.apply(

    lambda r: round(
        (
            r["Highly variable residues"]
            / r["Domain Length"]
            * 100
        ),
        2
    )
    if r["Domain Length"] > 0
    else 0.0,

    axis=1
)


# -------------------------
# Output folder
# -------------------------
folder = os.path.dirname(file) or "."


# -------------------------
# Save Processed Excel
# -------------------------
out_excel = os.path.join(
    folder,
    "Domain_Variability_Output.xlsx"
)

df.to_excel(
    out_excel,
    index=False
)

print(
    f"[✓] Saved updated Excel: {out_excel}"
)


# -------------------------
# Labels for figures
# -------------------------
labels = (
    df["Gene"].astype(str)
    + " - "
    + df["Domain/Motif"].astype(str)
)

x = range(len(df))

w = 0.5


# =========================================================
# Figure 1: Highly Variable Residue Counts
# =========================================================

plt.figure(figsize=(16, 7))

plt.bar(
    x,
    df["Highly variable residues"],
    width=w,
    label="Highly Variable"
)

plt.xticks(
    x,
    labels,
    rotation=90,
    fontsize=8
)

plt.ylabel("Highly Variable Residue Count")

plt.title(
    "Highly Variable Residues per Domain"
)

plt.tight_layout()

plt.legend()

plt.savefig(
    os.path.join(
        folder,
        "Figure1_Highly_Variable_Counts.png"
    ),
    dpi=600
)

plt.savefig(
    os.path.join(
        folder,
        "Figure1_Highly_Variable_Counts.pdf"
    )
)

plt.show()

plt.close()


# =========================================================
# Figure 2: Highly Variable Percentage
# =========================================================

plt.figure(figsize=(16, 7))

plt.bar(
    x,
    df["Highly Variable (%)"],
    width=w,
    label="Highly Variable %"
)

plt.xticks(
    x,
    labels,
    rotation=90,
    fontsize=8
)

plt.ylabel("Highly Variable Residues (%)")

plt.title(
    "Highly Variable Residue Percentage per Domain"
)

plt.tight_layout()

plt.legend()

plt.savefig(
    os.path.join(
        folder,
        "Figure2_Highly_Variable_Percentage.png"
    ),
    dpi=600
)

plt.savefig(
    os.path.join(
        folder,
        "Figure2_Highly_Variable_Percentage.pdf"
    )
)

plt.show()

plt.close()


# =========================================================
# Figure 3: Highly Variable Residue Composition
# =========================================================

plt.figure(figsize=(16, 7))

plt.bar(
    labels,
    df["Highly variable residues"],
    label="Highly Variable"
)

plt.xticks(
    rotation=90,
    fontsize=8
)

plt.ylabel("Highly Variable Residue Count")

plt.title(
    "Highly Variable Residues across Domains"
)

plt.tight_layout()

plt.legend()

plt.savefig(
    os.path.join(
        folder,
        "Figure3_Highly_Variable_Composition.png"
    ),
    dpi=600
)

plt.savefig(
    os.path.join(
        folder,
        "Figure3_Highly_Variable_Composition.pdf"
    )
)

plt.show()

plt.close()


# =========================================================
# Figure 4: Highly Variable Residues Heatmap
# =========================================================

var = df.pivot_table(
    index="Gene",
    columns="Domain/Motif",
    values="Highly variable residues",
    aggfunc="sum"
)


plt.figure(figsize=(12, 6))


if HAS_SNS:

    sns.heatmap(
        var,
        annot=True,
        cmap="magma",
        fmt=".0f",
        cbar_kws={
            'label': 'Highly Variable Residues'
        }
    )

else:

    plt.imshow(
        var.fillna(0),
        aspect="auto",
        cmap="magma"
    )

    plt.colorbar(
        label="Highly Variable Residues"
    )

    plt.xticks(
        range(len(var.columns)),
        var.columns,
        rotation=90
    )

    plt.yticks(
        range(len(var.index)),
        var.index
    )


plt.title(
    "Highly Variable Residues Heatmap"
)

plt.tight_layout()

plt.savefig(
    os.path.join(
        folder,
        "Figure4_Heatmap_Variable.png"
    ),
    dpi=600
)

plt.savefig(
    os.path.join(
        folder,
        "Figure4_Heatmap_Variable.pdf"
    )
)

plt.show()

plt.close()


# =========================================================
# Figure 5: Gene Summary
# =========================================================

summary = (
    df.groupby("Gene")
    [["Highly variable residues"]]
    .sum()
)


summary.to_excel(
    os.path.join(
        folder,
        "Summary_Table.xlsx"
    )
)

print(
    f"[✓] Saved summary table: "
    f"{os.path.join(folder, 'Summary_Table.xlsx')}"
)


plt.figure(figsize=(10, 5))

summary.plot(
    kind="bar",
    figsize=(10, 5)
)

plt.ylabel(
    "Highly Variable Residues"
)

plt.title(
    "Total Highly Variable Residues per Gene"
)

plt.tight_layout()

plt.savefig(
    os.path.join(
        folder,
        "Figure5_Gene_Summary.png"
    ),
    dpi=600
)

plt.savefig(
    os.path.join(
        folder,
        "Figure5_Gene_Summary.pdf"
    )
)

plt.show()

plt.close()


# =========================================================
# Statistics Report
# =========================================================

report_path = os.path.join(
    folder,
    "Statistics_Report.txt"
)


with open(
    report_path,
    "w"
) as f:

    f.write(
        "DOMAIN VARIABILITY ANALYSIS REPORT\n"
    )

    f.write(
        "=================================\n\n"
    )

    f.write(
        f"Total Highly Variable Residues: "
        f"{df['Highly variable residues'].sum()}\n"
    )

    f.write(
        f"Mean Highly Variable (%): "
        f"{df['Highly Variable (%)'].mean():.2f}\n\n"
    )


    # Most variable domain
    most_variable = df.loc[
        df["Highly Variable (%)"].idxmax()
    ]


    f.write(
        "Most Variable Domain:\n"
    )

    f.write(
        f"Gene: {most_variable['Gene']}\n"
    )

    f.write(
        f"Domain: {most_variable['Domain/Motif']}\n"
    )

    f.write(
        f"Highly Variable Residues: "
        f"{most_variable['Highly variable residues']}\n"
    )

    f.write(
        f"Highly Variable Percentage: "
        f"{most_variable['Highly Variable (%)']}%\n"
    )


print(
    f"[✓] Saved statistics report: "
    f"{report_path}"
)

print(
    "\n[✓] All tasks completed successfully."
)