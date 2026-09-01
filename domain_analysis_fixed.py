#!/usr/bin/env python
# =========================================================
# Domain Conservation Analysis Pipeline
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
# Select file (portable across environments)
# -------------------------
def get_file_path():
    """
    Get the Excel file path across environments (CLI, Tkinter dialog, or console input).
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

    # 3. Fallback to current directory check or manual entry
    default_filename = "domain wise- 26.7.26.xlsx"
    if os.path.isfile(default_filename):
        print(f"Using default file found in directory: {default_filename}")
        return default_filename

    print("Falling back to manual file path entry.")
    file = input("Enter the full path to your Excel file: ").strip().strip('"')
    if not file or not os.path.isfile(file):
        raise SystemExit("No valid file provided. Exiting.")
    return file


# -------------------------
# Robust Range Parsing
# -------------------------
def parse_domain_length(range_str):
    """
    Calculates total length from range strings like '1–332',
    '333–488 & 524–594', or '141–150, 208–221'.
    """
    if pd.isna(range_str):
        return 0
    
    clean_str = str(range_str).replace('–', '-').strip()
    total_length = 0
    
    # Split across commas or ampersands
    parts = re.split(r'[,&]', clean_str)
    for part in parts:
        part = part.strip()
        # Look for start-end patterns like 100-200
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

# Load Excel sheet
df = pd.read_excel(file)

# Drop unexpected index columns like 'Unnamed: 0' if it's an index column
if "Unnamed: 0" in df.columns:
    if df["Unnamed: 0"].isnull().all() or df["Unnamed: 0"].dtype == 'object':
        # In this dataset, Unnamed: 0 holds the Gene names (VP1, VP2, etc.)
        df.rename(columns={"Unnamed: 0": "Gene"}, inplace=True)

df.columns = df.columns.str.strip()

# Forward fill Gene names (VP1, VP2, VP3...) for empty rows
if "Gene" in df.columns:
    df["Gene"] = df["Gene"].ffill()

required = [
    "Gene",
    "Residue range",
    "Domain/Motif",
    "Conserved residues",
    "Highly variable residues"
]

missing = [c for c in required if c not in df.columns]
if missing:
    print("Detected columns:", list(df.columns))
    raise ValueError(f"Missing required column(s): {missing}")

# Handle swapped columns (Residue range <-> Domain/Motif) in rows 10-36
for idx, row in df.iterrows():
    val_res = str(row["Residue range"])
    val_dom = str(row["Domain/Motif"])
    
    # If Domain/Motif contains digits with dashes/commas and Residue range contains letters
    if re.search(r'\d', val_dom) and re.search(r'[a-zA-Z]', val_res):
        df.loc[idx, "Residue range"] = val_dom
        df.loc[idx, "Domain/Motif"] = val_res

# Calculate domain lengths and percentage metrics
df["Domain Length"] = df["Residue range"].apply(parse_domain_length)

# Prevent division by zero
df["Conserved (%)"] = df.apply(
    lambda r: round((r["Conserved residues"] / r["Domain Length"] * 100), 2)
    if r["Domain Length"] > 0 else 0.0, axis=1
)

df["Highly Variable (%)"] = df.apply(
    lambda r: round((r["Highly variable residues"] / r["Domain Length"] * 100), 2)
    if r["Domain Length"] > 0 else 0.0, axis=1
)

folder = os.path.dirname(file) or "."

# -------------------------
# Save Processed Excel
# -------------------------
out_excel = os.path.join(folder, "Domain_Analysis_Output.xlsx")
df.to_excel(out_excel, index=False)
print(f"[✓] Saved updated Excel: {out_excel}")

# Label generator for figures
labels = df["Gene"].astype(str) + " - " + df["Domain/Motif"].astype(str)
x = range(len(df))
w = 0.35

# -------------------------
# Figure 1: Grouped Counts
# -------------------------
plt.figure(figsize=(16, 7))
plt.bar([i - w / 2 for i in x], df["Conserved residues"], width=w, label="Conserved", color="#2b5c8f")
plt.bar([i + w / 2 for i in x], df["Highly variable residues"], width=w, label="Highly Variable", color="#d95f02")
plt.xticks(x, labels, rotation=90, fontsize=8)
plt.ylabel("Residue Count")
plt.title("Conserved vs Highly Variable Residue Counts per Domain")
plt.tight_layout()
plt.legend()
plt.savefig(os.path.join(folder, "Figure1_Grouped_Counts.png"), dpi=600)
plt.savefig(os.path.join(folder, "Figure1_Grouped_Counts.pdf"))
plt.show()
plt.close()

# -------------------------
# Figure 2: Grouped Percentages
# -------------------------
plt.figure(figsize=(16, 7))
plt.bar([i - w / 2 for i in x], df["Conserved (%)"], width=w, label="Conserved %", color="#2b5c8f")
plt.bar([i + w / 2 for i in x], df["Highly Variable (%)"], width=w, label="Highly Variable %", color="#d95f02")
plt.xticks(x, labels, rotation=90, fontsize=8)
plt.ylabel("Percentage (%)")
plt.title("Conserved vs Highly Variable Residue Percentage per Domain")
plt.tight_layout()
plt.legend()
plt.savefig(os.path.join(folder, "Figure2_Grouped_Percentage.png"), dpi=600)
plt.savefig(os.path.join(folder, "Figure2_Grouped_Percentage.pdf"))
plt.show()
plt.close()

# -------------------------
# Figure 3: Stacked Bar Chart
# -------------------------
plt.figure(figsize=(16, 7))
plt.bar(labels, df["Conserved residues"], label="Conserved", color="#2b5c8f")
plt.bar(labels, df["Highly variable residues"], bottom=df["Conserved residues"], label="Highly Variable", color="#d95f02")
plt.xticks(rotation=90, fontsize=8)
plt.ylabel("Residue Count")
plt.title("Stacked Composition of Residues per Domain")
plt.tight_layout()
plt.legend()
plt.savefig(os.path.join(folder, "Figure3_Stacked.png"), dpi=600)
plt.show()
plt.close()

# -------------------------
# Figure 4 & 5: Heatmaps
# -------------------------
counts = df.pivot_table(index="Gene", columns="Domain/Motif", values="Conserved residues", aggfunc="sum")

plt.figure(figsize=(12, 6))
if HAS_SNS:
    sns.heatmap(counts, annot=True, cmap="viridis", fmt=".0f", cbar_kws={'label': 'Conserved Residues'})
else:
    plt.imshow(counts.fillna(0), aspect="auto", cmap="viridis")
    plt.colorbar(label='Conserved Residues')
    plt.xticks(range(len(counts.columns)), counts.columns, rotation=90)
    plt.yticks(range(len(counts.index)), counts.index)
plt.title("Conserved Residues Heatmap")
plt.tight_layout()
plt.savefig(os.path.join(folder, "Figure4_Heatmap_Conserved.png"), dpi=600)
plt.show()
plt.close()

var = df.pivot_table(index="Gene", columns="Domain/Motif", values="Highly variable residues", aggfunc="sum")

plt.figure(figsize=(12, 6))
if HAS_SNS:
    sns.heatmap(var, annot=True, cmap="magma", fmt=".0f", cbar_kws={'label': 'Variable Residues'})
else:
    plt.imshow(var.fillna(0), aspect="auto", cmap="magma")
    plt.colorbar(label='Variable Residues')
    plt.xticks(range(len(var.columns)), var.columns, rotation=90)
    plt.yticks(range(len(var.index)), var.index)
plt.title("Highly Variable Residues Heatmap")
plt.tight_layout()
plt.savefig(os.path.join(folder, "Figure5_Heatmap_Variable.png"), dpi=600)
plt.show()
plt.close()

# -------------------------
# Figure 6: Gene Summary
# -------------------------
summary = df.groupby("Gene")[["Conserved residues", "Highly variable residues"]].sum()
summary.to_excel(os.path.join(folder, "Summary_Table.xlsx"))
print(f"[✓] Saved summary table: {os.path.join(folder, 'Summary_Table.xlsx')}")

summary.plot(kind="bar", figsize=(10, 5), color=["#2b5c8f", "#d95f02"])
plt.ylabel("Total Residues")
plt.title("Total Residues per Gene")
plt.tight_layout()
plt.savefig(os.path.join(folder, "Figure6_Gene_Summary.png"), dpi=600)
plt.show()
plt.close()

# -------------------------
# Statistics Report
# -------------------------
report_path = os.path.join(folder, "Statistics_Report.txt")
with open(report_path, "w") as f:
    f.write("DOMAIN ANALYSIS REPORT\n")
    f.write("======================\n\n")
    f.write(f"Total Conserved: {df['Conserved residues'].sum()}\n")
    f.write(f"Total Highly Variable: {df['Highly variable residues'].sum()}\n")
    f.write(f"Mean Conserved (%): {df['Conserved (%)'].mean():.2f}\n")
    f.write(f"Mean Highly Variable (%): {df['Highly Variable (%)'].mean():.2f}\n\n")
    
    most_conserved = df.loc[df["Conserved (%)"].idxmax()]
    f.write("Most Conserved Domain:\n")
    f.write(f"Gene: {most_conserved['Gene']}, Domain: {most_conserved['Domain/Motif']} ({most_conserved['Conserved (%)']}%\n\n")
    
    most_variable = df.loc[df["Highly Variable (%)"].idxmax()]
    f.write("Most Variable Domain:\n")
    f.write(f"Gene: {most_variable['Gene']}, Domain: {most_variable['Domain/Motif']} ({most_variable['Highly Variable (%)']}%\n")

print(f"[✓] Saved statistics report: {report_path}")
print("\n[✓] All tasks completed successfully.")