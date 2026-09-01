#!/usr/bin/env python3
"""
FINAL VARIABILITY ANALYSIS (Hotspots + Conserved)

INPUT:
- Excel:
    Columns: Position, Amino acid change, Count

OUTPUT:
- Excel:
    1_Merged_Data (all positions)
    2_Thresholds (V5, V95, V99)
    3_Hotspots_95
    4_Hotspots_99
    5_Conserved

- Plots:
    variability_curve.png
    variability_classification.png
"""

# ----------------------------
# Imports
# ----------------------------

import pandas as pd
from tkinter import Tk, filedialog, messagebox
import matplotlib.pyplot as plt
import os

# ----------------------------
# Variability calculation
# ----------------------------

def calculate_variability(df):

    df2 = df[['Position','Amino acid change','Count']].dropna()
    df2[['From','To']] = df2['Amino acid change'].str.split(' -> ', expand=True)

    result = []

    for pos, group in df2.groupby('Position'):
        unique_aa = group['To'].nunique()
        max_freq = group['Count'].max()   # ✅ CORRECT METHOD
        variability = unique_aa / max_freq if max_freq != 0 else 0

        result.append((pos, unique_aa, max_freq, variability))

    res_df = pd.DataFrame(result, columns=[
        'Position','Unique_AA','Max_Count','Variability'
    ])

    return res_df.sort_values(by='Position')

# ----------------------------
# Classification
# ----------------------------

def classify_positions(res_df):

    V95 = res_df['Variability'].quantile(0.95)
    V99 = res_df['Variability'].quantile(0.99)
    V5  = res_df['Variability'].quantile(0.05)

    def classify(v):
        if v >= V99:
            return "Extreme"
        elif v >= V95:
            return "Hotspot"
        elif v <= V5:
            return "Conserved"
        else:
            return "Moderate"

    res_df['Class'] = res_df['Variability'].apply(classify)

    return res_df, V5, V95, V99

# ----------------------------
# Main function
# ----------------------------

def run_analysis(input_file, output_file):

    # Load
    df = pd.read_excel(input_file)

    # Calculate variability
    res_df = calculate_variability(df)

    # Classify
    res_df, V5, V95, V99 = classify_positions(res_df)

    # Extract tables
    hotspots_95 = res_df[res_df['Variability'] >= V95]
    hotspots_99 = res_df[res_df['Variability'] >= V99]
    conserved   = res_df[res_df['Variability'] <= V5]

    thresholds_df = pd.DataFrame({
        "Metric": ["V5 (Conserved)", "V95 (Hotspot)", "V99 (Extreme)"],
        "Value": [V5, V95, V99]
    })

    # ----------------------------
    # Save Excel
    # ----------------------------

    with pd.ExcelWriter(output_file) as writer:
        res_df.to_excel(writer, sheet_name="1_All_Data", index=False)
        thresholds_df.to_excel(writer, sheet_name="2_Thresholds", index=False)
        hotspots_95.to_excel(writer, sheet_name="3_Hotspots_95", index=False)
        hotspots_99.to_excel(writer, sheet_name="4_Hotspots_99", index=False)
        conserved.to_excel(writer, sheet_name="5_Conserved", index=False)

    # ----------------------------
    # Plot option
    # ----------------------------

    root = Tk()
    root.withdraw()

    if messagebox.askyesno("Plots", "Generate plots?"):

        out_dir = os.path.dirname(output_file)

        # ----------------------------
        # Plot 1: Variability curve
        # ----------------------------

        plt.figure()
        plt.plot(res_df['Position'], res_df['Variability'])
        plt.axhline(V95, linestyle='--')
        plt.axhline(V99, linestyle='--')
        plt.axhline(V5, linestyle='--')

        plt.xlabel("Position")
        plt.ylabel("Variability")
        plt.title("Variability Profile")

        plt.savefig(os.path.join(out_dir, "variability_curve.png"), dpi=300)

        # ----------------------------
        # Plot 2: Classification scatter
        # ----------------------------

        plt.figure()

        for cls in ["Moderate","Hotspot","Extreme","Conserved"]:
            subset = res_df[res_df['Class'] == cls]
            plt.scatter(subset['Position'], subset['Variability'], label=cls)

        plt.xlabel("Position")
        plt.ylabel("Variability")
        plt.title("Variability Classification")
        plt.legend()

        plt.savefig(os.path.join(out_dir, "variability_classification.png"), dpi=300)

# ----------------------------
# GUI
# ----------------------------

root = Tk()
root.attributes("-topmost", True)
root.withdraw()

input_file = filedialog.askopenfilename(
    title="Select Excel file",
    filetypes=[("Excel files", "*.xlsx *.xls")]
)

output_file = filedialog.asksaveasfilename(
    title="Save output Excel",
    defaultextension=".xlsx",
    filetypes=[("Excel files", "*.xlsx")]
)

run_analysis(input_file, output_file)

print("VARIABILITY ANALYSIS COMPLETE ✅")