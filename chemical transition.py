#!/usr/bin/env python3
"""
GENERATE SIMPLIFIED AMINO ACID CHEMICAL PROPERTY TRANSITION MATRIX
- Merges Proline (P) into the Hydrophobic group for a 4x4 matrix view.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from tkinter import Tk, filedialog

# ----------------------------
# Constants & Mapping
# ----------------------------

# Simplified biochemical classification (Consolidating Proline into Hydrophobic)
amino_acid_groups = {
    'Hydrophobic': ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'G', 'P'], # Proline merged
    'Polar': ['S', 'T', 'N', 'Q', 'C'],
    'Basic': ['K', 'R', 'H'],
    'Acidic': ['D', 'E']
}
# Reverse lookup dictionary
aa_to_group = {aa: group for group, aas in amino_acid_groups.items() for aa in aas}

# ----------------------------
# Main function
# ----------------------------

def generate_chemical_matrix(input_file):
    # 1. Load data
    df = pd.read_excel(input_file)

    # 2. Extract 'From' and 'To' amino acids
    df[['From_AA', 'To_AA']] = df['Change'].str.split(' -> ', expand=True)

    # 3. Map amino acids to their chemical groups
    df['From_Group'] = df['From_AA'].map(aa_to_group)
    df['To_Group'] = df['To_AA'].map(aa_to_group)

    # 4. Sum the 'Count' for each chemical transition category
    # Ensure all groups are present even if count is zero for a clean 4x4 look
    all_groups = ['Hydrophobic', 'Polar', 'Basic', 'Acidic']
    matrix = df.groupby(['From_Group', 'To_Group'])['Count'].sum().unstack().fillna(0)
    matrix = matrix.reindex(index=all_groups, columns=all_groups, fill_value=0)

    # 5. Visualization: Create a professional 4x4 Heatmap
    plt.figure(figsize=(9, 7))
    sns.heatmap(matrix, annot=True, fmt='.0f', cmap='YlGnBu', 
                linewidths=.5, cbar_kws={'label': 'Substitution Frequency'})
    
    plt.title('Consolidated Amino Acid Chemical Property Transition Matrix', fontsize=14)
    plt.xlabel('To (Substituted Group)', fontsize=12)
    plt.ylabel('From (Original Group)', fontsize=12)
    plt.tight_layout()
    
    # Save outputs
    plt.savefig('simplified_chemical_matrix.png', dpi=300)
    matrix.to_excel('simplified_chemical_transition_data.xlsx')
    
    print(f"Analysis Successful ✅")
    print(f"4x4 Matrix saved to 'simplified_chemical_transition_data.xlsx'")
    print(f"Heatmap saved to 'simplified_chemical_matrix.png'")

# ----------------------------
# GUI Execution
# ----------------------------

root = Tk()
root.withdraw()
root.attributes("-topmost", True)

input_file = filedialog.askopenfilename(
    title="Select Mutation Input File",
    filetypes=[("Excel files", "*.xlsx *.xls")]
)

if input_file:
    generate_chemical_matrix(input_file)
    print("PROCESS COMPLETE ✅")
else:
    print("No file selected. Exiting.")