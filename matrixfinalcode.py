#!/usr/bin/env python3
"""
Enumerate all nucleotide-level possibilities for amino-acid substitutions
Input: Excel with Codon_Position | Change | Count
Output: Excel with all possible codon, nucleotide, and Ti/Tv paths (Count preserved)
"""

# ----------------------------
# Imports
# ----------------------------

import pandas as pd
from itertools import product
from tkinter import Tk, filedialog

# ----------------------------
# Genetic code
# ----------------------------

GENETIC_CODE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

AA_TO_CODONS = {}
for codon, aa in GENETIC_CODE.items():
    AA_TO_CODONS.setdefault(aa, []).append(codon)

TRANSITIONS = {('A','G'), ('G','A'), ('C','T'), ('T','C')}

# ----------------------------
# Helper functions
# ----------------------------

def classify_ti_tv(ref, alt):
    changes = []
    for r, a in zip(ref, alt):
        if r != a:
            if (r, a) in TRANSITIONS:
                changes.append("Ti")
            else:
                changes.append("Tv")
    if not changes:
        return "None"
    if all(c == "Ti" for c in changes):
        return "Transition"
    if all(c == "Tv" for c in changes):
        return "Transversion"
    return "Mixed"

def nucleotide_changes(ref, alt):
    return ", ".join(
        f"{r}->{a}(pos{idx+1})"
        for idx, (r, a) in enumerate(zip(ref, alt))
        if r != a
    )

# ----------------------------
# Main processing function
# ----------------------------

def generate_all_possibilities(input_excel, output_excel):

    df = pd.read_excel(input_excel)

    rows = []

    for _, row in df.iterrows():
        position = row.iloc[0]
        aa_change = row.iloc[1]
        count = row.iloc[2]

        ref_aa, alt_aa = aa_change.split(" -> ")

        ref_codons = AA_TO_CODONS.get(ref_aa, [])
        alt_codons = AA_TO_CODONS.get(alt_aa, [])

        for ref_codon, alt_codon in product(ref_codons, alt_codons):
            nt_diff = sum(r != a for r, a in zip(ref_codon, alt_codon))
            if nt_diff == 0:
                continue

            rows.append({
                "SNP": f"{position}_{ref_codon}>{alt_codon}",
                "Amino acid change": aa_change,
                "Codons": f"{ref_codon} → {alt_codon}",
                "Nucleotide": nucleotide_changes(ref_codon, alt_codon),
                "Nti or Ntv": classify_ti_tv(ref_codon, alt_codon),
                "Position": position,
                "Count": count
            })

    out_df = pd.DataFrame(rows)
    out_df.to_excel(output_excel, index=False)

# ----------------------------
# GUI Input
# ----------------------------

root = Tk()
root.attributes("-topmost", True)
root.withdraw()

input_excel = filedialog.askopenfilename(
    title="Select input Excel file (AA substitutions)",
    filetypes=[("Excel files", "*.xlsx *.xls")]
)

if not input_excel:
    raise SystemExit("No input file selected")

output_excel = filedialog.asksaveasfilename(
    title="Save output Excel file",
    defaultextension=".xlsx",
    filetypes=[("Excel files", "*.xlsx")]
)

if not output_excel:
    raise SystemExit("No output file selected")

# ----------------------------
# Run
# ----------------------------

generate_all_possibilities(input_excel, output_excel)
print("All nucleotide-level possibilities generated successfully (Count preserved)")
