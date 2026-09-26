# ============================================================
# AMINO ACID VARIABILITY ANALYSIS
# Based on the variability formula shown in the paper
#
# Variability =
# Number of different amino acids at a position
# ------------------------------------------------------------
# Frequency of the most common amino acid at that position
#
# Designed for the Excel structure:
# Sequence row
# Classification row
# Sequence row
# Classification row
# ...
# ============================================================

import pandas as pd
import numpy as np
from collections import Counter
import tkinter as tk
from tkinter import filedialog
import os


# ============================================================
# 1. SELECT INPUT EXCEL FILE
# ============================================================

root = tk.Tk()
root.withdraw()

input_file = filedialog.askopenfilename(
    title="Select your Excel dataset",
    filetypes=[
        ("Excel files", "*.xlsx *.xls"),
        ("All files", "*.*")
    ]
)

if not input_file:
    print("No file selected. Program stopped.")
    raise SystemExit

print("\nInput file selected:")
print(input_file)


# ============================================================
# 2. READ EXCEL FILE
# ============================================================

df = pd.read_excel(input_file)

print("\nDataset loaded.")
print("Rows:", df.shape[0])
print("Columns:", df.shape[1])


# ============================================================
# 3. IDENTIFY CODON COLUMNS
# ============================================================

id_column = df.columns[0]

codon_columns = [
    col for col in df.columns
    if str(col).startswith("Codon_")
]

print("\nNumber of codon positions:", len(codon_columns))


# ============================================================
# 4. GENETIC CODE
# ============================================================

genetic_code = {
    # Phenylalanine
    "TTT": "F", "TTC": "F",

    # Leucine
    "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",

    # Isoleucine
    "ATT": "I", "ATC": "I", "ATA": "I",

    # Methionine
    "ATG": "M",

    # Valine
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",

    # Serine
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "AGT": "S", "AGC": "S",

    # Proline
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",

    # Threonine
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",

    # Alanine
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",

    # Tyrosine
    "TAT": "Y", "TAC": "Y",

    # Histidine
    "CAT": "H", "CAC": "H",

    # Glutamine
    "CAA": "Q", "CAG": "Q",

    # Asparagine
    "AAT": "N", "AAC": "N",

    # Lysine
    "AAA": "K", "AAG": "K",

    # Aspartic acid
    "GAT": "D", "GAC": "D",

    # Glutamic acid
    "GAA": "E", "GAG": "E",

    # Cysteine
    "TGT": "C", "TGC": "C",

    # Tryptophan
    "TGG": "W",

    # Arginine
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGA": "R", "AGG": "R",

    # Glycine
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",

    # Stop codons
    "TAA": "*", "TAG": "*", "TGA": "*"
}


def translate_codon(codon):
    """
    Convert a codon into a one-letter amino-acid code.
    """
    if pd.isna(codon):
        return ""

    codon = str(codon).strip().upper()

    if codon in genetic_code:
        return genetic_code[codon]

    # Handle possible gap / missing / ambiguous codons
    if "-" in codon or len(codon) != 3:
        return "X"

    return genetic_code.get(codon, "X")


# ============================================================
# 5. FIND CONSENSUS ROW
# ============================================================

consensus_index = 0

consensus_codons = df.loc[consensus_index, codon_columns]

# Translate consensus codons into amino acids
consensus_aa = {}

for col in codon_columns:
    consensus_aa[col] = translate_codon(
        consensus_codons[col]
    )


# ============================================================
# 6. IDENTIFY ACTUAL SEQUENCE ROWS
# ============================================================

first_column = df[id_column].astype(str)

sequence_rows = []

for i in range(len(df)):

    value = str(df.iloc[i, 0])

    # Consensus is NOT an actual sequence
    if value == "Consensus":
        continue

    # Classification rows are NOT sequences
    if value.startswith("Classification"):
        continue

    # Everything else is an actual sequence
    sequence_rows.append(i)


print("\nNumber of actual sequences:", len(sequence_rows))


# ============================================================
# 7. CREATE SHEET 1
# ============================================================
#
# For each sequence and codon position:
#
# No Change:
#       consensus amino acid
#
# Synonymous:
#       amino acid translated from observed codon
#
# Non-synonymous:
#       amino acid translated from observed codon
#
# ============================================================

aa_data = []

for row_index in sequence_rows:

    sample_id = df.iloc[row_index, 0]

    row_result = {
        "Sample": sample_id
    }

    # Classification row is immediately after sequence row
    classification_row_index = row_index + 1

    if (
        classification_row_index < len(df)
        and str(df.iloc[classification_row_index, 0])
        .startswith("Classification")
    ):
        classification = df.iloc[classification_row_index]
    else:
        classification = None

    for col in codon_columns:

        # Observed codon
        observed_codon = df.loc[row_index, col]

        # Classification at this position
        if classification is not None:
            status = str(classification[col]).strip()
        else:
            status = ""

        # ----------------------------------------------------
        # NO CHANGE
        # Use consensus amino acid
        # ----------------------------------------------------

        if status == "No Change":

            aa = consensus_aa[col]

        # ----------------------------------------------------
        # SYNONYMOUS OR OTHER CHANGE
        # Translate observed codon
        # ----------------------------------------------------

        else:

            aa = translate_codon(observed_codon)

        row_result[col] = aa

    aa_data.append(row_result)


sheet1 = pd.DataFrame(aa_data)


# ============================================================
# 8. ADD CONSENSUS AS FIRST ROW
# ============================================================

consensus_output = {"Sample": "Consensus"}

for col in codon_columns:
    consensus_output[col] = consensus_aa[col]

consensus_df = pd.DataFrame([consensus_output])

sheet1 = pd.concat(
    [consensus_df, sheet1],
    ignore_index=True
)


print("\nSheet 1 created.")
print("Shape:", sheet1.shape)


# ============================================================
# 9. CREATE SHEET 2
# ============================================================
#
# For each amino-acid position:
#
# Most frequent amino acid
# Number of sequences having it
# Total number of sequences
# Frequency = count / total sequences
#
# ============================================================

# IMPORTANT:
# Do NOT include the consensus row in calculations
actual_aa = sheet1.iloc[1:].copy()

total_sequences = len(actual_aa)

frequency_results = []

for col in codon_columns:

    # Get amino acids at this position
    values = actual_aa[col].astype(str).str.strip()

    # Remove blank values and X if desired
    values = values[
        (values != "") &
        (values != "nan")
    ]

    counts = Counter(values)

    if len(counts) == 0:

        most_common_aa = ""
        most_common_count = 0
        frequency = np.nan

    else:

        # Most frequent amino acid
        most_common_aa, most_common_count = counts.most_common(1)[0]

        # Frequency
        frequency = most_common_count / total_sequences

    frequency_results.append({

        "Position": col.replace("Codon_", ""),

        "Codon_Position": col,

        "Consensus_AA": consensus_aa[col],

        "Most_Frequent_AA": most_common_aa,

        "Most_Frequent_Count": most_common_count,

        "Total_Sequences": total_sequences,

        "Most_Frequent_Frequency": frequency
    })


sheet2 = pd.DataFrame(frequency_results)


# ============================================================
# 10. CREATE SHEET 3
# ============================================================
#
# Variability =
#
# Number of different amino acids
# --------------------------------
# Frequency of most common amino acid
#
# ============================================================

variability_results = []

for col in codon_columns:

    values = actual_aa[col].astype(str).str.strip()

    values = values[
        (values != "") &
        (values != "nan")
    ]

    # Count different amino acids
    different_aas = sorted(values.unique())

    number_of_different_aas = len(different_aas)

    # Counts
    counts = Counter(values)

    if len(counts) == 0:

        most_common_aa = ""
        most_common_count = 0
        frequency = np.nan
        variability = np.nan

    else:

        most_common_aa, most_common_count = counts.most_common(1)[0]

        frequency = most_common_count / total_sequences

        # ----------------------------------------------------
        # PAPER'S VARIABILITY FORMULA
        # ----------------------------------------------------

        if frequency > 0:

            variability = (
                number_of_different_aas /
                frequency
            )

        else:

            variability = np.nan

    variability_results.append({

        "Position": col.replace("Codon_", ""),

        "Codon_Position": col,

        "Consensus_AA": consensus_aa[col],

        "Different_Amino_Acid_Count":
            number_of_different_aas,

        "Different_Amino_Acids":
            ", ".join(different_aas),

        "Most_Frequent_AA":
            most_common_aa,

        "Most_Frequent_Count":
            most_common_count,

        "Total_Sequences":
            total_sequences,

        "Most_Frequent_Frequency":
            frequency,

        "Variability":
            variability
    })


sheet3 = pd.DataFrame(variability_results)


# ============================================================
# 11. SAVE EVERYTHING INTO ONE EXCEL FILE
# ============================================================

output_directory = os.path.dirname(input_file)

output_file = os.path.join(
    output_directory,
    "Amino_Acid_Variability_Analysis.xlsx"
)


with pd.ExcelWriter(
    output_file,
    engine="openpyxl"
) as writer:

    sheet1.to_excel(
        writer,
        sheet_name="1_Amino_Acid_Data",
        index=False
    )

    sheet2.to_excel(
        writer,
        sheet_name="2_Most_Frequent_Frequency",
        index=False
    )

    sheet3.to_excel(
        writer,
        sheet_name="3_Variability",
        index=False
    )


# ============================================================
# 12. FINISHED
# ============================================================

print("\n" + "=" * 60)
print("ANALYSIS COMPLETED")
print("=" * 60)

print("\nTotal sequences:", total_sequences)
print("Total codon positions:", len(codon_columns))

print("\nOutput file:")
print(output_file)

print("\nSheets created:")
print("1. 1_Amino_Acid_Data")
print("2. 2_Most_Frequent_Frequency")
print("3. 3_Variability")

print("\nVariability formula used:")
print(
    "Variability = Number of different amino acids / "
    "Frequency of most common amino acid"
)

print("\nDone.")