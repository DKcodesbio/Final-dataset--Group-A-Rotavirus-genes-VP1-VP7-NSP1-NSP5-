import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import os
from tkinter import Tk, filedialog

# Integrated safety for sklearn (PCA)
try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
except ImportError:
    print("Optimization Hint: Run 'pip install scikit-learn' to enable PCA plots.")

# ---------------------------------------------------------
# 1. SETUP, MAPPING & FLEXIBLE COLUMNS
# ---------------------------------------------------------
categories = ['Hydrophobic', 'Polar', 'Basic', 'Acidic']
colors = sns.color_palette("husl", 4)
cat_colors = dict(zip(categories, colors))

def find_col(df, options):
    """Finds a column even if the name has different casing or extra spaces."""
    cols_clean = [str(c).strip().lower() for c in df.columns]
    for opt in options:
        if opt.lower() in cols_clean:
            return df.columns[cols_clean.index(opt.lower())]
    return None

# ---------------------------------------------------------
# 2. MAIN ANALYSIS ENGINE
# ---------------------------------------------------------
def generate_final_atlas():
    root = Tk(); root.withdraw(); root.attributes("-topmost", True)
    
    print("Step 1: Select the 'Merge Data' Excel (Stacked Matrices)")
    matrix_file = filedialog.askopenfilename(title="Select Merged Matrix Excel")
    
    print("Step 2: Select the FOLDER containing individual raw gene files")
    raw_folder = filedialog.askdirectory(title="Select Folder with individual Gene Excels")
    
    if not matrix_file or not raw_folder:
        print("Selection cancelled. Exiting."); return

    # --- PART A: PARSE MATRICES ---
    try:
        df_m = pd.read_excel(matrix_file)
        genes_data = {}
        for i in range(len(df_m)):
            gene_cell = df_m.iloc[i, 0]
            if pd.notna(gene_cell) and "From_Group" in str(df_m.iloc[i, 1]):
                gene_name = str(gene_cell).strip()
                data_block = df_m.iloc[i+1:i+5, 2:6].values.astype(float)
                genes_data[gene_name] = pd.DataFrame(data_block, index=categories, columns=categories)
    except Exception as e:
        print(f"Error parsing Matrix file: {e}"); return

    # --- PART B: PARSE RAW DATA (Manhattan & Load) ---
    raw_files = [f for f in os.listdir(raw_folder) if f.endswith(('.xlsx', '.xls'))]
    master_raw_list = []
    for f_name in raw_files:
        try:
            df_temp = pd.read_excel(os.path.join(raw_folder, f_name))
            gene_id = os.path.splitext(f_name)[0].split('_')[0]
            p_col = find_col(df_temp, ['Position', 'Codon_Position', 'Codon Position', 'Pos'])
            c_col = find_col(df_temp, ['Count', 'Mutation Count', 'Frequency'])
            if p_col and c_col:
                subset = df_temp[[p_col, c_col]].copy()
                subset.columns = ['Position', 'Count']
                subset['Gene'] = gene_id
                master_raw_list.append(subset)
        except: continue

    if not master_raw_list:
        print("Error: No raw data found. Check columns 'Position' and 'Count'."); return
    df_full_raw = pd.concat(master_raw_list, ignore_index=True)

    # ---------------------------------------------------------
    # 3. VISUALIZATION SUITE (THE 6 PLOTS)
    # ---------------------------------------------------------
    sns.set_theme(style="whitegrid")

    # 1. FACET ATLAS
    rows = math.ceil(len(genes_data) / 3)
    fig1, axes1 = plt.subplots(rows, 3, figsize=(15, rows*4))
    axes1 = axes1.flatten() if len(genes_data) > 1 else [axes1]
    for i, (name, m) in enumerate(genes_data.items()):
        sns.heatmap(m, annot=True, fmt='.0f', cmap='YlGnBu', ax=axes1[i], cbar=False)
        axes1[i].set_title(f"Segment: {name}", fontweight='bold')
    for j in range(i+1, len(axes1)): fig1.delaxes(axes1[j])
    plt.tight_layout(); plt.savefig('1_Facet_Atlas.png', dpi=300)

    # 2. GLOBAL GENOMIC MATRIX
    global_m = pd.DataFrame(0.0, index=categories, columns=categories)
    for m in genes_data.values(): global_m += m
    plt.figure(figsize=(9, 7))
    sns.heatmap(global_m, annot=True, fmt='.0f', cmap='YlGnBu', linewidths=1.5)
    plt.title('Global Genomic Chemical Transition Matrix', fontweight='bold', pad=20)
    plt.savefig('2_Global_Matrix.png', dpi=300)

    # 3. STACKED PERCENTAGE BAR PROFILE
    bar_data = []
    for gene, matrix in genes_data.items():
        total = matrix.values.sum()
        to_sums = matrix.sum(axis=0)
        for cat in categories:
            bar_data.append({'Gene': gene, 'Chemical Group': cat, 'Percentage': (to_sums[cat]/total)*100})
    plt.figure(figsize=(12, 6))
    sns.barplot(data=pd.DataFrame(bar_data), x='Gene', y='Percentage', hue='Chemical Group', palette='viridis')
    plt.title('Segmental Selection Profile (Target Group %)', fontweight='bold')
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left'); plt.tight_layout()
    plt.savefig('3_Stacked_Selection_Profile.png', dpi=300)

    # 4. PCA EVOLUTIONARY CLUSTERING
    try:
        pca_in = [genes_data[g].values.flatten() for g in genes_data.keys()]
        pca_res = PCA(n_components=2).fit_transform(StandardScaler().fit_transform(pca_in))
        plt.figure(figsize=(10, 8))
        sns.scatterplot(x=pca_res[:,0], y=pca_res[:,1], hue=list(genes_data.keys()), s=250, palette='tab20')
        for k, name in enumerate(genes_data.keys()):
            plt.annotate(name, (pca_res[k,0], pca_res[k,1]), xytext=(7,7), textcoords='offset points', fontweight='bold')
        plt.title('PCA: Segmental Mutational Similarity', fontweight='bold')
        plt.savefig('4_PCA_Clustering.png', dpi=300)
    except: pass

    # 5. DIRECTIONAL CHORD FLOW
    fig5, ax5 = plt.subplots(figsize=(10, 10)); ax5.set_xlim(-1.5, 1.5); ax5.set_ylim(-1.5, 1.5); ax5.axis('off')
    angles = np.linspace(0, 2*np.pi, 4, endpoint=False)
    pos = {cat: (np.cos(angles[idx]), np.sin(angles[idx])) for idx, cat in enumerate(categories)}
    for idx, cat in enumerate(categories):
        ax5.add_patch(plt.Circle(pos[cat], 0.15, color=colors[idx], alpha=0.7))
        ax5.text(pos[cat][0]*1.3, pos[cat][1]*1.3, cat, fontweight='bold', ha='center')
    max_v = global_m.values.max()
    for f in categories:
        for t in categories:
            if f != t and global_m.loc[f, t] > (max_v * 0.05):
                ax5.annotate("", xy=pos[t], xytext=pos[f], arrowprops=dict(arrowstyle="->", 
                             color=cat_colors[f], lw=global_m.loc[f, t]/max_v*15, alpha=0.4, connectionstyle="arc3,rad=.2"))
    plt.title("Genome-Wide Directional Chemical Flow", fontweight='bold'); plt.savefig('5_Chord_Flow.png', dpi=300)

    # 6. GENOME-WIDE MANHATTAN PLOT
    plt.figure(figsize=(15, 6))
    sns.scatterplot(data=df_full_raw, x='Position', y='Count', hue='Gene', palette='tab10', s=30, alpha=0.6, edgecolor='none')
    plt.title('Genome-Wide Manhattan Plot: Mutational Hotspots', fontsize=16, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left'); plt.grid(axis='y', alpha=0.3)
    plt.savefig('6_Manhattan_Plot.png', dpi=300, bbox_inches='tight')

    print("\nSUCCESS: All 6 Professional Plots Saved! ✅")

if __name__ == "__main__":
    generate_final_atlas()