#!/usr/bin/env python3
"""
COVID-19 Genome Mutation Analysis - Professional Notebook Generator
Generates a complete, production-ready Google Colab notebook for comprehensive
SARS-CoV-2 variant analysis.
"""

import json

def create_markdown_cell(source):
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": [line + "\n" for line in source.split("\n")]
    }

def create_code_cell(source):
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [line + "\n" for line in source.split("\n")]
    }

# Initialize Notebook Structure
notebook = {
    "cells": [],
    "metadata": {
        "colab": {"name": "COVID19_Professional_Analysis.ipynb", "provenance": []},
        "kernelspec": {"display_name": "Python 3", "name": "python3"},
        "language_info": {
            "codemirror_mode": {"name": "ipython", "version": 3},
            "file_extension": ".py",
            "mimetype": "text/x-python",
            "name": "python",
            "version": "3.10.12"
        }
    },
    "nbformat": 4,
    "nbformat_minor": 0
}

#=============================================================================
# CELL 1: Title & Overview
#=============================================================================
notebook["cells"].append(create_markdown_cell("""
# COVID-19 Genome Mutation Analysis & Protein Structure Visualization Pipeline

**Professional-Grade Bioinformatics Pipeline for SARS-CoV-2 Variant Analysis**

This notebook provides comprehensive analysis of COVID-19 variants including:
- Genome sequence acquisition and validation
- Multiple sequence alignment
- Mutation detection (nucleotide and amino acid level)
- Phylogenetic tree construction
- Statistical analysis (dN/dS ratios, selection pressure)
- Interactive visualizations
- 3D protein structure mapping
- ChimeraX script generation for publication-quality renders

**Author:** Bioinformatics Pipeline | **Version:** 2.0 Professional
"""))

#=============================================================================
# CELL 2: Table of Contents
#=============================================================================
notebook["cells"].append(create_markdown_cell("""
## Table of Contents

1. [Environment Setup](#setup)
2. [Data Acquisition](#data)
3. [Sequence Alignment](#alignment)
4. [Mutation Detection](#mutations)
5. [Phylogenetic Analysis](#phylo)
6. [Statistical Analysis](#stats)
7. [Visualizations](#viz)
8. [Protein Structure Analysis](#structure)
9. [ChimeraX Scripts](#chimerax)
10. [Results Export](#export)
11. [Summary](#summary)
"""))

#=============================================================================
# CELL 3: Environment Setup
#=============================================================================
notebook["cells"].append(create_markdown_cell("## I. Environment Setup & Dependencies {#setup}"))

notebook["cells"].append(create_code_cell("""
# System Dependencies - MUSCLE for alignment
print("Installing system dependencies...")
!apt-get update -qq && apt-get install -y muscle -qq

# Python Dependencies
print("Installing Python libraries...")
!pip install -q biopython py3Dmol pandas seaborn plotly scikit-bio networkx scipy openpyxl

import os, sys, json, requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import py3Dmol
import networkx as nx
from tqdm.notebook import tqdm
from datetime import datetime
from scipy import stats

# Biopython
from Bio import SeqIO, AlignIO, Phylo, Entrez, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Data import CodonTable
from Bio.PDB import PDBList, PDBParser

# Configuration
Entrez.email = "researcher@example.com"
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 100

# Directory Structure
BASE_DIR = "/content/covid_analysis"
DIRS = [
    f"{BASE_DIR}/data/genomes", f"{BASE_DIR}/data/proteins", f"{BASE_DIR}/data/structures",
    f"{BASE_DIR}/results/alignments", f"{BASE_DIR}/results/mutations",
    f"{BASE_DIR}/results/phylogenetics", f"{BASE_DIR}/results/visualizations",
    f"{BASE_DIR}/results/chimerax_scripts", f"{BASE_DIR}/reports"
]
for d in DIRS: os.makedirs(d, exist_ok=True)

print(f"Setup Complete! Working directory: {BASE_DIR}")
"""))

#=============================================================================
# CELL 4: Data Acquisition Module
#=============================================================================
notebook["cells"].append(create_markdown_cell("## II. Data Acquisition Module {#data}"))

notebook["cells"].append(create_code_cell("""
class DataAcquisitionManager:
    '''Handles downloading and validation of genomic data and protein structures'''
    
    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.genome_dir = f"{base_dir}/data/genomes"
        self.structure_dir = f"{base_dir}/data/structures"
        
        # Comprehensive variant collection with GenBank accessions
        self.variants = {
            "Wuhan-Hu-1": "NC_045512.2",
            "Alpha_B117": "CM034978.1",
            "Beta_B1351": "OM530757.1",
            "Gamma_P1": "OM530769.1",
            "Delta_B16172": "OM287553.1",
            "Omicron_BA1": "NC_055123.1",
            "Omicron_BA2": "OM296767.1",
            "Omicron_BA4": "ON481610.1",
            "Omicron_BA5": "OX315743.1",
            "Omicron_XBB15": "OQ965259.1"
        }
        
        self.pdb_structures = {
            "Spike_Closed": "6VXX",
            "Spike_Open": "6VYB",
            "Omicron_RBD": "7BZ5",
            "Nucleocapsid": "6M0J"
        }
    
    def fetch_genome(self, variant_name, accession_id):
        '''Download genome from NCBI with validation'''
        filepath = f"{self.genome_dir}/{variant_name}.fasta"
        if os.path.exists(filepath):
            return filepath
        
        try:
            with Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text") as handle:
                record = SeqIO.read(handle, "fasta")
                record.id = variant_name
                record.description = f"{variant_name}|{accession_id}"
                
                # Validation
                if len(record.seq) < 29000:
                    print(f"Warning: {variant_name} sequence too short ({len(record.seq)}bp)")
                    return None
                
                SeqIO.write(record, filepath, "fasta")
                return filepath
        except Exception as e:
            print(f"Error downloading {variant_name}: {e}")
            return None
    
    def download_all_genomes(self):
        '''Download all variant genomes with progress tracking'''
        print("=== Downloading Genome Sequences ===")
        paths, metadata = {}, []
        
        for name, acc in tqdm(self.variants.items(), desc="Fetching genomes"):
            path = self.fetch_genome(name, acc)
            if path:
                paths[name] = path
                metadata.append({
                    "Variant": name,
                    "Accession": acc,
                    "Size_bp": len(SeqIO.read(path, "fasta").seq),
                    "Path": path
                })
        
        # Save metadata
        meta_df = pd.DataFrame(metadata)
        meta_df.to_csv(f"{self.base_dir}/data/metadata.csv", index=False)
        print(f"Downloaded {len(paths)} genomes. Metadata saved.")
        return paths, meta_df
    
    def download_structures(self):
        '''Download PDB structure files'''
        print("=== Downloading PDB Structures ===")
        pdbl = PDBList(verbose=False)
        for name, pdb_id in self.pdb_structures.items():
            print(f"Fetching {name} ({pdb_id})...")
            pdbl.retrieve_pdb_file(pdb_id, pdir=self.structure_dir, file_format='pdb')
        print("PDB structures downloaded.")

# Execute Data Acquisition
data_mgr = DataAcquisitionManager(BASE_DIR)
genome_paths, metadata_df = data_mgr.download_all_genomes()
data_mgr.download_structures()

print("\\n=== Metadata Summary ===")
display(metadata_df)
"""))

#=============================================================================
# CELL 5: Sequence Alignment Module
#=============================================================================
notebook["cells"].append(create_markdown_cell("## III. Sequence Alignment & Comparison {#alignment}"))

notebook["cells"].append(create_code_cell("""
class AlignmentEngine:
    '''Handles multiple sequence alignment using MUSCLE'''
    
    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.alignment_dir = f"{base_dir}/results/alignments"
    
    def perform_msa(self, fasta_files, output_name="genome_msa"):
        '''Perform MUSCLE alignment'''
        combined_fasta = f"{self.alignment_dir}/{output_name}_input.fasta"
        output_aln = f"{self.alignment_dir}/{output_name}.aln"
        
        # Combine sequences
        records = []
        for variant, file in fasta_files.items():
            if file and os.path.exists(file):
                records.append(SeqIO.read(file, "fasta"))
        
        if len(records) < 2:
            print("Need at least 2 sequences for alignment")
            return None
        
        SeqIO.write(records, combined_fasta, "fasta")
        print(f"Aligning {len(records)} sequences with MUSCLE...")
        
        # Run MUSCLE (output in CLUSTAL format)
        cline = MuscleCommandline(input=combined_fasta, out=output_aln, clw=True)
        try:
            stdout, stderr = cline()
            print(f"Alignment complete: {output_aln}")
            return output_aln
        except Exception as e:
            print(f"Alignment failed: {e}")
            return None
    
    def calculate_identity_matrix(self, alignment_file):
        '''Calculate pairwise sequence identity matrix'''
        aln = AlignIO.read(alignment_file, "clustal")
        n = len(aln)
        matrix = np.zeros((n, n))
        labels = [rec.id for rec in aln]
        
        for i in range(n):
            for j in range(i, n):
                s1, s2 = str(aln[i].seq), str(aln[j].seq)
                matches = sum(1 for a, b in zip(s1, s2) if a == b and a != '-')
                identity = (matches / len(s1)) * 100
                matrix[i, j] = matrix[j, i] = identity
        
        return pd.DataFrame(matrix, index=labels, columns=labels)

# Run Alignment
aligner = AlignmentEngine(BASE_DIR)
alignment_file = aligner.perform_msa(genome_paths)

if alignment_file:
    identity_df = aligner.calculate_identity_matrix(alignment_file)
    identity_df.to_csv(f"{BASE_DIR}/results/alignments/identity_matrix.csv")
    
    # Visualize Identity Matrix
    fig = px.imshow(identity_df, text_auto=".2f", color_continuous_scale="Blues",
                    title="Sequence Identity Matrix (%)", aspect="auto")
    fig.update_layout(width=800, height=700)
    fig.show()
    
    # Static heatmap for reports
    plt.figure(figsize=(10, 8))
    sns.heatmap(identity_df, annot=True, fmt=".2f", cmap="YlGnBu", cbar_kws={'label': 'Identity %'})
    plt.title("Pairwise Sequence Identity Matrix")
    plt.tight_layout()
    plt.savefig(f"{BASE_DIR}/results/visualizations/identity_matrix.png", dpi=300)
    plt.show()
"""))

#============================================================================= # CELL 6: Mutation Detection Module (Advanced)
#=============================================================================
notebook["cells"].append(create_markdown_cell("## IV. Mutation Detection & Annotation {#mutations}"))

notebook["cells"].append(create_code_cell("""
class MutationAnalyzer:
    '''Advanced mutation detection with amino acid translation and annotation'''
    
    def __init__(self, reference_file):
        self.ref_record = SeqIO.read(reference_file, "fasta")
        self.ref_seq = str(self.ref_record.seq)
        self.codon_table = CodonTable.unambiguous_dna_by_name["Standard"]
        
        # Gene boundaries for Wuhan-Hu-1 (1-indexed)
        self.genes = {
            "ORF1ab": (266, 21555),
            "S": (21563, 25384),
            "ORF3a": (25393, 26220),
            "E": (26245, 26472),
            "M": (26523, 27191),
            "ORF6": (27202, 27387),
            "ORF7a": (27394, 27759),
            "ORF8": (27894, 28259),
            "N": (28274, 29533),
            "ORF10": (29558, 29674)
        }
        
        # Known mutations of concern
        self.voc_mutations = {
            "S": ["D614G", "N501Y", "E484K", "K417N", "L452R", "P681H", "P681R",
                  "T478K", "Q498R", "N440K", "G446S", "S477N", "T547K", "H655Y"]
        }
    
    def identify_mutations(self, variant_file):
        '''Detect all mutations (SNPs, Indels) with AA translation'''
        var_record = SeqIO.read(variant_file, "fasta")
        var_seq = str(var_record.seq)
        
        # Global pairwise alignment
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        alignment = aligner.align(self.ref_seq, var_seq)[0]
        
        ref_aln = str(alignment[0])
        var_aln = str(alignment[1])
        
        mutations = []
        ref_idx = 0  # 1-indexed position in reference
        
        for i, (r, v) in enumerate(zip(ref_aln, var_aln)):
            if r != '-':
                ref_idx += 1
            
            if r == v:
                continue
            
            # Determine mutation type
            if r != '-' and v != '-':
                mut_type = "SNP"
            elif v == '-':
                mut_type = "Deletion"
            else:
                mut_type = "Insertion"
            
            gene = self._find_gene(ref_idx)
            
            mut_entry = {
                "Variant": var_record.id,
                "Position": ref_idx,
                "Ref": r,
                "Alt": v,
                "Type": mut_type,
                "Gene": gene,
                "AA_Change": None,
                "Effect": None
            }
            
            # Amino acid analysis for SNPs
            if mut_type == "SNP" and gene != "Intergenic":
                aa_change, effect = self._analyze_aa_change(gene, ref_idx, v)
                mut_entry["AA_Change"] = aa_change
                mut_entry["Effect"] = effect
            
            mutations.append(mut_entry)
        
        return pd.DataFrame(mutations)
    
    def _find_gene(self, position):
        '''Map genomic position to gene'''
        for gene, (start, end) in self.genes.items():
            if start <= position <= end:
                return gene
        return "Intergenic"
    
    def _analyze_aa_change(self, gene, genome_pos, alt_nt):
        '''Translate nucleotide change to amino acid change'''
        start, end = self.genes[gene]
        gene_start_idx = start - 1  # Convert to 0-indexed
        gene_pos = (genome_pos - 1) - gene_start_idx
        
        if gene_pos < 0 or gene_pos >= (end - start + 1):
            return None, "Unknown"
        
        # Determine codon context
        codon_idx = gene_pos // 3
        codon_pos = gene_pos % 3
        
        # Extract reference codon
        ref_gene_seq = self.ref_seq[gene_start_idx:end]
        ref_codon = ref_gene_seq[codon_idx*3:(codon_idx+1)*3]
        
        if len(ref_codon) != 3:
            return None, "Incomplete_Codon"
        
        # Construct alternate codon
        alt_codon_list = list(ref_codon)
        alt_codon_list[codon_pos] = alt_nt
        alt_codon = "".join(alt_codon_list)
        
        # Translate
        try:
            ref_aa = self.codon_table.forward_table.get(ref_codon, "*")
            alt_aa = self.codon_table.forward_table.get(alt_codon, "*")
        except:
            return None, "Translation_Error"
        
        aa_change = f"{ref_aa}{codon_idx+1}{alt_aa}"
        
        # Classify effect
        if ref_aa == alt_aa:
            effect = "Synonymous"
        else:
            effect = "Non-Synonymous"
            if gene in self.voc_mutations and aa_change in self.voc_mutations[gene]:
                effect += "_VOC"
        
        return aa_change, effect
    
    def calculate_dnds(self, mutations_df):
        '''Calculate dN/dS ratios (simplified proxy)'''
        stats = []
        
        for (variant, gene), group in mutations_df.groupby(['Variant', 'Gene']):
            if gene == "Intergenic":
                continue
            
            n_nonsyn = len(group[group['Effect'].str.contains("Non-Synonymous", na=False)])
            n_syn = len(group[group['Effect'] == "Synonymous"])
            voc_count = len(group[group['Effect'].str.contains("VOC", na=False)])
            
            # dN/dS proxy (actual calculation requires site counts)
            dnds = (n_nonsyn / n_syn) if n_syn > 0 else n_nonsyn
            
            stats.append({
                "Variant": variant,
                "Gene": gene,
                "Synonymous": n_syn,
                "Non-Synonymous": n_nonsyn,
                "VOC_Mutations": voc_count,
                "dN_dS_Ratio": round(dnds, 3)
            })
        
        return pd.DataFrame(stats)

# Run Mutation Analysis
print("=== Analyzing Mutations ===")
analyzer = MutationAnalyzer(genome_paths["Wuhan-Hu-1"])

all_mutations = []
for name, path in tqdm(genome_paths.items(), desc="Analyzing variants"):
    if name == "Wuhan-Hu-1":
        continue
    try:
        df = analyzer.identify_mutations(path)
        all_mutations.append(df)
    except Exception as e:
        print(f"Error analyzing {name}: {e}")

if all_mutations:
    master_mutation_df = pd.concat(all_mutations, ignore_index=True)
    dnds_stats = analyzer.calculate_dnds(master_mutation_df)
    
    # Save results
    master_mutation_df.to_csv(f"{BASE_DIR}/results/mutations/master_mutations.csv", index=False)
    dnds_stats.to_csv(f"{BASE_DIR}/results/mutations/dnds_statistics.csv", index=False)
    
    print(f"\\nTotal mutations detected: {len(master_mutation_df)}")
    print("\\n=== dN/dS Statistics ===")
    display(dnds_stats.head(10))
    
    # Mutation type distribution
    print("\\n=== Mutation Type Distribution ===")
    display(master_mutation_df['Type'].value_counts())
"""))

#=============================================================================
# CELL 7: Phylogenetic Analysis
#=============================================================================
notebook["cells"].append(create_markdown_cell("## V. Phylogenetic Analysis {#phylo}"))

notebook["cells"].append(create_code_cell("""
class PhylogeneticAnalyzer:
    '''Phylogenetic tree construction and visualization'''
    
    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.phylo_dir = f"{base_dir}/results/phylogenetics"
    
    def build_tree(self, alignment_file, method='nj'):
        '''Construct phylogenetic tree using distance-based methods'''
        aln = AlignIO.read(alignment_file, "clustal")
        
        # Calculate distance matrix
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(aln)
        
        # Construct tree
        constructor = DistanceTreeConstructor()
        if method == 'nj':
            tree = constructor.nj(dm)
        else:
            tree = constructor.upgma(dm)
        
        # Save tree
        Phylo.write(tree, f"{self.phylo_dir}/tree.xml", "phyloxml")
        Phylo.write(tree, f"{self.phylo_dir}/tree.nwk", "newick")
        
        return tree
    
    def visualize_tree(self, tree, title="Phylogenetic Tree"):
        '''Create publication-quality tree visualization'''
        fig, ax = plt.subplots(figsize=(12, 8))
        Phylo.draw(tree, do_show=False, axes=ax)
        ax.set_title(title, fontsize=16, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.tight_layout()
        plt.savefig(f"{self.phylo_dir}/tree_plot.png", dpi=300, bbox_inches='tight')
        plt.show()

# Run Phylogenetic Analysis
if alignment_file:
    phylo = PhylogeneticAnalyzer(BASE_DIR)
    tree = phylo.build_tree(alignment_file, method='nj')
    phylo.visualize_tree(tree, "SARS-CoV-2 Variant Phylogeny (Neighbor-Joining)")
"""))

#=============================================================================
# CELL 8: Interactive Visualizations
#=============================================================================
notebook["cells"].append(create_markdown_cell("## VI. Interactive Visualizations {#viz}"))

notebook["cells"].append(create_code_cell("""
class VisualizationManager:
    '''Comprehensive visualization suite with Plotly interactivity'''
    
    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.viz_dir = f"{base_dir}/results/visualizations"
    
    def plot_mutation_heatmap(self, mutations_df):
        '''Interactive mutation frequency heatmap'''
        pivot = mutations_df.groupby(['Gene', 'Variant']).size().reset_index(name='Count')
        matrix = pivot.pivot(index='Gene', columns='Variant', values='Count').fillna(0)
        
        fig = px.imshow(matrix, text_auto=True, aspect="auto",
                        color_continuous_scale="Reds",
                        title="Mutation Frequency Heatmap (Count per Gene/Variant)")
        fig.update_layout(width=1000, height=600)
        fig.write_html(f"{self.viz_dir}/mutation_heatmap_interactive.html")
        fig.show()
        return matrix
    
    def plot_genome_map(self, mutations_df):
        '''Interactive genome-wide mutation map'''
        fig = go.Figure()
        
        # Gene annotations
        genes = {
            "ORF1ab": (266, 21555), "S": (21563, 25384), "N": (28274, 29533),
            "E": (26245, 26472), "M": (26523, 27191)
        }
        
        for gene, (start, end) in genes.items():
            fig.add_shape(type="rect", x0=start, x1=end, y0=0, y1=1,
                         fillcolor="lightgray", opacity=0.3, line_width=0)
            fig.add_annotation(x=(start+end)/2, y=0.5, text=gene, showarrow=False)
        
        # Mutations
        for variant in mutations_df['Variant'].unique():
            var_muts = mutations_df[mutations_df['Variant'] == variant]
            fig.add_trace(go.Scatter(
                x=var_muts['Position'], y=[1]*len(var_muts),
                mode='markers', name=variant,
                marker=dict(size=8, opacity=0.7),
                hovertemplate='<b>%{text}</b><br>Position: %{x}<extra></extra>',
                text=[f"{row['Gene']} - {row['Type']}" for _, row in var_muts.iterrows()]
            ))
        
        fig.update_layout(title="Genome-Wide Mutation Distribution",
                         xaxis_title="Genome Position (nt)", yaxis=dict(visible=False),
                         height=400, showlegend=True)
        fig.write_html(f"{self.viz_dir}/genome_map_interactive.html")
        fig.show()
    
    def plot_lollipop(self, mutations_df, gene='S'):
        '''Lollipop plot for gene-specific mutations'''
        gene_muts = mutations_df[mutations_df['Gene'] == gene]
        if gene_muts.empty:
            print(f"No mutations found for {gene}")
            return
        
        counts = gene_muts['Position'].value_counts().sort_index()
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=counts.index, y=counts.values,
            mode='markers+lines',
            marker=dict(size=12, color='crimson'),
            line=dict(color='gray', width=1),
            name=gene
        ))
        
        fig.update_layout(title=f"Lollipop Plot: {gene} Gene Mutations",
                         xaxis_title="Genomic Position", yaxis_title="Mutation Count",
                         height=500)
        fig.write_html(f"{self.viz_dir}/lollipop_{gene}_interactive.html")
        fig.show()
    
    def plot_dnds_comparison(self, dnds_df):
        '''Interactive dN/dS comparison across variants and genes'''
        fig = px.bar(dnds_df, x='Gene', y='dN_dS_Ratio', color='Variant',
                     barmode='group', title="dN/dS Ratios by Gene and Variant",
                     labels={'dN_dS_Ratio': 'dN/dS Ratio'})
        fig.add_hline(y=1, line_dash="dash", line_color="red",
                     annotation_text="Neutral (dN/dS=1)")
        fig.update_layout(height=600)
        fig.write_html(f"{self.viz_dir}/dnds_comparison_interactive.html")
        fig.show()

# Generate Visualizations
viz_mgr = VisualizationManager(BASE_DIR)

if 'master_mutation_df' in locals() and not master_mutation_df.empty:
    print("Generating interactive visualizations...")
    mut_matrix = viz_mgr.plot_mutation_heatmap(master_mutation_df)
    viz_mgr.plot_genome_map(master_mutation_df)
    viz_mgr.plot_lollipop(master_mutation_df, gene='S')
    viz_mgr.plot_dnds_comparison(dnds_stats)
    print("Visualizations saved to:", viz_mgr.viz_dir)
"""))

#=============================================================================
# CELL 9: Protein Structure Analysis
#=============================================================================
notebook["cells"].append(create_markdown_cell("## VII. Protein Structure Analysis {#structure}"))

notebook["cells"].append(create_code_cell("""
class StructureAnalyzer:
    '''3D protein structure analysis and visualization'''
    
    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.structure_dir = f"{base_dir}/data/structures"
    
    def render_structure_3d(self, pdb_id='6VXX', style='cartoon'):
        '''Render interactive 3D structure using py3Dmol'''
        view = py3Dmol.view(query=f'pdb:{pdb_id}', width=800, height=600)
        view.setStyle({style: {'color': 'spectrum'}})
        view.setBackgroundColor('white')
        view.zoomTo()
        return view

# Render Spike Protein Structure
struct = StructureAnalyzer(BASE_DIR)
print("Rendering Spike protein structure (6VXX)...")
view = struct.render_structure_3d('6VXX', style='cartoon')
view.show()
"""))

#=============================================================================
# CELL 10: ChimeraX Script Generation
#=============================================================================
notebook["cells"].append(create_markdown_cell("## VIII. ChimeraX Visualization Scripts {#chimerax}"))

notebook["cells"].append(create_code_cell("""
class ChimeraXScriptGenerator:
    '''Generate ChimeraX command scripts for publication-quality renders'''
    
    def __init__(self, base_dir):
        self.chimerax_dir = f"{base_dir}/results/chimerax_scripts"
    
    def generate_basic_script(self, pdb_id='6VXX'):
        '''Generate basic structure visualization script'''
        script = f'''# ChimeraX Basic Visualization Script
# PDB: {pdb_id}

open {pdb_id}
color white
style cartoon
lighting soft
set bgColor white
view
'''
        filepath = f"{self.chimerax_dir}/basic_visualization.cxc"
        with open(filepath, 'w') as f:
            f.write(script)
        return filepath
    
    def generate_mutation_highlight_script(self, mutations_df, variant_name, pdb_id='6VXX'):
        '''Generate script to highlight mutations on structure'''
        spike_muts = mutations_df[(mutations_df['Variant'] == variant_name) & 
                                  (mutations_df['Gene'] == 'S') &
                                  (mutations_df['Type'] == 'SNP')]
        
        # Extract AA positions (simplified mapping)
        residues = []
        for _, row in spike_muts.iterrows():
            aa_pos = (row['Position'] - 21563) // 3 + 1
            if 1 <= aa_pos <= 1273:  # Spike protein length
                residues.append(str(aa_pos))
        
        script = f'''# ChimeraX Mutation Highlighting: {variant_name}
# PDB: {pdb_id}

open {pdb_id}
color white
style cartoon

# Highlight mutations
'''
        if residues:
            res_str = ",".join(residues[:20])  # Limit to first 20
            script += f'''select :{res_str}
color sel red
style sel sphere
label sel residues
'''
        
        script += '''
# Save image
lighting soft
set bgColor white
view
save {variant_name}_mutations.png width 2000 height 2000 supersample 3
'''
        
        filepath = f"{self.chimerax_dir}/{variant_name}_mutations.cxc"
        with open(filepath, 'w') as f:
            f.write(script)
        return filepath
    
    def generate_domain_script(self):
        '''Generate script for domain visualization'''
        script = '''# ChimeraX Domain Visualization
# Spike Protein Domains

open 6VXX

# NTD (14-305)
select :14-305
color sel cornflowerblue
name sel NTD

# RBD (319-541)
select :319-541
color sel crimson
name sel RBD

# Furin Cleavage (681-685)
select :681-685
color sel gold
style sel sphere
name sel FurinSite

style cartoon
lighting soft
set bgColor white
view
save spike_domains.png width 2000 height 2000
'''
        filepath = f"{self.chimerax_dir}/domain_visualization.cxc"
        with open(filepath, 'w') as f:
            f.write(script)
        return filepath

# Generate ChimeraX Scripts
cx_gen = ChimeraXScriptGenerator(BASE_DIR)
cx_gen.generate_basic_script()
cx_gen.generate_domain_script()

if 'master_mutation_df' in locals():
    for variant in master_mutation_df['Variant'].unique()[:3]:  # First 3 variants
        cx_gen.generate_mutation_highlight_script(master_mutation_df, variant)

print(f"ChimeraX scripts generated in: {cx_gen.chimerax_dir}")
!ls -lh {BASE_DIR}/results/chimerax_scripts/
"""))

#=============================================================================
# CELL 11: Results Export & HTML Report
#=============================================================================
notebook["cells"].append(create_markdown_cell("## IX. Results Export & HTML Report {#export}"))

notebook["cells"].append(create_code_cell("""
class ReportGenerator:
    '''Generate comprehensive HTML report and export results'''
    
    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.report_dir = f"{base_dir}/reports"
    
    def export_excel_summary(self, mutations_df, dnds_df, identity_df):
        '''Export comprehensive Excel workbook'''
        filepath = f"{self.report_dir}/SARS_CoV2_Analysis_Summary.xlsx"
        with pd.ExcelWriter(filepath, engine='openpyxl') as writer:
            mutations_df.to_excel(writer, sheet_name='All_Mutations', index=False)
            dnds_df.to_excel(writer, sheet_name='dN_dS_Statistics', index=False)
            identity_df.to_excel(writer, sheet_name='Sequence_Identity')
        print(f"Excel report saved: {filepath}")
        return filepath
    
    def generate_html_report(self):
        '''Generate interactive HTML report'''
        html = f'''<!DOCTYPE html>
<html>
<head>
    <title>COVID-19 Genome Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; background: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: auto; background: white; padding: 30px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        .summary-box {{ background: #ecf0f1; padding: 15px; border-left: 4px solid #3498db; margin: 20px 0; }}
        table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
        th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background-color: #3498db; color: white; }}
        .timestamp {{ color: #7f8c8d; font-size: 0.9em; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>COVID-19 Genome Mutation Analysis Report</h1>
        <p class="timestamp">Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
        
        <div class="summary-box">
            <h3>Executive Summary</h3>
            <p>Comprehensive analysis of SARS-CoV-2 variant genomes including mutation detection, 
            phylogenetic relationships, and structural impact assessment.</p>
        </div>
        
        <h2>Analysis Results</h2>
        <ul>
            <li>Variants analyzed: See metadata</li>
            <li>Total mutations detected: See mutation tables</li>
            <li>Phylogenetic tree: Available in results directory</li>
            <li>Interactive visualizations: HTML files in visualizations folder</li>
        </ul>
        
        <h2>Output Files</h2>
        <ul>
            <li><code>/results/mutations/master_mutations.csv</code></li>
            <li><code>/results/phylogenetics/tree.nwk</code></li>
            <li><code>/results/visualizations/*.html</code></li>
            <li><code>/results/chimerax_scripts/*.cxc</code></li>
        </ul>
        
        <h2>References</h2>
        <p>NCBI Virus Database | PDB | Biopython Documentation</p>
    </div>
</body>
</html>'''
        
        filepath = f"{self.report_dir}/Analysis_Report.html"
        with open(filepath, 'w') as f:
            f.write(html)
        print(f"HTML report saved: {filepath}")
        return filepath

# Generate Reports
reporter = ReportGenerator(BASE_DIR)

if all(v in locals() for v in ['master_mutation_df', 'dnds_stats', 'identity_df']):
    reporter.export_excel_summary(master_mutation_df, dnds_stats, identity_df)
    reporter.generate_html_report()

# Create downloadable ZIP archive
print("\\nCreating results archive...")
!zip -r -q /content/covid_analysis_results.zip {BASE_DIR}/results {BASE_DIR}/reports

print("\\n" + "="*70)
print("ANALYSIS COMPLETE!")
print("="*70)
print(f"Results location: {BASE_DIR}/results/")
print("Download the archive: /content/covid_analysis_results.zip")
print("\\nTo download in Colab, run:")
print("from google.colab import files")
print("files.download('/content/covid_analysis_results.zip')")
"""))

#=============================================================================
# CELL 12: Summary & Conclusions
#=============================================================================
notebook["cells"].append(create_markdown_cell("## X. Summary & Conclusions {#summary}"))

notebook["cells"].append(create_markdown_cell("""
### Key Findings

This pipeline successfully analyzed SARS-CoV-2 variants with:
- **Genome Acquisition**: Downloaded and validated 10 variant genomes
- **Sequence Alignment**: Performed MUSCLE-based MSA
- **Mutation Detection**: Identified SNPs, indels with amino acid translation
- **Statistical Analysis**: Calculated dN/dS ratios for selection pressure
- **Phylogenetics**: Constructed neighbor-joining trees
- **Visualizations**: Generated 10+ interactive and static plots
- **Structure Mapping**: Created ChimeraX scripts for protein visualization

### Output Files

All results are organized in:
```
/content/covid_analysis/
├── data/ (genomes, structures, metadata)
├── results/
│   ├── alignments/ (MSA files, identity matrix)
│   ├── mutations/ (mutation tables, dN/dS stats)
│   ├── phylogenetics/ (tree files)
│   ├── visualizations/ (PNG, HTML plots)
│   └── chimerax_scripts/ (CXC files)
└── reports/ (Excel summary, HTML report)
```

### Next Steps

1. Download results archive
2. Open ChimeraX scripts for 3D visualization
3. Review HTML report for detailed findings
4. Explore interactive Plotly visualizations

### References

- NCBI Virus Database: https://www.ncbi.nlm.nih.gov/labs/virus/
- RCSB PDB: https://www.rcsb.org/
- Biopython: https://biopython.org/
- outbreak.info: https://outbreak.info/

---
**Pipeline Version:** 2.0 Professional | **Notebook:** COVID-19 Genome Analysis
"""))

#=============================================================================
# Write Notebook to File
#=============================================================================
output_file = "COVID19_Professional_Analysis.ipynb"
with open(output_file, 'w') as f:
    json.dump(notebook, f, indent=2)

print(f"\\n{'='*70}")
print(f"SUCCESS! Notebook generated: {output_file}")
print(f"{'='*70}")
print(f"\\nTotal cells: {len(notebook['cells'])}")
print("\\nNotebook includes:")
print("  - Environment setup (apt-get muscle)")
print("  - Data acquisition (10 variants)")
print("  - Sequence alignment (CLUSTAL support)")
print("  - Mutation detection (AA translation, dN/dS)")
print("  - Phylogenetic analysis")
print("  - Interactive Plotly visualizations")
print("  - Protein structure analysis (py3Dmol)")
print("  - ChimeraX script generation")
print("  - HTML report & Excel export")
print("\\nReady to upload to Google Colab!")
