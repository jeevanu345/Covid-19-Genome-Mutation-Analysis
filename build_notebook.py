import json
import os

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

notebook = {
    "cells": [],
    "metadata": {
        "colab": {
            "name": "COVID19_Genome_Mutation_Analysis_Professional.ipynb",
            "provenance": []
        },
        "kernelspec": {
            "display_name": "Python 3",
            "name": "python3"
        },
        "language_info": {
            "codemirror_mode": {
                "name": "ipython",
                "version": 3
            },
            "file_extension": ".py",
            "mimetype": "text/x-python",
            "name": "python",
            "version": "3.10.12"
        }
    },
    "nbformat": 4,
    "nbformat_minor": 0
}

# --- CELL 1: Header ---
notebook["cells"].append(create_markdown_cell("""
# COVID-19 Genome Mutation Analysis & Protein Structure Visualization Pipeline

**Project Title:** COVID-19 Genome Mutation Analysis & Protein Structure Visualization Pipeline

This notebook implements a professional-grade bioinformatics pipeline for analyzing SARS-CoV-2 variants.
"""))

# --- CELL 2: Setup ---
notebook["cells"].append(create_markdown_cell("## I. Environment Setup & Dependencies"))

notebook["cells"].append(create_code_cell("""
# Install System Dependencies (Crucial for MUSCLE)
!apt-get update -qq && apt-get install -y muscle -qq

# Install Python Libraries
!pip install -q biopython py3Dmol pandas seaborn plotly scikit-bio networkx

import os
import sys
import json
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
import py3Dmol
import networkx as nx
from tqdm.notebook import tqdm
from datetime import datetime

# BioPython Imports
from Bio import SeqIO, AlignIO, Phylo, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Data import CodonTable
from Bio.PDB import PDBList

# Configure Entrez
Entrez.email = "researcher@example.com"

# Setup Directory Structure
BASE_DIR = "/content/covid_analysis"
DIRS = [
    f"{BASE_DIR}/data/genomes",
    f"{BASE_DIR}/data/proteins",
    f"{BASE_DIR}/data/structures",
    f"{BASE_DIR}/results/alignments",
    f"{BASE_DIR}/results/mutations",
    f"{BASE_DIR}/results/phylogenetics",
    f"{BASE_DIR}/results/visualizations",
    f"{BASE_DIR}/results/chimerax_scripts",
    f"{BASE_DIR}/reports"
]

for d in DIRS:
    os.makedirs(d, exist_ok=True)

print(f"Environment Setup Complete. Working directory: {BASE_DIR}")
"""))

# --- CELL 3: Data Acquisition ---
notebook["cells"].append(create_markdown_cell("## II. Data Acquisition Module"))

notebook["cells"].append(create_code_cell("""
class DataAcquisitionManager:
    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.genome_dir = f"{base_dir}/data/genomes"
        self.structure_dir = f"{base_dir}/data/structures"
        
        # Comprehensive Variant List
        self.variants = {
            "Wuhan-Hu-1": "NC_045512.2",       # Reference
            "Alpha": "CM034978.1",             # B.1.1.7 
            "Beta": "OM530757.1",              # B.1.351
            "Gamma": "OM530769.1",             # P.1
            "Delta": "OM287553.1",             # B.1.617.2
            "Omicron_BA1": "NC_055123.1",      # BA.1
            "Omicron_BA2": "OM296767.1",       # BA.2
            "Omicron_BA5": "OX315743.1",       # BA.5
            "Omicron_XBB": "OQ965259.1",       # XBB.1.5
            "Epsilon": "OM287554.1"            # B.1.427
        }
        
        self.pdb_ids = {
            "Spike_Closed": "6VXX",
            "Spike_Open": "6VYB",
            "Omicron_RBD": "7BZ5",
            "Nucleocapsid": "6M0J"
        }

    def fetch_genome_from_ncbi(self, variant_name, accession_id):
        filepath = f"{self.genome_dir}/{variant_name}.fasta"
        if os.path.exists(filepath):
            return filepath
            
        print(f"Downloading {variant_name} ({accession_id})...")
        try:
            with Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text") as handle:
                seq_record = SeqIO.read(handle, "fasta")
                # Clean header
                seq_record.id = variant_name
                seq_record.description = f"{variant_name}|{accession_id}"
                
                # Validation: Minimal length check for complete genome
                if len(seq_record.seq) < 29000:
                    print(f"Warning: Sequence {variant_name} is too short ({len(seq_record.seq)} bp). Skipping.")
                    return None
                    
                SeqIO.write(seq_record, filepath, "fasta")
            return filepath
        except Exception as e:
            print(f"Failed to download {variant_name}: {e}")
            return None

    def download_all_variants(self):
        print("--- Starting Genome Download ---")
        paths = {}
        metadata = []
        
        for name, acc in tqdm(self.variants.items(), desc="Fetching Genomes"):
            path = self.fetch_genome_from_ncbi(name, acc)
            if path:
                paths[name] = path
                size = os.path.getsize(path)
                metadata.append({"Variant": name, "Accession": acc, "Size_Bytes": size, "Path": path})
        
        # Save Metadata
        pd.DataFrame(metadata).to_csv(f"{self.base_dir}/data/metadata.csv", index=False)
        print("Metadata saved to metadata.csv")
        return paths

    def download_structures(self):
        print("--- Starting PDB Structure Download ---")
        pdbl = PDBList(verbose=False)
        for name, pid in self.pdb_ids.items():
            print(f"Fetching {name} ({pid})...")
            pdbl.retrieve_pdb_file(pid, pdir=self.structure_dir, file_format='pdb')

# Instantiate and Run
data_manager = DataAcquisitionManager(BASE_DIR)
genome_paths = data_manager.download_all_variants()
data_manager.download_structures()
"""))

# Write to file
with open("COVID19_Genome_Mutation_Analysis_Professional.ipynb", "w") as f:
    json.dump(notebook, f, indent=2)

print("Notebook generated: COVID19_Genome_Mutation_Analysis_Professional.ipynb")
