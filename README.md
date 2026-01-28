# COVID-19 Genome Mutation Analysis Pipeline

**Professional-Grade Bioinformatics Pipeline for SARS-CoV-2 Variant Analysis**

Version 2.0 Professional

## Overview

This comprehensive Jupyter notebook pipeline provides end-to-end analysis of SARS-CoV-2 genomic variants, including mutation detection, phylogenetic analysis, and protein structure visualization. Designed for researchers, bioinformaticians, and computational biologists working with viral genomics.

## Key Features

- **Automated Genome Acquisition**: Download and validate variant sequences from NCBI GenBank
- **Multiple Sequence Alignment**: MUSCLE-based alignment with identity matrix calculation
- **Mutation Detection**: Nucleotide and amino acid-level mutation identification
- **Statistical Analysis**: dN/dS ratio calculations for selection pressure assessment
- **Phylogenetic Tree Construction**: Distance-based tree building (NJ/UPGMA)
- **Interactive Visualizations**: Plotly-based genome maps, heatmaps, and lollipop plots
- **3D Protein Structure Analysis**: py3Dmol integration for structure visualization
- **ChimeraX Script Generation**: Publication-quality rendering scripts
- **Automated Reporting**: Excel summaries and HTML reports

## Requirements

### System Dependencies

- Python 3.8+
- MUSCLE alignment tool
- Jupyter Notebook or Google Colab

### Python Libraries

```
biopython>=1.79
py3Dmol>=2.0
pandas>=1.3
numpy>=1.21
matplotlib>=3.4
seaborn>=0.11
plotly>=5.3
scikit-bio>=0.5
networkx>=2.6
scipy>=1.7
openpyxl>=3.0
```

## Installation

### Local Installation

```bash
# Clone repository (if applicable)
git clone <repository-url>
cd covid-analysis-pipeline

# Install system dependencies (Ubuntu/Debian)
sudo apt-get update
sudo apt-get install muscle

# Install Python dependencies
pip install -r requirements.txt

# Launch Jupyter
jupyter notebook COVID19_Professional_Analysis.ipynb
```

### Google Colab

1. Upload the notebook to Google Drive
2. Open with Google Colab
3. Run all cells sequentially (dependencies auto-install)

## Quick Start

The notebook is designed to run sequentially from top to bottom:

1. **Environment Setup**: Install dependencies and create directory structure
2. **Data Acquisition**: Download 10 SARS-CoV-2 variant genomes and PDB structures
3. **Sequence Alignment**: Perform multiple sequence alignment with MUSCLE
4. **Mutation Detection**: Identify SNPs, indels, and amino acid changes
5. **Phylogenetic Analysis**: Construct evolutionary trees
6. **Visualization**: Generate interactive plots and genome maps
7. **Structure Analysis**: Map mutations to 3D protein structures
8. **Report Generation**: Export results to Excel and HTML

## Analyzed Variants

The pipeline includes the following SARS-CoV-2 variants:

- **Wuhan-Hu-1** (Reference): NC_045512.2
- **Alpha (B.1.1.7)**: CM034978.1
- **Beta (B.1.351)**: OM530757.1
- **Gamma (P.1)**: OM530769.1
- **Delta (B.1.617.2)**: OM287553.1
- **Omicron BA.1**: NC_055123.1
- **Omicron BA.2**: OM296767.1
- **Omicron BA.4**: ON481610.1
- **Omicron BA.5**: OX315743.1
- **Omicron XBB.1.5**: OQ965259.1

## Output Structure

```
covid_analysis/
├── data/
│   ├── genomes/           # Downloaded FASTA sequences
│   ├── structures/        # PDB structure files
│   └── metadata.csv       # Variant metadata
├── results/
│   ├── alignments/
│   │   ├── genome_msa.aln
│   │   └── identity_matrix.csv
│   ├── mutations/
│   │   ├── master_mutations.csv
│   │   └── dnds_statistics.csv
│   ├── phylogenetics/
│   │   ├── tree.nwk
│   │   └── tree_plot.png
│   ├── visualizations/
│   │   ├── mutation_heatmap_interactive.html
│   │   ├── genome_map_interactive.html
│   │   └── [additional plots]
│   └── chimerax_scripts/
│       ├── basic_visualization.cxc
│       ├── domain_visualization.cxc
│       └── [variant-specific scripts]
└── reports/
    ├── SARS_CoV2_Analysis_Summary.xlsx
    └── Analysis_Report.html
```

## Key Outputs

### 1. Mutation Tables

**master_mutations.csv** contains:
- Variant name
- Genomic position
- Reference/alternate nucleotides
- Mutation type (SNP/insertion/deletion)
- Affected gene
- Amino acid change
- Functional effect classification

### 2. Statistical Analysis

**dnds_statistics.csv** includes:
- dN/dS ratios per gene and variant
- Synonymous vs. non-synonymous mutation counts
- Variants of concern (VOC) mutation tracking

### 3. Phylogenetic Trees

- Newick format (.nwk) for external tools
- PhyloXML format for advanced analysis
- Publication-quality PNG renders

### 4. Interactive Visualizations

- Mutation frequency heatmaps
- Genome-wide mutation distribution maps
- Gene-specific lollipop plots
- dN/dS comparison charts
- 3D protein structure viewers

### 5. ChimeraX Scripts

Ready-to-run command scripts for:
- Basic structure visualization
- Mutation highlighting on protein structures
- Domain-specific coloring
- Publication-quality rendering

## Usage Examples

### Running the Complete Pipeline

```python
# Execute all cells in sequence
# Or run specific sections:

# Data acquisition only
data_mgr = DataAcquisitionManager(BASE_DIR)
genome_paths, metadata_df = data_mgr.download_all_genomes()

# Mutation analysis for specific variant
analyzer = MutationAnalyzer(genome_paths["Wuhan-Hu-1"])
mutations = analyzer.identify_mutations(genome_paths["Omicron_BA5"])
```

### Customizing Variants

Edit the `variants` dictionary in the `DataAcquisitionManager` class:

```python
self.variants = {
    "Custom_Variant": "ACCESSION_ID",
    # Add more variants
}
```

### Visualizing Specific Genes

```python
viz_mgr = VisualizationManager(BASE_DIR)
viz_mgr.plot_lollipop(master_mutation_df, gene='N')  # Nucleocapsid gene
```

## Advanced Features

### dN/dS Ratio Interpretation

- **dN/dS < 1**: Purifying selection (conserved regions)
- **dN/dS = 1**: Neutral evolution
- **dN/dS > 1**: Positive selection (adaptive mutations)

### Variants of Concern (VOC) Tracking

The pipeline automatically flags known VOC mutations in the Spike protein:
- D614G, N501Y, E484K, K417N, L452R, P681H/R
- T478K, Q498R, N440K, G446S, and others

### ChimeraX Integration

Generated `.cxc` scripts can be opened directly in UCSF ChimeraX:

```bash
chimerax results/chimerax_scripts/domain_visualization.cxc
```

## Data Sources

- **Genomes**: NCBI GenBank Virus Database
- **Protein Structures**: RCSB Protein Data Bank (PDB)
- **Reference Sequence**: Wuhan-Hu-1 (NC_045512.2)

## Methodology

### Alignment Algorithm

MUSCLE (Multiple Sequence Comparison by Log-Expectation) with default parameters for high-accuracy nucleotide alignment.

### Phylogenetic Method

Neighbor-Joining (NJ) tree construction based on sequence identity distance matrices.

### Mutation Calling

Pairwise global alignment against reference with nucleotide-to-amino acid translation using standard genetic code.

## Troubleshooting

### Common Issues

**MUSCLE not found**
```bash
sudo apt-get install muscle
# Or on macOS: brew install muscle
```

**NCBI connection errors**
- Check internet connection
- Verify Entrez.email is set
- Try alternative accession IDs if sequences are deprecated

**Memory issues with large datasets**
- Reduce number of variants analyzed
- Process genes individually
- Use high-memory runtime in Colab

### Performance Optimization

- For faster processing, reduce the number of variants
- Skip structure downloads if only analyzing sequences
- Use cached alignment files to avoid re-running MUSCLE

## Citation

If you use this pipeline in your research, please cite:

```
COVID-19 Genome Mutation Analysis Pipeline v2.0
Available at: [repository URL]
```

Key dependencies to cite:
- Biopython: Cock et al. (2009) Bioinformatics
- MUSCLE: Edgar (2004) Nucleic Acids Research
- Plotly: Plotly Technologies Inc.

## References

1. NCBI Virus Database: https://www.ncbi.nlm.nih.gov/labs/virus/
2. RCSB PDB: https://www.rcsb.org/
3. Biopython Documentation: https://biopython.org/
4. MUSCLE Aligner: https://www.drive5.com/muscle/
5. outbreak.info: https://outbreak.info/
6. UCSF ChimeraX: https://www.cgl.ucsf.edu/chimerax/

## Contributing

Contributions are welcome. Please submit issues or pull requests for:
- Additional variant accessions
- Enhanced visualization features
- Performance improvements
- Bug fixes

## License

This pipeline is provided for research and educational purposes. Individual components may have their own licenses (Biopython: BSD, MUSCLE: public domain, etc.).

## Contact

For questions, issues, or collaboration:
- Open an issue in the repository
- Contact: researcher@example.com

## Changelog

### Version 2.0 Professional
- Complete rewrite with modular class structure
- Added ChimeraX script generation
- Interactive Plotly visualizations
- Comprehensive Excel and HTML reporting
- Enhanced mutation annotation with VOC tracking
- Automated directory structure creation
- Improved error handling and validation

### Version 1.0
- Initial release with basic mutation detection
- MUSCLE alignment integration
- Simple phylogenetic analysis

---

**Last Updated**: January 2025
**Platform Tested**: Google Colab, Jupyter Notebook, Python 3.10
**Status**: Production-ready
