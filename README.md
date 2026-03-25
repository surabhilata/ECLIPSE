# ECLIPSE —> Dark Proteome Exploration of ESKAPE Pathogens

**ECLIPSE** (ESKAPE Connectome Linkage and Inference for Proteome Sequence Exploration) is a modular, pathogen-agnostic computational pipeline for the systematic identification and prioritisation of functionally uncharacterised ("dark") protein families in bacterial panproteomes.

ECLIPSE embeds target-pathogen proteomes within the global sequence similarity network of the (https://uniprot3d.org/) (AFDB90v4, UniRef v.2022_03) and identifies connected components composed entirely of unannotated proteins. A taxonomic diversity framework and the Dark Proteome Prioritisation Score (DPPS) rank dark candidates by biological relevance for experimental follow-up.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Structure](#pipeline-structure)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input Files](#input-files)
- [Usage](#usage)
- [Output Files](#output-files)
- [Applying ECLIPSE to Other Pathogens](#applying-eclipse-to-other-pathogens)
- [Citation](#citation)
- [Contact](#contact)

---

## Overview

ECLIPSE is implemented as three sequential Jupyter notebooks:

| Notebook | Description |
|----------|-------------|
| `ECLIPSE_PartI.ipynb` | Atlas mapping, brightness classification, taxonomic diversity analysis |
| `ECLIPSE_PartII.ipynb` | Two-track stratification into Pseudomonas-specific and ESKAPE-enriched dark components |
| `ECLIPSE_DPPS_Scoring.ipynb` | DPPS scoring, weight sensitivity analysis, and visualisations |

Applied to the *P. aeruginosa* panproteome (3,460,657 proteins from 635 strains), ECLIPSE identified 120,985 proteins (4%) residing in completely dark connected components and prioritised four Tier I candidates using the DPPS framework.

---

## Pipeline Structure

```
ESKAPE pathogen FASTA files
        │
        ▼
┌─────────────────────────────────────┐
│  ECLIPSE Part I                     │
│  • MMseqs2 easy-search vs AFDB90v4  │
│  • Community + component mapping    │
│  • Brightness classification        │
│  • Taxonomic diversity analysis     │
│  Output: eclipse.csv                │
└──────────────────┬──────────────────┘
                   │
                   ▼
┌─────────────────────────────────────┐
│  ECLIPSE Part II                    │
│  • Track A: Pseudomonas-specific    │
│    (ESKAPE_proportion == 1.0)       │
│  • Track B: ESKAPE-enriched         │
│    (ESKAPE_proportion >= 0.5)       │
│  Output: two track CSV files        │
└──────────────────┬──────────────────┘
                   │
                   ▼
┌─────────────────────────────────────┐
│  ECLIPSE DPPS Scoring               │
│  • Median length filter (≥300 aa)   │
│  • MMseqs2 clustering (30%/80%)     │
│  • DPPS composite scoring (S1–S5)   │
│  • Monte Carlo sensitivity (n=500)  │
│  • Tier I–IV assignment             │
│  Output: scored CSV + figures       │
└─────────────────────────────────────┘
```

---

## Requirements

### Python packages

```
python >= 3.8
pandas >= 1.3
numpy >= 1.21
matplotlib >= 3.4
seaborn >= 0.11
tqdm
pyyaml
biopython
jupyter
```

### External tools

| Tool | Version | Purpose |
|------|---------|---------|
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | 17-b804f | Sequence search and clustering |

### Atlas data files (required, not included)

The following large Atlas data files must be downloaded separately from the (https://doi.org/10.5281/zenodo.19119408) before running Part I and Part II:

| File | Size | Used in |
|------|------|---------|
| `AFDB90v4_cc_data_uniprot_community_taxonomy_map_with_brightness.csv` | ~4 GB | Part I, Part II |
| `AFDB90v4_dust_uniprot_community_taxonomy_map_with_brightness.csv` | ~4 GB | Part I, Part II |
| `AFDB90v4_subgraphs_summary.csv` | ~500 MB | Part I |

> **Note:** These files are produced by Durairaj et al. (2023) and are available from the Protein Universe Atlas resource. See the Atlas paper for details (https://www.nature.com/articles/s41586-023-06622-3)
---

## Installation

```bash
No installation is required beyond having Jupyter and the standard Python scientific stack. Simply download or clone the repository, open each notebook in Jupyter, and run the cells from top to bottom in order.

```bash
git clone https://github.com/surabhilata/ECLIPSE.git
cd ECLIPSE
jupyter notebook
```

The only external tool that must be installed separately is **MMseqs2** (version 17-b804f), which is used before running Part I to search your proteome against the Atlas:

```bash
conda install -c conda-forge -c bioconda mmseqs2=17-b804f
```

---

## Input Files

### Part I inputs

| File | Description |
|------|-------------|
| `PA.m8` | MMseqs2 easy-search output — your pathogen proteome searched against AFDB90v4. Tab-separated, standard BLAST-6 format with columns: queryID, targetID, fident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bitscore |
| Atlas data files | See [Requirements](#requirements) above |

**To generate `PA.m8`**, run MMseqs2 easy-search against the Atlas AFDB90 target database:

```bash
After installation run in terminal: mmseqs easy-search PA_faa.tar.gz AFDBv4_90.fasta PA.m8 tmp --max-seqs 1
```

### Part II inputs

| File | Description |
|------|-------------|
| `eclipse.csv` | Output from Part I — proteins with mapped component IDs and brightness values |
| Atlas data files | Same as Part I — needed for P. aeruginosa proportion calculation |

### DPPS Scoring inputs

| File | Description |
|------|-------------|
| `eclipse.csv` | Part I output |
| `mapped_p_aer_dataset_pseudomonas_specific_dark_components.csv` | Part II Track A output |
| `mapped_p_aer_dataset_eskape_enriched_dark_components.csv` | Part II Track B output |

---

## Usage

Run the three notebooks in order. Each notebook must be run to completion before starting the next.

### Step 1 — Run Part I

```bash
jupyter notebook ECLIPSE_PartI.ipynb
```

Before running, update the file paths at the top of the notebook:

```python
# Update these paths for your system
atlas_search_results_info = ['./PA.m8']           # your MMseqs2 output

atlas_datafiles = [
    'AFDB90v4_cc_data_uniprot_community_taxonomy_map_with_brightness.csv',
    'AFDB90v4_dust_uniprot_community_taxonomy_map_with_brightness.csv'
]
```

**Runtime:** approximately 30-60 minutes depending on dataset size and available RAM. At least 32 GB RAM is recommended for loading the full Atlas data files.

### Step 2 — Run Part II

```bash
jupyter notebook ECLIPSE_PartII.ipynb
```

Update the file paths at the top of the notebook:

```python
# NOTE: you need to modify the paths here for your system
atlas_datafiles = [
    'AFDB90v4_cc_data_uniprot_community_taxonomy_map_with_brightness.csv',
    'AFDB90v4_dust_uniprot_community_taxonomy_map_with_brightness.csv'
]

mapped_p_aer_dataset = pd.read_csv('eclipse_search_results_component_dark.csv')
```

**Runtime:** approximately 30-60 minutes.

### Step 3 — Run DPPS Scoring

```bash
jupyter notebook ECLIPSE_DPPS_Scoring.ipynb
```

Update the configuration block at the top of the notebook:

```python
ECLIPSE_CSV        = './eclipse.csv'           # Part I output
COMPONENTS_CSV     = './AFDB90v4_subgraphs_summary.csv'
N_STRAINS          = 635                        # total strains in your dataset
```

**Runtime:** approximately 30–60 minutes including Monte Carlo sensitivity analysis.

---

## Output Files

### Part I outputs

| File | Description |
|------|-------------|
| `eclipse.csv` | All query proteins with mapped community IDs, component IDs, brightness values, and taxonomic diversity metrics (ESKAPE_proportion, ESKAPE_genus_evenness, ESKAPE_relative_evenness) |

### Part II outputs

| File | Description |
|------|-------------|
| `mapped_p_aer_dataset_pseudomonas_specific_dark_components.csv` | Track A — Pseudomonas-specific dark components (ESKAPE_proportion == 1.0, all genus == Pseudomonas) |
| `mapped_p_aer_dataset_eskape_enriched_dark_components.csv` | Track B — ESKAPE-enriched dark components (ESKAPE_proportion >= 0.5, excluding Track A) |

### DPPS Scoring outputs

| File | Description |
|------|-------------|
| `dpps_scored_pseudomonas_specific.csv` | Track A components with DPPS sub-scores, composite score, and tier assignment |
| `dpps_scored_eskape_enriched.csv` | Track B components with DPPS sub-scores (S1–S5), composite score, and tier assignment |
| `dpps_representatives_pseudomonas_specific.csv` | Track A non-redundant representative sequences with scores |
| `dpps_representatives_eskape_enriched.csv` | Track B non-redundant representative sequences with scores |
| `dpps_distribution.pdf / .png` | DPPS histogram by tier for both tracks |
| `dpps_scatter_s2_s3.pdf / .png` | Scatter plot of S2 vs S3 coloured by DPPS |
| `dpps_heatmap_tier1.pdf / .png` | Sub-score heatmap for Tier I candidates |
| `dpps_sensitivity.pdf / .png` | Monte Carlo weight sensitivity scatter plot |

### DPPS tier definitions

| Tier | DPPS range | Priority |
|------|-----------|----------|
| Tier I | ≥ 0.75 | Highest priority — experimental follow-up |
| Tier II | 0.50–0.75 | High priority |
| Tier III | 0.25–0.50 | Moderate priority |
| Tier IV | < 0.25 | Background |

---

## Applying ECLIPSE to Other Pathogens

ECLIPSE is designed to be pathogen-agnostic. To apply it to a different ESKAPE pathogen:

1. Download protein FASTA files for your target pathogen from [PATRIC](https://www.patricbrc.org) retaining only strains with complete genome annotation.

2. Run MMseqs2 easy-search against AFDB90v4 to generate your `.m8` search results file.

3. Update the file paths in each notebook to point to your new `.m8` file and output directory.

4. Update the strain count parameter in DPPS Scoring:
   ```python
   N_STRAINS = [your total strain count]  # e.g. 1200 for K. pneumoniae
   ```

5. Run Parts I, II, and DPPS Scoring in order.

No other changes to the code are required. The DPPS weights and tier thresholds were designed to be transferable across ESKAPE organisms without pathogen-specific calibration.

---

## Repository Structure

```
ECLIPSE/
├── README.md
├── LICENSE
├── requirements.txt
├── ECLIPSE_PartI.ipynb          # Atlas mapping and darkness classification
├── ECLIPSE_PartII.ipynb         # Two-track stratification
├── ECLIPSE_DPPS_Scoring.ipynb   # DPPS scoring and visualisations
└── example_data/
    └── example_input.m8         # Small example MMseqs2 output for testing
```

---

## Citation



Please also cite the Protein Universe Atlas:

> Durairaj J. et al. (2023). *Exploring the functional universe of protein families.* Nature, 626, 582–589. https://doi.org/10.1038/s41586-023-06716-0

And MMseqs2:

> Steinegger M. & Söding J. (2017). *MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.* Nature Biotechnology, 35, 1026–1028. https://doi.org/10.1038/nbt.3988

---

## Contact

**Surabhi Lata** — Department of Molecular Structural Biology, Helmholtz Centre for Infection Research (HZI), Braunschweig, Germany
email: surabhi.lata@helmholtz-hzi.de
       surabhilata94@gmail.com
For questions about the pipeline, please open a GitHub Issue or contact the authors directly.

---

## Licence

This project is licensed under the MIT Licence. See [LICENSE](LICENSE) for details.
