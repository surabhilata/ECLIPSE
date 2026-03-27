# ECLIPSE — Dark Proteome Exploration of ESKAPE Pathogens

**ECLIPSE** (ESKAPE Connectome Linkage and Inference for Proteome Sequence Exploration) is a modular, pathogen-agnostic computational pipeline for the systematic identification and prioritisation of functionally uncharacterised ("dark") protein families in bacterial panproteomes.

ECLIPSE maps target-pathogen proteomes onto the global sequence similarity network of the [Protein Universe Atlas](https://uniprot3d.org/) (AFDB90v4, UniRef v.2022_03), identifies protein families that are completely dark at the connected component level, characterises their taxonomic diversity across ESKAPE pathogens, and ranks them using the Dark Proteome Prioritisation Score (DPPS).

> **Associated publication:**
> Lata S. & Heinz D.W. (2025). *ECLIPSE: Exploring the dark proteome of ESKAPE pathogens through the sequence similarity network of the Protein Universe Atlas.* Bioinformatics, Oxford University Press. [DOI to be added upon acceptance]

---

## Table of Contents

- [Pipeline overview](#pipeline-overview)
- [Requirements](#requirements)
- [How to run](#how-to-run)
- [Input and output files](#input-and-output-files)
- [DPPS scoring system](#dpps-scoring-system)
- [Applying ECLIPSE to another ESKAPE pathogen](#applying-eclipse-to-another-eskape-pathogen)
- [Citation](#citation)
- [Contact](#contact)

---

## Pipeline overview

ECLIPSE runs as three sequential Jupyter notebooks. Each notebook must be completed before starting the next.

```
Your pathogen FASTA files
        |
        v  MMseqs2 easy-search (run before Part I)
        |
+-------------------------------------------+
|  ECLIPSE_PartI.ipynb                      |
|                                           |
|  - Load MMseqs2 search results (.m8)      |
|  - Map proteins -> Atlas communities      |
|  - Map communities -> connected components|
|  - Classify dark communities (0%)         |
|  - Classify dark components (0%)          |
|  - Compute ESKAPE taxonomic diversity:    |
|    - ESKAPE_proportion                    |
|    - ESKAPE_genus_evenness                |
|    - ESKAPE_relative_evenness             |
|                                           |
|  Output -> eclipse.csv                    |
+------------------+------------------------+
                   |
                   v
+-------------------------------------------+
|  ECLIPSE_PartII.ipynb                     |
|                                           |
|  - Load eclipse.csv + Atlas taxonomy      |
|  - Add target species proportion per      |
|    component from Atlas taxonomy files    |
|  - Track A -- Pathogen-specific:          |
|    ESKAPE_proportion == 1.0               |
|    all Atlas genus == target genus        |
|    e.g. 83 components for P. aeruginosa  |
|  - Track B -- ESKAPE-enriched:            |
|    ESKAPE_proportion >= 0.5               |
|    excluding Track A                      |
|    e.g. 215 components for P. aeruginosa |
|                                           |
|  Output -> two track CSV files            |
+------------------+------------------------+
                   |
                   v
+-------------------------------------------+
|  ECLIPSE_DPPS_score.ipynb                 |
|                                           |
|  - Median length filter (>= 300 aa)       |
|  - MMseqs2 easy-cluster (30% id, 80% cov) |
|  - One representative per cluster         |
|  - DPPS scoring:                          |
|    S1 darkness                            |
|    S2 target species proportion           |
|    S3 AMR-clade specificity               |
|    S4 strain coverage                     |
|    S5 ESKAPE enrichment (Track B only)    |
|  - Tier I-IV assignment                   |
|  - Monte Carlo weight sensitivity (n=500) |
|  - Figures and summary report             |
|                                           |
|  Output -> scored CSVs + figures          |
+-------------------------------------------+
```

---

## Requirements

### Python packages

All standard packages — install via pip if not already available:

```
pandas    numpy    matplotlib    seaborn    tqdm    pyyaml    biopython
```

### External tool

| Tool | Version | Used for |
|------|---------|----------|
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | 17-b804f | Sequence search (before Part I) and clustering (Part III) |

Install MMseqs2:

```bash
conda install -c conda-forge -c bioconda mmseqs2=17-b804f
```

### Atlas data files

The following large files from the Protein Universe Atlas must be downloaded separately and placed in the same directory as the notebooks. They are not included in this repository due to their size (~8 GB total).

| File | Used in |
|------|---------|
| `AFDB90v4_cc_data_uniprot_community_taxonomy_map_with_brightness.csv` | Part I, Part II |
| `AFDB90v4_dust_uniprot_community_taxonomy_map_with_brightness.csv` | Part I, Part II |
| `AFDB90v4_subgraphs_summary.csv` | Part I |

Download from Zenodo: https://doi.org/10.5281/zenodo.19119408

These files are produced by Durairaj et al. (2023) — see the original paper: https://www.nature.com/articles/s41586-023-06622-3

---

## How to run

Open each notebook in Jupyter and run all cells from top to bottom in order.

### Before Part I — run MMseqs2 search

Search your pathogen proteome against the Atlas AFDB90v4 database:

```bash
mmseqs easy-search your_proteome.fasta AFDBv4_90.fasta output.m8 tmp
```

This produces your `.m8` file — the MMseqs2 alignment results that Part I reads as input. For each query protein only the best match is retained. Output format is standard BLAST-6:

```
queryID  targetID  fident  alnlen  mismatch  gapopen  qstart  qend  tstart  tend  evalue  bitscore
```

> **QueryID format requirement:** QueryIDs must follow the format `STRAINNAME_proteinID` with at least one underscore — for example `PAO1_PA0001` or `KPNIH1_KPN_00001`. QueryIDs without an underscore cannot be attributed to a strain and are automatically excluded from strain coverage calculation with a printed warning.

### Part I

Update the file path at the top of the notebook:

```python
# Path to your MMseqs2 search output
atlas_search_results_info = ['./output.m8']
```

**Memory note:** at least 32 GB RAM recommended — the Atlas taxonomy files are ~4 GB each.

### Part II

Update the configuration block at the top of the notebook:

```python
# ── USER CONFIGURATION ──────────────────────────────────────────────────────
# Atlas data files -- paths to downloaded Atlas taxonomy files
atlas_datafiles = [
    'AFDB90v4_cc_data_uniprot_community_taxonomy_map_with_brightness.csv',
    'AFDB90v4_dust_uniprot_community_taxonomy_map_with_brightness.csv'
]

# Part I output file
mapped_target_dataset = pd.read_csv('eclipse_search_results_component_dark.csv')

# Target species and genus -- change for other ESKAPE pathogens
# Must match exactly how the species/genus is written in the Atlas taxonomy files
TARGET_SPECIES = 'Pseudomonas aeruginosa'   # e.g. 'Klebsiella pneumoniae'
TARGET_GENUS   = 'Pseudomonas'              # e.g. 'Klebsiella'
# ────────────────────────────────────────────────────────────────────────────
```

### Part III — DPPS Scoring

Update the configuration block at the top of the notebook:

```python
# ── USER CONFIGURATION ──────────────────────────────────────────────────────
# Input files (Part II outputs)
PS_CSV = './mapped_target_pathogen_specific_dark_components.csv'
ES_CSV = './mapped_target_eskape_enriched_dark_components.csv'

# Total number of strains in your input dataset
N_STRAINS   = 635     # 635 for P. aeruginosa -- change for other pathogens

# Minimum median component sequence length for structural characterisation
MIN_SEQ_LEN = 300     # amino acids
# ────────────────────────────────────────────────────────────────────────────
```

No Atlas files are needed for Part III — all mapping and proportion calculation was done in Parts I and II.

---

## Input and output files

### Part I

| File | Direction | Description |
|------|-----------|-------------|
| `output.m8` | Input | MMseqs2 search results — your proteome vs AFDB90v4 |
| Atlas CSV files | Input | Community and component brightness and taxonomy |
| `AFDB90v4_subgraphs_summary.csv` | Input | Component-level summary statistics |
| `eclipse.csv` | **Output** | All query proteins with community IDs, component IDs, brightness values, and ESKAPE diversity metrics |

### Part II

| File | Direction | Description |
|------|-----------|-------------|
| `eclipse.csv` | Input | Part I output |
| Atlas CSV files | Input | For target species proportion calculation |
| `mapped_target_pathogen_specific_dark_components.csv` | **Output** | Track A -- pathogen-specific dark components |
| `mapped_target_eskape_enriched_dark_components.csv` | **Output** | Track B -- ESKAPE-enriched dark components |

### Part III — DPPS Scoring

| File | Direction | Description |
|------|-----------|-------------|
| Track A and Track B CSVs | Input | Part II outputs |
| `dpps_components_track_a.csv` | **Output** | Track A components with DPPS sub-scores and tier |
| `dpps_components_track_b.csv` | **Output** | Track B components with DPPS sub-scores and tier |
| `dpps_representatives_track_a.csv` | **Output** | Track A representative sequences with scores |
| `dpps_representatives_track_b.csv` | **Output** | Track B representative sequences with scores |
| `dpps_tier1_track_a.csv` | **Output** | Tier I candidates Track A |
| `dpps_tier1_track_b.csv` | **Output** | Tier I candidates Track B |
| Figures (PDF + PNG) | **Output** | DPPS distributions, scatter plots, heatmaps, sensitivity analysis |

> **Note on file naming:** Output files in the notebooks follow *P. aeruginosa* naming conventions as used in the associated publication (e.g. `mapped_p_aer_dataset_pseudomonas_specific_dark_components.csv`). When applying ECLIPSE to other ESKAPE pathogens, rename output files accordingly to reflect your target organism.

---

## DPPS scoring system

Each dark connected component receives a composite DPPS calculated as a weighted sum of sub-scores. Weights differ between Track A (pathogen-specific) and Track B (ESKAPE-enriched) to reflect the distinct biological properties of each candidate set.

| Sub-score | Definition | Track A weight | Track B weight |
|-----------|-----------|----------------|----------------|
| S1 darkness | 1 - (brightness / 100) | 0.15 | 0.15 |
| S2 target species proportion | Fraction of Atlas members annotated as target species | 0.40 | 0.25 |
| S3 AMR-clade specificity | 1 - ESKAPE_relative_evenness | 0.25 | 0.20 |
| S4 strain coverage | Unique strains carrying component / total strains | 0.20 | 0.15 |
| S5 ESKAPE enrichment | ESKAPE_proportion x (1 - ESKAPE_genus_evenness) | -- | 0.25 |

All sub-scores are normalised to [0, 1] prior to weighting. Weights sum to 1.0 within each track.

### Priority tiers

| Tier | DPPS | Description |
|------|------|-------------|
| Tier I | >= 0.75 | Highest priority -- experimental follow-up recommended |
| Tier II | 0.50-0.75 | High priority |
| Tier III | 0.25-0.50 | Moderate priority |
| Tier IV | < 0.25 | Background |

Tier I stability scores are computed by Monte Carlo weight sensitivity analysis (500 Dirichlet-sampled weight vectors). Components with stability >= 0.80 are considered robustly prioritised independent of weight choice.

---

## Applying ECLIPSE to another ESKAPE pathogen

ECLIPSE is designed to be fully pathogen-agnostic. To apply it to any ESKAPE pathogen follow these steps.

**Step 1** — Download protein FASTA files for your target pathogen from [PATRIC](https://www.patricbrc.org), retaining only strains with complete genome annotation.

**Step 2** — Run MMseqs2 easy-search to generate your `.m8` file:

```bash
mmseqs easy-search your_proteome.fasta AFDBv4_90.fasta output.m8 tmp
```

**Step 3** — Update the input file path in Part I to point to your `.m8` file.

**Step 4** — In Part II update the two target configuration variables to match your organism. To check the exact species and genus strings used in the Atlas for your organism run this after loading the Atlas data:

```python
# Check exact species strings in Atlas for your organism
atlas_data[atlas_data['species'].str.contains('Klebsiella')]['species'].unique()

# Then set:
TARGET_SPECIES = 'Klebsiella pneumoniae'   # exact string from above
TARGET_GENUS   = 'Klebsiella'
```

**Step 5** — In the DPPS notebook update the strain count:

```python
N_STRAINS = [your total strain count]   # e.g. 1200 for K. pneumoniae
```

**Step 6** — Run Parts I, II, and III in order. No other changes to the code are required.

The DPPS weights and tier thresholds were designed to be transferable across ESKAPE organisms without pathogen-specific calibration.

---

## Citation

If you use ECLIPSE in your research please cite:

> Lata S. & Heinz D.W. (2025). *ECLIPSE: Exploring the dark proteome of ESKAPE pathogens through the sequence similarity network of the Protein Universe Atlas.* Bioinformatics, Oxford University Press. [DOI to be added upon acceptance]

Please also cite the Protein Universe Atlas:

> Durairaj J. et al. (2023). *Exploring the functional universe of protein families.* Nature, 626, 582-589. https://doi.org/10.1038/s41586-023-06716-0

And MMseqs2:

> Steinegger M. & Söding J. (2017). *MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.* Nature Biotechnology, 35, 1026-1028. https://doi.org/10.1038/nbt.3988

---

## Contact

**Surabhi Lata**
Department of Molecular Structural Biology
Helmholtz Centre for Infection Research (HZI), Braunschweig, Germany
surabhi.lata@helmholtz-hzi.de / surabhilata94@gmail.com

**Prof. Dirk W. Heinz** (corresponding author)
Department of Molecular Structural Biology
Helmholtz Centre for Infection Research (HZI), Braunschweig, Germany

For questions about the pipeline please open a GitHub Issue.

---

## Licence

MIT Licence — see [LICENSE](LICENSE) for details.
