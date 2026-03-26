# ECLIPSE — Dark Proteome Exploration of ESKAPE Pathogens

**ECLIPSE** (ESKAPE Connectome Linkage and Inference for Proteome Sequence Exploration) is a modular computational pipeline for the systematic identification and prioritisation of functionally uncharacterised ("dark") protein families in bacterial panproteomes.

ECLIPSE maps target-pathogen proteomes onto the global sequence similarity network of the (https://uniprot3d.org/) (AFDB90v4, UniRef v.2022_03), identifies protein families that are completely dark at the connected component level, characterises their taxonomic diversity across ESKAPE pathogens, and ranks them using the Dark Proteome Prioritisation Score (DPPS).

---

## Pipeline overview

ECLIPSE runs as three sequential Jupyter notebooks. Each notebook must be completed before starting the next.

```
Your pathogen FASTA files
        │
        ▼  MMseqs2 easy-search (run before Part I)
        │
┌───────────────────────────────────────────┐
│  ECLIPSE_PartI.ipynb                      │
│                                           │
│  • Load MMseqs2 search results (PA.m8)    │
│  • Map proteins → Atlas communities       │
│  • Map communities → connected components │
│  • Classify dark communities (0%)         │
│  • Classify dark components (0%)          │
│  • Compute ESKAPE taxonomic diversity:    │
│    - ESKAPE_proportion                    │
│    - ESKAPE_genus_evenness                │
│    - ESKAPE_relative_evenness             │
│                                           │
│  Output → eclipse.csv                     │
└──────────────────┬────────────────────────┘
                   │
                   ▼
┌───────────────────────────────────────────┐
│  ECLIPSE_PartII.ipynb                     │
│                                           │
│  • Load eclipse.csv + Atlas taxonomy      │
│  • Add P. aeruginosa proportion per       │
│    component from Atlas taxonomy files    │
│  • Track A — Pseudomonas-specific:        │
│    ESKAPE_proportion == 1.0               │
│    all Atlas genus == Pseudomonas         │
│    → 83 components                        │
│  • Track B — ESKAPE-enriched:             │
│    ESKAPE_proportion >= 0.5               │
│    excluding Track A                      │
│    → 215 components                       │
│                                           │
│  Output → two track CSV files             │
└──────────────────┬────────────────────────┘
                   │
                   ▼
┌───────────────────────────────────────────┐
│  ECLIPSE_DPPS_score.ipynb                 │
│                                           │
│  • Median length filter (>= 300 aa)       │
│  • MMseqs2 easy-cluster (30% id, 80% cov) │
│  • One representative per cluster         │
│  • DPPS scoring:                          │
│    S1 darkness · S2 PA proportion         │
│    S3 specificity · S4 PA strain coverage │
│    S5 ESKAPE enrichment (Track B only)    │
│  • Tier I–IV assignment                   │
│  • Monte Carlo weight sensitivity (n=500) │
│  • Figures and summary report             │
│                                           │
│  Output → scored CSVs + figures           │
└───────────────────────────────────────────┘
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

Download from the (https://doi.org/10.5281/zenodo.19119408)
Note: These files are produced by Durairaj et al. (2023) and are available from the paper (https://www.nature.com/articles/s41586-023-06622-3)

---

## How to run

Open each notebook in Jupyter and run all cells from top to bottom in order.

### Before Part I — run MMseqs2 search

Search your pathogen proteome against the Atlas AFDB90v4 database:

```bash
mmseqs easy-search PA_faa.tar.gz AFDBv4_90.fasta PA.m8 tmp
```

This produces `PA.m8` — the MMseqs2 alignment results file that Part I reads as input. For each query protein only the best match is retained. Output format is standard BLAST-6:

```
queryID  targetID  fident  alnlen  mismatch  gapopen  qstart  qend  tstart  tend  evalue  bitscore
```

### Part I

Update the file path at the top of the notebook:

```python
atlas_search_results_info = ['./PA.m8']   # path to your MMseqs2 output
```

**Memory note:** at least 32 GB RAM recommended — the Atlas taxonomy files are ~4 GB each.

### Part II

Update the file paths at the top of the notebook:

```python
# NOTE: update these paths for your system
atlas_datafiles = [
    'AFDB90v4_cc_data_uniprot_community_taxonomy_map_with_brightness.csv',
    'AFDB90v4_dust_uniprot_community_taxonomy_map_with_brightness.csv'
]
mapped_p_aer_dataset = pd.read_csv('eclipse_search_results_component_dark.csv')
```

### Part III — DPPS Scoring

Update the configuration block at the top of the notebook:

```python
PS_CSV       = './mapped_p_aer_dataset_pseudomonas_specific_dark_components.csv'
ES_CSV       = './mapped_p_aer_dataset_eskape_enriched_dark_components.csv'
N_PA_STRAINS = 635    # update to match your total strain count
MIN_SEQ_LEN  = 300    # minimum median component sequence length in aa
```

No Atlas files are needed for this notebook — all mapping was done in Parts I and II.

---

## Input and output files

### Part I

| File | Direction | Description |
|------|-----------|-------------|
| `PA.m8` | Input | MMseqs2 search results |
| Atlas CSV files | Input | Community and component brightness and taxonomy |
| `AFDB90v4_subgraphs_summary.csv` | Input | Component-level summary statistics |
| `eclipse.csv` | **Output** | All query proteins with community IDs, component IDs, brightness values, and ESKAPE diversity metrics |

### Part II

| File | Direction | Description |
|------|-----------|-------------|
| `eclipse.csv` | Input | Part I output |
| Atlas CSV files | Input | For P. aeruginosa proportion calculation |
| `mapped_p_aer_dataset_pseudomonas_specific_dark_components.csv` | **Output** | Track A — Pseudomonas-specific dark components |
| `mapped_p_aer_dataset_eskape_enriched_dark_components.csv` | **Output** | Track B — ESKAPE-enriched dark components |

### Part III — DPPS Scoring

| File | Direction | Description |
|------|-----------|-------------|
| Track A and Track B CSVs | Input | Part II outputs |
| `dpps_components_pseudomonas_specific.csv` | **Output** | Track A components with DPPS sub-scores and tier |
| `dpps_components_eskape_enriched.csv` | **Output** | Track B components with DPPS sub-scores and tier |
| `dpps_representatives_pseudomonas_specific.csv` | **Output** | Track A representative sequences with scores |
| `dpps_representatives_eskape_enriched.csv` | **Output** | Track B representative sequences with scores |
| `dpps_tier1_pseudomonas_specific.csv` | **Output** | Tier I candidates Track A |
| `dpps_tier1_eskape_enriched.csv` | **Output** | Tier I candidates Track B |
| Figures (PDF + PNG) | **Output** | DPPS distributions, scatter plots, heatmaps, sensitivity analysis |

---

## DPPS scoring system

Each dark component receives a composite score calculated as a weighted sum of sub-scores.

| Sub-score | Definition | Track A weight | Track B weight |
|-----------|-----------|----------------|----------------|
| S1 darkness | 1 − (brightness / 100) | 0.15 | 0.15 |
| S2 PA proportion | Fraction of Atlas members annotated as *P. aeruginosa* | 0.40 | 0.25 |
| S3 specificity | 1 − ESKAPE_relative_evenness | 0.25 | 0.20 |
| S4 PA strain coverage | Unique PA strains / total strains | 0.20 | 0.15 |
| S5 ESKAPE enrichment | ESKAPE_proportion × (1 − ESKAPE_genus_evenness) | — | 0.25 |

### Priority tiers

| Tier | DPPS | Description |
|------|------|-------------|
| Tier I | ≥ 0.75 | Highest priority — experimental follow-up recommended |
| Tier II | 0.50–0.75 | High priority |
| Tier III | 0.25–0.50 | Moderate priority |
| Tier IV | < 0.25 | Background |

Tier I stability scores are computed by Monte Carlo weight sensitivity analysis (500 Dirichlet-sampled weight vectors). Components with stability ≥ 0.80 are considered robustly prioritised independent of weight choice.

---

## Applying ECLIPSE to another ESKAPE pathogen

1. Download protein FASTA files for your target pathogen from (https://pseudomonas.com/) — retain only strains with complete genome annotation.
2. Run MMseqs2 easy-search to generate your `.m8` file.
3. Update all file paths in each notebook.
4. Update `N_PA_STRAINS` in the DPPS notebook to match your strain count.
5. Run Parts I, II, and III in order.

No other changes to the code are required.

---

## Citation

If you use ECLIPSE please cite:

Please also cite:

> Durairaj J. et al. (2023). *Exploring the functional universe of protein families.* Nature, 626, 582–589. https://doi.org/10.1038/s41586-023-06716-0

> Steinegger M. & Söding J. (2017). *MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets.* Nature Biotechnology, 35, 1026–1028. https://doi.org/10.1038/nbt.3988

---

## Contact

**Surabhi Lata**  
Department of Molecular Structural Biology
Helmholtz Centre for Infection Research (HZI), Braunschweig, Germany
surabhi.lata@helmholtz-hzi.de /
surabhilata94@gmail.com

For questions please open a GitHub Issue.

---

## Licence

MIT Licence — see [LICENSE](LICENSE) for details.
