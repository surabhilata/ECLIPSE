# ECLIPSE — Dark Proteome Exploration of ESKAPE Pathogens

**ECLIPSE** (ESKAPE Connectome Linkage and Inference for Proteome Sequence Exploration) is a modular, pathogen-agnostic computational pipeline for the systematic identification and prioritisation of functionally uncharacterised ("dark") protein families in bacterial panproteomes.

ECLIPSE maps target-pathogen proteomes onto the global sequence similarity network of the [Protein Universe Atlas](https://uniprot3d.org/) (AFDB90v4, UniRef v.2022_03), identifies protein families that are completely dark at the connected component level, characterises their taxonomic diversity across ESKAPE pathogens, and ranks them using the Dark Proteome Prioritisation Score (DPPS).


---

## Table of Contents

- [Overview](#overview)
- [Pipeline Architecture](#pipeline-architecture)
- [Prerequisites](#prerequisites)
- [Required Input Files](#required-input-files)
- [Notebook Guide](#notebook-guide)
  - [Part I — Atlas Mapping & Diversity Analysis](#part-i--atlas-mapping--diversity-analysis-eclipse_partiiipynb)
  - [Part II — Genus-specific & ESKAPE-enriched Component Extraction](#part-ii--genus-specific--eskape-enriched-component-extraction-eclipse_partiiiipynb)
  - [Part III — Clustering & DPPS Prioritization](#part-iii--clustering--dpps-prioritization-eclipse_dpps_scoreiipynb)
- [Adapting to Your Pathogen](#adapting-to-your-pathogen)
- [Output Files](#output-files)
- [Key Concepts](#key-concepts)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

---

## Overview

The ECLIPSE pipeline answers one core question:

> **Which proteins in my pathogen of interest are (1) functionally uncharacterized ("dark"), (2) evolutionarily restricted to AMR-relevant clades, and (3) conserved across many strains — making them promising unexplored drug targets?**

It does this in three sequential stages:

| Stage | Notebook | Purpose |
|-------|----------|---------|
| **I** | `ECLIPSE_PartI.ipynb` | Map proteome to Uniprot 3d Atlas → compute community/component darkness & taxonomic diversity |
| **II** | `ECLIPSE_PartII.ipynb` | Extract pathogen-specific and ESKAPE-enriched dark components |
| **III** | `ECLIPSE_DPPS_score.ipynb` | Cluster proteins → score components with DPPS → rank & visualize candidates |

---

## Pipeline Architecture

```
Your proteome FASTA files (one or many strains)
            │
            ▼
   MMseqs2 easy-search
   (vs. AFDBv4_90.fasta)
            │
            ▼
         PA.m8 (alignment results)
            │
            ▼
    ECLIPSE_PartI.ipynb
    ┌─────────────────────────────────────────────┐
    │ • Load Atlas community & component data      │
    │ • Map each protein → UniRef50 → community   │
    │ • Assign darkness (brightness) values        │
    │ • Compute ESKAPE taxonomic diversity metrics │
    │ • Save eclipse.csv / eclipse_seq.csv         │
    └─────────────────────────────────────────────┘
            │
            ▼
    ECLIPSE_PartII.ipynb
    ┌─────────────────────────────────────────────┐
    │ • Filter dark components (brightness ≤ 5%)  │
    │ • Track A: pathogen-genus-specific comps    │
    │ • Track B: ESKAPE-enriched components       │
    │ • Add species-level proportions             │
    │ • Save two track CSVs                       │
    └─────────────────────────────────────────────┘
            │
            ▼
    ECLIPSE_DPPS_score.ipynb
    ┌─────────────────────────────────────────────┐
    │ • Length filter (≥ 300 aa by median)        │
    │ • Export FASTA per track                    │
    │ • MMseqs2 easy-cluster (redundancy removal) │
    │ • Build component-level tables              │
    │ • Compute DPPS sub-scores (S1–S5)           │
    │ • Monte Carlo weight sensitivity (n=500)    │
    │ • Tier ranking (I–IV) + visualizations      │
    └─────────────────────────────────────────────┘
            │
            ▼
    Ranked candidate drug target lists
    (dpps_tier1_*.csv, dpps_all_representatives_combined.csv)
```

---

## Prerequisites

### Software

| Tool | Version | Purpose | Install |
|------|---------|---------|---------|
| Python | ≥ 3.8 | Core environment | [python.org](https://python.org) |
| Jupyter | any | Run notebooks | `pip install jupyter` |
| MMseqs2 | latest | Sequence search & clustering | [MMseqs2 GitHub](https://github.com/soedinglab/MMseqs2) |

### Python Packages

```bash
pip install pandas numpy matplotlib seaborn biopython tqdm pyyaml
```

Or using conda/mamba:

```bash
mamba create -n eclipse python=3.10
mamba activate eclipse
mamba install pip jupyter
pip install pandas numpy matplotlib seaborn biopython tqdm pyyaml
```

### MMseqs2

MMseqs2 must be accessible from the command line (`mmseqs` in PATH). Verify with:

```bash
mmseqs --version
```

---

## Required Input Files

You need to download the **ESM Metagenomic Atlas** files. These are large; download them once and reuse across analyses.

| File | Description | Where to Get |
|------|-------------|--------------|
| `AFDBv4_90.fasta` | Atlas UniRef50 representative sequences | [ESM Atlas](https://esmatlas.com/resources) |
| `AFDB90v4_cc_data_uniprot_community_taxonomy_map_with_brightness.csv` | Community-level Atlas data (connected components) | [ESM Atlas](https://esmatlas.com/resources) |
| `AFDB90v4_dust_uniprot_community_taxonomy_map_with_brightness.csv` | Community-level Atlas data (DUST regions) | [ESM Atlas](https://esmatlas.com/resources) |
| `AFDB90v4_subgraphs_summary.csv` | Component-level summary statistics | [ESM Atlas](https://esmatlas.com/resources) |

Your own input files:

| File | Description | Format |
|------|-------------|--------|
| `<organism>_faa.tar.gz` or individual `.faa` files | Protein FASTA files for your pathogen strains | FASTA (`.faa`) |

> **Naming convention for strain tracking:** Name each protein sequence as `STRAINNAME_PROTEINID` (e.g., `PAO1_PA0001`). The pipeline extracts the strain name from the prefix before the first underscore to compute strain coverage (S4 score). If your IDs do not follow this convention, strain coverage calculation will not work correctly — see [Adapting to Your Pathogen](#adapting-to-your-pathogen).

---

## Notebook Guide

### Part I — Atlas Mapping & Diversity Analysis (`ECLIPSE_PartI.ipynb`)

**What it does:** Maps every protein in your dataset to a community and component in the ESM Atlas, then computes darkness and taxonomic diversity metrics.

#### Step 0 — Run MMseqs2 (outside the notebook, in terminal)

Before opening the notebook, run the sequence search against the Atlas:

```bash
# If your FASTAs are in a tar.gz archive:
mmseqs easy-search <organism>_faa.tar.gz AFDBv4_90.fasta <organism>.m8 tmp --max-seqs 1

# If your FASTAs are in a directory of individual .faa files, merge them first:
cat ./faa/*.faa > all_proteins.faa
mmseqs easy-search all_proteins.faa AFDBv4_90.fasta <organism>.m8 tmp --max-seqs 1
```

> `--max-seqs 1` keeps only the single best Atlas match per query protein.

#### Step 1 — Configure paths (Cell 3)

```python
# Cell 3 — change the filename to match your MMseqs2 output
atlas_search_results_info = ['./<organism>.m8']
```

#### Step 2 — Load Atlas data (Cell 5)

Ensure the Atlas CSV files are in the same directory as the notebook (or update paths):

```python
atlas_datafiles = [
    'AFDB90v4_cc_data_uniprot_community_taxonomy_map_with_brightness.csv',
    'AFDB90v4_dust_uniprot_community_taxonomy_map_with_brightness.csv'
]
```

#### Step 3 — Update the ESKAPE genus list if needed (Cell 24–25)

By default the pipeline uses the full ESKAPE list:

```python
# Cell 25 — ESKAPE genera used for taxonomic labeling
ESKAPE_genera = [
    'Acinetobacter', 'Enterococcus', 'Escherichia',
    'Klebsiella', 'Pseudomonas', 'Staphylococcus', 'Enterobacter'
]
```

You can add or remove genera. For example, to include *Burkholderia*:

```python
ESKAPE_genera = [
    'Acinetobacter', 'Enterococcus', 'Escherichia',
    'Klebsiella', 'Pseudomonas', 'Staphylococcus', 'Enterobacter', 'Burkholderia'
]
```

#### Step 4 — Add sequences to results (Cells 39–44)

Update the directory path to point to your FASTA files:

```python
# Cell 40
directory_path = './faa'   # ← change to your FASTA directory
```

#### Outputs from Part I

| File | Contents |
|------|----------|
| `eclipse.csv` | All query proteins with Atlas mapping & diversity metrics |
| `eclipse_seq.csv` | Same as above with protein sequences added |
| `dark_communities.png` | Histogram of functional brightness at community level |
| `components_darkness.png` | Histogram of functional brightness at component level |
| `taxonomic_diversity_of_All_components.png` | Scatter plot of diversity metrics for all components |
| `taxonomic_diversity_of_dark_components.png` | Same, restricted to dark components |

---

### Part II — Genus-specific & ESKAPE-enriched Component Extraction (`ECLIPSE_PartII.ipynb`)

**What it does:** Takes the dark-component subset from Part I and separates it into two biologically meaningful tracks for downstream prioritization.

#### Step 1 — Load Atlas data (Cell 2)

Same Atlas CSV files as Part I. Update paths if needed.

#### Step 2 — Load the dark component dataset (Cell 3)

```python
# Cell 3 — point to your Part I output
mapped_dataset = pd.read_csv('eclipse_search_results_component_dark.csv', index_col='Unnamed: 0')
```

> **How to generate this file:** At the end of Part I, filter `eclipse_seq.csv` for dark components (e.g., `component_brightness <= 5`) and save it. This filtered file is the input to Part II.

#### Step 3 — Adapt genus-specific filtering (Cells 6–7)

**Track A — pathogen-genus-specific components:**

The default filters for *Pseudomonas*. Change the genus name for your organism:

```python
# Cell 7 — replace 'Pseudomonas' with your genus of interest
filtered_components = selected_atlas_data.groupby("componentIDs").filter(
    lambda x: (x["genus"] == "Pseudomonas").all()   # ← change here
)
```

For *Staphylococcus aureus*, use `"Staphylococcus"`. For *Klebsiella pneumoniae*, use `"Klebsiella"`, etc.

**Track B — ESKAPE-enriched components:**

```python
# Cell 12 — adjust ESKAPE_proportion threshold if desired (default 0.5)
mapped_dataset_eskape_enriched = mapped_dataset.loc[
    (mapped_dataset['ESKAPE_proportion'] >= 0.5) &
    (~mapped_dataset['componentID'].isin(genus_specific_components['componentID']))
]
```

#### Outputs from Part II

| File | Contents |
|------|----------|
| `mapped_<genus>_dataset_<genus>_specific_dark_components.csv` | Track A: genus-specific dark components |
| `mapped_<genus>_dataset_eskape_enriched_dark_components.csv` | Track B: ESKAPE-enriched dark components |

---

### Part III — Clustering & DPPS Prioritization (`ECLIPSE_DPPS_score.ipynb`)

**What it does:** Clusters redundant sequences, builds component-level scoring tables, computes the DPPS composite score, performs Monte Carlo weight sensitivity analysis, and produces ranked output tables and publication-ready figures.

#### Step 1 — Configure input files and parameters (Cell 3)

```python
# ── Input files (Part II outputs) ──────────────────────────────────────────
PS_CSV = './mapped_<genus>_dataset_<genus>_specific_dark_components.csv'   # Track A
ES_CSV = './mapped_<genus>_dataset_eskape_enriched_dark_components.csv'    # Track B

# ── Key parameter: total number of strains in YOUR dataset ─────────────────
N_STRAINS = 635   # ← MUST change this to match your dataset

# ── Sequence length filter ─────────────────────────────────────────────────
MIN_SEQ_LEN = 300   # minimum median component sequence length in amino acids

# ── MMseqs2 clustering parameters ─────────────────────────────────────────
MMSEQS_MIN_SEQ_ID = 0.3   # sequence identity threshold for clustering
MMSEQS_COVERAGE   = 0.8   # coverage threshold for clustering
```

> **`N_STRAINS` is the most critical parameter to update.** It determines the strain coverage sub-score (S4). Set it to the total number of genomes/strains in your input FASTA dataset.

#### Step 2 — Adjust DPPS weights if needed (Cell 3)

The default weights are tuned for *P. aeruginosa*. You can adjust them per track:

**Track A weights (pathogen-genus-specific):**

```python
WEIGHTS_PS = {
    'S1_darkness':           0.15,   # flat 1.0 for all dark components
    'S2b_pa_combined':       0.40,   # max(Atlas species proportion, strain fraction)
    'S3_specificity':        0.25,   # 1 - ESKAPE_relative_evenness
    'S4_pa_strain_coverage': 0.20,   # unique strains / N_STRAINS
}
# Weights must sum to 1.0
```

**Track B weights (ESKAPE-enriched):**

```python
WEIGHTS_ES = {
    'S1_darkness':           0.15,
    'S2_pa_proportion':      0.25,   # Atlas-derived species proportion
    'S3_specificity':        0.20,
    'S4_pa_strain_coverage': 0.15,
    'S5_eskape_enrich':      0.25,   # ESKAPE_proportion × (1 - ESKAPE_genus_evenness)
}
```

> If your pathogen has reliable species-level annotation in UniProt (unlike *P. aeruginosa*, which suffers from widespread "Pseudomonas sp." entries), you can replace `S2b_pa_combined` with the simpler `S2_pa_proportion` in Track A as well, and redistribute weights accordingly.

#### Step 3 — Rename species-specific column references

Throughout the notebook, `p_aeruginosa_proportion` refers to the fraction of Atlas members annotated as the target species. When you run Part II for a different organism, rename this column for clarity (or simply update the references in Part III):

```python
# In Part II, when computing species proportion, rename the output column:
components_with_prob.rename(columns={'p_aeruginosa_proportion': 's_aureus_proportion'}, inplace=True)

# Then update references in Part III accordingly:
comp['S2_pa_proportion'] = comp['s_aureus_proportion'].clip(0, 1)
```

#### DPPS Sub-score Definitions

| Sub-score | Formula | Track A weight | Track B weight |
|-----------|---------|----------------|----------------|
| S1 — darkness | `1 − (component_brightness / 100)` | 0.15 | 0.15 |
| S2b — combined species evidence *(Track A)* | `max(Atlas species proportion, strain fraction)` | 0.40 | — |
| S2 — species proportion *(Track B)* | Atlas fraction annotated as target species | — | 0.25 |
| S3 — taxonomic specificity | `1 − ESKAPE_relative_evenness` | 0.25 | 0.20 |
| S4 — strain coverage | `unique strains in component / N_STRAINS` | 0.20 | 0.15 |
| S5 — ESKAPE enrichment *(Track B)* | `ESKAPE_proportion × (1 − ESKAPE_genus_evenness)` | — | 0.25 |

#### Tier Thresholds

| Tier | DPPS range | Interpretation |
|------|-----------|----------------|
| I | ≥ 0.75 | High-priority candidates |
| II | 0.50 – 0.75 | Moderate priority |
| III | 0.25 – 0.50 | Low priority |
| IV | < 0.25 | Deprioritized |

#### Outputs from Part III

| File | Contents |
|------|----------|
| `dpps_components_pseudomonas_specific.csv` | Track A component-level scores |
| `dpps_components_eskape_enriched.csv` | Track B component-level scores |
| `dpps_representatives_pseudomonas_specific.csv` | Track A representative proteins with scores |
| `dpps_representatives_eskape_enriched.csv` | Track B representative proteins with scores |
| `dpps_tier1_pseudomonas_specific.csv` | Tier I candidates — Track A |
| `dpps_tier1_eskape_enriched.csv` | Tier I candidates — Track B |
| `dpps_all_representatives_combined.csv` | Both tracks merged, sorted by DPPS |
| `dpps_distribution_both_tracks.pdf` | DPPS histogram by tier |
| `dpps_pa_proportion_scatter.pdf` | Species proportion vs DPPS scatter |
| `dpps_heatmap_pseudomonas_specific.pdf` | Sub-score heatmap — top 30 Tier I (Track A) |
| `dpps_heatmap_eskape_enriched.pdf` | Sub-score heatmap — top 30 Tier I (Track B) |
| `pa_strain_coverage_distribution.pdf` | Strain fraction distribution |

---

## Adapting to Your Pathogen

Here is a concise checklist for running ECLIPSE on a new organism:

### 1. Prepare your FASTA files

- Collect annotated proteome FASTAs for as many strains as possible.
- Name each sequence: `STRAINNAME_PROTEINID` (e.g., `NCTC8325_SA0001`).
- Store FASTAs in a `./faa/` directory or package into `<organism>_faa.tar.gz`.

### 2. Run MMseqs2 (terminal)

```bash
mmseqs easy-search <organism>_faa.tar.gz AFDBv4_90.fasta <organism>.m8 tmp --max-seqs 1
```

### 3. Part I — update 3 things

| Cell | What to change | Example (S. aureus) |
|------|---------------|---------------------|
| Cell 3 | `atlas_search_results_info` | `['./SA.m8']` |
| Cell 25 | `genus_modifed` genus list (AMR genera) | keep default or extend |
| Cell 40 | `directory_path` | `'./faa'` |

### 4. Part II — update 3 things

| Cell | What to change | Example (S. aureus) |
|------|---------------|---------------------|
| Cell 3 | Input CSV filename | `'eclipse_search_results_component_dark.csv'` |
| Cell 5 | Species proportion — change `'Pseudomonas aeruginosa'` | `'Staphylococcus aureus'` |
| Cell 7 | Genus filter for Track A | `(x["genus"] == "Staphylococcus").all()` |

### 5. Part III — update 4 things

| Location | What to change | Example (S. aureus, 500 strains) |
|----------|---------------|----------------------------------|
| Cell 3 `PS_CSV` / `ES_CSV` | Input filenames | Match Part II output names |
| Cell 3 `N_STRAINS` | Total strain count in your dataset | `500` |
| Throughout | `p_aeruginosa_proportion` column name | `s_aureus_proportion` |
| Cell 3 weights | Adjust if annotation quality differs | Consider reducing S2b weight if species annotation is reliable |

---

## Output Files

### Key columns in `eclipse.csv` / `eclipse_seq.csv`

| Column | Description |
|--------|-------------|
| `queryID` | Protein identifier from your input FASTA |
| `targetID` | Best-matching UniRef50 representative in the Atlas |
| `fident` | Fraction of identical residues in the alignment |
| `communityID` | Atlas community to which the match belongs |
| `brightness` | Functional brightness of the community (0 = fully dark) |
| `componentID` | Atlas component (subgraph) to which the match belongs |
| `component_brightness` | Median brightness of the component |
| `ESKAPE_relative_evenness` | Shannon evenness across all genera (AMR genera merged); 0 = ESKAPE-exclusive |
| `ESKAPE_genus_evenness` | Shannon evenness within AMR genera only; 0 = single-genus |
| `ESKAPE_proportion` | Fraction of component proteins from AMR genera |

### Key columns in DPPS output files

| Column | Description |
|--------|-------------|
| `componentID` | Atlas component identifier |
| `DPPS` | Composite prioritization score (0–1, higher = better candidate) |
| `tier` | Priority tier (I = highest) |
| `S1_darkness` | Darkness sub-score |
| `S2b_pa_combined` / `S2_pa_proportion` | Species evidence sub-score |
| `S3_specificity` | Taxonomic specificity sub-score |
| `S4_pa_strain_coverage` | Strain coverage sub-score |
| `S5_eskape_enrich` | ESKAPE enrichment sub-score (Track B only) |
| `PA_strain_count` | Number of strains carrying this component |
| `PA_strain_fraction` | Fraction of total strains carrying this component |
| `tier1_stability` | Fraction of Monte Carlo replicates (n=500) where component reaches Tier I |
| `cluster_size` | Number of sequences in the MMseqs2 cluster represented by this protein |

---

## Key Concepts

**Functional brightness / darkness:** A metric from the ESM Metagenomic Atlas quantifying the fraction of proteins in a community that have functional annotations. A brightness of 0% means no protein in the community has a known function — a "dark" community. Components with median brightness ≤ 5% are considered dark in this pipeline.

**Community vs. component:** In the Atlas, proteins are grouped into *communities* based on structural similarity. Communities that are connected in the global similarity network form larger *components* (subgraphs). Darkness is assessed at both levels.

**ESKAPE_relative_evenness:** A normalized Shannon diversity index computed over all genera in a component, where all ESKAPE genera are merged into a single "AMR genus" label. A value of 0 means the component is exclusively occupied by ESKAPE proteins. A higher value indicates proteins from non-ESKAPE organisms are also present.

**ESKAPE_genus_evenness:** Same index but computed only within ESKAPE organisms. A value of 0 means all ESKAPE proteins in the component come from a single genus (maximally genus-specific). Higher values indicate multi-genus ESKAPE components.

**S2b (combined PA evidence):** Introduced for Track A because species names like "Pseudomonas sp." in UniProt suppress the Atlas-derived species proportion despite genuine conservation. S2b takes the maximum of the Atlas proportion and the query strain fraction, correcting for this annotation gap.

**Monte Carlo sensitivity analysis:** Weights are randomly sampled from a Dirichlet distribution (n=500 replicates) and DPPS is recomputed each time. `tier1_stability` reports how often a component reaches Tier I across all replicates — a value ≥ 0.8 indicates the ranking is robust to weight uncertainty.

---

## Troubleshooting

**`KeyError: 'LongLink'` in Part I Cell 43**
Some MMseqs2 versions truncate very long sequence names to "LongLink". Identify and drop these rows:
```python
atlas_search_results.drop(
    index=atlas_search_results.index[atlas_search_results['queryID'] == 'LongLink'],
    inplace=True
)
```

**Warning: `N_STRAINS` mismatch**
If `PA_strain_fraction` values exceed 1.0, your `N_STRAINS` is set too low. Count your strains:
```python
all_strains = set(df['queryID'].str.split('_').str[0])
print(f"Detected {len(all_strains)} unique strain prefixes")
```

**MMseqs2 clustering fails**
Ensure `mmseqs` is in PATH and the FASTA files were written correctly. Check:
```bash
wc -l track_ps_pseudomonas_specific.fasta  # should be 2× number of sequences
mmseqs easy-cluster track_ps_pseudomonas_specific.fasta cluster_ps tmp --min-seq-id 0.3 -c 0.8
```

**`values[0]` IndexError in Part I (mapping loops)**
This occurs when a `targetID` in the search results is not present in the Atlas data. This can happen if the Atlas CSV files are incomplete. Verify file integrity and that both `cc` and `dust` Atlas files were loaded.

**No Tier I components found**
Lower your `MMSEQS_MIN_SEQ_ID` to include more diverse sequences, or adjust `TIER_BINS` to broaden the Tier I threshold. Alternatively, re-examine whether your `N_STRAINS` value is correct — an incorrect count will suppress S4 scores.

---

## Citation

If you use ECLIPSE in your research, please cite:

> [Your citation here]

---

## License

[Your license here]

