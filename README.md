# ECLIPSE — Dark Proteome Exploration of ESKAPE Pathogens

**ECLIPSE** (ESKAPE Connectome Linkage and Inference for Proteome Sequence Exploration) is a modular, universal computational pipeline for the systematic identification and prioritisation of functionally uncharacterised ("dark") protein families in bacterial panproteomes.

ECLIPSE maps target-pathogen proteomes onto the global sequence similarity network of the [Protein Universe Atlas](https://uniprot3d.org/) (AFDB90v4, UniRef v.2022_03), identifies protein families that are completely dark at the connected component level, characterises their taxonomic diversity across ESKAPE pathogens, and ranks them using the Dark Proteome Prioritisation Score (DPPS).

> **Note:** This ECLIPSE repository has no overlap with the package "Eclipse: a Python package for alignment of two or more nontargeted LC-MS metabolomics datasets" (https://doi.org/10.1093/bioinformatics/btaf290).

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Architecture](#pipeline-architecture)
- [Prerequisites](#prerequisites)
- [Required Input Files](#required-input-files)
- [Universal Design — One Notebook Set for All ESKAPE Pathogens](#universal-design)
- [Notebook Guide](#notebook-guide)
  - [Part I — Atlas Mapping & Diversity Analysis](#part-i)
  - [Part II — Genus-specific & ESKAPE-enriched Component Extraction](#part-ii)
  - [Part III — Clustering & DPPS Prioritization](#part-iii)
- [Quick-Switch Guide: Running on a New Pathogen](#quick-switch-guide)
- [Strain Identity & S4 Scoring — How It Works](#strain-identity)
- [DPPS Formula Reference](#dpps-formula-reference)
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
| **I** | `ECLIPSE_PartI.ipynb` | Map proteome to Atlas → compute darkness & taxonomic diversity → build strain mapping |
| **II** | `ECLIPSE_PartII.ipynb` | Extract genus-specific and ESKAPE-enriched dark components |
| **III** | `ECLIPSE_DPPS_score.ipynb` | Cluster proteins → score with DPPS → rank & visualize candidates |

---

## Pipeline Architecture

```
Your proteome FASTA files (.faa, one file per strain)
            |
            v
   MMseqs2 easy-search
   (vs. AFDBv4_90.fasta)
            |
            v
      PATHOGEN.m8  (alignment results)
            |
            v
    ECLIPSE_PartI.ipynb
    +---------------------------------------------+
    | * Load Atlas community & component data      |
    | * Map each protein to UniRef50 -> community  |
    | * Assign darkness (brightness) values        |
    | * Compute ESKAPE taxonomic diversity metrics |
    | * Build queryID_to_strain.csv (strain map)   |
    | * Save eclipse.csv / eclipse_seq.csv         |
    +---------------------------------------------+
            |
            v
    ECLIPSE_PartII.ipynb
    +---------------------------------------------+
    | * Load dark component subset from Part I     |
    | * Track A: genus-specific components         |
    | * Track B: ESKAPE-enriched components        |
    | * Add species-level Atlas proportions        |
    | * Save two track CSVs for Part III           |
    +---------------------------------------------+
            |
            v
    ECLIPSE_DPPS_score.ipynb
    +---------------------------------------------+
    | * Length filter (>= 300 aa median)           |
    | * Export FASTA per track                     |
    | * MMseqs2 easy-cluster (redundancy removal)  |
    | * Auto-detect strain counting method         |
    | * Compute DPPS sub-scores (S1-S5)            |
    | * Monte Carlo weight sensitivity (n=500)     |
    | * Tier ranking (I-IV) + visualizations       |
    +---------------------------------------------+
            |
            v
    Ranked candidate drug target lists
    (dpps_tier1_*.csv, dpps_all_representatives_combined.csv)
```

---

## Prerequisites

### Software

| Tool | Version | Purpose | Install |
|------|---------|---------|---------|
| Python | >= 3.8 | Core environment | [python.org](https://python.org) |
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

Verify MMseqs2 is accessible:

```bash
mmseqs --version
```

---

## Required Input Files

Download these files once and reuse across all pathogen analyses. Available at: https://doi.org/10.5281/zenodo.19119408

> These files were produced by the Uniprot 3D Atlas (Durairaj et al., 2023).

| File | Description |
|------|-------------|
| `AFDBv4_90.fasta` | Atlas UniRef50 representative sequences (for MMseqs2 search) |
| `AFDB90v4_cc_data_uniprot_community_taxonomy_map_with_brightness.csv` | Community-level Atlas data (connected components) |
| `AFDB90v4_dust_uniprot_community_taxonomy_map_with_brightness.csv` | Community-level Atlas data (DUST regions) |
| `AFDB90v4_subgraphs_summary.csv` | Component-level summary statistics |

Your own input — one `.faa` file per genome/strain in a `./faa/` directory:

| File | Description |
|------|-------------|
| `./faa/*.faa` | Protein FASTA files, one per strain/genome assembly |

> **Filename convention for strains:** The pipeline extracts strain identity from the `.faa` **filename**, not from the sequence header. Name your files clearly. If downloaded from NCBI, the GCF/GCA accession is extracted automatically (e.g. `GCF_000013425.1_ASM1342v1_protein.faa` → strain `GCF_000013425.1`). Custom names like `PAO1.faa`, `MRSA252.faa` also work directly. See [Strain Identity & S4 Scoring](#strain-identity) for full details.

---

## Universal Design

**All three notebooks are designed to work for any ESKAPE pathogen** — *P. aeruginosa*, *S. aureus*, *K. pneumoniae*, *A. baumannii*, *E. faecium*, *E. cloacae*, or any other organism — **by changing only a small config block at the top of each notebook.** Nothing else needs to be touched.

### What you change per pathogen

| Notebook | Cell | Variables to change |
|----------|------|---------------------|
| Part I | Cell 3 | `M8_FILE` |
| Part II | Cell 3 | `GENUS`, `SPECIES`, `SPECIES_COL`, `INPUT_CSV`, `PS_OUT_CSV`, `ES_OUT_CSV` |
| Part III | Cell 2 | `GENUS`, `SPECIES`, `SPECIES_COL`, `SPECIES_ABBR`, `FAA_DIR`, `PS_CSV`, `ES_CSV` |

### Example config values for common pathogens

| Pathogen | `GENUS` | `SPECIES` | `SPECIES_COL` |
|----------|---------|-----------|---------------|
| *P. aeruginosa* | `"Pseudomonas"` | `"Pseudomonas aeruginosa"` | `"p_aeruginosa_proportion"` |
| *S. aureus* | `"Staphylococcus"` | `"Staphylococcus aureus"` | `"s_aureus_proportion"` |
| *K. pneumoniae* | `"Klebsiella"` | `"Klebsiella pneumoniae"` | `"k_pneumoniae_proportion"` |
| *A. baumannii* | `"Acinetobacter"` | `"Acinetobacter baumannii"` | `"a_baumannii_proportion"` |
| *E. faecium* | `"Enterococcus"` | `"Enterococcus faecium"` | `"e_faecium_proportion"` |

---

## Notebook Guide

---

### Part I — Atlas Mapping & Diversity Analysis {#part-i}

**Notebook:** `ECLIPSE_PartI.ipynb`

**What it does:** Maps every protein in your dataset to a community and component in the Atlas, computes darkness and ESKAPE taxonomic diversity metrics, and builds the universal strain-to-protein mapping file used by Part III.

#### Step 0 — Run MMseqs2 (in terminal, before opening the notebook)

```bash
# Merge all your .faa files into one:
cat ./faa/*.faa > all_proteins.faa

# Run the search against the Atlas:
mmseqs easy-search all_proteins.faa AFDBv4_90.fasta PATHOGEN.m8 tmp --max-seqs 1
```

> `--max-seqs 1` keeps only the single best Atlas match per query protein.

#### Step 1 — Set your .m8 file path (Cell 3) — THE ONLY CELL YOU NEED TO CHANGE

```python
# Cell 3 — UNIVERSAL CONFIG
M8_FILE = './PA.m8'   # <-- change to your MMseqs2 output file
# Examples:
# M8_FILE = './SA.m8'   # S. aureus
# M8_FILE = './KP.m8'   # K. pneumoniae
```

#### Step 2 — Atlas data files (Cell 5) — no change needed if files are in same directory

```python
atlas_datafiles = [
    'AFDB90v4_cc_data_uniprot_community_taxonomy_map_with_brightness.csv',
    'AFDB90v4_dust_uniprot_community_taxonomy_map_with_brightness.csv'
]
```

#### Step 3 — ESKAPE genus list (Cell 25) — no change needed for standard ESKAPE analysis

The full ESKAPE list is hardcoded. Add genera here if you want to expand the AMR list:

```python
# Cell 25 — replace list to add/remove AMR genera
genus_modifed = atlas_data_selected['genus'].replace(
    ['Acinetobacter', 'Enterococcus', 'Escherichia',
     'Klebsiella', 'Pseudomonas', 'Staphylococcus', 'Enterobacter'],
    'AMR genus'
)
```

#### Step 4 — FASTA directory (Cell 40) — update if your files are elsewhere

```python
# Cell 40
directory_path = './faa'   # <-- path to your .faa files
```

Cell 40 does three things automatically:
- Reads all `.faa` files and builds the `seq_df` table with sequences
- Extracts strain identity from each **filename** (not from sequence headers)
- Saves `queryID_to_strain.csv` — the universal strain mapping used by Part III for S4 scoring

Cell 41 filters `queryID_to_strain.csv` to only proteins that matched the Atlas, then re-saves it.

Cell 43 automatically removes any `LongLink` queryIDs (an occasional MMseqs2 artefact) without needing a hardcoded row index.

#### Outputs from Part I

| File | Contents |
|------|----------|
| `eclipse.csv` | All query proteins with Atlas mapping and diversity metrics |
| `eclipse_seq.csv` | Same as above with protein sequences and strain labels added |
| `queryID_to_strain.csv` | Universal protein-to-strain mapping (used by Part III) |
| `dark_communities.png` | Histogram of functional brightness at community level |
| `components_darkness.png` | Histogram of functional brightness at component level |
| `taxonomic_diversity_of_All_components.png` | Scatter plot of diversity metrics for all components |
| `taxonomic_diversity_of_dark_components.png` | Same, restricted to dark components |

---

### Part II — Genus-specific & ESKAPE-enriched Component Extraction {#part-ii}

**Notebook:** `ECLIPSE_PartII.ipynb`

**What it does:** Takes the dark-component subset from Part I and splits it into two biologically meaningful tracks for downstream prioritization. Everything is controlled by a single config block in Cell 3.

#### Step 1 — Set the config block (Cell 3) — THE ONLY CELL YOU NEED TO CHANGE

```python
# Cell 3 — UNIVERSAL CONFIG
GENUS        = "Pseudomonas"              # genus for genus-specific filter
SPECIES      = "Pseudomonas aeruginosa"   # species string for Atlas proportion
SPECIES_COL  = "p_aeruginosa_proportion"  # column name in output CSVs
INPUT_CSV    = "eclipse_search_results_component_dark .csv"  # your Part I dark subset
PS_OUT_CSV   = "./mapped_p_aer_dataset_pseudomonas_specific_dark_components.csv"
ES_OUT_CSV   = "./mapped_p_aer_dataset_eskape_enriched_dark_components.csv"
```

For *S. aureus* change to:

```python
GENUS        = "Staphylococcus"
SPECIES      = "Staphylococcus aureus"
SPECIES_COL  = "s_aureus_proportion"
INPUT_CSV    = "eclipse_search_results_component_dark.csv"
PS_OUT_CSV   = "./mapped_s_aur_dataset_staphylococcus_specific_dark_components.csv"
ES_OUT_CSV   = "./mapped_s_aur_dataset_eskape_enriched_dark_components.csv"
```

> **What is `INPUT_CSV`?** After running Part I, filter `eclipse_seq.csv` for dark components (e.g. `component_brightness <= 5`) and save it. This filtered file is the input to Part II.

#### What the rest of the notebook does automatically

**Cell 5** — Calculates the proportion of `SPECIES` proteins in each Atlas component. Stored in the `SPECIES_COL` column. No hardcoding.

**Cell 6** — Counts ESKAPE-specific components (`ESKAPE_proportion == 1.0`).

**Cell 7 — Track A (Genus-specific):** Filters for components where ALL Atlas proteins belong to `GENUS`. Uses the `GENUS` variable — no manual edits needed.

**Cell 12 — Track B (ESKAPE-enriched):** Selects components with `ESKAPE_proportion >= 0.5`, excluding Track A components. Threshold can be adjusted here if needed.

#### Variable names used throughout (all generic — no P. aeruginosa hardcoding)

| Variable | Meaning |
|----------|---------|
| `dataset` | All query proteins with dark component mapping |
| `dataset_genus_specific` | Track A: genus-specific dark components |
| `dataset_eskape_enriched` | Track B: ESKAPE-enriched dark components |
| `component_ids_genus_specific` | Component IDs passing the genus filter |
| `SPECIES_COL` | Column name for species proportion (set in config) |

#### Outputs from Part II

| File | Contents |
|------|----------|
| `PS_OUT_CSV` (set in config) | Track A: genus-specific dark components with all metrics |
| `ES_OUT_CSV` (set in config) | Track B: ESKAPE-enriched dark components with all metrics |

---

### Part III — Clustering & DPPS Prioritization {#part-iii}

**Notebook:** `ECLIPSE_DPPS_score.ipynb`

**What it does:** Clusters redundant sequences, computes the DPPS composite score for each component, performs Monte Carlo weight sensitivity analysis, and produces ranked output tables and figures.

#### Step 1 — Set the config block (Cell 2) — THE ONLY CELL YOU NEED TO CHANGE

```python
# Cell 2 — UNIVERSAL CONFIG
GENUS        = "Pseudomonas"
SPECIES      = "Pseudomonas aeruginosa"
SPECIES_COL  = "p_aeruginosa_proportion"
SPECIES_ABBR = "PA"                        # short label for plots
FAA_DIR      = "./faa"                     # folder of your .faa files

PS_CSV = './mapped_p_aer_dataset_pseudomonas_specific_dark_components.csv'  # Track A (Part II output)
ES_CSV = './mapped_p_aer_dataset_eskape_enriched_dark_components.csv'       # Track B (Part II output)
```

For *S. aureus* change to:

```python
GENUS        = "Staphylococcus"
SPECIES      = "Staphylococcus aureus"
SPECIES_COL  = "s_aureus_proportion"
SPECIES_ABBR = "SA"
FAA_DIR      = "./faa"

PS_CSV = './mapped_s_aur_dataset_staphylococcus_specific_dark_components.csv'
ES_CSV = './mapped_s_aur_dataset_eskape_enriched_dark_components.csv'
```

#### N_STRAINS — auto-counted from your FAA directory

`N_STRAINS` (the denominator for S4 strain coverage) is **automatically counted** from your `FAA_DIR` folder. You do not set it manually. The counter:

- Counts unique GCF/GCA accessions from filenames (so multi-file assemblies — chromosome + plasmid — count as one strain)
- Falls back to counting unique filename stems for custom-named files
- Prints the count at runtime so you can verify it

#### Strain counting — auto-detected method (Cell 14)

The notebook automatically detects which strain counting method is appropriate for your data:

**Format A — Custom-prefix queryIDs** (e.g. *P. aeruginosa* `PAO1_PA0001`):
- Detected when queryIDs contain short prefixes before an underscore (<=10 chars)
- Uses `queryID.split('_')[0]` to extract strain name
- This is the original published method for P. aeruginosa — results are identical

**Format B — NCBI accession queryIDs** (e.g. *S. aureus* `ABD20461.1`):
- Detected when queryIDs are plain accessions with no strain information
- Uses `queryID_to_strain.csv` generated by Part I Cell 40
- Required for any pathogen downloaded from NCBI with standard accession-based headers

No manual switching needed — the detection is automatic.

#### DPPS weights (Cell 2) — adjust if needed

**Track A weights (genus-specific):**

```python
WEIGHTS_PS = {
    'S1_darkness':           0.15,   # darkness of the component
    'S2b_pa_combined':       0.40,   # max(Atlas species proportion, strain fraction)
    'S3_specificity':        0.25,   # 1 - ESKAPE_relative_evenness
    'S4_pa_strain_coverage': 0.20,   # unique strains / N_STRAINS
}
# Must sum to 1.0
```

**Track B weights (ESKAPE-enriched):**

```python
WEIGHTS_ES = {
    'S1_darkness':           0.15,
    'S2_pa_proportion':      0.25,   # Atlas species proportion
    'S3_specificity':        0.20,
    'S4_pa_strain_coverage': 0.15,
    'S5_eskape_enrich':      0.25,   # ESKAPE_proportion x (1 - ESKAPE_genus_evenness)
}
```

> **Why S2b for Track A?** Genus-specific components often show zero Atlas species proportion despite being conserved across nearly all strains, because UniProt frequently uses non-canonical names (e.g. "Pseudomonas sp." instead of "Pseudomonas aeruginosa"). S2b = max(Atlas proportion, strain fraction) corrects for this annotation gap.

#### Outputs from Part III

| File | Contents |
|------|----------|
| `dpps_components_pseudomonas_specific.csv` | Track A component-level DPPS scores |
| `dpps_components_eskape_enriched.csv` | Track B component-level DPPS scores |
| `dpps_representatives_pseudomonas_specific.csv` | Track A representative proteins with scores |
| `dpps_representatives_eskape_enriched.csv` | Track B representative proteins with scores |
| `dpps_tier1_pseudomonas_specific.csv` | Tier I candidates — Track A |
| `dpps_tier1_eskape_enriched.csv` | Tier I candidates — Track B |
| `dpps_all_representatives_combined.csv` | Both tracks merged, sorted by DPPS |
| `dpps_distribution_both_tracks.pdf` | DPPS score histogram by tier |
| `dpps_pa_proportion_scatter.pdf` | Species proportion vs DPPS scatter plot |
| `dpps_heatmap_pseudomonas_specific.pdf` | Sub-score heatmap — top 30 Tier I Track A |
| `dpps_heatmap_eskape_enriched.pdf` | Sub-score heatmap — top 30 Tier I Track B |
| `pa_strain_coverage_distribution.pdf` | Strain fraction distribution |
| `dpps_sensitivity_both_tracks.pdf` | Monte Carlo weight sensitivity scatter |

---

## Quick-Switch Guide: Running on a New Pathogen

This is the complete checklist. Only these items change — nothing else.

### Terminal (before any notebook)

```bash
# Merge your .faa files and run MMseqs2
cat ./faa/*.faa > all_proteins.faa
mmseqs easy-search all_proteins.faa AFDBv4_90.fasta SA.m8 tmp --max-seqs 1
```

### Part I — Cell 3 only

```python
M8_FILE = './SA.m8'
```

### Part II — Cell 3 only

```python
GENUS        = "Staphylococcus"
SPECIES      = "Staphylococcus aureus"
SPECIES_COL  = "s_aureus_proportion"
INPUT_CSV    = "eclipse_search_results_component_dark.csv"
PS_OUT_CSV   = "./mapped_s_aur_dataset_staphylococcus_specific_dark_components.csv"
ES_OUT_CSV   = "./mapped_s_aur_dataset_eskape_enriched_dark_components.csv"
```

### Part III — Cell 2 only

```python
GENUS        = "Staphylococcus"
SPECIES      = "Staphylococcus aureus"
SPECIES_COL  = "s_aureus_proportion"
SPECIES_ABBR = "SA"
FAA_DIR      = "./faa"
PS_CSV       = './mapped_s_aur_dataset_staphylococcus_specific_dark_components.csv'
ES_CSV       = './mapped_s_aur_dataset_eskape_enriched_dark_components.csv'
```

**That is everything.** Run all cells in order. All scoring, plotting, and output naming adapts automatically.

---

## Strain Identity & S4 Scoring — How It Works {#strain-identity}

S4 (strain coverage) = `unique strains carrying a component / total strains (N_STRAINS)`

This requires knowing which strain each protein comes from. The pipeline handles this in two ways depending on your queryID format:

### Format A — Custom-prefix queryIDs (e.g. P. aeruginosa)

If your proteins are named `PAO1_PA0001`, `LESB58_PALES_01`, etc., strain identity is extracted by splitting on the first underscore:

```
PAO1_PA0001     ->  strain = PAO1
LESB58_PALES_01 ->  strain = LESB58
```

This is the original published method for *P. aeruginosa*. Auto-detected when queryID prefixes are short (<= 10 characters).

### Format B — NCBI accession queryIDs (e.g. S. aureus)

If your proteins are named `ABD20461.1`, `WP_000123456.1`, etc., there is no strain information in the queryID. Strain identity is instead extracted from the `.faa` **filename** in Part I Cell 40:

```
GCF_000013425.1_ASM1342v1_protein.faa  ->  strain = GCF_000013425.1
MRSA252.faa                             ->  strain = MRSA252
PAO1.faa                                ->  strain = PAO1
```

The mapping is saved as `queryID_to_strain.csv` and loaded by Part III automatically.

### Why filename-based extraction works universally

| Filename format | Extracted strain |
|----------------|-----------------|
| `GCF_000013425.1_ASM1342v1_protein.faa` | `GCF_000013425.1` |
| `GCF_000240185.2_Kpneumo_protein.faa` | `GCF_000240185.2` |
| `PAO1.faa` | `PAO1` |
| `MRSA252.faa` | `MRSA252` |
| `Ab307-0294.faa` | `Ab307-0294` |

GCF/GCA accession is detected by regex when present — this prevents multi-chromosome assemblies (chromosome.faa + plasmid.faa for the same strain) from being counted twice. N_STRAINS is also counted using the same logic.

---

## DPPS Formula Reference

| Sub-score | Formula | Track A weight | Track B weight |
|-----------|---------|----------------|----------------|
| S1 — darkness | `1 - (component_brightness / 100)` | 0.15 | 0.15 |
| S2b — combined species evidence *(Track A)* | `max(Atlas species proportion, strain fraction)` | 0.40 | — |
| S2 — species proportion *(Track B)* | Fraction of Atlas members annotated as target species | — | 0.25 |
| S3 — taxonomic specificity | `1 - ESKAPE_relative_evenness` | 0.25 | 0.20 |
| S4 — strain coverage | `unique strains carrying component / N_STRAINS` | 0.20 | 0.15 |
| S5 — ESKAPE enrichment *(Track B)* | `ESKAPE_proportion x (1 - ESKAPE_genus_evenness)` | — | 0.25 |

### Tier Thresholds

| Tier | DPPS range | Interpretation |
|------|-----------|----------------|
| I | >= 0.75 | High-priority candidates |
| II | 0.50 – 0.75 | Moderate priority |
| III | 0.25 – 0.50 | Low priority |
| IV | < 0.25 | Deprioritized |

### Monte Carlo Sensitivity Analysis

Weights are randomly sampled from a Dirichlet distribution (500 replicates). `tier1_stability` reports the fraction of replicates in which a component reaches Tier I. A value >= 0.8 means the ranking is robust to weight uncertainty.

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
| `strain` | Strain identity extracted from the .faa filename (added in Part I Cell 40) |

### Key columns in DPPS output files

| Column | Description |
|--------|-------------|
| `componentID` | Atlas component identifier |
| `DPPS` | Composite prioritization score (0-1, higher = better candidate) |
| `tier` | Priority tier (I = highest) |
| `S1_darkness` | Darkness sub-score |
| `S2b_pa_combined` / `S2_pa_proportion` | Species evidence sub-score |
| `S3_specificity` | Taxonomic specificity sub-score |
| `S4_pa_strain_coverage` | Strain coverage sub-score |
| `S5_eskape_enrich` | ESKAPE enrichment sub-score (Track B only) |
| `PA_strain_count` | Number of strains carrying this component |
| `PA_strain_fraction` | Fraction of total strains carrying this component |
| `tier1_stability` | Fraction of Monte Carlo replicates where component reaches Tier I |
| `cluster_size` | Number of sequences in the MMseqs2 cluster for this representative |

---

## Key Concepts

**Functional brightness / darkness:** A metric from the Protein Universe Atlas quantifying the fraction of proteins in a community that have functional annotations. A brightness of 0% means no protein in the community has a known function — a "dark" community. Components with median brightness <= 5% are considered dark in this pipeline.

**Community vs. component:** In the Atlas, proteins are grouped into communities based on structural similarity. Communities that are connected in the global similarity network form larger components (subgraphs). Darkness is assessed at both levels.

**ESKAPE_relative_evenness:** A normalized Shannon diversity index computed over all genera in a component, where all ESKAPE genera are merged into a single "AMR genus" label. A value of 0 means the component is exclusively occupied by ESKAPE proteins. Higher values indicate proteins from non-ESKAPE organisms are also present.

**ESKAPE_genus_evenness:** The same index computed only within ESKAPE organisms. A value of 0 means all ESKAPE proteins in the component come from a single genus (maximally genus-specific). Higher values indicate multi-genus ESKAPE components.

**Track A (genus-specific):** Components where all Atlas proteins belong to a single target genus (e.g. all Pseudomonas). These are the most taxon-restricted dark families — highest specificity for the target organism.

**Track B (ESKAPE-enriched):** Components where at least 50% of Atlas proteins are from any ESKAPE genus. These are broadly conserved across drug-resistant pathogens — potential broad-spectrum targets.

**S2b (combined species evidence):** Introduced for Track A because species names like "Pseudomonas sp." in UniProt suppress the Atlas-derived species proportion despite genuine conservation. S2b = max(Atlas proportion, strain fraction) corrects for this annotation gap without penalizing truly conserved components.

---

## Troubleshooting

**`FileNotFoundError: Cannot find ./queryID_to_strain.csv`**
This file is generated by Part I Cell 40. Run Part I completely before running Part III. If your working directory differs between notebooks, copy `queryID_to_strain.csv` to your Part III working directory.

**`UnicodeDecodeError` when reading FASTA files (Part I Cell 40)**
Some FASTA files use non-UTF-8 encoding. The updated Cell 40 automatically retries with Latin-1 encoding. If the error persists, check if the file is a valid FASTA with `head yourfile.faa`.

**LongLink entries in Part I Cell 43**
MMseqs2 occasionally truncates very long sequence names to "LongLink". Cell 43 now handles this automatically using a label-based filter — no hardcoded row index needed.

**`KeyError: SPECIES_COL` in Part III**
The column name set in `SPECIES_COL` in Part III Cell 2 must exactly match what was used in Part II Cell 3. Copy the value directly: if Part II used `"s_aureus_proportion"`, Part III must use the same string.

**`PA_strain_fraction` values > 1.0**
This means `N_STRAINS` was counted as fewer than the actual number of strains in the data. Check the auto-counter output printed when Cell 2 runs. If it looks wrong, verify that all your `.faa` files are present in `FAA_DIR` and are not nested in subdirectories.

**MMseqs2 clustering fails**
Ensure `mmseqs` is in PATH. Check that the FASTA files were written correctly:

```bash
wc -l track_ps_pseudomonas_specific.fasta   # should be 2x number of sequences
```

**No Tier I components found**
Possible causes: (1) `N_STRAINS` is too high, suppressing S4 scores — verify the auto-count; (2) too few strains in your dataset relative to `MIN_SEQ_LEN` filter removing most components; (3) your pathogen genuinely has no highly-scoring dark components — lower `TIER_BINS[3]` from 0.75 to inspect Tier II.

**`values[0]` IndexError in Part I mapping loops (Cells 9, 10, 16)**
This occurs when a `targetID` from your search results is not present in the Atlas CSV files. Verify that both Atlas files loaded correctly in Cell 5 and that your `AFDBv4_90.fasta` matches the Atlas CSV version.

---

## Citation

If you use ECLIPSE in your research, please cite this preprint (https://www.biorxiv.org/cgi/content/short/2026.03.30.715302v1).

Please also cite:

1. Durairaj, J. et al. (2023) Uncovering new families and folds in the natural protein universe. *Nature*, 622, 646–653.
2. Steinegger, M. and Söding, J. (2017) MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nat. Biotechnol.*, 35, 1026–1028.

---

## License

MIT License
