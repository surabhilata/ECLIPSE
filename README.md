# ECLIPSE — Dark Proteome Exploration of ESKAPE Pathogens

**ECLIPSE** (ESKAPE Connectome Linkage and Inference for Proteome Sequence Exploration) is a modular computational pipeline for the systematic identification and prioritisation of functionally ("dark") protein families in ESKAPE panproteomes.

ECLIPSE maps ESKAPE target-pathogen proteomes onto the global sequence similarity network of the [Protein Universe Atlas](https://uniprot3d.org/) (AFDB90v4, UniRef v.2022_03), identifies protein families that are completely dark at the connected-component level, characterises their taxonomic diversity across ESKAPE pathogens, and ranks them using the Dark Proteome Prioritisation Score (DPPS).

> **Note:** This ECLIPSE repository has no overlap with the package *"Eclipse: a Python package for alignment of two or more nontargeted LC-MS metabolomics datasets"* (<https://doi.org/10.1093/bioinformatics/btaf290>).

---

## Table of Contents

- [Installation](#installation)
- [Quick start — run ECLIPSE on test data](#quick-start--run-eclipse-on-test-data)
- [Running ECLIPSE on your own data](#running-eclipse-on-your-own-data)
- [Output files](#output-files)
- [Pipeline architecture and methods](#pipeline-architecture-and-methods)
- [DPPS formula reference](#dpps-formula-reference)
- [Strain identity & S4 scoring](#strain-identity--s4-scoring)
- [Notebook cell-by-cell reference](#notebook-cell-by-cell-reference)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

---

## Installation

### 1. Software requirements

| Tool    | Version | Purpose                      | Install                                                 |
| ------- | ------- | ---------------------------- | ------------------------------------------------------- |
| Python  | >= 3.8  | Core environment             | [python.org](https://python.org)                        |
| Jupyter | any     | Run notebook                 | `pip install jupyter`                                   |
| MMseqs2 | latest  | Sequence search & clustering | [MMseqs2 GitHub](https://github.com/soedinglab/MMseqs2) |

### 2. Python packages

```bash
pip install pandas numpy matplotlib seaborn biopython tqdm pyyaml
```

Or using conda/mamba:

```bash
mamba create -n eclipse python=3.10
mamba activate eclipse
mamba install -c bioconda mmseqs2
pip install jupyter pandas numpy matplotlib seaborn biopython tqdm pyyaml
```

Verify MMseqs2 is accessible:

```bash
mmseqs version
```

### 3. Clone this repository

```bash
git clone https://github.com/surabhilata/ECLIPSE.git
cd ECLIPSE
```

---

## Quick start — run ECLIPSE on test data

A small test dataset is provided to verify your installation and demonstrate the expected output. The dataset is a single *P. aeruginosa* PAO1 proteome, so the full pipeline was tested on CPU Apple M3 of clock frequency 2.40 GHz with 8 physical cores on Darwin 23.3.0 (x86_64) operating system. The test data runtime is 5.6 mins with peak memory of 6 gb. One can also launch the notebook through Virtual Studio Code. 

### 1. Unpack the test data

The test data is bundled as `test_data.zip` in this repository. Unzip it in the repository root:

```bash
unzip test_data.zip
```

This creates a `test_data/` folder containing:

| File / folder                    | Purpose                                                                            |
| -------------------------------- | ---------------------------------------------------------------------------------- |
| `test_data/faa/PAO1.faa`         | Single *P. aeruginosa* PAO1 proteome (5,079 proteins)                              |
| `test_data/PAO1.m8`              | Pre-computed MMseqs2 search result, so the test does not require `AFDBv4_90.fasta` |
| `test_data/expected_summary.txt` | Expected output values for verifying your run                                      |
| `test_data/output/`    | Reference Tier I CSVs for direct comparison                                                  |
| `ECLIPSE_test_data.ipynb/`   | Notebook to run                                                                          
| `eclipse_serach_results_component_dark.csv`  | Part 1 outputs file eclipse_seq.csv from which we select only component with brightness 0 and name it as eclipse_serach_results_component_dark.csv |

### 2. Download the Atlas reference data

The Atlas reference files are required for any ECLIPSE run (test or your own data). They are too large to host in the repository and are available from Zenodo: <https://doi.org/10.5281/zenodo.19119408>.

Download these files and place them in the unzipped downloaded test_data directory:

| File                                                                  | Purpose                                                |
| --------------------------------------------------------------------- | ------------------------------------------------------ |
| `AFDB90v4_cc_data_uniprot_community_taxonomy_map_with_brightness.csv` | Community-level Atlas data (connected components)      |
| `AFDB90v4_dust_uniprot_community_taxonomy_map_with_brightness.csv`    | Community-level Atlas data (DUST regions)              |
| `AFDB90v4_subgraphs_summary.csv`                                      | Component-level summary statistics                     |
| `AFDBv4_90.fasta` (optional)                                          | Only needed if you want to re-run MMseqs2 from scratch |

> The Atlas CSVs only need to be downloaded once and are reused across all ECLIPSE runs.

### 3. (Optional) Re-run MMseqs2 search

A pre-computed `PAO1.m8` is provided so you can skip this step. If you want to verify it from scratch (requires `AFDBv4_90.fasta`):

```bash
mmseqs easy-search test_data/faa/PAO1.faa AFDBv4_90.fasta test_data/PAO1.m8 tmp --max-seqs 1
```

### 4. Open and run the notebook

```bash
## open jupyter notebook after activating the eclipse environment as stated above and then launch jupyter notebook from terminal. Go to folder test_data ( cd test_data).  Use the "ECLIPSE_test_data.ipynb" present in test_data folder. User can also launch this notebook through virtual studio code. 
```

In the **Configuration cell** (Section 1 of the notebook), the default values are already set for the test data:

```python
GENUS        = "Pseudomonas"
SPECIES      = "Pseudomonas aeruginosa"
SPECIES_COL  = "p_aeruginosa_proportion"
SPECIES_ABBR = "PA"
FAA_DIR      = "./faa"
M8_FILE      = "./PAO1.m8"
MIN_SEQ_LEN  = 300
```

Run all cells (`Cell → Run All`). The full pipeline (Part I → Part II → Part III) executes end-to-end in one notebook.

### 5. Expected output

A successful run on this single-strain test produces:

- **Part I** — 4,040 proteins mapped to Atlas across 881 components; 134 dark components (~3 %)
- **Part II** — 56 *Pseudomonas*-specific dark components (Track A) and 302 ESKAPE-enriched dark components (Track B)
- **Part III** — **24 Tier I** *Pseudomonas*-specific components (mean DPPS 0.888) and **9 Tier I** ESKAPE-enriched components (mean DPPS 0.800), after the ≥ 300 aa median-length filter and MMseqs2 redundancy removal
- Full list of CSV and figure outputs in the working directory

A reference summary is provided in `test_data/expected_summary.txt`, and reference Tier I CSVs for direct comparison are in `test_data/output/`. Small numerical variation (± 1 component) is possible due to library versions and Monte Carlo sampling.

> **Note on the single-strain test:** Because the test uses only PAO1, the S4 strain-coverage sub-score saturates at 1.0 for every component by construction. Multi-strain behaviour of S4 is demonstrated in the manuscript on the full *P. aeruginosa* panproteome and is not exercised by this minimal test set.

Approximate runtime: on a laptop of CPU Apple M3 of clock frequency 2.40 GHz with 8 physical cores on Darwin 23.3.0 (x86_64) operating system. The test data runtime is 5.6 mins with peak memory of 6 gb
---

## Running ECLIPSE on your own data

Once the test run succeeds, switching to your own pathogen requires changes in **one configuration cell only**.

### 1. Prepare your input

Create a directory with one `.faa` file per strain/genome assembly:

```
./faa/
   GCF_000013425.1_ASM1342v1_protein.faa
   GCF_000017085.1_ASM1708v1_protein.faa
   ...
```

> **Filename convention:** The pipeline extracts strain identity from the `.faa` **filename**, not from the sequence header. NCBI GCF/GCA accessions are detected automatically; custom names like `PAO1.faa` or `MRSA252.faa` also work. See [Strain identity & S4 scoring](#strain-identity--s4-scoring) for details.

### 2. Run MMseqs2

```bash
cat ./faa/*.faa > all_proteins.faa
mmseqs easy-search all_proteins.faa AFDBv4_90.fasta YOUR_PATHOGEN.m8 tmp --max-seqs 1
```

`--max-seqs 1` keeps only the single best Atlas match per query protein.

### 3. Edit the configuration cell

Open `ECLIPSE_complete.ipynb` and edit the Configuration cell (Section 1). Example values for common ESKAPE pathogens:

| Pathogen        | `GENUS`            | `SPECIES`                   | `SPECIES_COL`               | `SPECIES_ABBR` |
| --------------- | ------------------ | --------------------------- | --------------------------- | -------------- |
| *P. aeruginosa* | `"Pseudomonas"`    | `"Pseudomonas aeruginosa"`  | `"p_aeruginosa_proportion"` | `"PA"`         |
| *S. aureus*     | `"Staphylococcus"` | `"Staphylococcus aureus"`   | `"s_aureus_proportion"`     | `"SA"`         |
| *K. pneumoniae* | `"Klebsiella"`     | `"Klebsiella pneumoniae"`   | `"k_pneumoniae_proportion"` | `"KP"`         |
| *A. baumannii*  | `"Acinetobacter"`  | `"Acinetobacter baumannii"` | `"a_baumannii_proportion"`  | `"AB"`         |
| *E. faecium*    | `"Enterococcus"`   | `"Enterococcus faecium"`    | `"e_faecium_proportion"`    | `"EF"`         |

Example for *S. aureus*:

```python
GENUS        = "Staphylococcus"
SPECIES      = "Staphylococcus aureus"
SPECIES_COL  = "s_aureus_proportion"
SPECIES_ABBR = "SA"
FAA_DIR      = "./faa"
M8_FILE      = "./SA.m8"
MIN_SEQ_LEN  = 300
```

### 4. Run all cells

`Cell → Run All`. The full pipeline (mapping → stratification → scoring) executes in one pass. `N_STRAINS` is counted automatically from `FAA_DIR`. Strain counting method is auto-detected. No other manual edits are required.

---

## Output files

The notebook writes the following files to the working directory:

### Part I outputs

| File                                         | Contents                                                     |
| -------------------------------------------- | ------------------------------------------------------------ |
| `eclipse.csv`                                | All query proteins with Atlas mapping and diversity metrics  |
| `eclipse_seq.csv`                            | Same, with protein sequences and strain labels added         |
| `queryID_to_strain.csv`                      | Protein-to-strain mapping (used internally by Part III)      |
| `dark_communities.png`                       | Histogram of functional brightness at community level        |
| `components_darkness.png`                    | Histogram of functional brightness at component level        |
| `taxonomic_diversity_of_All_components.png`  | Scatter of diversity metrics for all components              |
| `taxonomic_diversity_of_dark_components.png` | Same, restricted to dark components                          |

### Part II outputs

| File                                                              | Contents                                                  |
| ----------------------------------------------------------------- | --------------------------------------------------------- |
| `mapped_*_dataset_<genus>_specific_dark_components.csv`           | Track A: genus-specific dark components with all metrics  |
| `mapped_*_dataset_eskape_enriched_dark_components.csv`            | Track B: ESKAPE-enriched dark components with all metrics |

### Part III outputs (DPPS)

| File                                            | Contents                                    |
| ----------------------------------------------- | ------------------------------------------- |
| `dpps_components_<genus>_specific.csv`          | Track A component-level DPPS scores         |
| `dpps_components_eskape_enriched.csv`           | Track B component-level DPPS scores         |
| `dpps_representatives_<genus>_specific.csv`     | Track A representative proteins with scores |
| `dpps_representatives_eskape_enriched.csv`      | Track B representative proteins with scores |
| `dpps_tier1_<genus>_specific.csv`               | Tier I candidates — Track A                 |
| `dpps_tier1_eskape_enriched.csv`                | Tier I candidates — Track B                 |
| `dpps_all_representatives_combined.csv`         | Both tracks merged, sorted by DPPS          |
| `dpps_distribution_both_tracks.pdf`             | DPPS score histogram by tier                |
| `dpps_pa_proportion_scatter.pdf`                | Species proportion vs DPPS scatter plot     |
| `dpps_heatmap_<genus>_specific.pdf`             | Sub-score heatmap — top 30 Tier I Track A   |
| `dpps_heatmap_eskape_enriched.pdf`              | Sub-score heatmap — top 30 Tier I Track B   |
| `pa_strain_coverage_distribution.pdf`           | Strain fraction distribution                |
| `dpps_sensitivity_both_tracks.pdf`              | Monte Carlo weight sensitivity scatter      |

---

## Pipeline architecture and methods

ECLIPSE runs in three sequential stages, all contained in `ECLIPSE_complete.ipynb`:

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
    Part I — Darkness estimation
    +---------------------------------------------+
    | * Load Atlas community & component data     |
    | * Map each protein to UniRef50 -> community |
    | * Assign darkness (brightness) values       |
    | * Compute ESKAPE taxonomic diversity metrics|
    | * Build queryID_to_strain.csv (strain map)  |
    | * Save eclipse.csv / eclipse_seq.csv        |
    +---------------------------------------------+
            |
            v
    Part II — Two-track stratification
    +---------------------------------------------+
    | * Track A: genus-specific dark components   |
    | * Track B: ESKAPE-enriched dark components  |
    | * Add species-level Atlas proportions       |
    +---------------------------------------------+
            |
            v
    Part III — DPPS scoring & sensitivity analysis
    +---------------------------------------------+
    | * Length filter (>= 300 aa median)          |
    | * MMseqs2 easy-cluster (redundancy removal) |
    | * Auto-detect strain counting method        |
    | * Compute DPPS sub-scores (S1-S5)           |
    | * Monte Carlo weight sensitivity (n=500)    |
    | * Tier ranking (I-IV) + visualizations      |
    +---------------------------------------------+
            |
            v
    Ranked candidate drug target lists
```

The pipeline answers the question:

> *Which proteins in my pathogen of interest are (1) functionally uncharacterized ("dark"), (2) evolutionarily restricted to AMR-relevant clades, and (3) conserved across many strains — making them promising unexplored drug targets?*

---

## DPPS formula reference

| Sub-score                                   | Formula                                               | Track A weight | Track B weight |
| ------------------------------------------- | ----------------------------------------------------- | -------------- | -------------- |
| S1 — darkness                               | `1 - (component_brightness / 100)`                    | 0.15           | 0.15           |
| S2b — combined species evidence *(Track A)* | `max(Atlas species proportion, strain fraction)`      | 0.40           | —              |
| S2 — species proportion *(Track B)*         | Fraction of Atlas members annotated as target species | —              | 0.25           |
| S3 — taxonomic specificity                  | `1 - ESKAPE_relative_evenness`                        | 0.25           | 0.20           |
| S4 — strain coverage                        | `unique strains carrying component / N_STRAINS`       | 0.20           | 0.15           |
| S5 — ESKAPE enrichment *(Track B)*          | `ESKAPE_proportion x (1 - ESKAPE_genus_evenness)`     | —              | 0.25           |

### Tier thresholds

| Tier | DPPS range  | Interpretation           |
| ---- | ----------- | ------------------------ |
| I    | >= 0.75     | High-priority candidates |
| II   | 0.50 – 0.75 | Moderate priority        |
| III  | 0.25 – 0.50 | Low priority             |
| IV   | < 0.25      | Deprioritized            |

### Monte Carlo sensitivity analysis

Weights are randomly sampled from a Dirichlet distribution (α = 1, 500 replicates). `tier1_stability` reports the fraction of replicates in which a component reaches Tier I. A value ≥ 0.8 indicates the ranking is robust to weight uncertainty.

### Why S2b for Track A?

Genus-specific components often show zero Atlas species proportion despite being conserved across nearly all strains, because UniProt frequently uses non-canonical names (e.g. "Pseudomonas sp." instead of "Pseudomonas aeruginosa"). S2b = max(Atlas proportion, strain fraction) corrects for this annotation gap without penalizing truly conserved components.

---

## Strain identity & S4 scoring

S4 (strain coverage) = `unique strains carrying a component / total strains (N_STRAINS)`. This requires knowing which strain each protein comes from. The pipeline handles this in two ways depending on your queryID format, and auto-detects which to apply:

### Format A — Custom-prefix queryIDs (e.g. *P. aeruginosa*)

If your proteins are named `PAO1_PA0001`, `LESB58_PALES_01`, etc., strain identity is extracted by splitting on the first underscore:

```
PAO1_PA0001     ->  strain = PAO1
LESB58_PALES_01 ->  strain = LESB58
```

Auto-detected when queryID prefixes are short (≤ 10 characters).

### Format B — NCBI accession queryIDs (e.g. *S. aureus*)

If your proteins are named `ABD20461.1`, `WP_000123456.1`, etc., there is no strain information in the queryID. Strain identity is instead extracted from the `.faa` **filename**:

| Filename                                | Extracted strain  |
| --------------------------------------- | ----------------- |
| `GCF_000013425.1_ASM1342v1_protein.faa` | `GCF_000013425.1` |
| `GCF_000240185.2_Kpneumo_protein.faa`   | `GCF_000240185.2` |
| `PAO1.faa`                              | `PAO1`            |
| `MRSA252.faa`                           | `MRSA252`         |
| `Ab307-0294.faa`                        | `Ab307-0294`      |

GCF/GCA accession is detected by regex when present — this prevents multi-chromosome assemblies (chromosome + plasmid for the same strain) from being counted twice. `N_STRAINS` is also counted using the same logic.

---

## Notebook cell-by-cell reference

The `ECLIPSE_complete.ipynb` notebook is divided into four sections, mirroring the pipeline:

1. **Configuration** — the single config cell users edit.
2. **Part I — Darkness estimation** — Atlas mapping, brightness assignment, ESKAPE diversity metrics, strain mapping.
3. **Part II — Two-track stratification** — genus-specific and ESKAPE-enriched extraction.
4. **Part III — DPPS scoring & sensitivity analysis** — length filter, clustering, sub-score computation, Monte Carlo sensitivity, tier ranking, visualisation.

Every code cell contains an in-line comment describing what it does. Users who want a deeper understanding of intermediate variables (e.g. `atlas_search_results_maped_with_componentID`, `ESKAPE_relative_evenness`) can read the comments directly in the notebook — they do not need to be consulted for a standard run.

### Key concepts referenced in the notebook

**Functional brightness / darkness:** Fraction of proteins in an Atlas community that have functional annotations. Brightness of 0% means no protein in the community has a known function. Components with median brightness 0% are considered dark.

**Community vs. component:** In the Atlas, proteins are grouped into communities based on structural similarity. Communities that are connected in the global similarity network form larger components (subgraphs). Darkness is assessed at both levels.

**ESKAPE_relative_evenness:** Normalised Shannon diversity index across all genera in a component, with all ESKAPE genera merged into a single "AMR genus" label. A value of 0 means the component is exclusively occupied by ESKAPE proteins.

**ESKAPE_genus_evenness:** The same index computed only within ESKAPE organisms. A value of 0 means all ESKAPE proteins in the component come from a single genus.

**Track A (genus-specific):** Components where all Atlas proteins belong to a single target genus — most taxon-restricted dark families.

**Track B (ESKAPE-enriched):** Components where ≥ 50% of Atlas proteins are from any ESKAPE genus — potential broad-spectrum targets.

### Key columns in output files

| Column                                   | Description                                                           |
| ---------------------------------------- | --------------------------------------------------------------------- |
| `queryID`                                | Protein identifier from your input FASTA                              |
| `targetID`                               | Best-matching UniRef50 representative in the Atlas                    |
| `componentID`                            | Atlas component (subgraph) to which the match belongs                 |
| `component_brightness`                   | Median brightness of the component                                    |
| `ESKAPE_relative_evenness`               | Shannon evenness across all genera (AMR merged); 0 = ESKAPE-exclusive |
| `ESKAPE_genus_evenness`                  | Shannon evenness within AMR genera only; 0 = single-genus             |
| `ESKAPE_proportion`                      | Fraction of component proteins from AMR genera                        |
| `strain`                                 | Strain identity extracted from the .faa filename                      |
| `DPPS`                                   | Composite prioritisation score (0–1)                                  |
| `tier`                                   | Priority tier (I = highest)                                           |
| `S1_darkness`–`S5_eskape_enrich`         | Individual sub-scores                                                 |
| `PA_strain_count` / `PA_strain_fraction` | Strain coverage (raw count and normalised)                            |
| `tier1_stability`                        | Fraction of Monte Carlo replicates where component reaches Tier I     |
| `cluster_size`                           | Number of sequences in the MMseqs2 cluster for this representative    |

---

## Troubleshooting

**`FileNotFoundError: Cannot find ./queryID_to_strain.csv`** — Generated in Part I. Make sure Part I cells have run before Part III. If you split execution across notebooks, copy `queryID_to_strain.csv` to the Part III working directory.

**`UnicodeDecodeError` when reading FASTA files** — Some FASTA files use non-UTF-8 encoding. The notebook automatically retries with Latin-1. If the error persists, check the file with `head yourfile.faa`.

**`KeyError: SPECIES_COL` in Part III** — The column name in the Configuration cell must exactly match what was used upstream. The single-notebook layout prevents this in most cases; if you ran sections independently, verify the string is identical.

**`PA_strain_fraction` values > 1.0** — `N_STRAINS` was counted as fewer than the actual number of strains. Verify the auto-count printed by the Configuration cell, and check that all `.faa` files are in `FAA_DIR` without nested subdirectories.

**MMseqs2 clustering fails** — Ensure `mmseqs` is in `PATH`. Verify FASTA was written:
```bash
wc -l track_ps_<genus>_specific.fasta   # should be 2x number of sequences
```

**No Tier I components found** — Possible causes: (1) `N_STRAINS` is too high, suppressing S4 — verify the auto-count; (2) too few strains relative to the `MIN_SEQ_LEN` filter removing most components; (3) the pathogen genuinely has no high-scoring dark components — inspect Tier II by lowering the threshold.

**`values[0]` IndexError in mapping loops** — A `targetID` in your search result is missing from the Atlas CSV files. Verify the Atlas CSVs loaded correctly and that `AFDBv4_90.fasta` matches the Atlas CSV version.

---

## Citation

If you use ECLIPSE in your research, please cite this preprint: <https://www.biorxiv.org/cgi/content/short/2026.03.30.715302v1>.

Please also cite:

1. Durairaj, J. et al. (2023) Uncovering new families and folds in the natural protein universe. *Nature*, 622, 646–653.
2. Steinegger, M. and Söding, J. (2017) MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nat. Biotechnol.*, 35, 1026–1028.

---

## License

MIT License — see [LICENSE](LICENSE).
