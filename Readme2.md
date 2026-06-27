# ECLIPSE

**E**SKAPE **C**onnectome **L**inkage and **I**nference for **P**roteome **S**equence **E**xploration — a network-based pipeline for identifying and prioritising functionally dark protein families in bacterial panproteomes.

> If you use ECLIPSE, please cite: Lata & Heinz (2026), *Bioinformatics*.

---

## 1. Installation

```bash
git clone https://github.com/surabhilata/ECLIPSE.git
cd ECLIPSE
conda env create -f environment.yml
conda activate eclipse
```

Requires Python ≥ 3.9, MMseqs2 ≥ 14.7e284, and ~16 GB RAM.

---

## 2. Quick-start: run ECLIPSE on test data

A small test dataset (`test_data/`) is provided so you can verify your installation in under 5 minutes.

```bash
# Step 1 — MMseqs2 search (pre-computed test .m8 already provided)
# Skip to Step 2 if using the bundled test data.

# Step 2 — Open the notebook
jupyter notebook ECLIPSE_complete.ipynb
```

In the `CONFIG` cell at the top, the defaults already point at `test_data/`. Run all cells.

**Expected output on the test dataset:**
| Metric | Expected value |
|---|---|
| Sequences mapped to Atlas | ~4,900 |
| Dark components identified | ~25 |
| *Pseudomonas*-specific dark components | ~3 |
| ESKAPE-enriched dark components | ~7 |
| Tier I components (Track A + B combined) | 1–2 |
| Runtime on a 16 GB-RAM laptop | < 5 minutes |

If your output matches these values within ~10%, the pipeline is working correctly.

---

## 3. Running ECLIPSE on your own data

### Inputs you need

1. **Protein FASTA files**, one per strain, in a directory (e.g. `./faa/`).
2. **MMseqs2 `easy-search` result (`.m8` file)** of your panproteome against the Atlas AFDB_90 database. Run this once:

```bash
   mmseqs easy-search your_panproteome.fasta atlas_afdb90 your_panproteome.m8 tmp \
     --max-seqs 1 -e 1e-4
```

   (Atlas AFDB_90 database download: https://uniprot3d.org/atlas)

### Edit the CONFIG cell

Open `ECLIPSE_complete.ipynb` and edit the top `CONFIG` cell:

```python
GENUS        = "Pseudomonas"            # your target genus
SPECIES      = "Pseudomonas aeruginosa" # your target species
SPECIES_COL  = "p_aeruginosa_proportion"
SPECIES_ABBR = "PA"
FAA_DIR      = "./faa"                  # path to FASTA directory
M8_FILE      = "./PA.m8"                # path to your .m8 file
MIN_SEQ_LEN  = 300                      # length filter; user-configurable
```

Then run all cells in order.

### Expected runtime and memory at panproteome scale

| Dataset | Sequences | Runtime | Peak memory |
|---|---|---|---|
| Test data | ~5,000 | < 5 min | < 2 GB |
| Full *P. aeruginosa* panproteome (635 strains) | 3.46M | ~34 min | ~5.5 GB |

---

## 4. Outputs

| File | Contents |
|---|---|
| `eclipse_seq.csv` | Per-protein dark/bright classification with Atlas community and component IDs |
| `mapped_*_pseudomonas_specific_dark_components.csv` | Track A candidates |
| `mapped_*_eskape_enriched_dark_components.csv` | Track B candidates |
| `dpps_components_*.csv` | DPPS scores, tier assignments, Monte Carlo stability |
| `figures/` | Output plots (brightness histograms, taxonomic diversity, DPPS distributions, sub-score heatmaps) |

---

## 5. Adapting ECLIPSE to other ESKAPE pathogens

To run ECLIPSE on a different ESKAPE pathogen (e.g. *Klebsiella pneumoniae*), only the `CONFIG` block needs editing — change `GENUS`, `SPECIES`, `SPECIES_COL`, `SPECIES_ABBR`, `FAA_DIR`, and `M8_FILE`. No other code changes are required. The pipeline is genus-agnostic by construction.

---

## 6. Technical reference: notebook structure

For users who want to understand or modify individual analytical steps, the notebook is organised into four sections:

- **Section 1 — Configuration** (cells 1–3): user parameters and imports
- **Section 2 — Part I: Darkness estimation** (cells 4–43): Atlas mapping, community/component brightness, taxonomic diversity
- **Section 3 — Part II: Two-track stratification** (cells 44–55): *Pseudomonas*-specific and ESKAPE-enriched candidate selection
- **Section 4 — Part III: DPPS scoring** (cells 56–74): length filter, sequence clustering, sub-score computation, tier assignment, sensitivity analysis, runtime summary

Inline comments throughout each cell describe what that step does. The full methodological detail is in the manuscript (Lata & Heinz, 2026).

---

## 7. Citation
