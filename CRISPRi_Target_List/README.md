# *R. palustris* CGA009 — CRISPRi targets for PHB on p-coumarate

A minimal, reproducible release of the genome-scale modeling pipeline that nominates single-gene
**CRISPRi down-regulation targets** to increase PHB/PHA production in *Rhodopseudomonas palustris* CGA009
(model **iAN1128**; Navid 2019) growing on **p-coumarate**, for a DBTL round-0 screen.

**The whole release leads to one figure:**
[`Results/nb60_targets/figures/nb71_downselected.png`](Results/nb60_targets/figures/nb71_downselected.png)
— *Round-0 CRISPRi targets, after literature downselection (41 kept)* — produced by notebook `03`.

> Everything here is a **prioritized hypothesis set for a round-0 screen — not a validated design.**

---

## Contents

```
Release_github/
  Model/                     the model, before and after the fix
    12859_2019_2844_MOESM3_ESM.xml   published supplementary SBML (Navid/Alsiyabi 2019) — UNLOADABLE as-is
    CGA009_model_biomass_fix.xml     the working model the pipeline loads (published file + the one fix)
  Notebooks/
    01_model_fix.ipynb               how the published SBML was repaired into the working model
    02_reproduce_navid_2019.ipynb    validation: reproduces the published paper's results (12/12)
    03_targets_and_figure.ipynb      the full target-selection pipeline → the headline figure
  lib/
    nb60_lib.py                      shared config: medium, IDs, CRISPRi repression, filters, viability
    nb71_lib.py                      curated AI-prediction table + the figure renderer
    FluxRETAP/                       vendored JBEI FluxRETAP library (unmodified)
  builders/
    build_01_model_fix.py            regenerates notebook 01
    build_02_reproduce_navid.py      regenerates notebook 02
    build_03_targets_and_figure.py   regenerates notebook 03 (self-contained; every cell embedded)
  Results/
    nb60_targets/                    method CSVs, selection, downselected list, cut list, caches, figures/
    paper_figs/                      paper figures extracted for the side-by-side in notebook 02
    nb05_*.png                       notebook-02 reproduction figures
  requirements.txt
```

### The two model files
| File | Role |
|---|---|
| `Model/12859_2019_2844_MOESM3_ESM.xml` | The **published** supplementary SBML. **cobra cannot load it** — it references a biomass boundary species (`XC54_b`) it never declares. This is the *before*. |
| `Model/CGA009_model_biomass_fix.xml` | The **working model** the pipeline loads: the published file plus the one-line `XC54_b` fix (see notebook `01`). This is the *after*. |

---

## The three notebooks

| # | Notebook | What it does | Key result |
|---|---|---|---|
| **01** | `01_model_fix.ipynb` | Repairs the published SBML (adds the missing `XC54_b` biomass boundary species) and proves the repaired model is byte-identical to the shipped working model. | growth = 0.0824 /h (Δ = 0 vs shipped) |
| **02** | `02_reproduce_navid_2019.ipynb` | Validates the working model by reproducing the published Navid/Alsiyabi 2019 results. | **12/12 checks reproduce the paper** |
| **03** | `03_targets_and_figure.ipynb` | Runs the complete target-selection pipeline end-to-end and renders the headline figure. | 51 combined targets → **41 kept after downselection** |

### What notebook 03 runs (all in one notebook, each part feeds the next)
| Part | Method / step | Writes |
|---|---|---|
| 1 | **FVSEOF** — flux-response down-regulation scan (51 targets) | `nb60_fvseof_targets.csv` |
| 2 | **FluxRETAP** — JBEI 1/overlap flux-range shift (53 targets) | `nb61_fluxretap_targets.csv` |
| 3 | **CASOP** — yield-stratified importance via flux sampling (11 targets) | `nb63_casop_targets.csv` |
| 4 | **Consensus + round-0 selection** (30 genes) | `nb64_master_target_table.csv`, `nb64_round0_selection.csv` |
| 5 | **Merge curated AI predictions + literature downselection → the figure** | `nb71_downselected_targets.csv`, `nb71_cut_list.csv`, `figures/nb71_downselected.{png,svg}` |

The figure's four columns are the three FBA methods (FVSEOF / FluxRETAP / CASOP) plus a curated
**AI-prediction** column (the literature table embedded in `lib/nb71_lib.py`). Part 5 cuts 10 redundant/
counter-indicated metabolic targets and keeps the 5 model-blind regulators; the reasoning for every cut
is printed and saved to `nb71_cut_list.csv`.

---

## Setup — download, install, Run All

There is exactly **one** setup step. The notebooks use the standard `python3` kernel that ships with every
Jupyter install, so no custom kernel registration is needed.

```bash
pip install -r requirements.txt      # Python 3.12; installs the exact pinned versions used here
```

Then open the notebooks in `Notebooks/` (Jupyter Lab / Notebook / VS Code) and **Run All**, in order
`01 → 02 → 03`. Each is self-contained; notebook `03` runs the whole target pipeline top-to-bottom and
renders the figure. No commercial solver is needed — **GLPK** (bundled with cobra/optlang via swiglpk) is
the default and is what this release was reproduced with. `lib/nb60_lib.py` sets
`cobra.Configuration().processes = 1` (Windows multiprocess FVA is broken).

> The version pins in `requirements.txt` are load-bearing for *exact* reproduction: the cached scans in
> `Results/nb60_targets/` are pandas pickles (the CASOP one is gzip-compressed, `*.pkl.gz`, loaded
> transparently), and GLPK/cobra determinism is what makes the numbers byte-stable. A very different
> environment may shift low-order digits.

The notebooks in `Notebooks/` are **already executed** (open them to see the outputs without running
anything). To re-run headless from the command line instead of the Jupyter UI:

```bash
cd Notebooks
for nb in 01_model_fix 02_reproduce_navid_2019 03_targets_and_figure ; do
  python -m nbconvert --to notebook --execute --inplace \
         --ExecutePreprocessor.timeout=2400 $nb.ipynb
done
```

To regenerate the blank notebooks from their builders first:
```bash
python builders/build_01_model_fix.py
python builders/build_02_reproduce_navid.py
python builders/build_03_targets_and_figure.py
```

### Determinism
- **FVSEOF** (Part 1) recomputes its FVA scan every run and is fully deterministic (GLPK).
- **FluxRETAP** (Part 2) and **CASOP** (Part 3) load cached scans from `Results/nb60_targets/` (`*.pkl` /
  `*.pkl.gz`)
  by default (`FORCE_RERUN = False` in the notebook). **Keep CASOP cached** — its flux sampler is not
  seeded, so recomputing can shift the target list slightly.
- The **AI-prediction** column is a fixed curated table, not a live model call.

## Expected key numbers
| Quantity | Value |
|---|---|
| Model-fix growth (acetate point, XR57=1.96) | **0.0824 /h** |
| Paper reproduction (notebook 02) | **12 / 12 checks** |
| Coumarate μmax (growth state) | **0.0773 /h** (doubling 8.97 h) |
| FVSEOF / FluxRETAP / CASOP targets | **51 / 53 / 11** |
| Round-0 selection (Part 4) | **30 genes** |
| Combined FBA + AI list (Part 5) | **51 targets** |
| **After downselection (the figure)** | **41 kept** (10 cut, 5 regulators retained) |

## Citations
- Alsiyabi, Immethun & Saha 2019, *BMC Bioinformatics* 20:233 — the iAN1128 model and the p-coumarate /
  photon (36.6) conditions.
- McKinlay & Harwood 2010, *PNAS* — CO₂ fixation as the dominant electron sink on reduced substrates.
- FVSEOF: Choi et al. 2010 *AEM*; Park et al. 2012 *BMC Syst Biol*. FluxRETAP: Tibocha-Bonilla et al.
  (JBEI). CASOP: Hädicke & Klamt 2010 *J Biotechnol*.
