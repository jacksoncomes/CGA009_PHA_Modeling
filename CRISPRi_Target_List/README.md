# *R. palustris* CGA009 — CRISPRi targets for PHB on p-coumarate

A minimal, reproducible release of the genome-scale modeling pipeline that nominates single-gene
**CRISPRi down-regulation targets** to increase PHB/PHA production in *Rhodopseudomonas palustris* CGA009
(model **iAN1128**; Navid/Alsiyabi 2019) growing on **p-coumarate**, for a DBTL round-0 screen.

**The whole release leads to one figure:**
[`Results/nb60_targets/figures/nb71_downselected.png`](Results/nb60_targets/figures/nb71_downselected.png)
— *Round-0 CRISPRi targets, after literature downselection (40 kept)* — produced by notebook `03`.

> Everything here is a **prioritized hypothesis set for a round-0 screen — not a validated design.**

---

## Contents

```
CRISPRi_Target_List/
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
    nb60_targets/                    method CSVs, selection, downselected list, cut list, FluxRETAP cache, figures/
    paper_figs/                      paper figures extracted for the side-by-side in notebook 02
    nb05_*.png                       notebook-02 reproduction figures
  requirements.txt
```

> The CASOP flux-sampling cache (`_casop_coumarate_samples.pkl.gz`, ~70 MB) is **not** shipped: notebook 03
> regenerates it deterministically from the fixed seed (see [Determinism](#determinism)). The only committed
> scan cache is the small FluxRETAP one.

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
| **03** | `03_targets_and_figure.ipynb` | Runs the complete target-selection pipeline end-to-end and renders the headline figure. | 50 combined targets → **40 kept after downselection** |

### What notebook 03 runs (all in one notebook, each part feeds the next)
| Part | Method / step | Writes |
|---|---|---|
| 1 | **FVSEOF** — flux-response down-regulation scan (51 targets) | `nb60_fvseof_targets.csv` |
| 2 | **FluxRETAP** — JBEI 1/overlap flux-range shift (53 targets) | `nb61_fluxretap_targets.csv` |
| 3 | **CASOP** — yield-stratified importance via **seeded** flux sampling (`SEED=42`, `N=15000`); scores 12 genes, contributes its **rank-1** gene `RPA1578` | `nb63_casop_targets.csv` |
| 4 | **Consensus + round-0 selection** (27 genes) | `nb64_master_target_table.csv`, `nb64_round0_selection.csv` |
| 5 | **Merge curated AI predictions + literature downselection → the figure** | `nb71_downselected_targets.csv`, `nb71_cut_list.csv`, `figures/nb71_downselected.{png,svg}` |

The figure's four columns are the three FBA methods (FVSEOF / FluxRETAP / CASOP) plus a curated
**AI-prediction** column (the literature table embedded in `lib/nb71_lib.py`). Scoring is global — every
displayed row shows its CASOP score whether or not CASOP *selected* it — while CASOP contributes only its
single rank-1 gene to the target set. Part 5 cuts 10 redundant/counter-indicated metabolic targets and keeps
the 5 model-blind regulators; the reasoning for every cut is printed and saved to `nb71_cut_list.csv`.

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

> The version pins in `requirements.txt` are load-bearing for *exact* reproduction. The FluxRETAP scan is a
> small committed cache (`_fluxretap_coumarate_ranges.pkl`); the CASOP scan is **seeded** and regenerated on
> demand rather than shipped (see [Determinism](#determinism)). GLPK/cobra determinism is what makes the
> numbers byte-stable — a very different environment may shift low-order digits. The first run of notebook 03
> takes a few extra minutes because CASOP draws its 15,000 samples (then caches them locally, git-ignored).

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
- **FluxRETAP** (Part 2) loads a small committed FVA cache (`_fluxretap_coumarate_ranges.pkl`) as a
  byte-stability anchor; recomputation is numerically identical but reorders exact ties, and the 40-gene
  target set is unaffected.
- **CASOP** (Part 3) is **seeded** (`SEED = 42`, `N_SAMPLES = 15000`): the OptGP sampler is deterministic, so
  its scores are byte-identical run to run and across the `FORCE_RERUN` True/False paths. Its ~70 MB sample
  cache (`_casop_coumarate_samples.pkl.gz`) is **not committed** — notebook 03 regenerates it from the seed
  whenever it is absent, so a fresh clone reproduces the figure from code + seed + model alone. CASOP
  contributes only its **rank-1** gene (`RPA1578`): its score *magnitude* drifts between seeds but its rank
  does not, so rank-1 selection is seed-invariant (verified on seeds 1/2/42).
- The **AI-prediction** column is a fixed curated table, not a live model call.

## Expected key numbers
| Quantity | Value |
|---|---|
| Model-fix growth (acetate point, XR57=1.96) | **0.0824 /h** |
| Paper reproduction (notebook 02) | **12 / 12 checks** |
| Coumarate μmax (growth state) | **0.0773 /h** (doubling 8.97 h) |
| FVSEOF / FluxRETAP / CASOP scored | **51 / 53 / 12** |
| CASOP selects (rank-1) | **RPA1578** |
| Round-0 selection (Part 4) | **27 genes** |
| Combined FBA + AI list (Part 5) | **50 targets** |
| **After downselection (the figure)** | **40 kept** (10 cut, 5 regulators retained) |

## Citations
- Navid/Alsiyabi 2019, *BMC Bioinformatics* 20:233 — the iAN1128 model and the p-coumarate /
  photon (36.6) conditions.
- McKinlay & Harwood 2010, *PNAS* — CO₂ fixation as the dominant electron sink on reduced substrates.
- FVSEOF: Choi et al. 2010 *AEM*; Park et al. 2012 *BMC Syst Biol*. FluxRETAP: Tibocha-Bonilla et al.
  (JBEI). CASOP: Hädicke & Klamt 2010 *J Biotechnol*.
