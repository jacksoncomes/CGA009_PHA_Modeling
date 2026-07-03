"""
nb60_lib.py — shared configuration & helpers for the CLEAN NB60+ CRISPRi/PHB target pipeline.
================================================================================================
Rhodopseudomonas palustris CGA009 (Alsiyabi, Immethun & Saha 2019, "iRpa940") on p-coumarate.

This is an INDEPENDENT reimplementation written for the NB60+ redo. It deliberately does NOT import
the prior-effort modules (crispri_shared.py, phb_combo_lib.py, coumarate_methods_faithful.py) — every
modeling choice below was re-derived and re-verified against the SBML model on 2026-06-29. Importing
the vendored JBEI FluxRETAP library is allowed and is done only inside the FluxRETAP notebook.

Verified against Model/CGA009_model_biomass_fix.xml (2026-06-29):
  * growth state (coumarate, NH3 replete): mu_max = 0.0773 1/h, doubling 8.97 h  (paper: ~9.4 h)
  * coumarate uptake self-limits to 0.327 mmol/gDW/h under photon 36.6  (paper: ~0.33)
  * couAB route R01616 (4-coumarate:CoA ligase) -> XR225 -> acetyl-CoA carries the uptake flux
  * phaC R04254 is the SOLE non-demand producer of PHB metabolite C06143
  * coumarate (XR242) is the only OPEN organic-carbon UPTAKE valve after the medium is set;
    XR98 (Manganese) is an ESSENTIAL trace-metal valve (NOT carbon) and is left open;
    XR350 (sulfoacetate) is an organic-C valve carrying ZERO flux -> closed here for a clean medium.

Literature anchors (cite in the notebooks):
  Alsiyabi, Immethun & Saha 2019 BMC Bioinformatics — the model + acetate(1.96)/coumarate/photon(36.6).
  McKinlay & Harwood 2010 PNAS — CO2 fixation as the dominant electron sink (coumarate is more reduced
    than acetate, so PHB/NADPH redox nodes matter more on coumarate).
  Pan et al. 2012 J Bacteriol; Hirakawa et al. 2020 Mol Cell Proteomics — couAB p-coumarate catabolism.
"""
from __future__ import annotations
import io, contextlib
from pathlib import Path

import numpy as np
import pandas as pd
import cobra
from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba, flux_variability_analysis

# On Windows cobra's default multiprocess FVA spawns subprocesses that re-import the calling module
# (breaks notebooks and FluxRETAP's internal FVA). Force single-process solving everywhere.
cobra.Configuration().processes = 1


# ── paths ─────────────────────────────────────────────────────────────────────
def _find_root() -> Path:
    here = Path(__file__).resolve()
    for p in [here.parent] + list(here.parents):
        if (p / "Model" / "CGA009_model_biomass_fix.xml").exists():
            return p
    return Path.cwd()

ROOT = _find_root()
MODEL_PATH = ROOT / "Model" / "CGA009_model_biomass_fix.xml"
OUT = ROOT / "Results" / "nb60_targets"            # all NB60+ outputs land here
FIG = OUT / "figures"
for _d in (OUT, FIG):
    _d.mkdir(parents=True, exist_ok=True)


# ── verified reaction / metabolite IDs ─────────────────────────────────────────
BIOMASS       = "XR229"
PHB_MET       = "C06143"            # Poly-beta-hydroxybutyrate
PHB_DEMAND    = "DM_C06143"         # runtime PHB sink (sole readout of excess PHB)
PHB_SYNTHASE  = "R04254"            # phaC : sole producer of C06143
PHB_REDUCT    = "R01977"            # phaB
COUAB_LIGASE  = "R01616"            # 4-coumarate:CoA ligase (couAB route entry)

PHOTON   = "XR55"
NH3      = "XR90"
CO2_OUT  = "XR73"
CO2_IN   = "XR72"
ACETATE_UP    = "XR57"              # disk-pinned to 1.96 (acetate point)
COUMARATE_UP  = "XR242"
COUMARATE_SEC = "XR243"

# Every organic-carbon UPTAKE valve in the model. The medium closes ALL of these, then re-opens only
# the chosen source. Re-derived by inspection; note XR350 (sulfoacetate) — missed by the prior pool —
# is included so coumarate is unambiguously the only open organic-C valve. (XR98 = Manganese is an
# ESSENTIAL inorganic trace-metal valve, NOT carbon, and is intentionally NOT in this list.)
ORGANIC_C_VALVES = ["XR57", "XR242", "XR80", "XR339", "XR62", "XR96", "XR94", "XR350"]
#                   acetate coumarate succ  pyruvate benz  lactate fumar sulfoacetate

# ── feeding / knockdown constants (literature- and validation-grounded) ────────
PHOTON_UB     = 36.6      # light-limited; Alsiyabi 2019
NH3_REPLETE   = 10.0      # replete N for the growth state
COUMARATE_UB  = 1000.0    # leave coumarate uptake unconstrained; photon (36.6) self-limits it to ~0.33
ACETATE_FIX   = 1.96      # acetate point (fixed), for optional acetate comparisons only

CRISPRI_REMAINING = 0.10  # 90% knockdown -> 10% residual flux (the experimentally relevant operating pt)
MIN_GROWTH_FRAC   = 0.10  # viability floor: KD must keep growth >= 10% of WT mu_max
MIN_FLUX          = 1e-6

PHB_PATHWAY = {PHB_SYNTHASE: "phaC PHB synthase", PHB_REDUCT: "phaB reductase",
               PHB_DEMAND: "PHB demand sink"}
# The PHB-pathway reactions are POSITIVE CONTROLS, not CRISPRi-down targets: repressing phaC/phaB
# lowers PHB. They can surface in a flux-response scan (e.g. phaB's FVA band is driven down as enforced
# PHB rises, a sign artifact), so they are excluded from every method's knockdown output.
PHB_PATHWAY_EXCLUDE = {PHB_SYNTHASE, PHB_REDUCT, PHB_DEMAND}
# Feedstock catabolic-entry reactions: repressing these reduces carbon SUPPLY (starves the cell), so they
# are not valid PHB-redirection CRISPRi targets. They can surface in yield-based scoring (high PHB-yield
# correlates with low carbon throughput). R01616 = 4-coumarate:CoA ligase (the couAB catabolic entry).
FEEDSTOCK_EXCLUDE = {COUAB_LIGASE}


# ── model construction ─────────────────────────────────────────────────────────
def read_model():
    with contextlib.redirect_stderr(io.StringIO()):
        return read_sbml_model(str(MODEL_PATH))


def add_phb_demand(m):
    if PHB_DEMAND not in [r.id for r in m.reactions]:
        dm = cobra.Reaction(PHB_DEMAND, name="PHB demand (sink)")
        dm.add_metabolites({m.metabolites.get_by_id(PHB_MET): -1.0})
        dm.bounds = (0.0, 1000.0)
        m.add_reactions([dm])
    return m


def set_medium(m, carbon="coumarate", nh3_ub=NH3_REPLETE):
    """Open exactly one organic-carbon valve; set photon/NH3/CO2. Anaerobic (O2 valves already 0)."""
    ids = {r.id for r in m.reactions}
    for rid in ORGANIC_C_VALVES:
        if rid in ids:
            m.reactions.get_by_id(rid).bounds = (0.0, 0.0)
    if carbon == "coumarate":
        m.reactions.get_by_id(COUMARATE_UP).bounds = (0.0, COUMARATE_UB)   # unconstrained uptake
        m.reactions.get_by_id(COUMARATE_SEC).bounds = (0.0, 0.0)           # no secretion leak
    elif carbon == "acetate":
        m.reactions.get_by_id(ACETATE_UP).bounds = (ACETATE_FIX, ACETATE_FIX)
    else:
        raise ValueError(f"unknown carbon source {carbon!r}")
    m.reactions.get_by_id(PHOTON).bounds = (-PHOTON_UB, PHOTON_UB)
    m.reactions.get_by_id(NH3).bounds = (0.0, nh3_ub)
    m.reactions.get_by_id(CO2_OUT).bounds = (-1000.0, 1000.0)
    m.reactions.get_by_id(CO2_IN).bounds = (-1000.0, 1000.0)
    return m


def build_growth_state(carbon="coumarate"):
    """NH3-replete reference physiological state; objective = biomass."""
    m = read_model()
    set_medium(m, carbon, NH3_REPLETE)
    add_phb_demand(m)
    m.objective = BIOMASS
    m.objective.direction = "max"
    return m


def build_production_state(carbon="coumarate"):
    """N-starved PHA-accumulation state; objective = PHB demand."""
    m = read_model()
    set_medium(m, carbon, 0.0)
    add_phb_demand(m)
    m.objective = PHB_DEMAND
    m.objective.direction = "max"
    return m


def open_organic_carbon_audit(m):
    """Sanity check: list organic-carbon uptake valves currently open (should be only coumarate)."""
    organic = []
    for r in m.reactions:
        if not r.id.startswith("XR"):
            continue
        subs = [mt for mt, c in r.metabolites.items() if c < 0]
        pros = [mt for mt, c in r.metabolites.items() if c > 0]
        if len(r.metabolites) == 2 and subs and pros and \
           subs[0].id.endswith("_b") and pros[0].id.endswith("XT") and r.upper_bound > 0:
            nm = (subs[0].name or "").lower()
            inorganic = ("water", "photon", "co2", "o2", "oxygen", "hydrogen", "proton", "ammon",
                         "nh3", "phosphate", "sulfate", "potassium", "sodium", "calcium", "magnesium",
                         "mangan", "mn2", "molybd", "zinc", "iron", "fe2", "fe3", "cu2", "cobalt",
                         "chloride", "biotin", "nitrate", "h2o", "thiosulf")
            if not any(k in nm for k in inorganic):
                organic.append((r.id, subs[0].name, r.bounds))
    return organic


# ── solver helpers ─────────────────────────────────────────────────────────────
def safe_slim(m):
    try:
        v = m.slim_optimize()
        return float(v) if v == v and v is not None else 0.0
    except Exception:
        return 0.0


def safe_fva(m, rxns, frac=0.0):
    try:
        return flux_variability_analysis(m, reaction_list=rxns, fraction_of_optimum=frac, processes=1)
    except Exception as e:
        print("  FVA error:", e)
        return None


def apply_repression(rxn, ref_flux, rep_remaining=CRISPRI_REMAINING, eps=MIN_FLUX):
    """Emulate CRISPRi: cap |flux| at max(eps, rep_remaining*|ref_flux|), sign-preserving so a
    reversible reaction cannot flip direction under knockdown."""
    cap = max(eps, rep_remaining * abs(ref_flux))
    if ref_flux >= 0:
        rxn.lower_bound = min(rxn.lower_bound, cap)
        rxn.upper_bound = min(rxn.upper_bound, cap)
    else:
        rxn.upper_bound = max(rxn.upper_bound, -cap)
        rxn.lower_bound = max(rxn.lower_bound, -cap)


# ── reaction metadata ──────────────────────────────────────────────────────────
def metadata():
    m = read_model()
    add_phb_demand(m)
    return dict(
        genes={r.id: "|".join(sorted(g.id for g in r.genes)) for r in m.reactions},
        gpr={r.id: r.gene_reaction_rule for r in m.reactions},
        name={r.id: r.name for r in m.reactions},
        subsys={r.id: (r.subsystem or "") for r in m.reactions},
    )


def gene_candidates(meta):
    """Reactions with a real (non-gapfilled) gene — the universe scanned by each method."""
    return [r for r, g in meta["genes"].items() if g and g.lower() != "gapfilled"]


def drop_phb_pathway(df, reaction_col="reaction", verbose=True):
    """Remove PHB-pathway positive-control reactions from a knockdown target list (see
    PHB_PATHWAY_EXCLUDE). Applied to every method so phaC/phaB never appear as 'repress-me' targets."""
    if df is None or df.empty:
        return df
    n0 = len(df)
    hit = sorted(set(df[reaction_col]) & PHB_PATHWAY_EXCLUDE)
    d = df[~df[reaction_col].isin(PHB_PATHWAY_EXCLUDE)].reset_index(drop=True)
    if verbose and hit:
        print(f"    dropped {n0 - len(d)} PHB-pathway positive-control rxn(s): {hit}")
    return d


def drop_feedstock(df, reaction_col="reaction", verbose=True):
    """Remove feedstock catabolic-entry reactions (see FEEDSTOCK_EXCLUDE): repressing them cuts carbon
    supply rather than redirecting it to PHB, so they are not actionable knockdown targets."""
    if df is None or df.empty:
        return df
    n0 = len(df)
    hit = sorted(set(df[reaction_col]) & FEEDSTOCK_EXCLUDE)
    d = df[~df[reaction_col].isin(FEEDSTOCK_EXCLUDE)].reset_index(drop=True)
    if verbose and hit:
        print(f"    dropped {n0 - len(d)} feedstock catabolic-entry rxn(s): {hit}")
    return d


# ── shared CRISPRi gene filter (applied to every method) ───────────────────────
def crispri_gene_filter(df, meta, reaction_col="reaction", collapse_by_gene=True, verbose=True):
    """Keep only CRISPRi-actionable single targets:
        (a) has a gene, (b) not 'gapfilled', (c) not transport/exchange (XR*/EX_),
        (d) single-enzyme (GPR has no 'or'/isozymes);
    then optionally collapse duplicate reactions sharing a gene set (keep best-scoring)."""
    if df is None or df.empty:
        return df
    d = df.copy(); n0 = len(d)
    genes = d[reaction_col].map(lambda r: meta["genes"].get(r, ""))
    gpr   = d[reaction_col].map(lambda r: meta["gpr"].get(r, ""))
    keep = (genes.str.len() > 0) & (genes.str.lower() != "gapfilled") \
         & (~d[reaction_col].astype(str).str.match(r"^(XR|EX_)")) \
         & (~gpr.str.lower().str.contains(r"\bor\b", regex=True, na=False))
    d = d[keep].copy()
    if collapse_by_gene and not d.empty:
        d["_gs"] = d[reaction_col].map(lambda r: meta["genes"].get(r, ""))
        score_cols = [c for c in d.columns if c.endswith("_score") or c == "score" or c.endswith("_gap")]
        if score_cols:
            d = d.sort_values(score_cols[0], ascending=False)
        d = d.drop_duplicates("_gs", keep="first").drop(columns="_gs")
    if verbose:
        print(f"    gene filter: {n0} -> {len(d)}")
    return d.reset_index(drop=True)


# ── WT reference + viability (90%-KD growth >= 10% WT) ──────────────────────────
def wt_growth_reference(carbon="coumarate"):
    """pFBA growth-state flux distribution + mu_max (the knockdown reference)."""
    m = build_growth_state(carbon)
    sol = pfba(m)
    if sol.status != "optimal":
        raise RuntimeError(f"WT growth state not optimal on {carbon}")
    return sol.fluxes, float(sol.fluxes[BIOMASS])


def kd_growth_frac(rxn_id, growth_fluxes, mu_max, carbon="coumarate"):
    """Growth retained when this single reaction is held to 10% of its WT growth-state flux."""
    ref = float(growth_fluxes.get(rxn_id, 0.0))
    if abs(ref) < MIN_FLUX:
        return 1.0                       # carries ~no WT flux -> KD can't reduce growth
    m = build_growth_state(carbon)
    with m:
        apply_repression(m.reactions.get_by_id(rxn_id), ref, CRISPRI_REMAINING)
        return safe_slim(m) / mu_max if mu_max > 0 else 0.0


def viability_filter(df, growth_fluxes, mu_max, carbon="coumarate"):
    """Keep targets whose 90%-KD growth stays >= MIN_GROWTH_FRAC * mu_max. Essential genes survive as
    valid PARTIAL-KD targets (a full-KO gate would wrongly drop them)."""
    if df is None or df.empty:
        return df
    df = df.copy()
    df["kd_growth_frac"] = [kd_growth_frac(r, growth_fluxes, mu_max, carbon) for r in df["reaction"]]
    return df[df.kd_growth_frac >= MIN_GROWTH_FRAC].reset_index(drop=True)
