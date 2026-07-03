"""
nb71_lib.py — shared engine for NB71 (downselection of the combined FBA+AI CRISPRi target list).
Reproduces the NB70 combined table (FBA finalized-30 + 27 locus-reconciled AI targets, gene-level FBA method
cells) and exposes it as a dataframe + a heatmap renderer, so NB71 can just: load -> draw full -> cut -> draw
downselected. Depends only on nb60_lib (read-only) and the NB60-64 result CSVs. Does NOT modify any notebook.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as _LSC
sys.path.insert(0, str(Path(__file__).resolve().parent))
import nb60_lib as L

def _seq(name, cols):  # build a light->dark sequential colormap from a list of hex stops
    return _LSC.from_list(name, cols)

# ── the 27 authoritative AI targets (symbol, ai_locus, function, ai_module, growth, score, model_rxn,
#    curated model subunit ("" -> derive representative from GPR), mechanism, citation) ─────────────────
NAN = float("nan")
AI_ROWS = [
 ("phaZ", "RPA1786", "PHB depolymerase",                       "phb",   "all anaerobic phototrophic",            1.00, "",       "",        "KD blocks PHB remobilization -> raises net granule accumulation (titer, not flux)", "Kobayashi 2019 (R. sphaeroides phaZ del -> 2.9x PHB); Jendrossek 2002"),
 ("glnA", "RPA4182", "Glutamine synthetase I (GS)",            "nreg",  "N-limited phototrophic, p-coumarate",   0.92, "R00253", "RPA2967", "partial repression lowers NH4+ assimilation -> high C/N -> PHB overflow storage", "Ranaivoarisoa/spectrum 2023 (GlnA1 primary NH4 assimilation); McKinlay & Harwood 2010"),
 ("gltA", "RPA4656", "Citrate synthase (TCA entry)",           "accoa", "phototrophic + p-coumarate",            0.88, "R00351", "RPA2907", "KD of first committed TCA step raises cytoplasmic acetyl-CoA for PhaA (partial only)", "Lin 2021 (H. mediterranei +76% PHBV); Halomonas TD01 +8%"),
 ("accA", "RPA3510", "Acetyl-CoA carboxylase alpha (AccAD)",   "accoa", "phototrophic + p-coumarate",            0.82, "R04386", "RPA0071", "blocks first committed/rate-limiting FAS step -> frees acetyl-CoA for PHB", "Davis 2000; Nomura 2004 Metab Eng"),
 ("fabH", "RPA2368", "beta-Ketoacyl-ACP synthase III (FAS init)", "accoa", "phototrophic + p-coumarate",         0.78, "R04960", "RPA0426", "FAS initiation; downstream of accA (redundant, malonyl-CoA still made)", "Heath & Rock 1996; Zheng 2004"),
 ("glnB", "RPA0454", "PII nitrogen sensor (GlnB1)",            "nreg",  "N-limited phototrophic",                0.75, "",       "",        "PII sensor; one of 3 redundant paralogs (GlnB/GlnK1/GlnK2) that buffer each other", "Connolly 2006; Arcondeguy 2001"),
 ("ackA", "RPA3344", "Acetate kinase",                         "accoa", "phototrophic + p-coumarate",            0.68, "R00315", "RPA4566", "2nd step of pta-ackA overflow; downstream of pta; dAckA accumulates toxic acetyl-P", "Kim 2010 (ackA instability/acetyl-P); Castano-Cerezo 2009"),
 ("pta",  "RPA3343", "Phosphotransacetylase",                  "accoa", "phototrophic + p-coumarate",            0.65, "R00921", "RPA4567", "raises acetyl-CoA/CoA ratio; 'key factor' for PHB (upstream overflow valve)", "PTA-as-key-factor review; Alsiyabi 2022 (TFA)"),
 ("icd",  "RPA2989", "Isocitrate dehydrogenase",               "tca",   "anaerobic phototrophic",               0.62, "R01899", "RPA3834", "partial KD slows oxidative TCA + lowers 2-OG (secondary N-starve mimic)", "Sunnarborg 1990; Stephanopoulos & Vallino 1991"),
 ("ntrB", "RPA1711", "Histidine kinase NtrB (N two-component)", "nreg", "N-limited phototrophic",                0.58, "",       "",        "NtrB kinase; signals through the same NtrC regulon N-starvation triggers (indirect)", "Ninfa & Magasanik 1986"),
 ("gltB", "RPA3150", "Glutamate synthase large subunit (GOGAT)", "nreg", "N-limited phototrophic + p-coumarate", 0.55, "R00114", "RPA0891", "second N-assimilation enzyme; redundant with glnA (GDH can bypass GOGAT)", "van Heeswijk 2013"),
 ("sucA", "RPA1538", "2-Oxoglutarate dehydrogenase E1 (OGDH)", "tca",   "phototrophic + p-coumarate",            0.52, "R03316", "RPA0189", "blocks 2-OG->succinyl-CoA; diverts to glutamate (N mimic) + cuts TCA drain", "McKinlay 2014; Carlson & Srienc 2004"),
 ("phaR", "RPA1795", "PHB granule transcriptional repressor",  "phb",   "all anaerobic phototrophic",            0.50, "",       "",        "derepresses phaC/phaA/phaB; but also derepresses phaZ depolymerase (mixed)", "Yamada 2002; Stubbe & Tian 2003"),
 ("glnK", "RPA0584", "PII-like nitrogen sensor (GlnK)",        "nreg",  "N-limited phototrophic",                0.48, "",       "",        "GlnB paralog; GlnK1 predominates under N-excess; redundant PII sensor", "Detsch & Stuelke 2003; Connolly 2006"),
 ("sdhA", "RPA0957", "Succinate dehydrogenase flavoprotein",   "tca",   "anaerobic phototrophic",               0.45, "R00408", "RPA0216", "KD slows succinate->fumarate; shunts carbon to PHB as alt electron sink", "Sauer 1999"),
 ("fabF", "RPA2790", "beta-Ketoacyl-ACP synthase II (FAS elong)", "accoa", "phototrophic + p-coumarate",         0.42, "R04960", "RPA2019", "FAS elongation; downstream of accA (redundant)", "Magnuson 1993; Lim 2011"),
 ("acnA", "RPA0724", "Aconitase A (TCA)",                      "tca",   "anaerobic phototrophic",               0.38, "R01900", "RPA0202", "citrate->isocitrate; slows TCA cycling; synergistic with gltA KD", "Jordan 1999"),
 ("regB", "RPA3332", "Redox-sensing histidine kinase RegB",    "redox", "anaerobic phototrophic",               0.35, "",       "",        "global redox two-component regulator (>50 targets); pleiotropic, high-risk", "Elsen 2004; Laguri 2003"),
 ("prpC", "RPA2576", "2-Methylcitrate synthase",               "accoa", "phototrophic + p-coumarate",            0.32, "R00931", "RPA2394", "diverts propionyl-CoA (from aromatic catabolism) to methylcitrate; KD retains CoA", "Horswill & Escalante-Semerena 1999"),
 ("ppc",  "RPA1702", "PEP carboxylase (anaplerosis)",          "tca",   "phototrophic + p-coumarate",            0.28, "R00345", "RPA1772", "diverts PEP to OAA; partial KD slightly raises acetyl-CoA (modest)", "Sauer & Eikmanns 2005"),
 ("fumC", "RPA3151", "Fumarate hydratase C (TCA)",             "tca",   "anaerobic phototrophic",               0.25, "R01082", "RPA1329", "fumarate->malate; cuts downstream TCA flux + NADH", "iAN1128 FBA"),
 ("hbaA", "RPA4616", "4-Hydroxybenzoyl-CoA reductase A",       "other", "phototrophic + p-coumarate",            0.22, "R05316", "RPA0670", "aromatic ring reduction on the MAIN p-coumarate catabolic route -> KD cuts carbon SUPPLY", "Egland 1997; Gibson & Harwood 2002"),
 ("nuoF", "RPA2018", "NADH:ubiquinone oxidoreductase F (Cx I)", "redox", "phototrophic (excess reductant)",      0.18, "XR14",   "",        "Complex I; reducing it raises NADH/NAD+ toward PHB but risks lethality (obligate CBB sink)", "Brandt 2006; Wang 2010 (CBB obligate)"),
 ("nifA", "RPA4630", "Nitrogenase transcriptional activator",  "nreg",  "N2-fixing phototrophic",               0.15, "",       "",        "nitrogenase activator; relevant only under diazotrophic (N2-fixing) growth", "Dixon & Kahn 2004"),
 ("mdh",  "RPA3406", "Malate dehydrogenase",                   "tca",   "anaerobic phototrophic",               0.12, "R00342", "RPA0192", "malate->OAA; slightly reduces OAA+acetyl-CoA condensation (less citrate)", "Gourdon 2000"),
 ("pntA", "RPA3355", "NAD(P)+ transhydrogenase alpha",         "redox", "anaerobic phototrophic",               0.10, "XR13",   "RPA4180", "tunes NADPH/NADH for PhaB; direction-dependent, uncertain; redox risk", "Sauer 2004; Chin 2011"),
 ("atpA", "RPA4709", "ATP synthase alpha (F1F0)",              "nuc",   "anaerobic phototrophic",               0.08, "R00086", "RPA0175", "energy-stress lever; pleiotropic, lethality risk (ATP burst); late-DBTL", "Koebmann 2002"),
]
COLNAMES = ["symbol", "ai_locus", "name", "ai_module", "growth_condition", "ai_pred",
            "model_rxn", "model_locus_cur", "mechanism", "citation"]

MODULES = [("redox", ("Redox: ferredoxin /", "NAD(P)H balance")), ("backbone", ("Central carbon backbone", "(glycolysis / PEP)")),
           ("accoa", ("Acetyl-CoA / FAS", "conservation & PHB precursor")), ("ppp", ("Pentose phosphate", "(NADPH & precursors)")),
           ("tca", ("TCA cycle throttling", "")), ("nreg", ("N-starvation mimicry", "(AI: PII / GS / GOGAT / Ntr)")),
           ("phb", ("PHB maintenance", "(AI: depolymerase / PhaR)")), ("aa", ("Amino-acid biosynthesis", "drains")),
           ("nuc", ("Nucleotide / energy", "charge")), ("other", ("Other / various", ""))]
MOD_COLOR = {"redox": "#3182bd", "backbone": "#2ca25f", "accoa": "#756bb1", "ppp": "#cc8b00",
             "tca": "#d9694a", "nreg": "#1f9e89", "phb": "#d62a7a", "aa": "#2a9d8f", "nuc": "#b5539c", "other": "#8c8c8c"}

# ── COLOUR PALETTES ────────────────────────────────────────────────────────────────────────────────────
# Each palette = a sequential cmap for the FBA (cool) columns, one for the AI (warm) column, and a cohesive
# 10-colour categorical set for the module side-bars. Colour-theory notes per palette below. Switch with the
# `palette=` arg to heatmap() (default = ACTIVE_PALETTE). "original" is LOCKED as the reference version.
PALETTES = {
  # LOCKED reference — the version we can always come back to (pure Blues/Oranges complements + rainbow bars).
  "original": dict(fba=plt.cm.Blues, ai=plt.cm.Oranges, mods=dict(MOD_COLOR)),

  # Muted split-complementary: teal (FBA) vs terracotta (AI) are lowered in chroma so they harmonise instead
  # of vibrating; module bars are one tonally-uniform "dusty" family (mid saturation/lightness) with nreg/redox
  # echoing the teal and tca echoing the terracotta, so the side-bars belong to the same world as the cells.
  "teal-clay": dict(
     fba=_seq("fba_teal", ["#f3f8f7", "#c3ddd8", "#7cb8b0", "#3d938b", "#1b6f67"]),
     ai =_seq("ai_clay",  ["#fbf3ec", "#f2cfb3", "#e0a377", "#cc7a4c", "#a4512a"]),
     mods=dict(redox="#5b86a6", backbone="#7a9e63", accoa="#8f6f9c", ppp="#c1a052", tca="#cc7a5c",
               nreg="#3d938b", phb="#c286a4", aa="#94a05e", nuc="#a5789e", other="#9b968c")),

  # Cool/analogous with a single warm accent: slate-indigo (FBA) vs soft gold/amber (AI); jewel-muted bars.
  "slate-amber": dict(
     fba=_seq("fba_slate", ["#f1f3f7", "#c6cfe0", "#8e9cc4", "#5a6ea0", "#38497a"]),
     ai =_seq("ai_amber",  ["#fbf6e9", "#efdca6", "#dcb35d", "#bf8d2f", "#93691c"]),
     mods=dict(redox="#4f6b9e", backbone="#4a8c6f", accoa="#7d5ba6", ppp="#c39a3d", tca="#c15e52",
               nreg="#3f8f88", phb="#a95585", aa="#6f8f52", nuc="#8a5f9e", other="#8c8c8c")),

  # University of Washington brand: Husky Purple (#4b2e83, FBA) vs Husky Gold (#b7a57a -> saturated, AI).
  # Ramps run pale -> FULL-saturation UW colour at the top of the scale; module bars are a cohesive
  # purple/gold-family set (tints & shades of the two brand colours + neutral grey).
  "uw": dict(
     fba=_seq("fba_uw",   ["#f4f1f9", "#cabbe3", "#9a7fc6", "#6a47a0", "#4b2e83"]),
     ai =_seq("ai_uw",    ["#faf6ea", "#e8d6a3", "#cdb265", "#b7a57a", "#8a6f28"]),
     mods=dict(redox="#4b2e83", backbone="#7d5ba6", accoa="#a58fc9", ppp="#b7a57a", tca="#85754d",
               nreg="#8a5a9c", phb="#c9b78c", aa="#615a86", nuc="#a37f9e", other="#8c8c8c")),

  # Desaturated "editorial": steel blue-grey (FBA) vs rust (AI), earthy bars — lowest-contrast / calmest.
  "steel-rust": dict(
     fba=_seq("fba_steel", ["#f2f4f5", "#c8d2d6", "#90a4ac", "#5b7681", "#3a4f57"]),
     ai =_seq("ai_rust",   ["#faf1ec", "#eec6ae", "#d99a72", "#c06a42", "#9a4926"]),
     mods=dict(redox="#5f7d8a", backbone="#6d9068", accoa="#8a6f93", ppp="#bd9a56", tca="#c07156",
               nreg="#4c9088", phb="#b57e97", aa="#88976a", nuc="#94749a", other="#98938b")),
}
ACTIVE_PALETTE = "uw"            # <- default; NB71 can override per-call via heatmap(..., palette="...")

SRC_COLOR = {"FBA": "#222", "AI": "#7d3c98", "FBA+AI": "#1a7d3c", "AI+FBA": "#5a9e2f"}
COLS = ["FVSEOF", "FluxRETAP", "CASOP"]
AA_SUBS = {"glycine, serine and threonine metabolism", "valine, leucine and isoleucine biosynthesis",
           "valine, leucine and isoleucine degradation", "lysine biosynthesis", "lysine degradation",
           "histidine metabolism", "methionine metabolism", "cysteine and methionine metabolism",
           "phenylalanine, tyrosine and tryptophan biosynthesis", "arginine and proline metabolism",
           "alanine, aspartate and glutamate metabolism", "glutamate metabolism", "beta-alanine metabolism",
           "tyrosine metabolism", "tryptophan metabolism", "urea cycle and metabolism of amino groups",
           "aminoacyl trna biosynthesis", "selenocompound metabolism"}

_STATE = {}

def get_data():
    """Build the combined 51-row table + gene-level FBA method maps. Returns (rows, ai)."""
    meta = L.metadata()
    NAME = lambda r: meta["name"].get(r, ""); SUB = lambda r: meta["subsys"].get(r, "")
    master = (pd.read_csv(L.OUT / "nb64_master_target_table.csv")
              .sort_values("consensus_score", ascending=False).drop_duplicates("gene").reset_index(drop=True))
    SEL30 = set(pd.read_csv(L.OUT / "nb64_round0_selection.csv").gene)
    FINAL30 = (SEL30 - {"RPA0164"}) | {"RPA4394"}

    METHODS = {"FVSEOF": ("nb60_fvseof_targets.csv", "q_slope"),
               "FluxRETAP": ("nb61_fluxretap_targets.csv", "flux_diff"),
               "CASOP": ("nb63_casop_targets.csv", "casop_ko_score")}
    gmraw, gmmag, gene_methods = {}, {}, {}
    for m, (f, col) in METHODS.items():
        d = pd.read_csv(L.OUT / f); best = {}
        for _, r in d.iterrows():
            for g in str(r.genes).split("|"):
                if g and (g not in best or abs(r[col]) > abs(best[g])): best[g] = r[col]
                if g: gene_methods.setdefault(g, set()).add(m)
        gmraw[m] = best
        av = {g: abs(v) for g, v in best.items()}; lo, hi = min(av.values()), max(av.values())
        gmmag[m] = {g: (a - lo) / (hi - lo) if hi > lo else 1.0 for g, a in av.items()}
    FBA_FLAGGED = set(gene_methods)

    def module_of(r):
        nm = (NAME(r) or "").lower(); ss = (SUB(r) or "").lower()
        if "ferredoxin" in nm: return "redox"
        if any(k in nm for k in ["acetate phosphotransferase", "phosphate acetyltransferase", "citrate oxaloacetate-lyase",
                "malate glyoxylate-lyase", "hydroxybutanoyl", "carbon-dioxide ligase", "carboxyl-carrier"]): return "accoa"
        if ss == "pentose phosphate pathway" or any(k in nm for k in ["ribulose-5-phosphate", "ribose-5-phosphate",
                "6-phospho-d-gluconate", "glucono", "transketolase", "xylulose"]): return "ppp"
        if ss == "glycolysis / gluconeogenesis" or any(k in nm for k in ["2-phospho-d-glycerate hydro-lyase",
                "pyruvate,orthophosphate", "phosphoenolpyruvate", "glucose-6-phosphate ketol-isomerase", "3-phospho-d-glyceroyl"]): return "backbone"
        if ss in ("citrate cycle (tca cycle)", "glyoxylate and dicarboxylate metabolism") or any(k in nm for k in
                ["isocitrate", "(s)-malate", "fumarate", "succinate", "2-oxoglutarate", "citrate hydro-lyase"]): return "tca"
        if ss in AA_SUBS or "trna" in nm or "amidotransferase" in nm: return "aa"
        if ss in ("purine metabolism", "pyrimidine metabolism") or any(k in nm for k in
                ["amp phosphotransferase", "nucleoside-diphosphate", "pyrophosphate phosphohydrolase"]): return "nuc"
        return "other"

    m = L.read_model()
    rxn_genes = {r.id: sorted(g.id for g in r.genes) for r in m.reactions}
    rxn_name = {r.id: r.name for r in m.reactions}
    ai = pd.DataFrame(AI_ROWS, columns=COLNAMES); ai["ai_pred"] = ai["ai_pred"].astype(float)
    def reconcile(row):
        rx, cur = row.model_rxn, row.model_locus_cur
        if not rx:
            return pd.Series({"canonical_gene": row.ai_locus.split("/")[0], "in_model": False, "n_gpr": 0, "model_rxn_name": ""})
        gpr = rxn_genes[rx]; rep = cur if cur else gpr[0]
        return pd.Series({"canonical_gene": rep, "in_model": True, "n_gpr": len(gpr), "model_rxn_name": rxn_name.get(rx, "")})
    ai = pd.concat([ai, ai.apply(reconcile, axis=1)], axis=1)
    ai["model_locus_disp"] = [g + (f" (+{n-1} isoz)" if n > 1 else "") for g, n in zip(ai.canonical_gene, ai.n_gpr)]

    fba = master[master.gene.isin(FINAL30)].copy()
    fba["module"] = fba.reaction.map(module_of)
    fba["source"] = "FBA"; fba["in_model"] = True; fba["ai_flag"] = False
    fba["ai_pred"] = np.nan; fba["ai_symbol"] = ""; fba["growth_condition"] = ""; fba["mechanism"] = ""; fba["citation"] = ""
    fba["blue_gene"] = fba["gene"]
    fba["label"] = [f"{r.gene} · {str(r['name'])[:30]}" for _, r in fba.iterrows()]
    fba = fba.set_index("gene", drop=False)
    new_rows = []
    for _, a in ai.iterrows():
        gpr = rxn_genes.get(a.model_rxn, []) if a.in_model else []
        hit = set(gpr) & set(fba.index)
        if hit:
            cg = sorted(hit)[0]
            fba.loc[cg, "source"] = "FBA+AI"; fba.loc[cg, "ai_pred"] = a.ai_pred; fba.loc[cg, "ai_flag"] = True
            fba.loc[cg, "ai_symbol"] = a.symbol; fba.loc[cg, "mechanism"] = a.mechanism
            fba.loc[cg, "citation"] = a.citation; fba.loc[cg, "growth_condition"] = a.growth_condition
            fba.loc[cg, "label"] = f"{cg} · {str(fba.loc[cg, 'name'])[:24]}  «{a.symbol}»"
        else:
            flagged = [g for g in gpr if g in FBA_FLAGGED]
            blue_gene = flagged[0] if flagged else a.canonical_gene
            src = "AI+FBA" if flagged else "AI"
            key = a.model_rxn if a.in_model and a.model_rxn else f"AI_{a.symbol}"
            nm = a.model_rxn_name if a.in_model and a.model_rxn_name else a["name"]
            star = "" if a.in_model else " *"
            new_rows.append({"reaction": key, "gene": a.canonical_gene, "blue_gene": blue_gene, "name": nm,
                             "module": a.ai_module, "source": src, "in_model": a.in_model, "ai_flag": True,
                             "ai_pred": a.ai_pred, "ai_symbol": a.symbol, "growth_condition": a.growth_condition,
                             "mechanism": a.mechanism, "citation": a.citation,
                             "label": f"{a.symbol} ({a.model_locus_disp}{star}) · {str(a['name'])[:20]}"})
    keep = ["reaction", "gene", "blue_gene", "name", "module", "source", "in_model", "ai_flag", "ai_pred",
            "ai_symbol", "growth_condition", "mechanism", "citation", "label"]
    rows = pd.concat([fba[keep].reset_index(drop=True), pd.DataFrame(new_rows)[keep]], ignore_index=True)
    assert rows.gene.nunique() == len(rows), "double-counted gene!"
    rows["n_fba_methods"] = rows.blue_gene.map(lambda g: len(gene_methods.get(g, set())))
    rows["fba_methods"] = rows.blue_gene.map(lambda g: "|".join(sorted(gene_methods.get(g, set()))))
    _STATE.update(gmraw=gmraw, gmmag=gmmag)
    return rows, ai


def heatmap(df, title, fname, ax=None, palette=None):
    """Render the FBA+AI heatmap for the given row subset. Gene-level cool (FBA) cells + warm (AI) column.
    `palette` selects a colour scheme from PALETTES (default ACTIVE_PALETTE)."""
    pal = PALETTES[palette or ACTIVE_PALETTE]
    MODc = pal["mods"]
    gmraw, gmmag = _STATE["gmraw"], _STATE["gmmag"]
    df = df.copy()
    df["_bri"] = [max([gmmag[m].get(r.blue_gene, 0.0) for m in COLS] + [r.ai_pred if r.ai_pred == r.ai_pred else 0.0])
                  for _, r in df.iterrows()]
    blocks = [(k, lab, df[df.module == k].sort_values("_bri", ascending=False)) for k, lab in MODULES if (df.module == k).any()]
    rows = pd.concat([s for _, _, s in blocks]).reset_index(drop=True)
    n = len(rows); ncol = len(COLS) + 1
    V = np.full((n, ncol), np.nan); Cc = np.full((n, ncol), np.nan)
    for i, r in rows.iterrows():
        bg = r["blue_gene"]
        for j, mth in enumerate(COLS):
            if bg in gmraw[mth]: V[i, j] = gmraw[mth][bg]; Cc[i, j] = gmmag[mth][bg]
        if r["ai_pred"] == r["ai_pred"]: V[i, len(COLS)] = r["ai_pred"]; Cc[i, len(COLS)] = r["ai_pred"]
        elif r["ai_flag"]: Cc[i, len(COLS)] = 0.10
    own = ax is None
    if own: fig, ax = plt.subplots(figsize=(9.8, max(4.5, 0.33 * n + 2.0)))
    else: fig = ax.figure
    coolcm = pal["fba"].copy(); coolcm.set_bad((0, 0, 0, 0))
    ax.imshow(np.ma.masked_invalid(np.where(np.arange(ncol) < len(COLS), Cc, np.nan)), cmap=coolcm, vmin=0, vmax=1, aspect="auto")
    warmcm = pal["ai"].copy(); warmcm.set_bad((0, 0, 0, 0))
    so = np.full((n, ncol), np.nan); so[:, len(COLS)] = Cc[:, len(COLS)]
    ax.imshow(np.ma.masked_invalid(so), cmap=warmcm, vmin=0, vmax=1, aspect="auto")
    for i in range(n):
        for j in range(ncol):
            if not np.isnan(V[i, j]):
                ax.text(j, i, f"{V[i, j]:.2f}", ha="center", va="center", fontsize=8, color="white" if Cc[i, j] > 0.55 else "#222")
    ax.set_xticks(range(ncol)); ax.set_xticklabels(COLS + ["AI prediction\n(norm. score)"], fontsize=9)
    ax.set_yticks(range(n)); ax.set_yticklabels(rows.label, fontsize=6.8)
    for lbl, src in zip(ax.get_yticklabels(), rows.source):
        lbl.set_color(SRC_COLOR[src]); lbl.set_fontweight("bold" if src == "FBA+AI" else "normal")
    for x in range(1, ncol): ax.axvline(x - 0.5, color="white", lw=2.5)
    ax.axvline(len(COLS) - 0.5, color="#888", lw=2.5); ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    xbar = ncol - 0.5 + 0.16; xtxt = ncol - 0.5 + 0.42; y0 = 0
    for key, (t1, t2), sub in blocks:
        y1 = y0 + len(sub) - 1
        ax.plot([xbar, xbar], [y0 - 0.32, y1 + 0.32], color=MODc[key], lw=6.5, solid_capstyle="round", clip_on=False)
        ym = (y0 + y1) / 2.0
        ax.text(xtxt, ym - (0.18 if t2 else 0), t1, color=MODc[key], fontsize=9.4, fontweight="bold", va="center", ha="left", clip_on=False)
        if t2: ax.text(xtxt, ym + 0.42, t2, color="#6b6b6b", fontsize=8.0, va="center", ha="left", clip_on=False)
        ax.axhline(y1 + 0.5, color="#ededed", lw=1.0); y0 = y1 + 1
    ax.set_title(title, fontsize=12, fontweight="bold", pad=16)
    ax.text(0.0, 1.0 + 0.5 / max(n, 8),
            "blue = FBA within-method strength (gene-level) · orange = AI score (0-1) · black=FBA · "
            "green/bold=FBA+AI (top-30) · olive=AI+FBA (below top-30) · purple=AI-only · * = not in iAN1128",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6.6, style="italic", color="#777")
    if own:
        for ext in ("png", "svg"): fig.savefig(L.FIG / f"{fname}.{ext}", bbox_inches="tight", dpi=160, transparent=True)
        plt.show()
    print(f"rendered {fname} ({n} genes)")
    return fig
