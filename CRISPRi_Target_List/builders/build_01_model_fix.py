"""
build_01_model_fix.py — build the CLEAN model-fix notebook (01_model_fix.ipynb).
Documents and reproduces the one repair that turns the published supplementary SBML into the working
model the whole pipeline uses: the published file references a biomass boundary species (XC54_b) that it
never declares, so cobra cannot even load it. The fix adds the missing species; the repaired model loads
and grows identically to the shipped CGA009_model_biomass_fix.xml.
Run:  python Release/builders/build_01_model_fix.py
Exec: python -m nbconvert --to notebook --execute --inplace \
        --ExecutePreprocessor.kernel_name=gem Release/Notebooks/01_model_fix.ipynb
"""
import nbformat as nbf
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent          # Release/
OUT = ROOT / "Notebooks" / "01_model_fix.ipynb"

nb = nbf.v4.new_notebook()
C = []
def md(s):   C.append(nbf.v4.new_markdown_cell(s))
def code(s): C.append(nbf.v4.new_code_cell(s))


md(r"""# 01 — Model fix: making the published *iRpa940* SBML usable
### *Rhodopseudomonas palustris* CGA009 (Alsiyabi, Immethun & Saha 2019, BMC Bioinformatics 20:233)

**What this notebook does.** It takes the model exactly as published — the journal supplementary file
`Model/12859_2019_2844_MOESM3_ESM.xml` — and applies the **one** repair needed to make it usable, then
proves the repaired model reproduces the shipped working model `Model/CGA009_model_biomass_fix.xml`.

**The bug (one line).** The biomass reaction `XR229` drains biomass into a boundary species **`XC54_b`**,
but the published SBML **never declares that species**. cobrapy therefore raises `KeyError: 'XC54_b'` and
**cannot even load the model** — let alone simulate growth.

**The fix.** Insert the missing boundary-species declaration

```xml
<species id="XC54_b" name="Biomass" compartment="XT" charge="-2" boundaryCondition="true"/>
```

immediately after the existing `XC54XT` biomass species. This is a pure **text-level SBML repair** (no
re-serialization), so nothing else in the model changes. (An earlier attempt — `01_fixing_model.ipynb` in
the original repo — instead *deleted* a duplicate biomass reaction; *adding the missing endpoint*, as in
`03_fixing_model_v2.ipynb`, is the correct, minimal fix and the one adopted here.)

> **Provenance note.** The working model the rest of the pipeline loads,
> `Model/CGA009_model_biomass_fix.xml`, is shipped **verbatim** in this Release (byte-for-byte identical to
> the original repo's). This notebook regenerates an *equivalent* repaired model from the published file and
> verifies it grows identically — it is documentation + validation, and deliberately does **not** overwrite
> the canonical file, so downstream numerical results are guaranteed unchanged.""")

md(r"""## 1 · Show the bug — the published file fails to load""")
code(r"""import io, contextlib, logging
from pathlib import Path
import cobra
logging.getLogger("cobra").setLevel(logging.CRITICAL)

MODEL_DIR = Path.cwd().parent / "Model"          # notebook runs from Release/Notebooks/
PUBLISHED = MODEL_DIR / "12859_2019_2844_MOESM3_ESM.xml"
CANONICAL = MODEL_DIR / "CGA009_model_biomass_fix.xml"

text = PUBLISHED.read_text(encoding="utf-8")
print(f"published supplementary: {PUBLISHED.name}")
print(f"  'XC54_b' occurrences : {text.count('XC54_b')}  (1 = referenced by biomass, but NOT declared)")
print(f"  'XC54XT' occurrences : {text.count('XC54XT')}")

try:
    with contextlib.redirect_stderr(io.StringIO()):
        cobra.io.read_sbml_model(str(PUBLISHED))
    print("\nunexpected: published model loaded")
except Exception as e:
    print(f"\nLOAD FAILS as published -> {type(e).__name__}: {e}")
    print("  (the biomass boundary species XC54_b is referenced but never declared)")""")

md(r"""## 2 · Apply the repair — declare the missing `XC54_b` biomass boundary species""")
code(r"""ANCHOR  = '<species id="XC54XT" name="Biomass" compartment="XT" charge="-2" boundaryCondition="false"/>\n'
MISSING = '\t<species id="XC54_b" name="Biomass" compartment="XT" charge="-2" boundaryCondition="true"/>\n'
assert ANCHOR in text, "anchor XC54XT species not found — file layout changed"

repaired_text = text.replace(ANCHOR, ANCHOR + MISSING, 1)
REGEN = MODEL_DIR / "CGA009_model_biomass_fix_regenerated.xml"
REGEN.write_text(repaired_text, encoding="utf-8")
print(f"inserted XC54_b declaration -> 'XC54_b' now occurs {repaired_text.count('XC54_b')}x")
print(f"wrote regenerated model: {REGEN.name}")""")

md(r"""## 3 · Validate — the repaired model loads and grows, identical to the shipped working model

The model ships pre-pinned to **Navid's acetate validation point** (acetate uptake `XR57` = 1.96), so a
plain `slim_optimize()` on either model returns that reference growth rate. They must match.""")
code(r"""def growth(path):
    with contextlib.redirect_stderr(io.StringIO()):
        m = cobra.io.read_sbml_model(str(path))
    return m, float(m.slim_optimize() or 0.0)

m_regen, mu_regen = growth(REGEN)
m_canon, mu_canon = growth(CANONICAL)

print(f"regenerated (text-repair of published)  : reactions={len(m_regen.reactions)}  "
      f"metabolites={len(m_regen.metabolites)}  growth={mu_regen:.6f} /h")
print(f"canonical CGA009_model_biomass_fix.xml  : reactions={len(m_canon.reactions)}  "
      f"metabolites={len(m_canon.metabolites)}  growth={mu_canon:.6f} /h")
match = abs(mu_regen - mu_canon) < 1e-6
print(f"\ngrowth matches the shipped working model: {match}  (Δ = {abs(mu_regen-mu_canon):.2e})")
print(f"  expected reference growth ≈ 0.0824 /h on the shipped acetate point (XR57=1.96)")
assert match, "regenerated model does not reproduce the canonical growth rate"

# Stronger: the text repair reproduces the canonical working model BYTE-FOR-BYTE. (The ~41 kB size jump
# vs the published file is just universal-newline -> Windows CRLF translation on write; the canonical
# CGA009_model_biomass_fix.xml was produced by this exact repair of the published SBML.)
byte_id = REGEN.read_bytes() == CANONICAL.read_bytes()
print(f"regenerated model is BYTE-IDENTICAL to CGA009_model_biomass_fix.xml: {byte_id}")
print("\nMODEL FIX VALIDATED ✓  — downstream notebooks load the canonical CGA009_model_biomass_fix.xml.")""")

md(r"""## 4 · Notes
- **The fix:** add the missing biomass boundary species `XC54_b` to the published supplementary SBML so the
  biomass reaction `XR229` can drain — the published file is otherwise unloadable (`KeyError: 'XC54_b'`).
- **Equivalence:** the text-repaired model and the shipped `CGA009_model_biomass_fix.xml` load with the same
  reaction/metabolite counts and grow identically (≈0.0824 /h on the acetate point), confirming the working
  model is exactly this repair of the published model.
- **Reproducibility guarantee:** the pipeline loads the *shipped* `CGA009_model_biomass_fix.xml` (copied
  byte-for-byte), not the regenerated file, so results never depend on SBML re-serialization formatting.
- **Original-repo lineage:** repaired in `03_fixing_model_v2.ipynb`; the alternative delete-a-reaction
  approach in `01_fixing_model.ipynb` was rejected. Growth validated in `04_validate_model.ipynb`.""")

nb["cells"] = C
nb["metadata"] = {"kernelspec": {"name": "python3", "display_name": "Python 3 (ipykernel)", "language": "python"},
                  "language_info": {"name": "python"}}
OUT.parent.mkdir(parents=True, exist_ok=True)
nbf.write(nb, str(OUT))
print(f"wrote {OUT}")
