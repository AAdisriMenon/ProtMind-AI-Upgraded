# ProtMind AI

**An Explainable Protein Mutation Intelligence & Therapeutic Discovery Platform**

ProtMind AI is a Streamlit web application that takes a single protein mutation (a UniProt ID + an
amino-acid substitution, e.g. `P04637` / `R175H`) and runs it through a ten-stage, explainable
analysis pipeline — combining live biological database queries, established biochemical scales, and
structural visualization into one coherent, shareable report.

🔗 **Live App:** https://protmind-ai-bqnbemvxmfpxvm2wr7ebpk.streamlit.app/
📦 **Repository:** https://github.com/AAdisriMenon/ProtMind-AI

---

## Pipeline Overview

```text
1. User Input                → validate protein ID, mutation, and sequence
2. Data Retrieval            → live UniProt + AlphaFold + known-variant lookup
3. Preprocessing             → physicochemical property comparison, sequence alignment
4. Mutation Analysis         → pathogenicity estimate (BLOSUM62 + ClinVar/UniProt record) + STRING PPI network
5. Structure Prediction      → interactive 3D AlphaFold structure viewer
6. Thermodynamic Stability   → ΔΔG estimate (Kyte-Doolittle + Zamyatnin + BLOSUM62)
7. Stereochemical Validation → real Phi/Psi backbone angle calculation + Ramachandran plot
8. Explainable AI            → composite risk score + component breakdown, built from steps 3–7
9. Therapeutic Targeting     → live ChEMBL query for existing drugs/mechanisms against the target
10. Clinical Report          → consolidated, downloadable PDF using the same session-wide numbers
```

All ten stages are implemented. The project's ongoing work is about **deepening the modeling layer**
(see *Methodology & Honesty Notes* below), not adding missing pages.

---

## Methodology & Honesty Notes

This project is intentionally built in two layers, and understanding the difference matters for how
its output should be interpreted:

**Real, live data (no simulation):**
- Protein identity and metadata — live UniProt REST API
- 3D structure — live AlphaFold API + interactive py3Dmol rendering
- Interaction partners — live STRING database API
- Known clinical/variant records — live UniProt/ClinVar variation API
- Backbone torsion angles — computed directly from real AlphaFold atomic coordinates (Biopython)
- Existing drug mechanisms — live ChEMBL API

**Computed estimates (transparent formulas, not full physics simulation):**
- **Pathogenicity estimate** — when no documented clinical record exists for the exact substitution
  (the common case for novel mutations), a BLOSUM62 substitution-likelihood score is used as a
  stand-in. This is a real, standard bioinformatics matrix, but it is *not* equivalent to running
  SIFT, PolyPhen-2, or AlphaMissense themselves — a natural next upgrade is calling the Ensembl VEP
  API for a genuine multi-tool consensus prediction.
- **Thermodynamic stability (ΔΔG)** — estimated from a weighted combination of the Kyte-Doolittle
  hydrophobicity scale, Zamyatnin residue volumes, and BLOSUM62 score, plus known structural-breaker
  residues (Pro/Gly/Cys). This is a legitimate lightweight approximation used ahead of full
  physics-based tools, but it is not a molecular dynamics simulation (e.g., FoldX, Rosetta ddG).
- **Composite risk score (Step 8)** — a documented, weighted combination of the above real/estimated
  values (40% stability, 35% pathogenicity, 15% network exposure, 10% backbone geometry). Every input
  is already computed elsewhere in the pipeline; nothing is a newly fabricated number.

Every page in the app that shows an estimate rather than a database fact says so explicitly.

---

## Tech Stack

| Layer | Tools |
|---|---|
| App framework | Streamlit (multi-page) |
| Bioinformatics | Biopython (alignment, structure parsing, BLOSUM62) |
| Structure visualization | py3Dmol, stmol |
| Charts | Plotly (gauges, radar, Ramachandran), Matplotlib (report figures) |
| Reporting | fpdf |
| External data | UniProt REST API, AlphaFold API, STRING API, ChEMBL API |

## Running Locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Future Scope

- Swap the BLOSUM62 pathogenicity fallback for a live Ensembl VEP / AlphaMissense API call
- Swap the empirical ΔΔG formula for an established structure-based tool (FoldX API, Rosetta)
- Map protein-level mutations to genomic coordinates to enable real gnomAD population-frequency lookup
- Batch mode for analyzing multiple mutations across a protein at once
- Cross-species conservation comparison
