import streamlit as st
import plotly.graph_objects as go
from backend import (
    apply_custom_css, extract_features, predict_pathogenicity,
    predict_thermodynamic_stability, fetch_string_ppi, compute_composite_risk
)

st.set_page_config(page_title="Step 8 | ProtMind AI", page_icon="🧠")
apply_custom_css()

st.markdown("<h2>🧠 Step 8: Explainable Intelligence Summary</h2>", unsafe_allow_html=True)

if not st.session_state.get('payload'):
    st.warning("⚠️ No active protein data found. Please complete Step 1: User Input first.")
    st.stop()

payload = st.session_state['payload']
results = st.session_state.get('results', {})

# Reuse results already computed in earlier steps where available (single source of truth);
# fall back to recomputing only what's missing so this page still works if visited out of order.
if 'features' not in results:
    st.info("ℹ️ Steps 3, 4, 6 and 7 haven't all been run yet — computing the missing pieces now. "
            "For the most consistent numbers, visit them in order first.")

features = results.get('features') or extract_features(payload['sequence'], payload['mutation'])
pathogenicity = results.get('pathogenicity') or predict_pathogenicity(payload['uniprot_id'], payload['mutation'])
stability = results.get('stability') or predict_thermodynamic_stability(
    features['Wildtype AA'], features['Mutant AA'], features['WT Property'], features['Mutant Property']
)
ppi = results.get('ppi')
if ppi is None:
    ppi = fetch_string_ppi(payload['uniprot_id'])
angles = results.get('angles')  # None means Step 7 hasn't run — handled as "unknown", not "valid"

composite = compute_composite_risk(stability, pathogenicity, len(ppi), angles is not None)
st.session_state.setdefault('results', {})
st.session_state['results']['composite'] = composite

col1, col2 = st.columns([1, 1])

with col1:
    st.markdown("<h3>Composite Risk Score</h3>", unsafe_allow_html=True)
    gauge = go.Figure(go.Indicator(
        mode="gauge+number",
        value=composite["composite_score"],
        number={'suffix': "/100"},
        title={'text': composite["tier"]},
        gauge={
            'axis': {'range': [0, 100]},
            'bar': {'color': "#00F5D4"},
            'steps': [
                {'range': [0, 35], 'color': "#1b3a2f"},
                {'range': [35, 65], 'color': "#4a3a12"},
                {'range': [65, 100], 'color': "#4a1414"},
            ],
        }
    ))
    gauge.update_layout(height=280, margin=dict(l=20, r=20, t=60, b=10),
                         paper_bgcolor="rgba(0,0,0,0)", font={'color': "white"})
    st.plotly_chart(gauge, use_container_width=True)

with col2:
    st.markdown("<h3>Risk Component Breakdown</h3>", unsafe_allow_html=True)
    comp = composite["components"]
    radar = go.Figure(go.Scatterpolar(
        r=list(comp.values()) + [list(comp.values())[0]],
        theta=list(comp.keys()) + [list(comp.keys())[0]],
        fill='toself',
        line_color="#9B5DE5"
    ))
    radar.update_layout(
        polar=dict(radialaxis=dict(visible=True, range=[0, 100], color="white"),
                   angularaxis=dict(color="white")),
        showlegend=False, height=280, margin=dict(l=40, r=40, t=30, b=30),
        paper_bgcolor="rgba(0,0,0,0)", font={'color': "white"}
    )
    st.plotly_chart(radar, use_container_width=True)

st.divider()
st.markdown("### 📋 Automated Pathogenicity Narrative")

geometry_note = (
    "Backbone geometry has not been validated yet (Step 7 not run)." if angles is None else
    "Backbone torsion angles fall within expected conformational space." if angles else
    "⚠️ Backbone torsion angles could not be resolved at this position."
)

st.info(
    f"The mutation **{payload['mutation']}** in target **{payload['uniprot_id']}** produced a composite risk "
    f"score of **{composite['composite_score']}/100 ({composite['tier']})**, combining thermodynamic, "
    f"evolutionary, and network-level evidence:"
)

st.markdown(f"""
* **Evolutionary Context:** The substitution shows a BLOSUM62 score of **{stability['blosum62']}** and a
  property shift of *{features['Property Shift']}* — {"a substantial jump in chemical class, which is generally more disruptive to folding." if features['WT Property'] != features['Mutant Property'] else "a change within the same chemical class, which is generally better tolerated."}
* **Structural Physics:** Estimated ΔΔG is **{stability['ddG']} kcal/mol** ({stability['status']}), driven by a
  hydrophobicity shift of {stability['hydrophobicity_delta']} (Kyte-Doolittle) and a volume shift of
  {stability['volume_delta']} Å³ (Zamyatnin). {geometry_note}
* **Systemic Impact:** {len(ppi)} candidate interaction partner(s) were found in the STRING database, so a
  destabilizing change at this residue could plausibly propagate to {"several connected pathways." if len(ppi) >= 3 else "a limited part of the interaction network." if ppi else "no currently known partners — network impact looks limited."}
* **Evidence Basis:** {pathogenicity['basis']}
""")

if composite["tier"] == "High Concern":
    st.error("🔴 **Overall Signal:** Multiple independent lines of evidence point toward a disruptive mutation. "
              "Recommend prioritizing this variant for structural/experimental follow-up.")
elif composite["tier"] == "Moderate Concern":
    st.warning("🟠 **Overall Signal:** Mixed evidence — some destabilizing signal, but not conclusive across all "
               "dimensions. Treat as a candidate for further review rather than a confirmed pathogenic call.")
else:
    st.success("🟢 **Overall Signal:** Current evidence suggests this mutation is likely well-tolerated.")

st.caption("This summary combines empirical/statistical estimates from public databases and established "
           "biochemical scales. It is a research and educational aid, not a clinical diagnostic.")
