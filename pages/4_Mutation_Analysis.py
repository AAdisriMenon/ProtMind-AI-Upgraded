import streamlit as st
import time
import pandas as pd
import plotly.graph_objects as go
from backend import apply_custom_css, fetch_string_ppi, predict_pathogenicity

st.set_page_config(page_title="Step 4 | ProtMind AI", page_icon="🧬")
apply_custom_css()

st.markdown("<h2>🧬 Step 4: Mutation Analysis & Functional Proteomics</h2>", unsafe_allow_html=True)

if not st.session_state.get('payload'):
    st.warning("⚠️ No active protein data found. Please complete Step 1: User Input first.")
    st.stop()

payload = st.session_state['payload']

with st.spinner(f"🔬 Analyzing pathogenicity and querying STRING database for {payload['uniprot_id']}..."):
    ppi_data = fetch_string_ppi(payload['uniprot_id'])
    pathogenicity = predict_pathogenicity(payload['uniprot_id'], payload['mutation'])
    time.sleep(1.0)

# Persist real results for later stages (Explainable AI summary, Clinical Report)
st.session_state.setdefault('results', {})
st.session_state['results']['pathogenicity'] = pathogenicity
st.session_state['results']['ppi'] = ppi_data

col1, col2 = st.columns(2)

with col1:
    st.markdown("<h3>⚠️ Pathogenicity Estimate</h3>", unsafe_allow_html=True)

    if pathogenicity.get("has_clinical_record"):
        st.success("✅ Matched to a documented variant record in a public database.")
    else:
        st.info("ℹ️ No documented clinical record for this exact substitution — showing a modeled estimate.")

    risk_pct = round(pathogenicity["normalized_risk"] * 100, 1)
    gauge = go.Figure(go.Indicator(
        mode="gauge+number",
        value=risk_pct,
        number={'suffix': "%"},
        title={'text': "Estimated Deleteriousness"},
        gauge={
            'axis': {'range': [0, 100]},
            'bar': {'color': "#9B5DE5"},
            'steps': [
                {'range': [0, 35], 'color': "#1b3a2f"},
                {'range': [35, 65], 'color': "#4a3a12"},
                {'range': [65, 100], 'color': "#4a1414"},
            ],
        }
    ))
    gauge.update_layout(height=260, margin=dict(l=20, r=20, t=50, b=10),
                         paper_bgcolor="rgba(0,0,0,0)", font={'color': "white"})
    st.plotly_chart(gauge, use_container_width=True)

    st.markdown(f"**BLOSUM62 substitution score:** {pathogenicity['BLOSUM62 Score']['score']} "
                f"({pathogenicity['BLOSUM62 Score']['prediction']})")
    st.caption(pathogenicity["basis"])

with col2:
    st.markdown("<h3>🕸️ Protein-Protein Interactions (PPI)</h3>", unsafe_allow_html=True)
    st.write("Top interacting partners (live STRING database query) that may be disrupted by this mutation:")

    if ppi_data:
        df = pd.DataFrame(ppi_data)
        df['Interaction Score'] = df['Interaction Score'].apply(lambda x: f"{x * 100:.1f}%")
        st.dataframe(df, use_container_width=True, hide_index=True)
    else:
        st.warning(f"No significant interactions found in the STRING database for {payload['uniprot_id']}.")

st.divider()
n_partners = len(ppi_data)
if pathogenicity["normalized_risk"] > 0.5 and n_partners > 0:
    st.warning(f"⚠️ Elevated substitution risk combined with {n_partners} known interaction partner(s) — "
               "this mutation may warrant closer structural review in the next steps.")
else:
    st.success("✅ Mutation analysis complete — see Step 8 for the full explainable summary.")
