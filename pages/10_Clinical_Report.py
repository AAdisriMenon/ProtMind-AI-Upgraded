import streamlit as st
from backend import (
    apply_custom_css, generate_pdf_report, fetch_alphafold_url, fetch_pdb_content,
    calculate_phi_psi, extract_features, predict_pathogenicity, predict_thermodynamic_stability,
    fetch_string_ppi, compute_composite_risk
)

st.set_page_config(page_title="Step 10 | ProtMind AI", page_icon="📄")
apply_custom_css()

st.markdown("<h2>📄 Step 10: Automated Clinical Report</h2>", unsafe_allow_html=True)

if not st.session_state.get('payload'):
    st.warning("⚠️ No active protein data found. Please complete Step 1: User Input first.")
    st.stop()

payload = st.session_state['payload']
results = st.session_state.get('results', {})

st.markdown("### Export Full Pipeline Analytics")
st.write("Generate a comprehensive PDF document summarizing sequence metrics, structural disruption profiles, "
         "and therapeutic targeting insights — built from the same numbers computed in Steps 3–9.")

if len(results) < 4:
    st.info("ℹ️ Not every earlier step has been visited this session, so a few values below will be computed "
            "fresh rather than reused. For a fully consistent report, run Steps 3, 4, 6, 7 and 9 first.")

# Reuse everything already computed this session; only recompute what's genuinely missing.
features = results.get('features') or extract_features(payload['sequence'], payload['mutation'])
pathogenicity = results.get('pathogenicity') or predict_pathogenicity(payload['uniprot_id'], payload['mutation'])
stability = results.get('stability') or predict_thermodynamic_stability(
    features['Wildtype AA'], features['Mutant AA'], features['WT Property'], features['Mutant Property']
)
ppi = results.get('ppi')
if ppi is None:
    ppi = fetch_string_ppi(payload['uniprot_id'])

angles = results.get('angles')
if angles is None:
    mut_pos = int(payload['mutation'][1:-1]) if payload['mutation'][1:-1].isdigit() else 1
    af_pdb_url = fetch_alphafold_url(payload['uniprot_id'])
    pdb_content = fetch_pdb_content(af_pdb_url)
    angles = calculate_phi_psi(pdb_content, mut_pos)

composite = results.get('composite') or compute_composite_risk(stability, pathogenicity, len(ppi), angles is not None)

report_results = {
    "features": features, "pathogenicity": pathogenicity, "stability": stability,
    "ppi": ppi, "composite": composite,
}

with st.spinner("📑 Assembling report and generating structural plots..."):
    pdf_bytes = generate_pdf_report(payload=payload, results=report_results, angles=angles)

st.markdown("<br>", unsafe_allow_html=True)

col1, col2, col3 = st.columns(3)
col1.metric("Composite Risk", f"{composite['composite_score']}/100", composite['tier'])
col2.metric("Estimated ΔΔG", f"{stability['ddG']} kcal/mol", stability['status'])
col3.metric("Interaction Partners", len(ppi))

st.download_button(
    label="📥 Download Clinical Report (PDF)",
    data=pdf_bytes,
    file_name=f"{payload['uniprot_id']}_{payload['mutation']}_Clinical_Report.pdf",
    mime="application/pdf"
)

st.success("✅ Report generated with embedded structural and biophysical profiles, consistent with the values "
           "shown throughout the pipeline.")
