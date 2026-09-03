import streamlit as st
import time
from backend import apply_custom_css, fetch_uniprot_data, fetch_alphafold_url, fetch_known_variant_annotation

st.set_page_config(page_title="Step 2 | ProtMind AI", page_icon="🔍")
apply_custom_css()

st.markdown("<h2>🔍 Step 2: Data Retrieval Results</h2>", unsafe_allow_html=True)

if not st.session_state.get('payload'):
    st.warning("⚠️ No active protein data found. Please complete Step 1: User Input first.")
    st.stop()

payload = st.session_state['payload']
mutation = payload['mutation']
wt_aa, mut_aa, position = mutation[0], mutation[-1], int(mutation[1:-1])

with st.spinner(f"🔬 Querying live databases for {payload['uniprot_id']}..."):
    uniprot_metadata = fetch_uniprot_data(payload['uniprot_id'])
    af_pdb_url = fetch_alphafold_url(payload['uniprot_id'])
    known_variant = fetch_known_variant_annotation(payload['uniprot_id'], position, wt_aa, mut_aa)
    time.sleep(0.5)

tab1, tab2, tab3, tab4 = st.tabs(["🔵 UniProt", "🏥 Clinical Record", "🧬 Structure", "ℹ️ Notes"])

with tab1:
    st.markdown("<h3>Protein Knowledgebase</h3>", unsafe_allow_html=True)
    st.json(uniprot_metadata)

with tab2:
    st.markdown("<h3>Known Variant / Clinical Significance</h3>", unsafe_allow_html=True)
    st.metric(label="Target Mutation", value=mutation)
    if known_variant:
        st.success(f"✅ This exact substitution has a documented record: **{known_variant['clinical_significance']}**")
        if known_variant["sources"]:
            st.caption("Sources: " + ", ".join(known_variant["sources"]))
    else:
        st.info("ℹ️ No documented record for this exact substitution was found in the public UniProt/ClinVar "
                "variation database. This is common and expected for novel or hypothetical mutations — it "
                "does not mean the mutation is harmless, only that it hasn't been catalogued yet. "
                "Step 4 provides a computed substitution-likelihood estimate instead.")

with tab3:
    st.markdown("<h3>Predicted 3D Structure</h3>", unsafe_allow_html=True)
    if af_pdb_url:
        st.success("AlphaFold structure located!")
        st.markdown(f"**[Download .PDB File]({af_pdb_url})**")
    else:
        st.warning("No AlphaFold structure found.")

with tab4:
    st.markdown("<h3>About This Step</h3>", unsafe_allow_html=True)
    st.write(
        "This stage only surfaces data that genuinely exists for this protein/mutation — protein identity "
        "and structure come from live UniProt and AlphaFold queries, and the clinical record above is a real "
        "database lookup. Population-frequency data (e.g., gnomAD allele frequencies) is intentionally not "
        "shown here, because computing it correctly requires mapping this protein-level mutation to a genomic "
        "coordinate first — a good candidate for a future upgrade to this pipeline."
    )
