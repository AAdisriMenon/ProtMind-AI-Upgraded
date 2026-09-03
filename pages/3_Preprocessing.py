import streamlit as st
import time
from backend import apply_custom_css, extract_features, calculate_alignment_score

st.set_page_config(page_title="Step 3 | ProtMind AI", page_icon="⚙️")
apply_custom_css()

st.markdown("<h2>⚙️ Step 3: Preprocessing & Feature Extraction</h2>", unsafe_allow_html=True)

if not st.session_state.get('payload'):
    st.warning("⚠️ No active protein data found. Please complete Step 1: User Input first.")
    st.stop()

payload = st.session_state['payload']

with st.spinner("⏳ Extracting physicochemical features and running MSA..."):
    features = extract_features(payload['sequence'], payload['mutation'])
    
    mutant_list = list(payload['sequence'])
    if features["Positional Match"]:
        mutant_list[features["Position"] - 1] = features["Mutant AA"]
    mutant_sequence = "".join(mutant_list)
    
    alignment_score = calculate_alignment_score(payload['sequence'], mutant_sequence)
    time.sleep(1)

# Persist real computed results so later pipeline stages (Explainable AI, Clinical Report)
# can build on the SAME numbers instead of recomputing or hardcoding them.
st.session_state.setdefault('results', {})
st.session_state['results']['features'] = features
st.session_state['results']['alignment_score'] = alignment_score

col_feat, col_align = st.columns(2)

with col_feat:
    st.markdown("<h3>Extracted Features</h3>", unsafe_allow_html=True)
    if not features["Positional Match"]:
        st.error(f"⚠️ Wildtype mismatch at position {features['Position']}.")
    st.json(features)
    
with col_align:
    st.markdown("<h3>Evolutionary Analysis (MSA)</h3>", unsafe_allow_html=True)
    st.metric(label="Sequence Identity", value=f"{alignment_score:.2f}%")
    st.progress(alignment_score / 100.0)
    
    if features["WT Property"] != features["Mutant Property"]:
        st.warning(f"**Significant Shift:** {features['Property Shift']}")
    else:
        st.info(f"**Conserved Property:** {features['WT Property']}")