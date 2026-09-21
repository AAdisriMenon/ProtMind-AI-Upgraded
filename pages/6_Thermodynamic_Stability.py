import streamlit as st
import time
from backend import apply_custom_css, extract_features, predict_thermodynamic_stability

st.set_page_config(page_title="Step 6 | ProtMind AI", page_icon="🔥")
apply_custom_css()

st.markdown("<h2>🔥 Step 6: Thermodynamic Stability (ΔG)</h2>", unsafe_allow_html=True)

if not st.session_state.get('payload'):
    st.warning("⚠️ No active protein data found. Please complete Step 1: User Input first.")
    st.stop()

payload = st.session_state['payload']

with st.spinner("⚛️ Calculating biophysical potential energy and ΔG shifts..."):
    features = extract_features(payload['sequence'], payload['mutation'])
    
    stability_data = predict_thermodynamic_stability(
        wt_aa=features["Wildtype AA"],
        mut_aa=features["Mutant AA"],
        wt_property=features["WT Property"],
        mut_property=features["Mutant Property"]
    )
    time.sleep(1)

st.session_state.setdefault('results', {})
st.session_state['results']['stability'] = stability_data

col1, col2 = st.columns([1, 1])

# Safely extract the score using 'ddG' (or 'delta_g' depending on your backend key)
delta_g_val = stability_data.get('ddG', 0.0)

with col1:
    st.markdown("<h3>⚙️ Gibbs Free Energy (ΔG)</h3>", unsafe_allow_html=True)
    st.info("ΔG measures the conformational energy state. **Negative values** indicate a stable, favorable state; **positive values** indicate high structural energy.")
    
    # Format readout cleanly with proper sign (+ or -)
    formatted_val = f"+{delta_g_val:.2f}" if delta_g_val > 0 else f"{delta_g_val:.2f}"
    
    st.metric(
        label="Calculated ΔG (kcal/mol)", 
        value=formatted_val,
        delta="Unstable" if delta_g_val > 0 else "Stable",
        delta_color="inverse"
    )

with col2:
    st.markdown("<h3>📉 Biophysical Verdict</h3>", unsafe_allow_html=True)
    st.write("Based on the calculated energy shift, the structural impact is:")
    
    st.markdown(f"## {stability_data['color']} {stability_data['status']}")
    
    if stability_data['status'] == "Highly Destabilizing":
        st.error(stability_data['alert'])
    elif stability_data['status'] == "Mildly Destabilizing":
        st.warning(stability_data['alert'])
    else:
        st.success(stability_data['alert'])

st.divider()
st.markdown("### 🧪 Calculation Parameters")
st.write(f"- **Wildtype Core:** {features['WT Property']} ({features['Wildtype AA']})")
st.write(f"- **Mutant Core:** {features['Mutant Property']} ({features['Mutant AA']})")
st.write(f"- **Hydrophobicity Shift (Kyte-Doolittle):** {stability_data['hydrophobicity_delta']}")
st.write(f"- **Residue Volume Shift (Å³, Zamyatnin):** {stability_data['volume_delta']}")
st.write(f"- **BLOSUM62 Substitution Score:** {stability_data['blosum62']}")
st.caption(stability_data["method"])
