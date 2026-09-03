import streamlit as st
import pandas as pd
from backend import apply_custom_css, fetch_chembl_drugs

st.set_page_config(page_title="Step 9 | ProtMind AI", page_icon="💊")
apply_custom_css()

st.markdown("<h2>💊 Step 9: Therapeutic Targeting</h2>", unsafe_allow_html=True)

if not st.session_state.get('payload'):
    st.warning("⚠️ No active protein data found. Please complete Step 1: User Input first.")
    st.stop()

payload = st.session_state['payload']
composite = st.session_state.get('results', {}).get('composite')

st.info(f"Querying the **ChEMBL** database for molecules with a documented mechanism of action against "
        f"**{payload['uniprot_id']}** — i.e., existing drug-discovery leads that already target this protein.")

with st.spinner(f"💊 Querying ChEMBL for known ligands of {payload['uniprot_id']}..."):
    drugs = fetch_chembl_drugs(payload['uniprot_id'])

st.session_state.setdefault('results', {})
st.session_state['results']['drugs'] = drugs

if composite:
    if composite["tier"] in ("High Concern", "Moderate Concern"):
        st.warning(f"This mutation was flagged as **{composite['tier']}** in Step 8 — existing ligands for this "
                   "target may be worth evaluating as a starting point for follow-up research.")
    else:
        st.caption("This mutation showed low composite risk in Step 8; targeting information is shown for "
                   "completeness.")

st.divider()

if drugs:
    st.markdown(f"### 🎯 {len(drugs)} Known Ligand(s)/Mechanism(s) Found")
    df = pd.DataFrame(drugs)
    st.dataframe(df, use_container_width=True, hide_index=True)

    action_counts = df["Action Type"].value_counts()
    st.markdown("**Mechanism types represented:**")
    st.bar_chart(action_counts)

    st.caption("Source: ChEMBL `mechanism` endpoint, filtered to targets matching this UniProt accession. "
               "Presence here means a compound has a documented mechanism against this protein target in "
               "general — not that it specifically corrects the effect of this exact mutation.")
else:
    st.warning(f"No documented drug mechanisms were found in ChEMBL for **{payload['uniprot_id']}**. "
               "This may mean the target has no approved/investigational ligands yet, or that ChEMBL does not "
               "list this particular accession — it does not necessarily mean the protein is undruggable.")

st.divider()
st.success("✅ Therapeutic targeting check complete. Proceed to **Step 10** for the compiled clinical report.")
