from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from backend import extract_features, predict_thermodynamic_stability

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class MutationRequest(BaseModel):
    sequence: str
    mutation: str

@app.post("/api/analyze")
def analyze_mutation(data: MutationRequest):
    features = extract_features(data.sequence, data.mutation)
    stability = predict_thermodynamic_stability(
        features["Wildtype AA"], features["Mutant AA"], 
        features["WT Property"], features["Mutant Property"]
    )
    return {
        "features": features,
        "stability": stability
    }
