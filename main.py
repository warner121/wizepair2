import json
import logging
import pandas as pd
from typing import List, Dict, Any
from fastapi import FastAPI, HTTPException, Body, Query
from pydantic import BaseModel
from src.wizepair2.mmp import MMP, Reactor, Desalinator

# Initialize FastAPI app
app = FastAPI(title="WizePair2", description="API for Matched Molecular Pair Analysis")


# Define Pydantic models for request validation
class MMPInput(BaseModel):
    smiles1: str
    smiles2: str


class ReactorInput(BaseModel):
    smirks: str
    smiles: str


class BQRemoteInput(BaseModel):
    calls: List[List[str]]


# MMP endpoint
@app.post("/wizepair2/api/v1.0/mmp", response_model=List[Any])
async def mmp_endpoint(data: List[MMPInput] = Body(...), strictness: int = Query(...)):
    """
    Perform Matched Molecular Pair analysis on input SMILES pairs.
    """
    try:
        df = pd.DataFrame([item.dict() for item in data])
        df = df.apply(lambda x: MMP(x.smiles1, x.smiles2, strictness).execute(), axis=1)
        return df.tolist()
    except Exception as e:
        logging.error(f"Error in MMP endpoint: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))


# Helper function for reactor calls
def strip_and_react(smirks: str, smiles: str) -> List[Any]:
    smiles = Desalinator(smiles).getSmiles()
    return Reactor(smirks).generate_products(smiles)


# Reactor endpoint
@app.post("/wizepair2/api/v1.0/reactor", response_model=List[Any])
async def reactor_endpoint(data: List[ReactorInput] = Body(...)):
    """
    Perform chemical reactions based on SMIRKS patterns and SMILES.
    """
    try:
        df = pd.DataFrame([item.dict() for item in data])
        df = df.apply(lambda x: strip_and_react(x.smirks, x.smiles), axis=1)
        return df.tolist()
    except Exception as e:
        logging.error(f"Error in reactor endpoint: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))


# BQ Remote endpoint
@app.post("/", response_model=Dict[str, List[str]])
async def bqremote_endpoint(data: BQRemoteInput = Body(...)):
    """
    BigQuery remote function endpoint for chemical reactions.
    """
    try:
        return_value = []
        calls = data.calls

        for call in calls:
            try:
                call_data = json.loads(call[0])
                productlist = json.dumps(strip_and_react(call_data['smirks'], call_data['smiles']))
            except ValueError as e:
                logging.error(f"Error processing call: {str(e)}")
                logging.error(f"Call data: {call[0]}")
                productlist = '[]'

            return_value.append(productlist)

        return {"replies": return_value}
    except Exception as e:
        logging.error(f"Error in BQ remote endpoint: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8000)