import logging
import pandas as pd

from flask import jsonify, request

from classes.mmp import MMP

from app import app

# define routes
@app.route('/wizepair2/api/v1.0/mmp', methods=['POST'])
def mmp():
    
    # ensure json and read to data frame
    strictness = int(request.args['strictness'])
    if not request.json: abort(400)
    df = pd.DataFrame(request.get_json())
    
    # perform mmpa
    df = df.sample(frac=1).apply(lambda x: MMP(x.smiles1, x.smiles2, strictness, correspondence=2).execute(), axis=1)
        
    # format response and return
    return jsonify(df.tolist())
