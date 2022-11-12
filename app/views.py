import json
import logging
import pandas as pd

from flask import jsonify, request, abort

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
    df = df.sample(frac=1).apply(lambda x: MMP(x.smiles1, x.smiles2, strictness).execute(), axis=1)

    # format response and return
    return jsonify(df.tolist())

# define routes
@app.route('/', methods=['POST'])
def bqremote():
    
    # ensure json and read to data frame
    return_value = []
    request_json = request.get_json()
    calls = request_json['calls']
    # df = pd.json_normalize(pd.Series(calls).explode().apply(json.loads))
    for call in calls:
        call = json.loads(call[0])
        try: 
            mmp = json.dumps(MMP(call['smiles1'], call['smiles2'], call['strictness']).execute())
        except ValueError: 
            print(json.loads(call[0])) 
            mmp = '{}' 
        return_value.append(mmp)
    
    # format response and return
    return jsonify({ "replies" :  return_value})

