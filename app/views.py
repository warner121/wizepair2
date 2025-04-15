import json
import logging
import pandas as pd

from flask import jsonify, request, abort
from src.wizepair2.mmp import MMP, Reactor, Desalinator
from app import app

# define routes
@app.route('/wizepair2/api/v1.0/mmp', methods=['POST'])
def mmp():

    # ensure json and read to data frame
    strictness = int(request.args['strictness'])
    if not request.json: abort(400)
    df = pd.DataFrame(request.get_json())

    # perform mmpa
    df = df.apply(lambda x: MMP(x.smiles1, x.smiles2, strictness).execute(), axis=1)

    # format response and return
    return jsonify(df.tolist())

# define function for reactor calls
def strip_and_react(smirks, smiles):
    smiles = Desalinator(smiles).getSmiles()
    return Reactor(smirks).generate_products(smiles)
    
# define routes
@app.route('/wizepair2/api/v1.0/reactor', methods=['POST'])
def reactor():

    # ensure json and read to data frame
    if not request.json: abort(400)
    df = pd.DataFrame(request.get_json())

    # perform reactions
    df = df.apply(lambda x: strip_and_react(x.smirks, x.smiles), axis=1)

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
            productlist = json.dumps(strip_and_react(call['smirks'], call['smiles']))
        except ValueError: 
            print(json.loads(call[0])) 
            productlist = '[]' 
        return_value.append(productlist)
    
    # format response and return
    return jsonify({ "replies" :  return_value})

