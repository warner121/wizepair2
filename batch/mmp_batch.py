import sys
sys.path.append("./")

import pandas as pd
from classes.mmp import MMP

def batch(): 

    # obtain in/out files from command line
    infile = sys.argv[1]
    outfile = sys.argv[2]

    # read mmp data into frame and retain reference
    df = pd.read_json(infile, compression='gzip', lines=True)
    s = df.mmp_search_uuid

    # process mmpa
    df = df.apply(lambda x: 
            MMP(x.canonical_smiles_1, x.canonical_smiles_2, strictness=5).execute(), 
            axis=1).rename('response')

    # rejoin with references
    df = df.to_frame().join(s).set_index('mmp_search_uuid')
    df = df.response.explode().reset_index()

    # write to out file
    df.to_json(outfile, orient='records', compression='gzip', lines=True)
    return

if __name__ == "__main__":
    batch()
