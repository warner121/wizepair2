# importing the sys module
import sys
import pandas as pd
from classes.mmp import MMP

def batch(): 
    infile = sys.argv[1]
    outfile = sys.argv[2]
    df = pd.read_json(infile, compression='gzip', lines=True)[1:10]
    s = df.wizepair2_uuid
    print(df.columns)
    df = df.apply(lambda x: 
            MMP(x.canonical_smiles_1, x.canonical_smiles_2, strictness=5).execute(), 
            axis=1).rename('response')
    df = df.to_frame().join(s).set_index('wizepair2_uuid')
    df = df.response.explode().reset_index()
    df.to_json(outfile, orient='records', compression='gzip', lines=True)
    return

if __name__ == "__main__":
    batch()
