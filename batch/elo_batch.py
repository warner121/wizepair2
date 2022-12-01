import sys
import re
import glob
import logging
import json
import pandas as pd
from skelo.model.elo import EloEstimator

def elo(df, return_ratings=False):
    
    # create a table where winner / loser is defined
    if df.empty: return None
    
    # fit model
    model = EloEstimator(
        key1_field="fragment1",
        key2_field="fragment2",
        timestamp_field="publication_date_greatest",
        initial_time=pd.to_datetime('1970-01-01'),
        default_k=20
    ).fit(df, df.label)
    if return_ratings: return model.rating_model.to_frame()

    # calculate proba
    df['proba'] = model.transform(df, output_type='prob', strict_past_data=True)
    return df[['assay_id', 'wizepair2_uuid', 'chessleague_uuid', 'proba', 'label']]

def batch(): 

    # obtain in/out files from command line
    infile = sys.argv[1]
    outfile = sys.argv[2]

    # read elo input data
    logging.info(f'infile = {infile}, outfile = {outfile}')
    logging.info(json.dumps(glob.glob(infile)))
    infiles = re.sub('[0-9]{12}', '*', infile)
    logging.info(f'infiles = {infiles}')
    logging.info(json.dumps(glob.glob(infiles)))
    infiles = pd.Series(glob.glob(infiles)).sample(frac=1, replace=True)
    df = pd.concat(infiles.apply(pd.read_csv, compression='gzip').tolist())

    # ensure data is in chronological order
    df.publication_date_greatest = pd.to_datetime(df.publication_date_greatest)
    df.sort_values('publication_date_greatest', inplace=True)
    df.reset_index(inplace=True, drop=True)

    # add label
    df['label'] = df.standard_change=='increase'
    
    # run elo scoring
    df = df.groupby(['chessleague_uuid']).apply(elo, return_ratings=True)
    
    # write to out file
    df[df.valid_to.isna()].to_csv(outfile, compression='gzip')
    return

if __name__ == "__main__":
    batch()
