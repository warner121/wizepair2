import sys
import re
import glob
import logging
import json
import numpy as np
import pandas as pd
from skelo.model.elo import EloEstimator

# Imports the Cloud Logging client library
import google.cloud.logging

# Instantiates a client
client = google.cloud.logging.Client()

# Retrieves a Cloud Logging handler based on the environment
# you're running in and integrates the handler with the
# Python logging module. By default this captures all logs
# at INFO level and higher
client.setup_logging()

logger = client.logger("batch_task_logs")

def elo(df, return_ratings=False):
    
    # create a table where winner / loser is defined
    if df.empty: return None
    
    # fit model
    model = EloEstimator(
        key1_field="fragment1",
        key2_field="fragment2",
        timestamp_field="publication_date_greatest",
        initial_time=pd.to_datetime('1970-01-01'),
        default_k=10
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
    #logger.log_text(f'infile = {infile}, outfile = {outfile}')
    #logger.log_text(json.dumps(glob.glob(infile)))
    infiles = re.sub('[0-9]{12}', '*', infile)
    #logger.log_text(f'infiles = {infiles}')
    #logger.log_text(json.dumps(glob.glob(infiles)))
    infiles = pd.Series(glob.glob(infiles)).sample(frac=0.2, replace=True)
    #logger.log_text(json.dumps(infiles.tolist()))
    df = pd.concat(infiles.apply(pd.read_csv, compression='gzip').tolist())
    
    # suffle publication dates
    publication_date_shuffled = df.publication_date_greatest.unique()
    np.random.shuffle(publication_date_shuffled)
    df_shuffled = pd.DataFrame(
        zip(df.publication_date_greatest.unique(), publication_date_shuffled), 
        columns=['publication_date_greatest', 'publication_date_shuffled'])
    df = df.merge(df_shuffled, on='publication_date_greatest').drop('publication_date_greatest', axis=1)
    df.rename({'publication_date_shuffled': 'publication_date_greatest'}, axis=1, inplace=True)

    # ensure data is in chronological order
    df.publication_date_greatest = pd.to_datetime(df.publication_date_greatest)
    df.sort_values('publication_date_greatest', inplace=True)
    df.reset_index(inplace=True, drop=True)

    # add label
    df['label'] = df.standard_change=='increase'
    
    # run elo scoring
    df = df.groupby(['chessleague_uuid']).apply(elo, return_ratings=True)
    df.reset_index(inplace=True)
    
    # write to out file
    df.sort_values('valid_from', inplace=True)
    df.drop_duplicates(['chessleague_uuid', 'key'], keep='last', inplace=True)
    #df[df.valid_to.isna()].to_csv(outfile, compression='gzip')
    logger.log_text(f'{outfile} write beginning')
    df.to_csv(outfile, compression='gzip', index=False)
    logger.log_text(f'{outfile} write complete')
    return

if __name__ == "__main__":
    batch()
