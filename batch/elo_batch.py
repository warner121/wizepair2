import sys
import re
import glob
import json
import numpy as np
import pandas as pd
import concurrent.futures
import itertools
from skelo.model.elo import EloEstimator
import google.cloud.logging

# Setup Google Cloud logging
client = google.cloud.logging.Client()
client.setup_logging()
logger = client.logger("batch_task_logs")

def elo(df, return_ratings=False):
    """
    Compute ELO ratings or probabilities for a given DataFrame slice.
    If return_ratings is True, return a DataFrame of ratings over time.
    Otherwise, return the input DataFrame enriched with a 'proba' column.
    """
    if df.empty:
        return None

    # Fit the Elo model
    model = (
        EloEstimator(
            key1_field="fragment1",
            key2_field="fragment2",
            timestamp_field="publication_date_greatest",
            initial_time=pd.to_datetime('1970-01-01'),
            default_k=10
        )
        .fit(df, df.label)
    )

    if return_ratings:
        # Obtain the time-series of ratings
        ratings_df = model.rating_model.to_frame()
        # Reset index to turn MultiIndex into columns
        ratings_df = ratings_df.reset_index()
        # Tag with the chessleague identifier
        ratings_df['chessleague_uuid'] = df['chessleague_uuid'].iloc[0]
        return ratings_df

    # Otherwise, compute win probabilities for each record
    df['proba'] = model.transform(df, output_type='prob', strict_past_data=True)
    return df[['assay_id', 'wizepair2_uuid', 'chessleague_uuid', 'proba', 'label']]


def batch():
    # Parse input/output paths
    infile = sys.argv[1]
    outfile = sys.argv[2]
    logger.log_text(f'infile = {infile}, outfile = {outfile}')

    # Expand pattern for input files
    infiles_pattern = re.sub(r'[0-9]{12}', '*', infile)
    logger.log_text(f'infiles pattern = {infiles_pattern}')

    # Sample 10% of matching files (with replacement)
    file_list = pd.Series(glob.glob(infiles_pattern))
    sample_files = file_list.sample(frac=0.2, replace=True).tolist()
    logger.log_text(json.dumps(sample_files))

    # Read and concatenate data
    df = pd.concat(
        [pd.read_csv(f, compression='gzip') for f in sample_files],
        ignore_index=True
    )

    # Shuffle publication dates to anonymize ordering
    unique_dates = df.publication_date_greatest.unique()
    np.random.shuffle(unique_dates)
    df_dates = pd.DataFrame({
        'publication_date_greatest': df.publication_date_greatest.unique(),
        'publication_date_shuffled': unique_dates
    })
    df = (
        df
        .merge(df_dates, on='publication_date_greatest', how='left')
        .drop(columns=['publication_date_greatest'])
        .rename(columns={'publication_date_shuffled': 'publication_date_greatest'})
    )

    # Sort chronologically
    df['publication_date_greatest'] = pd.to_datetime(df['publication_date_greatest'])
    df.sort_values('publication_date_greatest', inplace=True)
    df.reset_index(drop=True, inplace=True)

    # Create binary label for increases
    df['label'] = df.standard_change == 'increase'

    # Parallel ELO rating computation
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = list(executor.map(
            elo,
            (group for _, group in df.groupby('chessleague_uuid')),
            itertools.repeat(True),
        ))

    # Combine results, dropping any None returns
    df_ratings = pd.concat([r for r in results if r is not None], ignore_index=True)

    # Final sorting and deduplication
    df_ratings.sort_values('valid_from', inplace=True)
    df_final = (
        df_ratings
        .drop_duplicates(['chessleague_uuid', 'key'], keep='last')
        .reset_index(drop=True)
    )

    # Write output
    logger.log_text(f'Writing output to {outfile}')
    df_final.to_csv(outfile, compression='gzip', index=False)
    logger.log_text(f'Completed writing {outfile}')


if __name__ == "__main__":
    batch()
