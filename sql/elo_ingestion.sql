LOAD DATA OVERWRITE `cloudrun.elo_ratings_temp`
FROM FILES(
  skip_leading_rows=1,
  format='CSV',
  compression='GZIP',
  uris = ['gs://wizepair2_batch/elo_output/elo_input-*.csv.gz']
  );

create or replace table `wizepair2.cloudrun.elo_ratings` as
select
  ert.chessleague_uuid,
  ert.key,
  avg(ert.rating) as rating_avg,
  stddev(ert.rating) as rating_stddev,
  min(ert.valid_from) as valid_from_min,
  max(ert.valid_from) as valid_from_max
from `cloudrun.elo_ratings_temp` ert
group by
  ert.chessleague_uuid,
  ert.key;

drop table `cloudrun.elo_ratings_temp`;