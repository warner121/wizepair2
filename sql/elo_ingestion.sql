LOAD DATA OVERWRITE `cloudrun.elo_ratings_temp`
FROM FILES(
  skip_leading_rows=1,
  format='CSV',
  compression='GZIP',
  uris = ['gs://wizepair2_batch/elo_output/elo-*.csv.gz']
  );

create or replace table `wizepair2.cloudrun.elo_ratings_temp` as
select
  ert.chessleague_uuid,
  ert.key,
  min(ert.valid_from) as valid_from_min,
  max(ert.valid_from) as valid_from_max,
  --JSON_QUERY_ARRAY(replace(to_json_string(array_agg(ert.rating)), ']', concat(repeat(',1500', 10-count(*)), ']'))) as ratings
  array_agg(ert.rating) as ratings
from `cloudrun.elo_ratings_temp` ert
group by
  ert.chessleague_uuid,
  ert.key;

create or replace table `wizepair2.cloudrun.elo_ratings` as
select
  ert.chessleague_uuid,
  ert.key,
  ert.valid_from_min,
  ert.valid_from_max,
  array_agg(rating) as ratings,
  avg(cast(rating as float64)) as rating_avg,
  stddev(cast(rating as float64)) as rating_stddev
from `cloudrun.elo_ratings_temp` ert, ert.ratings as rating
group by
  ert.chessleague_uuid,
  ert.key,
  ert.valid_from_min,
  ert.valid_from_max;

drop table `cloudrun.elo_ratings_temp`;