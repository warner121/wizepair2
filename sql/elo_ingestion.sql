LOAD DATA OVERWRITE `cloudrun.elo_ratings_temp`
FROM FILES(
  skip_leading_rows=1,
  format='CSV',
  compression='GZIP',
  uris = ['gs://wizepair2_batch/elo_output/elo_input-*.csv.gz']
  );

create or replace table `wizepair2.cloudrun.elo_ratings_temp` as
select
  ert.chessleague_uuid,
  ert.key,
  min(ert.valid_from) as valid_from_min,
  max(ert.valid_from) as valid_from_max,
  JSON_QUERY_ARRAY(replace(to_json_string(array_agg(ert.rating)), ']', concat(repeat(',1500', 10-count(*)), ']'))) as ratings
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
  avg(cast(rating as float64)) as rating_avg,
  stddev(cast(rating as float64)) as rating_stddev
from `cloudrun.elo_ratings_temp` ert, ert.ratings as rating
group by
  ert.chessleague_uuid,
  ert.key,
  ert.valid_from_min,
  ert.valid_from_max;

drop table `cloudrun.elo_ratings_temp`;

create or replace table `wizepair2.cloudrun.elo_training_agg` as
select
  chessleague_uuid,
  radius,
  pref_name,
  standard_type,
  fragment1, 
  fragment2,
  count(*) as deltas_count,
  count(distinct et.assay_id) as assay_count,
  count(distinct et.doc_id_greatest) as doc_count,
  count(distinct et.wizepair2_uuid) as wizepair2_count,
  avg(
    case 
      when standard_change='increase' then 1.0
      when standard_change='no-change' then 0.0
      when standard_change='decrease' then -1.0 
    end) as deltas_avg,
  stddev(
    case 
      when standard_change='increase' then 1.0
      when standard_change='no-change' then 0.0
      when standard_change='decrease' then -1.0 
    end) as deltas_stddev
from `cloudrun.elo_training` et
group by
  chessleague_uuid,
  radius,
  pref_name,
  standard_type,
  fragment1, 
  fragment2