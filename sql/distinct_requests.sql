#standardSQL
create or replace table `wizepair2.cloudrun.mmp_searches` 
partition by RANGE_BUCKET(rand_partition, GENERATE_ARRAY(1, 4000, 1)) as
select 
  *,
  cast(ceiling(rand()*4000) as int64) as rand_partition
from (
  select distinct
    mmp_search_uuid,
    canonical_smiles_1, 
    canonical_smiles_2
  from `wizepair2.cloudrun.mmp_deltas`);

  EXPORT DATA OPTIONS(
  uri='gs://wizepair2_batch/mmp_input/mmp-*.json.gz',
  format='JSON',
  overwrite=true,
  compression='GZIP') AS select * from `cloudrun.mmp_searches`