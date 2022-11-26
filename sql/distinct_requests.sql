#standardSQL
create or replace table `wizepair2.cloudrun.mmp_requests` 
partition by RANGE_BUCKET(rand_partition, GENERATE_ARRAY(1, 4000, 1)) as
select 
  *,
  cast(ceiling(rand()*4000) as int64) as rand_partition
from (
  select distinct
    wizepair2_uuid,
    canonical_smiles_1, 
    canonical_smiles_2
  from `wizepair2.cloudrun.mmp_deltas`)