create or replace view `cloudrun.elo_training_base` as 
SELECT
  md.wizepair2_uuid,
  md.assay_id,
  md.standard_type,
  md.standard_change,
  md.publication_date_greatest,
  mr.response.fragment1,
  mr.response.fragment2,
  mr.response.radius
FROM `cloudrun.mmp_deltas` md
join `cloudrun.mmp_responses` mr using (wizepair2_uuid) 
where 
  mr.response.valid and
  md.standard_change is not null and
  mr.response.fragment1 is not null and 
  mr.response.fragment2 is not null