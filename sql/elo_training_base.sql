create or replace view `cloudrun.elo_training_base` as 
SELECT
  md.mmp_delta_uuid,
  md.assay_id,
  md.standard_type,
  md.standard_change,
  md.publication_date_greatest,
  mr.response.fragment1,
  mr.response.fragment2,
  mr.response.radius
FROM `cloudrun.mmp_deltas` md
join `cloudrun.mmp_results` mr using (mmp_search_uuid) 
where 
  md.standard_change is not null and
  mr.response.valid and
  mr.response.fragment1 is not null and 
  mr.response.fragment2 is not null