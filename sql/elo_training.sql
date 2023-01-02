#standardSQL
create or replace table `cloudrun.elo_training`
partition by RANGE_BUCKET(rand_partition, GENERATE_ARRAY(1, 1000, 1)) as
select distinct * from (
  SELECT
    etb.*,
    'ANY' as pref_name,
    to_hex(md5(concat('ANY', ':', etb.standard_type, ':', etb.radius))) as chessleague_uuid,
    cast(ceiling(rand()*1000) as int64) as rand_partition
  FROM `cloudrun.elo_training_base` etb
  join `patents-public-data.ebi_chembl.assays_29` a using (assay_id)
  union all  
  SELECT
    etb.*,
    td.pref_name as pref_name,
    to_hex(md5(concat(td.pref_name, ':', etb.standard_type, ':', etb.radius))) as chessleague_uuid,
    cast(ceiling(rand()*1000) as int64) as rand_partition
  FROM `cloudrun.elo_training_base` etb
  join `patents-public-data.ebi_chembl.assays_29` a using (assay_id)
  join `patents-public-data.ebi_chembl.target_dictionary_29` td on a.tid = td.tid
  union all
  SELECT
    etb.*,
    td.pref_name as pref_name,
    to_hex(md5(concat(td.pref_name, ':', etb.standard_type, ':', etb.radius))) as chessleague_uuid,
    cast(ceiling(rand()*1000) as int64) as rand_partition
  FROM `cloudrun.elo_training_base` etb
  join `patents-public-data.ebi_chembl.assays_29` a using (assay_id)
  join `patents-public-data.ebi_chembl.target_relations_29` tr1 on a.tid = tr1.tid and tr1.relationship = 'SUBSET OF'
  join `patents-public-data.ebi_chembl.target_dictionary_29` td on tr1.related_tid = td.tid
  union all
  SELECT
    etb.*,
    td.pref_name as pref_name,
    to_hex(md5(concat(td.pref_name, ':', etb.standard_type, ':', etb.radius))) as chessleague_uuid,
    cast(ceiling(rand()*1000) as int64) as rand_partition
  FROM `cloudrun.elo_training_base` etb
  join `patents-public-data.ebi_chembl.assays_29` a using (assay_id)
  join `patents-public-data.ebi_chembl.target_relations_29` tr1 on a.tid = tr1.tid and tr1.relationship = 'SUBSET OF'
  join `patents-public-data.ebi_chembl.target_relations_29` tr2 on tr1.related_tid = tr2.tid and tr2.relationship = 'SUBSET OF'
  join `patents-public-data.ebi_chembl.target_dictionary_29` td on tr2.related_tid = td.tid
  union all
  SELECT
    etb.*,
    td.pref_name as pref_name,
    to_hex(md5(concat(td.pref_name, ':', etb.standard_type, ':', etb.radius))) as chessleague_uuid,
    cast(ceiling(rand()*1000) as int64) as rand_partition
  FROM `cloudrun.elo_training_base` etb
  join `patents-public-data.ebi_chembl.assays_29` a using (assay_id)
  join `patents-public-data.ebi_chembl.target_relations_29` tr1 on a.tid = tr1.tid and tr1.relationship = 'SUBSET OF'
  join `patents-public-data.ebi_chembl.target_relations_29` tr2 on tr1.related_tid = tr2.tid and tr2.relationship = 'SUBSET OF'
  join `patents-public-data.ebi_chembl.target_relations_29` tr3 on tr2.related_tid = tr3.tid and tr3.relationship = 'SUBSET OF'
  join `patents-public-data.ebi_chembl.target_dictionary_29` td on tr3.related_tid = td.tid);

EXPORT DATA OPTIONS(
  uri='gs://wizepair2_batch/elo_input/elo-*.csv.gz',
  format='CSV',
  overwrite=true,
  header=true,
  field_delimiter=',',
  compression='GZIP') AS select * from `cloudrun.elo_training`