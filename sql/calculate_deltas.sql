#standardSQL
create or replace table `wizepair2.cloudrun.mmp_deltas` as
select 
  *, 
  greatest(heavy_atoms_1, heavy_atoms_2) as heavy_atoms_greatest,
  greatest(publication_date_1, publication_date_2) as publication_date_greatest,
  case 
    when 
      standard_value_1 > standard_value_2 and 
      standard_relation_1 not in ('<', '<<') and 
      standard_relation_2 not in ('>', '>>')
    then 'decrease'
    when
      standard_value_1 < standard_value_2 and 
      standard_relation_1 not in ('>', '>>') and 
      standard_relation_2 not in ('<', '<<') 
    then 'increase'
    when
      standard_value_1 = standard_value_2 and 
      standard_relation_1 in ('=', '~') and 
      standard_relation_2 in ('=', '~') 
    then 'no-change'
    else null
  end as standard_change,
  to_json_string(struct(canonical_smiles_1 as smiles1, canonical_smiles_2 as smiles2, 5 as strictness)) as request,
  to_hex(md5(to_json_string(struct(canonical_smiles_1 as smiles1, canonical_smiles_2 as smiles2, 5 as strictness)))) as wizepair2_uuid
from (
  select 
    act.assay_id,
    act.standard_type,
    cast(act.standard_value as numeric) as standard_value_1,
    act.standard_relation as standard_relation_1,
    cast(act.pchembl_value as numeric) as pchembl_value_1,
    count(*) over (partition by act.assay_id) as count_activities_1,
    act.molregno as molregno_1,
    com.canonical_smiles as canonical_smiles_1,
    cast(cmp.heavy_atoms as int64) as heavy_atoms_1,
  date(
    coalesce(cast(d.year as int64), 1970), 
    coalesce(cast(ceiling(percent_rank() over (
      partition by d.journal, d.year order by SAFE_CAST(d.first_page as int64)) * 11) + 1 as int64), 1),
    coalesce(cast(ceiling(percent_rank() over (
      partition by d.journal, d.year order by SAFE_CAST(d.first_page as int64)) * 27) + 1 as int64), 1)) as publication_date_1
  FROM `patents-public-data.ebi_chembl.activities_29` act
  join `patents-public-data.ebi_chembl.compound_structures_29` com using (molregno)
  join `patents-public-data.ebi_chembl.compound_properties_29` cmp using (molregno)
  left join `patents-public-data.ebi_chembl.docs_29` d using (doc_id)
  where standard_type in (select distinct standard_type from`patents-public-data.ebi_chembl.activities_29` where pchembl_value is not null)
  ) a1
join (
  select 
    act.assay_id,
    act.standard_type,
    cast(act.standard_value as numeric) as standard_value_2,
    act.standard_relation as standard_relation_2,
    cast(act.pchembl_value as numeric) as pchembl_value_2,
    count(*) over (partition by act.assay_id) as count_activities_2,    
    act.molregno as molregno_2,
    com.canonical_smiles as canonical_smiles_2, 
    cast(cmp.heavy_atoms as int64) as heavy_atoms_2,
  date(
    coalesce(cast(d.year as int64), 1970), 
    coalesce(cast(ceiling(percent_rank() over (
      partition by d.journal, d.year order by SAFE_CAST(d.first_page as int64)) * 11) + 1 as int64), 1),
    coalesce(cast(ceiling(percent_rank() over (
      partition by d.journal, d.year order by SAFE_CAST(d.first_page as int64)) * 27) + 1 as int64), 1)) as publication_date_2
  FROM `patents-public-data.ebi_chembl.activities_29` act
  join `patents-public-data.ebi_chembl.compound_structures_29` com using (molregno)
  join `patents-public-data.ebi_chembl.compound_properties_29` cmp using (molregno)
  left join `patents-public-data.ebi_chembl.docs_29` d using (doc_id)
  where standard_type in (select distinct standard_type from`patents-public-data.ebi_chembl.activities_29` where pchembl_value is not null)
  ) a2 using (assay_id, standard_type)
where 
  a1.molregno_1 != a2.molregno_2 and
  a1.count_activities_1 < 50 and 
  a2.count_activities_2 < 50 and 
  a1.heavy_atoms_1 between 2 and 20 and
  a2.heavy_atoms_2 between 2 and 20 and
  a1.standard_value_1 is not null and 
  a2.standard_value_2 is not null