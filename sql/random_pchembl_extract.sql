#standardSQL
with activities as (
select 
  act.pchembl_value,
  com.canonical_smiles, 
    cast(cmp.heavy_atoms as int64) as heavy_atoms,
  mol.chembl_id as molecule_chembl_id,
  ass.chembl_id as assay_chembl_id,
  count(*) over (partition by act.assay_id) as count_activities
FROM `patents-public-data.ebi_chembl.activities_29` act
join `patents-public-data.ebi_chembl.compound_structures_29` com using (molregno)
join `patents-public-data.ebi_chembl.compound_properties_29` cmp using (molregno)
join `patents-public-data.ebi_chembl.molecule_dictionary_29` mol using (molregno)
join `patents-public-data.ebi_chembl.assays_29` ass using (assay_id)
where pchembl_value is not null)

select *
from activities a1
join activities a2 using (assay_chembl_id)
where
  a1.molecule_chembl_id != a2.molecule_chembl_id and
  a1.count_activities < 50 and 
  a2.count_activities < 50 and 
  a1.heavy_atoms between 2 and 50 and
  a2.heavy_atoms between 2 and 50 and
  rand() < 0.0001