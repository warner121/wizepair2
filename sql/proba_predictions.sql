with 

elo_rating_radii as (
  select distinct 
    chessleague_uuid,
    radius,
    pref_name,
    standard_type
  from `cloudrun.elo_training_agg`
  where pref_name in ('HERG', 'ANY') 
    and standard_type = 'IC50' 
),

elo_ratings_enhanced as (
  select * 
  from `cloudrun.elo_ratings` et
  left join elo_rating_radii using (chessleague_uuid)
),

reactor_products as (
  select distinct
    mr.response.fragment1,
    mr.response.fragment2,
    mr.response.smirks,
    ere1.radius,
    ere1.pref_name,
    ere1.standard_type,
    to_json_string(ere1.ratings) as ratings1,
    to_json_string(ere2.ratings) as ratings2,
    ere1.rating_avg as rating1,
    ere2.rating_avg as rating2,
    ere3.rating_avg as rating3,
    ere4.rating_avg as rating4,
    1.0 / (1 + pow(10, (ere2.rating_avg - ere1.rating_avg) / 400.0)) as proba1,
    1.0 / (1 + pow(10, (ere4.rating_avg - ere3.rating_avg) / 400.0)) as proba2,
    eta1.deltas_count as deltas_count1,
    eta1.deltas_avg as deltas_avg_1,
    to_json_string(eta1.mmp_delta_array) as mmp_delta_array1,
    eta2.deltas_count as deltas_count2,
    eta2.deltas_avg as deltas_avg_2,
    to_json_string(eta2.mmp_delta_array) as mmp_delta_array2
  from `cloudrun.mmp_results` mr
  join elo_ratings_enhanced ere1 on 
    ere1.key = mr.response.fragment1 and
    ere1.radius = mr.response.radius and
    ere1.pref_name = 'HERG' and
    ere1.standard_type = 'IC50'
  join elo_ratings_enhanced ere2 on 
    ere2.key = mr.response.fragment2 and
    ere2.radius = mr.response.radius and 
    ere2.pref_name = 'HERG' and
    ere2.standard_type = 'IC50'
  join elo_ratings_enhanced ere3 on 
    ere3.key = mr.response.fragment1 and
    ere3.radius = mr.response.radius and 
    ere3.pref_name = 'ANY' and
    ere3.standard_type = 'IC50'
  join elo_ratings_enhanced ere4 on 
    ere4.key = mr.response.fragment2 and
    ere4.radius = mr.response.radius and 
    ere4.pref_name = 'ANY' and
    ere4.standard_type = 'IC50'
  left join `cloudrun.elo_training_agg` eta1 on
    mr.response.fragment1 = eta1.fragment1 and
    mr.response.fragment2 = eta1.fragment2 and
    ere1.chessleague_uuid = eta1.chessleague_uuid
  left join `cloudrun.elo_training_agg` eta2 on
    mr.response.fragment1 = eta2.fragment1 and
    mr.response.fragment2 = eta2.fragment2 and
    ere3.chessleague_uuid = eta2.chessleague_uuid
  order by proba1 * (1-(2*abs(proba2-0.5))) desc
)

select 
  *,
  'c1ccccc1CCc2ccc(CN)cc2' as smiles,
  `wizepair2.cloudrun.reactor`(to_json_string(struct(smirks, 'c1ccccc1CCc2ccc(CN)cc2' as smiles))) as products
from reactor_products