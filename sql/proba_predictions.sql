with 

elo_rating_radii as (
  select distinct 
    chessleague_uuid,
    radius,
    standard_type,
    pref_name
  from `cloudrun.elo_training`
),

elo_ratings_augmented as (
  select * 
  from `cloudrun.elo_ratings` et
  join elo_rating_radii using (chessleague_uuid)
)

select distinct
  mr.response.fragment1, 
  mr.response.fragment2, 
  mr.response.smirks, 
  era1.radius,
  era1.pref_name,
  era1.standard_type,
  era1.rating as rating1,
  era2.rating as rating2,
  1.0 / (1 + pow(10, (era1.rating - era2.rating) / 400.0)) as proba
from `cloudrun.mmp_responses` mr
join elo_ratings_augmented era1 on 
  mr.response.fragment1 = era1.key and 
  mr.response.radius = era1.radius
join elo_ratings_augmented era2 on 
  mr.response.fragment2 = era2.key and 
  mr.response.radius = era2.radius
where 
  era1.chessleague_uuid = '001114ce71206ed6d15807e3eb08f0f1' and
  era2.chessleague_uuid = '001114ce71206ed6d15807e3eb08f0f1'
order by proba desc