with 

elo_rating_radii as (
  select distinct 
    chessleague_uuid,
    radius,
    pref_name,
    standard_type
  from `cloudrun.elo_training`
  where pref_name in ('HERG', 'ANY') 
    and standard_type = 'IC50' 
    and radius = 2
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
  era1.rating_avg as rating1,
  era2.rating_avg as rating2,
  era3.rating_avg as rating3,
  era4.rating_avg as rating4,
  1.0 / (1 + pow(10, (era1.rating_avg - era2.rating_avg) / 400.0)) as proba1,
  1.0 / (1 + pow(10, (era3.rating_avg - era4.rating_avg) / 400.0)) as proba2
from `cloudrun.mmp_responses` mr
join elo_ratings_augmented era1 on 
  era1.key = mr.response.fragment1 and
  era1.radius = mr.response.radius and
  era1.pref_name = 'HERG' and
  era1.standard_type = 'IC50'
join elo_ratings_augmented era2 on 
  era2.key = mr.response.fragment2 and
  era2.radius = mr.response.radius and 
  era2.pref_name = 'HERG' and
  era2.standard_type = 'IC50'
join elo_ratings_augmented era3 on 
  era3.key = mr.response.fragment1 and
  era3.radius = mr.response.radius and 
  era3.pref_name = 'ANY' and
  era3.standard_type = 'IC50'
join elo_ratings_augmented era4 on 
  era4.key = mr.response.fragment2 and
  era4.radius = mr.response.radius and 
  era4.pref_name = 'ANY' and
  era4.standard_type = 'IC50'
order by proba1 * (1-proba2) desc