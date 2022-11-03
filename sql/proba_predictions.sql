select distinct
  mr.response.fragment1, 
  mr.response.fragment2, 
  mr.response.radius,
  er1.rating,
  er2.rating,
  1.0 / (1 + pow(10, (er2.rating - er1.rating) / 400.0)) as proba
from `cloudrun.mmp_responses` mr
join `cloudrun.elo_ratings` er1 on 
  to_hex(md5(mr.response.fragment1)) = er1.key and 
  mr.response.radius = er1.radius
join `cloudrun.elo_ratings` er2 on 
  to_hex(md5(mr.response.fragment2)) = er2.key and 
  mr.response.radius = er2.radius
where 
  er1.tid = 107 and 
  er2.tid = 107 and 
  er1.standard_type = 'EC50' and
  er2.standard_type = 'EC50'
order by proba desc
