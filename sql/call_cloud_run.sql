SELECT `wizepair2.cloudrun.mmp`(to_json_string(struct(smiles1, smiles2, strictness)))
from (select "c1ccccc1" as smiles1, "c1ccccn1" as smiles2, 5 as strictness)