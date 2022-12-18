SELECT `wizepair2.cloudrun.mmp`(to_json_string(struct(smirks, smiles)))
from (select "[#6:4](:[#6:3](:[#6:2]-[H])-[H])-[H]>>[#6:4](:[#7:3]:[#6:2]-[H])-[H]" as smirks, "c1ccccc1" as smiles)