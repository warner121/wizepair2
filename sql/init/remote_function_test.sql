SELECT
  `wizepair2.cloudrun.reactor`(TO_JSON_STRING(STRUCT(smirks,
        smiles)))
FROM (
  SELECT
    "[#6:4](:[#6:3](:[#6:2]-[H])-[H])-[H]>>[#6:4](:[#7:3]:[#6:2]-[H])-[H]" AS smirks,
    "c1ccccc1" AS smiles );