import pytest
import pandas as pd
from wizepair2.mmp import MMP

@pytest.mark.parametrize("smiles1, smiles2, expected", [
    ('Cc1cccnc1', 'Cc1ccccn1', (0.7142857142857143, 4, 1, 1)),
    ('Cc1oc(C)cc1', 'Cc1ccc(C)cc1', (0.75, 4, 1, 0)),
    ('c1([N+](=O)[O-])ccccc1', 'c1(C(=O)OC)ccccc1', (0.7, 4, 1, 0)),
    ('N1CCC1', 'N1CCNCC1', (0.5, 4, 1, 0)),
    ('CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12',
     'CCCc1nc(C)c2c(=O)nc(-c3cc(S(=O)(=O)N4CCN(CC)CC4)ccc3OCC)[nH]n12', (14/17, None, None, None)),
    ('c1c(F)c(F)ccc1Oc2ncc(F)cc2C(=O)NC3CCC(NC(=O)c4cc(C)ccc4(O))CC3',
     'c1c(F)cc2C(=O)N(C3CCC(NC(=O)c5cc(C)ccc5(O))CC3)C(=O)N(c4cc(F)c(F)cc4)c2n1', (0.868421052631579, 4, 1, 0)),
])
def test_mmp_transformations(smiles1, smiles2, expected):
    percentmcs, valid, radius, biproducts = expected
    df = pd.json_normalize(MMP(smiles1, smiles2, strictness=7).execute())
    if percentmcs is not None:
        assert df.percentmcs.mean() == pytest.approx(percentmcs)
    if valid is not None:
        assert df.valid.sum() == valid
    if radius is not None:
        assert df[df.valid].radius.min() == radius
    if biproducts is not None:
        assert df[df.valid].biproducts.sum() == biproducts


@pytest.mark.parametrize("smiles1, smiles2, expected", [
    ('Nc1ccccc1NC(=O)c1ccc(-c2ncc(CN3CCC3)cc2F)cc1',
     'CCN1CCN(Cc2cnc(-c3ccc(C(=O)Nc4ccccc4N)cc3)c(C)c2)CC1', (0.8125, 2, 3, 0)),
])
def test_hdac(smiles1, smiles2, expected):
    percentmcs, valid, radius, biproducts = expected
    df = pd.json_normalize(MMP(smiles1, smiles2, strictness=7).execute())
    assert df.percentmcs.mean() == pytest.approx(percentmcs)
    assert df.valid.sum() == valid
    assert df[df.valid].radius.min() == radius
    assert df[df.valid].biproducts.sum() == biproducts


@pytest.mark.parametrize("smiles1, smiles2, expected", [
    ('CC(C)NC[C@@H](O)c1ccc(O)c(O)c1', 'CC(C)NC[C@H](O)c1ccc(O)c(O)c1', (0.9375, 4, 1, 0)),
    ('CNC[C@H](O)c1cccc(O)c1', 'CNCC(=O)c1ccc(O)c(O)c1', (0.6923076923076923, 3, 2, 0)),
    ('CC(C)NC[C@H](O)c1ccc(NS(C)(=O)=O)c(O)c1', 'CC(C)NC[C@H](O)c1ccc(O)c(CS(C)(=O)=O)c1', (0.9, 4, 1, 0)),
    ('CC(C)c1cc(C(O)CN)ccc1O', 'CCc1ccc(C(O)CN)cc1O', (9/14, 4, 1, 0)),
    ('CNC[C@@H](SC)c1ccc(O)c(O)c1', 'CC[C@H](NC(C)C)[C@H](O)c1ccc(O)c(O)c1', (0.5789473684210527, 4, 1, 0)),
])
def test_beta2(smiles1, smiles2, expected):
    percentmcs, valid, radius, biproducts = expected
    df = pd.json_normalize(MMP(smiles1, smiles2, strictness=7).execute())
    assert df.percentmcs.mean() == pytest.approx(percentmcs)
    assert df.valid.sum() == valid
    assert df[df.valid].radius.min() == radius
    assert df[df.valid].biproducts.sum() == biproducts


@pytest.mark.parametrize("smiles1, smiles2, expected", [
    ('C[C@]12C[C@H](O)[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]43C)[C@@H]1CC[C@]2(O)C(=O)CO',
     'CCCC1O[C@@H]2C[C@H]3[C@@H]4CCC5=CC(=O)C=C[C@]5(C)[C@H]4[C@@H](O)C[C@]3(C)[C@]2(C(=O)CO)O1', (0.6944444444444444, 2, 3, 0)),
    ('COc1ccc(F)cc1C(C)(C)CC(O)(Cn1cnc2ccccc21)C(F)(F)F',
     'C[C@]12C[C@H](O)[C@H]3[C@@H](CCC4=CC(=O)C=C[C@@]43C)[C@@H]1CC[C@]2(O)C(=O)CO', (1/15, 4, 1, 0)),
])
def test_nr3c1(smiles1, smiles2, expected):
    percentmcs, valid, radius, biproducts = expected
    df = pd.json_normalize(MMP(smiles1, smiles2, strictness=7).execute())
    assert df.percentmcs.mean() == pytest.approx(percentmcs)
    assert df.valid.sum() == valid
    assert df[df.valid].radius.min() == radius
    assert df[df.valid].biproducts.sum() == biproducts
