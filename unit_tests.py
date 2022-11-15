import unittest
import pandas as pd

from classes.mmp import MMP

#@unittest.skip("showing class skipping")
class TestMMP(unittest.TestCase):

    def test_3_methylpyridine_to_2_methylpyridine(self):
        response = MMP('Cc1cccnc1', 'Cc1ccccn1', strictness=5).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.7142857142857143)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)
        self.assertEqual(df[df.valid].biproducts.sum(), 1)
        
    def test_2_5_dimethylfuran_to_1_4_dimethylbenzene(self):
        response = MMP('Cc1oc(C)cc1', 'Cc1ccc(C)cc1', strictness=5).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.75)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)
        self.assertEqual(df[df.valid].biproducts.sum(), 0)
        
    def test_nitro_to_ester(self):
        response = MMP('c1([N+](=O)[O-])ccccc1', 'c1(C(=O)OC)ccccc1', strictness=5).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.7)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)
        self.assertEqual(df[df.valid].biproducts.sum(), 0)
        
    def test_azetidine_to_piperazine(self):
        response = MMP('N1CCC1', 'N1CCNCC1', strictness=5).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.5)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)
        self.assertEqual(df[df.valid].biproducts.sum(), 0)

    def test_sildenafil_to_vardenafil(self):
        response = MMP(
            'CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12', 
            'CCCc1nc(C)c2c(=O)nc(-c3cc(S(=O)(=O)N4CCN(CC)CC4)ccc3OCC)[nH]n12', 
            strictness=5, 
            ).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.8)
        #self.assertEqual(df.valid.sum(), 4)
        #self.assertEqual(df[df.valid].radius.min(), 1)
        #self.assertEqual(df[df.valid].biproducts.sum(), 0)
        
    def test_pfizer_to_azpde4(self):
        response = MMP(
            'c1c(F)c(F)ccc1Oc2ncc(F)cc2C(=O)NC3CCC(NC(=O)c4cc(C)ccc4(O))CC3', 
            'c1c(F)cc2C(=O)N(C3CCC(NC(=O)c5cc(C)ccc5(O))CC3)C(=O)N(c4cc(F)c(F)cc4)c2n1', 
            strictness=5,
            ).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.868421052631579)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)
        self.assertEqual(df[df.valid].biproducts.sum(), 0)
        
#@unittest.skip("showing class skipping")
class TestHDAC(unittest.TestCase):
    
    def test_hdac_example(self):
        response = MMP(
            'Nc1ccccc1NC(=O)c1ccc(-c2ncc(CN3CCC3)cc2F)cc1', 
            'CCN1CCN(Cc2cnc(-c3ccc(C(=O)Nc4ccccc4N)cc3)c(C)c2)CC1', 
            strictness=5, 
            ).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.8125)
        self.assertEqual(df.valid.sum(), 2)
        self.assertEqual(df[df.valid].radius.min(), 3)
        self.assertEqual(df[df.valid].biproducts.sum(), 0)

#@unittest.skip("showing class skipping")
class TestBeta2(unittest.TestCase):
    
    def test_beta2_example1(self):
        response = MMP(
            'CC(C)NC[C@@H](O)c1ccc(O)c(O)c1', 
            'CC(C)NC[C@H](O)c1ccc(O)c(O)c1', 
            strictness=5, 
            ).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.9375)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)
        self.assertEqual(df[df.valid].biproducts.sum(), 0)
    
    def test_beta2_example2(self):
        response = MMP(
            'CNC[C@H](O)c1cccc(O)c1', 
            'CNCC(=O)c1ccc(O)c(O)c1', 
            strictness=5, 
            ).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.6923076923076923)
        self.assertEqual(df.valid.sum(), 3)
        self.assertEqual(df[df.valid].radius.min(), 2)
        self.assertEqual(df[df.valid].biproducts.sum(), 0)
        
    def test_beta2_example3(self):
        response = MMP(
            'CC(C)NC[C@H](O)c1ccc(NS(C)(=O)=O)c(O)c1', 
            'CC(C)NC[C@H](O)c1ccc(O)c(CS(C)(=O)=O)c1', 
            strictness=5, 
            ).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.9)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)    
        self.assertEqual(df[df.valid].biproducts.sum(), 0)
        
    def test_beta2_example4(self):
        response = MMP(
            'CC(C)c1cc(C(O)CN)ccc1O', 
            'CCc1ccc(C(O)CN)cc1O', 
            strictness=5, 
            ).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.7857142857142857)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)    
        self.assertEqual(df[df.valid].biproducts.sum(), 0)
        
    def test_beta2_example5(self):
        response = MMP(
            'CNC[C@@H](SC)c1ccc(O)c(O)c1', 
            'CC[C@H](NC(C)C)[C@H](O)c1ccc(O)c(O)c1', 
            strictness=5, 
            ).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.5789473684210527)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)
        self.assertEqual(df[df.valid].biproducts.sum(), 0)

#@unittest.skip("showing class skipping")
class TestNR3C1(unittest.TestCase):
    
    def test_NR3C1_example1(self):
        response = MMP(
            'C[C@]12C[C@H](O)[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]43C)[C@@H]1CC[C@]2(O)C(=O)CO', 
            'CCCC1O[C@@H]2C[C@H]3[C@@H]4CCC5=CC(=O)C=C[C@]5(C)[C@H]4[C@@H](O)C[C@]3(C)[C@]2(C(=O)CO)O1', 
            strictness=5, 
            ).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.6944444444444444)
        self.assertEqual(df.valid.sum(), 2)
        self.assertEqual(df[df.valid].radius.min(), 3)
        self.assertEqual(df[df.valid].biproducts.sum(), 0)
        
    def test_NR3C1_example2(self):
        response = MMP(
            'COc1ccc(F)cc1C(C)(C)CC(O)(Cn1cnc2ccccc21)C(F)(F)F', 
            'C[C@]12C[C@H](O)[C@H]3[C@@H](CCC4=CC(=O)C=C[C@@]43C)[C@@H]1CC[C@]2(O)C(=O)CO', 
            strictness=5, 
            ).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 1/30)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)
        self.assertEqual(df[df.valid].biproducts.sum(), 0)
        
#@unittest.skip("showing class skipping")
class TestCanonicalization(unittest.TestCase):
    
    def test_canonicalization(self):
        mmp1 = MMP('CC(=O)CCc1ccc2ccccc2c1', 'CC(=O)CCc1ccc2cc(Cl)ccc2c1', 
                   strictness=5).execute()
        mmp2 = MMP('S=c1[nH]ccn1Cc1ccccc1Cl', 'S=c1[nH]ccn1Cc1cc(Cl)ccc1Cl', 
                   strictness=5).execute()
        mmp3 = MMP('NC1=N[C@@H](CCc2ccccc2)CO1', 'NC1=N[C@@H](CCc2cccc(Cl)c2)CO1', 
                   strictness=5).execute()
        df1 = pd.json_normalize(mmp1)
        df2 = pd.json_normalize(mmp2)
        df3 = pd.json_normalize(mmp3)
        self.assertListEqual(df1[df1.radius<=2].smirks.tolist(), df2[df2.radius<=2].smirks.tolist())
        self.assertListEqual(df2[df2.radius<=2].smirks.tolist(), df3[df3.radius<=2].smirks.tolist())
        self.assertListEqual(df1[df1.radius<=2].smirks.tolist(), df3[df3.radius<=2].smirks.tolist())

    def test_unknown_salts(self):
        mmp = MMP('C[n+]1cccc2[nH]c3ccccc3c21.O=S(=O)([O-])C(F)(F)F', 'C[n+]1cccc2[nH]c3ccccc3c21.[Cl-]', 
                   strictness=5).execute()
        df = pd.json_normalize(mmp)
        self.assertEqual(df.percentmcs.mean(), 1)

if __name__ == '__main__':
    unittest.main()
