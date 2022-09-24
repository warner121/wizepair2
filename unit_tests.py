import unittest
import pandas as pd

from mmpa.mmp import MMP

#@unittest.skip("showing class skipping")
class TestMMP(unittest.TestCase):

    def test_3_methylpyridine_to_2_methylpyridine(self):
        response = MMP('Cc1cccnc1', 'Cc1ccccn1', strictness=5, correspondence=2).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.7142857142857143)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)
        
    def test_2_5_dimethylfuran_to_1_4_dimethylbenzene(self):
        response = MMP('Cc1oc(C)cc1', 'Cc1ccc(C)cc1', strictness=5, correspondence=2).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.75)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)
        
    def test_nitro_to_ester(self):
        response = MMP('c1([N+](=O)[O-])ccccc1', 'c1(C(=O)OC)ccccc1', strictness=5, correspondence=2).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.7)
        #self.assertEqual(df.valid.sum(), 4)
        #self.assertEqual(df[df.valid].radius.min(), 1)
        
    def test_azetidine_to_piperazine(self):
        response = MMP('N1CCC1', 'N1CCNCC1', strictness=5, correspondence=2).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.5)
        self.assertEqual(df.valid.sum(), 1)
        self.assertEqual(df[df.valid].radius.min(), 1)

    def test_sildenafil_to_vardenafil(self):
        response = MMP(
            'CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12', 
            'CCCc1nc(C)c2c(=O)nc(-c3cc(S(=O)(=O)N4CCN(CC)CC4)ccc3OCC)[nH]n12', 
            strictness=5, 
            correspondence=2).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.8)
        #self.assertEqual(df.valid.sum(), 4)
        #self.assertEqual(df[df.valid].radius.min(), 1)

#@unittest.skip("showing class skipping")
class TestHDAC(unittest.TestCase):
    
    def test_hdac_example(self):
        response = MMP(
            'Nc1ccccc1NC(=O)c1ccc(-c2ncc(CN3CCC3)cc2F)cc1', 
            'CCN1CCN(Cc2cnc(-c3ccc(C(=O)Nc4ccccc4N)cc3)c(C)c2)CC1', 
            strictness=5, 
            correspondence=2).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.8125)
        self.assertEqual(df.valid.sum(), 2)
        self.assertEqual(df[df.valid].radius.min(), 3)

#@unittest.skip("showing class skipping")
class TestBeta2(unittest.TestCase):
    
    def test_beta2_example1(self):
        response = MMP(
            'CC(C)NC[C@@H](O)c1ccc(O)c(O)c1', 
            'CC(C)NC[C@H](O)c1ccc(O)c(O)c1', 
            strictness=5, 
            correspondence=1).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.9375)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)
    
    def test_beta2_example2(self):
        response = MMP(
            'CNC[C@H](O)c1cccc(O)c1', 
            'CNCC(=O)c1ccc(O)c(O)c1', 
            strictness=5, 
            correspondence=1).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.6923076923076923)
        self.assertEqual(df.valid.sum(), 3)
        self.assertEqual(df[df.valid].radius.min(), 2)
        
    def test_beta2_example3(self):
        response = MMP(
            'CC(C)NC[C@H](O)c1ccc(NS(C)(=O)=O)c(O)c1', 
            'CC(C)NC[C@H](O)c1ccc(O)c(CS(C)(=O)=O)c1', 
            strictness=5, 
            correspondence=1).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.9)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)    
        
    def test_beta2_example4(self):
        response = MMP(
            'CC(C)c1cc(C(O)CN)ccc1O', 
            'CCc1ccc(C(O)CN)cc1O', 
            strictness=5, 
            correspondence=1).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.7857142857142857)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)    
        
    def test_beta2_example5(self):
        response = MMP(
            'CNC[C@@H](SC)c1ccc(O)c(O)c1', 
            'CC[C@H](NC(C)C)[C@H](O)c1ccc(O)c(O)c1', 
            strictness=5, 
            correspondence=1).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.5789473684210527)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)

#@unittest.skip("showing class skipping")
class TestNR3C1(unittest.TestCase):
    
    def test_NR3C1_example1(self):
        response = MMP(
            'C[C@]12C[C@H](O)[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]43C)[C@@H]1CC[C@]2(O)C(=O)CO', 
            'CCCC1O[C@@H]2C[C@H]3[C@@H]4CCC5=CC(=O)C=C[C@]5(C)[C@H]4[C@@H](O)C[C@]3(C)[C@]2(C(=O)CO)O1', 
            strictness=5, 
            correspondence=1).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 0.6944444444444444)
        self.assertEqual(df.valid.sum(), 2)
        self.assertEqual(df[df.valid].radius.min(), 3)
        
    def test_NR3C1_example2(self):
        response = MMP(
            'COc1ccc(F)cc1C(C)(C)CC(O)(Cn1cnc2ccccc21)C(F)(F)F', 
            'C[C@]12C[C@H](O)[C@H]3[C@@H](CCC4=CC(=O)C=C[C@@]43C)[C@@H]1CC[C@]2(O)C(=O)CO', 
            strictness=5, 
            correspondence=1).execute()
        df = pd.json_normalize(response)
        self.assertEqual(df.percentmcs.mean(), 1/30)
        self.assertEqual(df.valid.sum(), 4)
        self.assertEqual(df[df.valid].radius.min(), 1)
            
if __name__ == '__main__':
    unittest.main()