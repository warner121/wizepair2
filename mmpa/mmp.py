import json
import networkx

from rdkit import Chem
from func_timeout import func_timeout, FunctionTimedOut
from networkx.algorithms.clique import find_cliques

class MMP():

    @staticmethod
    def __molFromSmiles(smiles: str):
        
        # pasre smiles
        try: mol = Chem.MolFromSmiles(smiles)
        except: return None
            
        # add hydrogen where defining isomer
        isomerics = []
        for atom in mol.GetAtoms():
            if not atom.HasProp('_CIPCode'): continue
            isomerics.append(atom.GetIdx())
        mol = Chem.AddHs(mol, onlyOnAtoms=isomerics, explicitOnly=True)
                      
        # clear mappings and initialise radii (assume all atoms are RECS)
        for atom in mol.GetAtoms(): 
            atom.SetProp('molAtomRadius','0')
            atom.ClearProp('molAtomMapNumber')
            
        # return
        return mol
    
    def __init__(self, smiles_x: str, smiles_y: str):
        
        # initialise molecules for comparison
        self._mol1 = self.__molFromSmiles(smiles_x)
        self._mol2 = self.__molFromSmiles(smiles_y)
        
        # intialise correspondance graph
        self._graph = networkx.Graph()
       
    def createCorrespondence(self, penalty=3.0):
        '''
        Build the correspondance matrix from which to determine the MCS. 
        
        Each atomic pairing is assigned a providional score, from 0 (identical) to 1 (at least [penalty] differences detected). 
        This pairwise score is attached to the nodes alongside the atomic indices.
        '''
        
        # create local ref for input molecules
        mol1 = self._mol1
        mol2 = self._mol2
            
        # iterate over all potential atom-atom pairings
        for atom1 in mol1.GetAtoms():
            for atom2 in mol2.GetAtoms():
                
                # store the CIP codes somewhere that doesn't throw errors on comparison when missing
                try: atom1._CIPCode = atom1.GetProp('_CIPCode')
                except KeyError: atom1._CIPCode = None
                try: atom2._CIPCode = atom2.GetProp('_CIPCode')
                except KeyError: atom2._CIPCode = None            
    
                # set penalties - [penalty] strikes and you're out!
                __tempscore = 0
                if atom1.GetAtomicNum() != atom2.GetAtomicNum(): __tempscore += 1
                if atom1.GetImplicitValence() != atom2.GetImplicitValence(): __tempscore += 1
                if atom1.GetExplicitValence() != atom2.GetExplicitValence(): __tempscore += 1
                if atom1.GetFormalCharge() != atom2.GetFormalCharge(): __tempscore += 1
                if atom1.GetIsAromatic() != atom2.GetIsAromatic(): __tempscore += 1
                if atom1.GetDegree() != atom2.GetDegree(): __tempscore += 1
                if atom1.IsInRing() != atom2.IsInRing(): __tempscore += 1
                if atom1._CIPCode != atom2._CIPCode: __tempscore += 1
                
                # set upper limit on penalty to 1 and append to node
                __tempscore = min(1, __tempscore/penalty)              
                mapping = (atom1.GetIdx(), atom2.GetIdx(), __tempscore)
                if __tempscore < 1: self._graph.add_node(mapping) 
        
        # calculate distance matrices
        __dmat1 = Chem.GetDistanceMatrix(mol1)
        __dmat2 = Chem.GetDistanceMatrix(mol2)  

        # create correspondance graph edges
        for map1 in self._graph.nodes():
            for map2 in self._graph.nodes():

                # test if criteria are met for correspondance
                correspondance = __dmat1[map1[0]][map2[0]] == __dmat2[map1[1]][map2[1]]
                if correspondance: self._graph.add_edge(map1, map2)       
                    
    def findCliques(self, timeout=60):
                    
        # define function with no arguments for use with timeout function
        def findCliquesNoArgs(): return list(find_cliques(self._graph))
        
        # try finding cliques within [timeout] seconds
        try:
            cliques = func_timeout(timeout, findCliquesNoArgs)
        except FunctionTimedOut:
            logging.warning(json.dumps({"message": "failed to find cliques in {} seconds".format(timeout)}))
            return None
            
        # score (largest cliques first)
        bestscore = -1E800
        cliques.sort(key=len, reverse=True)
        for clique in cliques:   
            
            # initialize scores
            score = len(clique)
            mcs = list(clique)
            if score <= bestscore: continue # score cannot be greater than clique size
    
            # iterate over the mappings identified from MCSS
            for pair in clique:
                score = score - pair[2] # lookup predefined score penalty
                if pair[2]: mcs.remove(pair) # postive score marks discrepancy, mcs will contain only identical pairings

            # store results if best so far
            if score > bestscore:
                bestscore = score
                bestmcs = mcs
                bestclique = clique  
                
        # store results
        self._clique = bestclique
        self._mcs = bestmcs
    
    def __setAtomMapNumbers(self):        
        '''
        Use the indices of the best scoring clique to define the atom mappings.
        '''

        # iterate over the remaining mappings identified from MCSS
        for pair in self._clique:
            
            # increment the index to prevent atom mappings of 0
            mapIdx = self._clique.index(pair) + 1
            
            # map first atom and set radius to 99 (atom part of MCS)
            atom1 = self._mol1.GetAtomWithIdx(pair[0])
            atom1.SetProp('molAtomMapNumber','%d'%mapIdx)
            
            # map second atom and set radius to 99 (atom part of MCS)
            atom2 = self._mol2.GetAtomWithIdx(pair[1])
            atom2.SetProp('molAtomMapNumber','%d'%mapIdx)
    
    def __setAtomRadii(self):
        '''
        Use the atomic radii to denote which atoms are part of the MCS. By elimination, those atoms with radii of 0 will form the RECS.
        '''

        # iterate over the remaining mappings identified from MCSS
        for pair in self._mcs:
            
            # map first atom and set radius to 99 (atom part of MCS)
            atom1 = self._mol1.GetAtomWithIdx(pair[0])
            atom1.SetProp('molAtomRadius','99')
            
            # map second atom and set radius to 99 (atom part of MCS)
            atom2 = self._mol2.GetAtomWithIdx(pair[1])
            atom2.SetProp('molAtomRadius','99')

    def eliminateMCS(self):
        
        # mark up atom mappings and MCS/RECS split
        self.__setAtomMapNumbers()
        self.__setAtomRadii()
        
        def eliminate(mol):
        
            # tag atoms within 4 bonds of attachment
            toRemove = set(range(mol.GetNumAtoms()))
            for atom in mol.GetAtoms():
                if atom.GetProp('molAtomRadius') ==  '0':
                    for idx in Chem.FindAtomEnvironmentOfRadiusN(mol, 4, atom.GetIdx()):
                        envBond = mol.GetBondWithIdx(idx)
                        toRemove.discard(envBond.GetBeginAtom().GetIdx())
                        toRemove.discard(envBond.GetEndAtom().GetIdx())
                        
            # remove environment from core
            toRemove = list(toRemove)
            toRemove.sort(reverse=True)
            frag = Chem.EditableMol(mol)
            for atom in toRemove: frag.RemoveAtom(atom)
            frag = frag.GetMol()
            return frag
      
        self._frag1 = eliminate(self._mol1)
        self._frag2 = eliminate(self._mol2)
        
    def getFragment1(self):
        
        frag = Chem.EditableMol(self._frag1).GetMol()
        for atom in frag.GetAtoms(): atom.ClearProp('molAtomMapNumber')
        return Chem.MolToSmiles(frag)
    
    def getFragment2(self):
        
        frag = Chem.EditableMol(self._frag2).GetMol()
        for atom in frag.GetAtoms(): atom.ClearProp('molAtomMapNumber')
        return Chem.MolToSmiles(frag)        

    def getSmirks(self):
        '''
        Generate SMIRKS, ensuring reaction denotes a simple 1:1 mapping of fragments.
        
        Hydrogens are added to ensure the transformation may only be applied in the appropriate molecular environment.
        '''
        
        smirks = '{}>>{}'.format(Chem.MolToSmarts(Chem.AddHs(self._frag1)), Chem.MolToSmarts(Chem.AddHs(self._frag2)))
        rxn = Chem.rdChemReactions.ReactionFromSmarts(smirks)
        if rxn.GetNumReactantTemplates() != 1 or rxn.GetNumProductTemplates() != 1: return None
        return smirks