import json
import logging
import networkx
import numpy as np

from rdkit import Chem
from func_timeout import func_timeout, FunctionTimedOut
from networkx.algorithms.clique import find_cliques

class CorrespondenceGraph(networkx.Graph):
    '''
    Build the correspondance matrix of putative atom pairings, from which to determine the maximal clique (comprising the MCS). 
    '''
    
    def __init__(self, mol1, mol2, fuzziness):
        '''
        Instantiate the correspondance graph.
        
        Each atomic pairing is assigned a providional score, from 0 (identical) to 1 (fewer than [fuzziness] differences detected). 
        This pairwise score is included in a tuple representing the each node, alongside the atomic indices.
        '''
        
        # inherit from networkx.Graph
        super().__init__(self)
        
        # store fuzziness for scoring
        self._fuzz = fuzziness
        
        # calculate distance matrices
        dmat1 = Chem.GetDistanceMatrix(mol1)
        dmat2 = Chem.GetDistanceMatrix(mol2) 
        
        # extract propery in such a way error is not thrown on comparison
        def getCIPCode(atom):
            try: return atom.GetProp('_CIPCode')
            except KeyError: return None
            
        # create description of how central in molecule atom
        def getPeripherality(atom, dmat):
            peripherality = dmat[atom.GetIdx()]
            return np.mean(peripherality / np.max(peripherality))            
            
        # iterate over all potential atom-atom pairings
        for atom1 in mol1.GetAtoms():
            for atom2 in mol2.GetAtoms():

                # set penalties - [penalty] strikes and you're out!
                mismatches = 0
                if atom1.GetAtomicNum() != atom2.GetAtomicNum(): mismatches += 1
                if atom1.GetImplicitValence() != atom2.GetImplicitValence(): mismatches += 1
                if atom1.GetExplicitValence() != atom2.GetExplicitValence(): mismatches += 1
                if atom1.GetFormalCharge() != atom2.GetFormalCharge(): mismatches += 1
                if atom1.GetIsAromatic() != atom2.GetIsAromatic(): mismatches += 1
                if atom1.GetDegree() != atom2.GetDegree(): mismatches += 1
                if atom1.IsInRing() != atom2.IsInRing(): mismatches += 1
                if getCIPCode(atom1) != getCIPCode(atom2): mismatches += 1

                # apply jitter in the event of a tie
                peripherality1 = getPeripherality(atom1, dmat1)
                peripherality2 = getPeripherality(atom2, dmat2)
                mismatches += 0.001 * np.linalg.norm(peripherality1 - peripherality2)
                
                # set upper limit on score deductions to 1 and append to node
                mismatches = min(1, mismatches/self._fuzz)              
                mapping = (atom1.GetIdx(), atom2.GetIdx(), mismatches)
                if mismatches < 1: self.add_node(mapping)  

        # create correspondance graph edges
        for map1 in self.nodes():
            for map2 in self.nodes():

                # test if criteria are met for correspondance
                correspondance = dmat1[map1[0]][map2[0]] == dmat2[map1[1]][map2[1]]
                if correspondance: self.add_edge(map1, map2)       
                    
    def execute(self, timeout=60):
        '''
        Execute the maximal clique search.
        
        Returns:
            bestclique - the best scoring clique (including partially matching atoms, and used to derive reaction mappings)
            bestmcs - the subset of bestclique containing exact chemical matches only (to be discarded to produce the RECS)
        '''
                    
        # define function with no arguments (for use with timeout function)
        def findCliquesNoArgs(): return list(find_cliques(self))
        
        # try finding cliques within [timeout] seconds
        try:
            cliques = func_timeout(timeout, findCliquesNoArgs)
        except FunctionTimedOut:
            logging.warning(json.dumps({"message": "failed to find cliques in {} seconds".format(timeout)}))
            return list(), list()
            
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
                if pair[2] >= (1.0/self._fuzz): mcs.remove(pair) # score > threshold marks atom discrepancy, mcs will contain only identical pairings

            # store results if best so far
            if score > bestscore:
                bestscore = score
                bestmcs = mcs
                bestclique = clique  
                
        # return results
        return bestclique, bestmcs
        
class MMP():

    @staticmethod
    def __molFromSmiles(smiles: str):
        
        # parse smiles
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
    
    def __init__(self, smiles_x: str, smiles_y: str, fuzziness=4):
        '''
        Initialise the matched molecular pair.
        
        smiles_x: First molecule to compare.
        smiles_y: Second molecule to compare.
        fuzziness: Integer (1-8) to indicate how tolerant the algortitm should to be to atom-wise differences in the construction of the correspondance graph. 
            1 (fastest) atoms chemically identical to be considered part of mcss. 
            8 (slowest) purely topological comparison of structures. 
        '''
        
        # canonicalise smiles
        self._smiles1 = Chem.MolToSmiles(Chem.MolFromSmiles(smiles_x))
        self._smiles2 = Chem.MolToSmiles(Chem.MolFromSmiles(smiles_y))
         
        # initialise molecules for comparison
        self._mol1 = self.__molFromSmiles(self._smiles1)
        self._mol2 = self.__molFromSmiles(self._smiles2)
        
        # intialise correspondance graph
        self._graph = CorrespondenceGraph(self._mol1, self._mol2, fuzziness)
        
        # dummy vars
        self._clique = None
        self._mcs = None
        
    def __search(self):
        
        # find the MCS
        self._clique, self._mcs = self._graph.execute()

        # determine the % of largest molecule covered by MCS
        self._percentmcs = len(self._mcs) / max(self._mol1.GetNumAtoms(), self._mol2.GetNumAtoms())        

    def __setAtomMapNumbers(self):        
        '''
        Use the indices of the best scoring clique to define the atom mappings.
        '''

        # iterate over the mappings identified from MCSS
        for pair in self._clique:
            
            # increment the index to prevent atom mappings of 0
            mapIdx = self._clique.index(pair) + 1
            
            # map first atom and set radius to 99 (atom part of MCS)
            atom1 = self._mol1.GetAtomWithIdx(pair[0])
            atom1.SetProp('molAtomMapNumber', '%d'%mapIdx)
            
            # map second atom and set radius to 99 (atom part of MCS)
            atom2 = self._mol2.GetAtomWithIdx(pair[1])
            atom2.SetProp('molAtomMapNumber', '%d'%mapIdx)
    
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

    def execute(self, radii=4):
                
        # search, mark up atom mappings and MCS/RECS split
        self.__search()
        self.__setAtomMapNumbers()
        self.__setAtomRadii()
        
        # define function for elimination of MCS
        def eliminate(mol, radius):
        
            # tag atoms within 4 bonds of attachment
            toRemove = set(range(mol.GetNumAtoms()))
            for atom in mol.GetAtoms():
                if atom.GetProp('molAtomRadius') ==  '0':
                    for idx in Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atom.GetIdx()):
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
      
        # loop from 4 down to 1 bond radius to find smallest valid transformation
        responselist = list()
        for radius in reversed(range(radii+1)):
            
            # return list of valid transformations
            if radius == 0: return responselist
            
            # initialise response object
            response = {'smiles1': self._smiles1,
                        'smiles2': self._smiles2,
                        'percentmcs': self._percentmcs,
                        'radius': radius,
                        'valid': False}
            
            # Define reaction as SMIRKS while mappings still present
            frag1 = eliminate(self._mol1, radius)
            frag2 = eliminate(self._mol2, radius)   
            smirks = '{}>>{}'.format(Chem.MolToSmarts(Chem.AddHs(frag1)), Chem.MolToSmarts(Chem.AddHs(frag2)))
            response['smirks'] = smirks
            
            # verify 1:1 reaction
            rxn = Chem.rdChemReactions.ReactionFromSmarts(smirks)
            if rxn.GetNumReactantTemplates() != 1 or rxn.GetNumProductTemplates() != 1: 
                logging.info(json.dumps({'radius': radius, "message": "no 1:1 reaction could be generated"}))
                responselist.append(response)
                continue

            # verify derived reaction produces original 'product'
            productset = rxn.RunReactants((Chem.AddHs(Chem.MolFromSmiles(self._smiles1)),))
            productlist = list()
            for product in productset:
                productlist.append('.'.join([Chem.MolToSmiles(Chem.RemoveHs(productpart)) for productpart in product]))
            if self._smiles2 not in productlist:
                logging.info(json.dumps({'radius': radius, "message": "second molecule not found amongst products enumerated from first"}))
                responselist.append(response)
                continue

            # remove mappings to yield clean fragments
            for atom in frag1.GetAtoms(): atom.ClearProp('molAtomMapNumber')
            frag1 = Chem.MolToSmiles(frag1, allHsExplicit=True)
            for atom in frag2.GetAtoms(): atom.ClearProp('molAtomMapNumber')
            frag2 = Chem.MolToSmiles(frag2, allHsExplicit=True)

            # return key response elements
            response['valid'] = True
            response['fragment1'] = frag1
            response['fragment2'] = frag2
            responselist.append(response)