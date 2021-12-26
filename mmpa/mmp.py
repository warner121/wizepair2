import json
import logging
import networkx
import numpy as np

from rdkit import Chem
from func_timeout import func_timeout, FunctionTimedOut
from networkx.algorithms.clique import enumerate_all_cliques, find_cliques, find_cliques_recursive, max_weight_clique

class CorrespondenceGraph(networkx.Graph):
    '''
    Build the correspondence matrix of putative atom pairings, from which to determine the maximal clique (comprising the MCS). 
    '''
    
    def __init__(self):

        # inherit from networkx.Graph
        super().__init__(self)

    def build(self, mol1, mol2, strictness, correspondence):
        '''
        Build the correspondence graph.
        
        Each atomic pairing is assigned a providional score, from 0 (most different) to 8000 (identical). 
        This pairwise score is appended to the tuple representing each node, following the atomic indices.
        '''
        
        # store strictness for scoring
        self._strict = strictness * 1000
        self._corr = correspondence
        
        # calculate distance matrices
        self._dmat1 = Chem.GetDistanceMatrix(mol1)
        self._dmat2 = Chem.GetDistanceMatrix(mol2)
        
        # extract propery in such a way error is not thrown on comparison
        def getCIPCode(atom):
            try: return atom.GetProp('_CIPCode')
            except KeyError: return None
            
        # create description of how central in molecule atom
        def getPeripherality(atom, dmat):
            peripherality = dmat[atom.GetIdx()]
            return np.mean(peripherality / np.max(peripherality))            
                        
        # build numpy matrix for atom/atom combos
        self._atommaps = []

        # iterate over all potential atom-atom pairings
        for atom1 in mol1.GetAtoms():
            for atom2 in mol2.GetAtoms():

                # set penalties - [penalty] strikes and you're out!
                # TODO: move atomic properies to graph 
                score = 0
                if atom1.GetAtomicNum() == atom2.GetAtomicNum(): score += 1000
                if atom1.GetImplicitValence() == atom2.GetImplicitValence(): score += 1000
                if atom1.GetExplicitValence() == atom2.GetExplicitValence(): score += 1000
                if atom1.GetFormalCharge() == atom2.GetFormalCharge(): score += 1000
                if atom1.GetIsAromatic() == atom2.GetIsAromatic(): score += 1000
                if atom1.GetDegree() == atom2.GetDegree(): score += 1000
                if atom1.IsInRing() == atom2.IsInRing(): score += 1000
                if getCIPCode(atom1) == getCIPCode(atom2): score += 1000

                # apply jitter in the event of a tie
                peripherality1 = getPeripherality(atom1, self._dmat1)
                peripherality2 = getPeripherality(atom2, self._dmat2)
                score += 999 - int(np.floor(999 * abs(peripherality1 - peripherality2)))
                
                # accept node with greater than specified match level
                if score >= self._strict: 
                    newmap = (atom1.GetIdx(), atom2.GetIdx())
                    self._atommaps.append(newmap)
                    self.add_node(self._atommaps.index(newmap), weight=score)
                    
        # build numpy matrix for correspondance combos
        self._corrmat = np.zeros((len(self.nodes), len(self.nodes)), dtype='int64')
        
        # create correspondence graph edges
        for node1 in self.nodes():
            map1 = self._atommaps[node1]
            for node2 in self.nodes():
                map2 = self._atommaps[node2]
                                
                # ensure any given atom is not mapped twice in a clique
                if map1[0] == map2[0] or map1[1] == map2[1]: continue

                # test if criteria are met for correspondence
                correspondence = abs(self._dmat1[map1[0]][map2[0]] - self._dmat2[map1[1]][map2[1]])
                score = int(np.floor(1000/((1+correspondence)**2)))

                if correspondence < self._corr: 
                    self.add_edge(node1, node2, weight=score)
                    self._corrmat[node1][node2] = score

    def solve(self, solver, timeout=60):
        '''
        Enumerate cliques and rescore:
            solver=find_cliques - analyse maximal cliques only (default)
            solver=find_cliques_recursive - (recursively) analyse maximal cliques only
            solver=enumerate_all_cliques - analyse all cliques (not recommended for performance reasons)
        
        Returns:
            bestclique - the best scoring clique (including partially matching atoms, and used to derive reaction mappings)
            bestmcs - the subset of bestclique containing exact chemical matches only (to be discarded to produce the RECS)
        '''
                    
        # define function with no arguments (for use with timeout function)
        def findCliquesNoArgs(): return list(solver(self))
        
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
            
            # lookup scores from matrices
            score = np.sum(self._corrmat[clique].T[clique])/2
            score += sum([self.nodes[mapping]['weight'] for mapping in clique])
            
            # store results if best so far
            if score > bestscore:
                bestscore = score
                bestclique = clique  

        # work out which atoms comprising the best clique neighbour atoms which remain unpaired
        bestmcs = [x for x in bestclique if self.nodes[x]['weight'] >= 8000] 
        mcs1, mcs2 = zip(*[self._atommaps[x] for x in bestmcs])
        dmat1 = np.delete(self._dmat1, list(mcs1), axis=1)
        dmat2 = np.delete(self._dmat2, list(mcs2), axis=1)
        if dmat1.size != 0: idx1 = np.where(dmat1.min(axis=1) > 1)[0]
        else: idx1 = mcs1
        if dmat2.size != 0: idx2 = np.where(dmat2.min(axis=1) > 1)[0]
        else: idx2 = mcs2
            
        # eliminate those nodes where the constituent atoms do not neighbour outsiders or have atomic differences
        bestmcs = [x for x in bestmcs if self._atommaps[x][0] in idx1 and self._atommaps[x][1] in idx2]
        
        # replace indices with actual mappings for downstream
        bestmcs = [self._atommaps[x] for x in bestmcs]
        bestclique = [self._atommaps[x] for x in bestclique]

        # return results
        return bestclique, bestmcs
                            
    def solve_weighted(self, timeout=60):
        '''
        For testing only, please do not use.
        
        Execute the maximum weight clique search.
        
        Returns:
            clique - the best scoring clique (including partially matching atoms, and used to derive reaction mappings)
            mcs - the subset of bestclique containing exact chemical matches only (to be discarded to produce the RECS)
        '''
                    
        # define function with no arguments (for use with timeout function)
        def findCliquesNoArgs(): return max_weight_clique(self)
        
        # try finding cliques within [timeout] seconds
        try:
            clique, maxweight = func_timeout(timeout, findCliquesNoArgs)
        except FunctionTimedOut:
            logging.warning(json.dumps({"message": "failed to find cliques in {} seconds".format(timeout)}))
            return list(), list()

        # lookup scores from matrices
        score = np.sum(self._corrmat[clique].T[clique])/2
        score += sum([self.nodes[mapping]['weight'] for mapping in clique])

        # work out which atoms comprising the best clique neighbour atoms which remain unpaired
        mcs = [x for x in clique if self.nodes[x]['weight'] >= 8000] 
        mcs1, mcs2 = zip(*[self._atommaps[x] for x in mcs])
        dmat1 = np.delete(self._dmat1, list(mcs1), axis=1)
        dmat2 = np.delete(self._dmat2, list(mcs2), axis=1)
        if dmat1.size != 0: idx1 = np.where(dmat1.min(axis=1) > 1)[0]
        else: idx1 = mcs1
        if dmat2.size != 0: idx2 = np.where(dmat2.min(axis=1) > 1)[0]
        else: idx2 = mcs2
            
        # eliminate those nodes where the constituent atoms do not neighbour outsiders or have atomic differences
        mcs = [x for x in mcs if self._atommaps[x][0] in idx1 and self._atommaps[x][1] in idx2]
        
        # replace indices with actual mappings for downstream
        mcs = [self._atommaps[x] for x in mcs]
        clique = [self._atommaps[x] for x in clique]
    
        # return results
        return clique, mcs
        
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
    
    def __init__(self, smiles_x: str, smiles_y: str, strictness=4, correspondence=1):
        '''
        Initialise the matched molecular pair.
        
        smiles_x: First molecule to compare.
        smiles_y: Second molecule to compare.
        strictness: Integer (1-8) to indicate how tolerant the algortithm should to be to atom-wise chemical differences. 
            1 (slowest) all atom types match.
            8 (fastest) atoms chemically identical to be considered part of mcss.  
        correspondence: Integer (1-8) to indicate how tolerant the algortithm should to be to topological differences. 
            1 (fastest) standard MCS using exact correspondance matrix only.
            4 (slowest) atoms are allowed to 'drift' up to [correspondance] bonds away from neighbouring counterparts.  
        '''
        
        if strictness-1 not in range(8): return
        if correspondence-1 not in range(4): return
        
        # canonicalise smiles
        self._smiles1 = Chem.MolToSmiles(Chem.MolFromSmiles(smiles_x))
        self._smiles2 = Chem.MolToSmiles(Chem.MolFromSmiles(smiles_y))
         
        # initialise molecules for comparison
        self._mol1 = self.__molFromSmiles(self._smiles1)
        self._mol2 = self.__molFromSmiles(self._smiles2)
        
        # intialise correspondence graph
        self._graph = CorrespondenceGraph()
        self._graph.build(self._mol1, self._mol2, strictness, correspondence)
        
        # dummy vars
        self._clique = None
        self._mcs = None

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

    def execute(self, radii=4, solver=find_cliques):
        '''
        solver = find_cliques, find_cliques_recursive, enumerate_all_cliques, max_weight_clique
        '''

        # find the MCS
        if solver == max_weight_clique: self._clique, self._mcs = self._graph.solve_weighted()
        else: self._clique, self._mcs = self._graph.solve(solver=solver)

        # determine the % of largest molecule covered by MCS
        self._percentmcs = len(self._mcs) / max(self._mol1.GetNumAtoms(), self._mol2.GetNumAtoms())        

        # search, mark up atom mappings and MCS/RECS split
        self.__setAtomMapNumbers()
        self.__setAtomRadii()
                
        # define function for elimination of MCS
        def eliminate(mol, radius):
            
            # environment fails if radius > max distance
            radius = int(min(radius, np.max(Chem.GetDistanceMatrix(mol))))
        
            # tag atoms within 4 bonds of attachment
            toRemove = set(range(mol.GetNumAtoms()))
            for atom in mol.GetAtoms():
                if atom.GetProp('molAtomRadius') == '0':
                    for idx in Chem.FindAtomEnvironmentOfRadiusN(mol, radius-1, atom.GetIdx()):
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
            try:
                for product in productset:
                    productlist.append('.'.join([Chem.MolToSmiles(Chem.RemoveHs(productpart)) for productpart in product]))
            except Chem.KekulizeException:
                print(self._smiles1, self._smiles2, smirks)
            except Chem.AtomValenceException:
                print(self._smiles1, self._smiles2, smirks)
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