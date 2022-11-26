import json
import re
import logging
import networkx as nx
import numpy as np
import timeit

from hashlib import md5
from rdkit import Chem, RDLogger
from rdkit.Chem import SaltRemover
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
from func_timeout import func_timeout, FunctionTimedOut
from networkx.algorithms.clique import enumerate_all_cliques, find_cliques, find_cliques_recursive, max_weight_clique
from multiprocessing import Pool

# disable C++ logger for production
RDLogger.DisableLog('rdApp.*')

class CorrespondenceGraph(nx.Graph):
    '''
    Build the correspondence matrix of putative atom pairings, from which to determine the maximal clique (comprising the MCS). 
    '''
    
    def __init__(self):

        # inherit from nx.Graph
        super().__init__(self)

    def build(
        self, 
        mol1, 
        mol2, 
        strictness, 
        sfunc=np.full(8, 10)):
        '''
        Build the correspondence graph.
        
        Each atomic pairing is assigned a providional score, from 0 (most different) to 8000 (identical). 
        This pairwise score is appended to the tuple representing each node, following the atomic indices.
        '''
        # ensure scoring function suitable
        assert len(sfunc) == 8
        assert (sfunc > 1).all()
        
        # store strictness for scoring
        self._mol1 = mol1
        self._mol2 = mol2
        self._strict = (sfunc.mean() * strictness) ** 2 # np.sort(sfunc)[:5].sum() ** 2
        
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
                        
        # create lookup for atomic index mappings
        self._idxmaps = []

        # iterate over all potential atom-atom pairings
        for atom1 in self._mol1.GetAtoms():
            for atom2 in self._mol2.GetAtoms():

                # score putative nodes based on atom:atom similarity (0-88 + up to 10 point centrality bonus)
                score = np.zeros(8)
                if atom1.GetAtomicNum() == atom2.GetAtomicNum(): score[0] = 1
                if atom1.GetImplicitValence() == atom2.GetImplicitValence(): score[1] = 1
                if atom1.GetExplicitValence() == atom2.GetExplicitValence(): score[2] = 1
                if atom1.GetFormalCharge() == atom2.GetFormalCharge(): score[3] = 1
                if atom1.GetIsAromatic() == atom2.GetIsAromatic(): score[4] = 1
                if atom1.GetDegree() == atom2.GetDegree(): score[5] = 1
                if atom1.IsInRing() == atom2.IsInRing(): score[6] = 1
                if getCIPCode(atom1) == getCIPCode(atom2): score[7] = 1
                score = (score * sfunc).sum()

                # apply jitter in the event of a tie
                peripherality1 = getPeripherality(atom1, self._dmat1)
                peripherality2 = getPeripherality(atom2, self._dmat2)
                score += sfunc.min() * (1 - abs(peripherality1 - peripherality2)) * 0.999
                score = int(np.floor(score**2))
                
                # accept node with greater than specified match level
                if score >= self._strict: 
                    newmap = (atom1.GetIdx(), atom2.GetIdx())
                    self._idxmaps.append(newmap)
                    self.add_node(self._idxmaps.index(newmap), weight=score)
                    
        # build numpy matrices for weights
        self._nodeweights = np.array([self.nodes[x]['weight'] for x in self.nodes()], dtype='int64') / (sfunc.sum() ** 2)
        #self._edgeweights = np.zeros((len(self.nodes), len(self.nodes)), dtype='int64')
        
        # create correspondence graph edges
        for node1 in self.nodes():
            map1 = self._idxmaps[node1]
            for node2 in self.nodes():
                map2 = self._idxmaps[node2]
                
                # only build 1/2 matrix
                if node1 > node2: continue
                                
                # ensure any given atom is not mapped twice in a clique
                if map1[0] == map2[0] or map1[1] == map2[1]: continue

                # test if criteria are met for correspondence
                #correspondence = abs(self._dmat1[map1[0]][map2[0]] - self._dmat2[map1[1]][map2[1]])
                #score = int(np.floor(1000/((1+correspondence)**2)))

                # check comparative distance between mapped atoms is within tolerance
                if (2/3) <= self._dmat1[map1[0]][map2[0]] / self._dmat2[map1[1]][map2[1]] <= (3/2):
                    self.add_edge(node1, node2, weight=0)
                    #self._edgeweights[node1][node2] = 0
                    
        # get weighted degrees
        self._embedding = [np.log(self._nodeweights[node] * val) for (node, val) in self.degree()]
        boundaries = np.sort(np.concatenate((np.linspace(3, 7, 15), [np.inf, np.NINF]), axis=0))
        self._embedding = np.histogram(self._embedding, bins=boundaries)[0].tolist()
        
        # predict solution time using simple linear model
        self._predsolversecs = np.sum(self._embedding * np.array([
            0.00136326, 0.00136326, 0.00188707, 0.00241089, 0.0029347 ,
            0.00345852, 0.00398233, 0.00450614, 0.01288632, 0.01488049,
            0.02397065, 0.02899514, 0.07183043, 0.41705952, 0.7622886 ,
            0.7622886 ]))
           
    def score_clique(self, clique):
            
        # lookup scores from matrices
        score = np.sum(self._nodeweights[clique])
        #score += np.sum(self._edgeweights[clique].T[clique])   
        return score

    def filter_mcs(self, clique, clipper=2):
        
        # remove any atom pairings with less than perfect score, excluding bonus i.e. (11*8)**2 = 7744
        mcs = [x for x in clique if self._nodeweights[x] >= 1] 
        
        # replace integer node numbers with atomic index tuples
        mcs = [self._idxmaps[x] for x in mcs]
        clique = [self._idxmaps[x] for x in clique]
        if not len(mcs): return clique, mcs
        
        # split the tuples and homogenise the distance matrices
        idx1, idx2 = zip(*mcs)
        dmat1 = self._dmat1[list(idx1)].T[list(idx1)]
        dmat2 = self._dmat2[list(idx2)].T[list(idx2)]       
        
        # take the difference in the reduced distance matries
        dmatdiff = np.clip(dmat1, 0, clipper) - np.clip(dmat2, 0, clipper)
        idx1 = set([idx1[x] for x in np.where(dmatdiff != 0)[0]])
        idx2 = set([idx2[x] for x in np.where(dmatdiff != 0)[1]])

        # retain those nodes where the constituent atoms have not drifted (or explicit i.e. chiral H)
        mcs = [x for x in mcs
               if (x[0] not in idx1 or self._mol1.GetAtomWithIdx(x[0]).GetAtomicNum() == 1)
               and (x[1] not in idx2 or self._mol2.GetAtomWithIdx(x[1]).GetAtomicNum() == 1)]
        
        # return
        return clique, mcs
                            
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

        # lookup scores from matrices (purely for comparison with return from max_weight_cliques)
        score = self.score_clique(clique)

        # remap to indices and remove atomic/drift based discrepancies from mcs
        clique, mcs = self.filter_mcs(clique)
        
        # return results
        return clique, mcs

class Reactor():
    
    def __init__(self, smirks):        
        '''
        Instantiate MMP 'Reactor'.
        
        smirks: SMIRKS encoded reaction.
        '''
        
        self._rxn = Chem.rdChemReactions.ReactionFromSmarts(smirks)

    def assert_one2one(self):
        '''
        Assert 1:1 relationship between reactants and products.
        '''
            
        try: 
            assert self._rxn.GetNumReactantTemplates() == 1
            assert self._rxn.GetNumProductTemplates() == 1
            return True
        except AssertionError:
            logging.info(json.dumps({"message": "no 1:1 reaction could be generated"}))
            return False
     
    def generate_products(self, smiles):
        '''
        Return products as list of SMILES.
        
        smiles: SMILES to serve as seed or reactant.
        '''
            
        reactant = Chem.AddHs(Chem.MolFromSmiles(smiles))
        products = self._rxn.RunReactants((reactant,))
        productset = set()
        for product in products:
            try:
                productparts = [Chem.MolToSmiles(Chem.RemoveHs(productpart)) for productpart in product]
                productset.add('.'.join(productparts))
            except (Chem.AtomValenceException, Chem.AtomKekulizeException, Chem.KekulizeException):
                logging.info(json.dumps({"message": "MolSanitizeException raised on product enumeration"}))
        return list(productset)
    
class SMIRKSEncoder():
    
    def __init__(self):
        
        return None
        
    def encode(self, frag1, frag2, smiles1, smiles2):
        '''
        Returns tuple of (
            smirks: str, 
            valid: bool, 
            error: str,
            biproducts: int)
        '''

        def getMolAtomMapNumber(atom):
            try: return atom.GetProp('molAtomMapNumber')
            except KeyError: return None

        # save backup in case this fails validation
        smarts1 = Chem.MolToSmarts(frag1)
        smarts2 = Chem.MolToSmarts(frag2)
        backup = '{}>>{}'.format(smarts1, smarts2)

        # extract index/mapping lookup of atom all (including None) mappings
        mappings1 = list(enumerate([getMolAtomMapNumber(atom) for atom in frag1.GetAtoms()]))
        mappings2 = list(enumerate([getMolAtomMapNumber(atom) for atom in frag2.GetAtoms()]))

        # strip of mappings and call MolToSmiles to define _smilesAtomOutputOrder
        for atom in frag1.GetAtoms(): atom.ClearProp('molAtomMapNumber')
        for atom in frag2.GetAtoms(): atom.ClearProp('molAtomMapNumber')
        junk1 = Chem.MolToSmiles(frag1)
        junk2 = Chem.MolToSmiles(frag2)

        # put the original mappings back
        for mapping in mappings1:
            if not mapping[1]: continue
            frag1.GetAtomWithIdx(mapping[0]).SetIntProp('molAtomMapNumber', int(mapping[1]))
        for mapping in mappings2:
            if not mapping[1]: continue
            frag2.GetAtomWithIdx(mapping[0]).SetIntProp('molAtomMapNumber', int(mapping[1]))

        # renumber according to mapping-free output order
        frag1 = Chem.RenumberAtoms(frag1, frag1.GetPropsAsDict(True,True)["_smilesAtomOutputOrder"])
        frag2 = Chem.RenumberAtoms(frag2, frag2.GetPropsAsDict(True,True)["_smilesAtomOutputOrder"])
        smarts1 = Chem.MolToSmarts(frag1)
        smarts2 = Chem.MolToSmarts(frag2)

        # arbitrarily define dominant fragment and 
        hashstr1 = re.sub('(?<=:)[0-9]+(?=])', 'X', smarts1)
        hashstr2 = re.sub('(?<=:)[0-9]+(?=])', 'X', smarts2)
        if hashstr1 > hashstr2: lookup = re.findall('(?<=:)[0-9]+(?=])', smarts1)
        else: lookup = re.findall('(?<=:)[0-9]+(?=])', smarts2)

        # replace with sequentially ordered mappings from 1
        smirks = '{}>>{}'.format(smarts1, smarts2)
        for idx, swap in enumerate(lookup):
            smirks = re.sub(':{}]'.format(swap), ':X{}]'.format(idx+1), smirks)
            backup = re.sub(':{}]'.format(swap), ':X{}]'.format(idx+1), backup)
        smirks = re.sub(':X', ':', smirks)
        backup = re.sub(':X', ':', backup)

        # finally add missing +0 charges where mismatched
        chargemaps1 = set(re.findall('(?<=[+-]:)[0-9]+(?=])', smarts1))
        chargemaps2 = set(re.findall('(?<=[+-]:)[0-9]+(?=])', smarts2))
        for chargemap in chargemaps1.symmetric_difference(chargemaps2):
            idx = lookup.index(chargemap)
            smirks = re.sub('(?<=[0-9]):{}]'.format(idx+1), '+0:{}]'.format(idx+1), smirks)
            backup = re.sub('(?<=[0-9]):{}]'.format(idx+1), '+0:{}]'.format(idx+1), backup)

        # verify 1:1 reaction
        reactor = Reactor(smirks)
        if not reactor.assert_one2one():
            return smirks, False, 'not one2one reaction', None

        # verify derived reaction produces original 'product'
        for test in [smirks, backup]:
            reactor = Reactor(test)
            productlist = reactor.generate_products(smiles1)
            if smiles2 in productlist:
                return test, True, None, len(productlist) - 1
        return backup, False, 'second molecule not found amongst products enumerated from first', len(productlist)
    
class MMP():

    @staticmethod
    def __saltRemover(smiles: str):
        
        # parse smiles
        try: mol = Chem.MolFromSmiles(smiles)
        except: return None
    
        # remove salts
        remover = SaltRemover.SaltRemover()
        mol, salts = remover.StripMolWithDeleted(mol)
        
        # retain largest fragment
        largest = Chem.Mol()
        for submol in Chem.GetMolFrags(mol, asMols=True): 
            if submol.GetNumHeavyAtoms() <= largest.GetNumHeavyAtoms(): continue
            largest = submol
        mol = largest
        
        # return
        return mol
    
    @staticmethod
    def __prepareMol(mol: Chem.Mol):

        # add hydrogen where defining isomer
        isomerics = []
        for atom in mol.GetAtoms():
            if not atom.HasProp('_CIPCode'): continue
            isomerics.append(atom.GetIdx())
        mol = Chem.AddHs(mol, onlyOnAtoms=isomerics, explicitOnly=True)
                      
        # clear mappings and initialise radii (assume all atoms are RECS)
        for atom in mol.GetAtoms(): 
            atom.SetIntProp('molAtomRadius', 0)
            atom.ClearProp('molAtomMapNumber')
            
        # return
        return mol
    
    def __init__(self, smiles_x: str, smiles_y: str, strictness=4):
        '''
        Initialise the matched molecular pair.
        
        smiles_x: First molecule to compare.
        smiles_y: Second molecule to compare.
        strictness: Integer (1-8) to indicate how tolerant the algortithm should to be to atom-wise chemical differences. 
            1 (slowest) all atom types match.
            8 (fastest) atoms chemically identical to be considered part of mcss.  
        '''
        
        if strictness-1 not in range(8): return
        
        # initialise molecules for comparison
        self._mol1 = self.__saltRemover(smiles_x)
        self._mol2 = self.__saltRemover(smiles_y)

        # canonicalise salt-free smiles
        self._smiles1 = Chem.MolToSmiles(self._mol1)
        self._smiles2 = Chem.MolToSmiles(self._mol2)
                 
        # add chiral Hs etc. ahead of MCSS
        self._mol1 = self.__prepareMol(self._mol1)
        self._mol2 = self.__prepareMol(self._mol2)

        # intialise correspondence graph
        self._graph = CorrespondenceGraph()
        self._graph.build(self._mol1, self._mol2, strictness)
        
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
            atom1.SetIntProp('molAtomMapNumber', mapIdx)
            
            # map second atom and set radius to 99 (atom part of MCS)
            atom2 = self._mol2.GetAtomWithIdx(pair[1])
            atom2.SetIntProp('molAtomMapNumber', mapIdx)

    def __setAtomRadii(self):
        '''
        Use the atomic radii to denote which atoms are part of the MCS. By elimination, those atoms with radii of 0 will form the RECS.
        '''

        # iterate over the remaining mappings identified from MCSS
        for pair in self._mcs:
            
            # map first atom and set radius to 99 (atom part of MCS)
            atom1 = self._mol1.GetAtomWithIdx(pair[0])
            atom1.SetIntProp('molAtomRadius',99)
            
            # map second atom and set radius to 99 (atom part of MCS)
            atom2 = self._mol2.GetAtomWithIdx(pair[1])
            atom2.SetIntProp('molAtomRadius',99)

    def execute(self, radii=4):
        '''
        solver = max_weight_clique
        '''
        
        # predict timeout
        if self._graph._predsolversecs > np.inf: 
            return [{
                'embedding': self._graph._embedding,
                'predsolversecs': self._graph._predsolversecs,
                'error': 'timeout expected - skipping'
            }]

        # find the MCS
        self._solversecs = timeit.default_timer()
        self._clique, self._mcs = self._graph.solve_weighted()
        self._solversecs = timeit.default_timer() - self._solversecs

        # determine the % of largest molecule covered by MCS
        maxnumatoms = max(self._mol1.GetNumAtoms(), self._mol2.GetNumAtoms())
        if not maxnumatoms: 
            return [{
                'valid': False,
                'error': 'neither mol has any non-salt atoms'
            }]
        self._percentmcs = len(self._mcs) / maxnumatoms 
        if not (0 < self._percentmcs < 1): 
            return [{
                'valid': False,
                'percentmcs': self._percentmcs,
                'error': 'either none, or complete common substructure'
            }]
        
        # search, mark up atom mappings and MCS/RECS split
        self.__setAtomMapNumbers()
        self.__setAtomRadii()
                
        # define function for elimination of MCS
        def eliminate(mol, radius):
            
            # preserve atoms within <radius> bonds of attachment
            toRemove = set(range(mol.GetNumAtoms()))
            for atom in mol.GetAtoms():
                if atom.GetIntProp('molAtomRadius') > 0: continue
                for x in reversed(range(radius+1)):
                    env = list(Chem.FindAtomEnvironmentOfRadiusN(mol, x, atom.GetIdx()))
                    if not len(env): continue
                    break
                for idx in env:
                    envBond = mol.GetBondWithIdx(idx)
                    toRemove.discard(envBond.GetBeginAtom().GetIdx())
                    toRemove.discard(envBond.GetEndAtom().GetIdx())
                if radius == 0:
                    toRemove.discard(atom.GetIdx())
                        
            # remove environment from core
            toRemove = list(toRemove)
            toRemove.sort(reverse=True)
            frag = Chem.EditableMol(mol)
            for atom in toRemove: frag.RemoveAtom(atom)
            frag = frag.GetMol()
            frag = Chem.AddHs(frag)
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
                        'valid': False,
                        'solversecs': self._solversecs,
                        'embedding': self._graph._embedding,
                        'predsolversecs': self._graph._predsolversecs,
                        'error': None}
            
            # Define reaction as SMIRKS while mappings still present
            frag1 = eliminate(self._mol1, radius)
            frag2 = eliminate(self._mol2, radius)   
            encoder = SMIRKSEncoder()
            smirks, valid, error, biproducts = encoder.encode(frag1, frag2, self._smiles1, self._smiles2)

            # remove mappings to yield clean fragments
            for atom in frag1.GetAtoms(): atom.ClearProp('molAtomMapNumber')
            for atom in frag2.GetAtoms(): atom.ClearProp('molAtomMapNumber')
            junk1 = Chem.MolToSmiles(frag1)
            junk2 = Chem.MolToSmiles(frag2)

            # renumber according to mapping-free output order
            frag1 = Chem.RenumberAtoms(frag1, frag1.GetPropsAsDict(True,True)["_smilesAtomOutputOrder"])
            frag2 = Chem.RenumberAtoms(frag2, frag2.GetPropsAsDict(True,True)["_smilesAtomOutputOrder"])
            frag1 = md5(Chem.MolToSmarts(frag1).encode()).hexdigest()
            frag2 = md5(Chem.MolToSmarts(frag2).encode()).hexdigest()
            
            # return key response elements
            response['smirks'] = smirks
            response['biproducts'] = biproducts
            response['valid'] = valid
            response['error'] = error
            response['fragment1'] = frag1
            response['fragment2'] = frag2
            responselist.append(response)