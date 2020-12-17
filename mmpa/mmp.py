import csv
#import pprint
import sys
import os
import networkx
#import matplotlib.pyplot as plt

from optparse import OptionParser
from rdkit import Chem
from rdkit.Chem.rdChemReactions import ChemicalReaction, ReactionToRxnBlock
from networkx.algorithms.clique import find_cliques
import timeit
    
class MMP(networkx.Graph):

    def setMol1(self, mol=None):
        
        # TODO(warner121@hotmail.com): be good to check here, I guess?
        if mol: self._mol1 = mol
            
        # add hydrogen where defining isomer
        isomerics = []
        for atom in self._mol1.GetAtoms():
            if not atom.HasProp('_CIPCode'): continue
            isomerics.append(atom.GetIdx())
        self._mol1 = Chem.AddHs(self._mol1, onlyOnAtoms=isomerics, explicitOnly=True)
                      
        # clear mappings and initialise radii (assume all atoms are RECS)
        for atom in self._mol1.GetAtoms(): 
            atom.SetProp('molAtomRadius','0')
            atom.ClearProp('molAtomMapNumber')
  
    def setMol2(self, mol=None):
        
        # TODO(warner121@hotmail.com): be good to check here, I guess?
        if mol: self._mol2 = mol
            
        # add hydrogen where defining isomer
        isomerics = []
        for atom in self._mol2.GetAtoms():
            if not atom.HasProp('_CIPCode'): continue
            isomerics.append(atom.GetIdx())
        self._mol2 = Chem.AddHs(self._mol2, onlyOnAtoms=isomerics, explicitOnly=True)
            
        # clear mappings and initialise radii (assume all atoms are RECS)
        for atom in self._mol2.GetAtoms(): 
            atom.SetProp('molAtomRadius','0')
            atom.ClearProp('molAtomMapNumber')
                                          
    def createCorrespondence(self, penalty=3.0):
        
        mol1 = self._mol1
        mol2 = self._mol2
            
        for atom1 in mol1.GetAtoms():
            for atom2 in mol2.GetAtoms():
                
                # store the CIP codes somewhere that doesn't throw errors on comparison when missing
                try: atom1._CIPCode = atom1.GetProp('_CIPCode')
                except KeyError: atom1._CIPCode = None
                try: atom2._CIPCode = atom2.GetProp('_CIPCode')
                except KeyError: atom2._CIPCode = None            
    
                # set penalties - 3 strikes and you're out!
                __tempscore = 0
                if atom1.GetImplicitValence() != atom2.GetImplicitValence(): __tempscore += 1
                if atom1.GetExplicitValence() != atom2.GetExplicitValence(): __tempscore += 1
                if atom1.GetAtomicNum() != atom2.GetAtomicNum(): __tempscore += 1
                if atom1.GetDegree() != atom2.GetDegree(): __tempscore += 1
                if atom1.IsInRing() != atom2.IsInRing(): __tempscore += 1
                if atom1._CIPCode != atom2._CIPCode: __tempscore += 1
                
                # set upper limit on penalty to 1
                __tempscore = min(1, __tempscore/penalty)              
                mapping = (atom1.GetIdx(), atom2.GetIdx(), __tempscore)
                if __tempscore < 1: self.add_node(mapping) 
        
        # calculate distance matrices
        __dmat1 = Chem.GetDistanceMatrix(mol1)
        __dmat2 = Chem.GetDistanceMatrix(mol2)  

        # create correspondance graph edges
        for map1 in self.nodes():
            for map2 in self.nodes():

                # test if criteria are met for correspondance
                correspondance = __dmat1[map1[0]][map2[0]] == __dmat2[map1[1]][map2[1]]
                if correspondance: self.add_edge(map1, map2)        
    
    def __setMappings(self, clique):        

        # iterate over the remaining mappings identified from MCSS
        for pair in clique:
            
            # increment the index to prevent atom mappings of 0
            mapIdx = clique.index(pair) + 1
            
            # map first atom and set radius to 99 (atom part of MCS)
            atom1 = self._mol1.GetAtomWithIdx(pair[0])
            atom1.SetProp('molAtomMapNumber','%d'%mapIdx)
            
            # map second atom and set radius to 99 (atom part of MCS)
            atom2 = self._mol2.GetAtomWithIdx(pair[1])
            atom2.SetProp('molAtomMapNumber','%d'%mapIdx)
    
    def __flagMCS(self, mcs):

        # iterate over the remaining mappings identified from MCSS
        for pair in mcs:
            
            # map first atom and set radius to 99 (atom part of MCS)
            atom1 = self._mol1.GetAtomWithIdx(pair[0])
            atom1.SetProp('molAtomRadius','99')
            
            # map second atom and set radius to 99 (atom part of MCS)
            atom2 = self._mol2.GetAtomWithIdx(pair[1])
            atom2.SetProp('molAtomRadius','99')
                                            
    def scoreCliques(self, cliques):
        
        # initialise
        bestscore = -1E800
        bestmcs = None    
        bestclique = None            
        
        # score largest cliques first
        cliques.sort(key=len, reverse=True)
        for clique in cliques:   
            
            # initialize scores
            score = len(clique)
            mcs = list(clique)
            if score <= bestscore: continue # score cannot be greater than clique size
    
            # iterate over the mappings identified from MCSS
            for pair in clique:
                score = score - pair[2] # lookup predefined score penalty
                if pair[2]: mcs.remove(pair) # flag atoms as RECS in addition to those not mapped

            # store results if best so far
            if score > bestscore:
                bestscore = score
                bestmcs = mcs
                bestclique = clique  

        # map all atoms in clique and flag MCS
        self.__setMappings(bestclique)
        self.__flagMCS(bestmcs)
        
    def eliminateMCS(self):
        
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
#            frag.Debug()  
#            Chem.AssignStereochemistry(frag, cleanIt=True, force=True)
            return frag
      
        self.frag1 = eliminate(self._mol1)
        self.frag2 = eliminate(self._mol2)
        
    def getFragmentA(self):
        
        frag = Chem.EditableMol(self.frag1).GetMol()
        for atom in frag.GetAtoms(): atom.ClearProp('molAtomMapNumber')
        return Chem.MolToSmiles(frag)
    
    def getFragmentB(self):
        
        frag = Chem.EditableMol(self.frag2).GetMol()
        for atom in frag.GetAtoms(): atom.ClearProp('molAtomMapNumber')
        return Chem.MolToSmiles(frag)        
