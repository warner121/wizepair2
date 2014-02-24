'''
Created on 24 Nov 2012

@author: Dan
'''

__all__ = []
__version__ = 0.1
__date__ = '2012-11-19'
__updated__ = '2013-07-06'

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
            
        # clear mappings and initialise radii (assume all atoms are RECS)
        for atom in self._mol1.GetAtoms(): 
            atom.SetProp('molAtomRadius','0')
            atom.ClearProp('molAtomMapNumber')
            
    def setMol2(self, mol=None):
        
        # TODO(warner121@hotmail.com): be good to check here, I guess?
        if mol: self._mol2 = mol
            
        # clear mappings and initialise radii (assume all atoms are RECS)
        for atom in self._mol2.GetAtoms(): 
            atom.SetProp('molAtomRadius','0')
            atom.ClearProp('molAtomMapNumber')
                        
    def createCorrespondence(self, penalty=3.0):
            
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
        for map1 in self.nodes_iter():
            for map2 in self.nodes_iter():

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
                    for idx in Chem.FindAtomEnvironmentOfRadiusN(mol, 3, atom.GetIdx()):
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

if __name__ == '__main__':
    
        '''Command line options.'''
        
        program_name = os.path.basename(sys.argv[0])
        program_version = "v%f" %__version__
        program_build_date = "%s" % __updated__
     
        program_version_string = '%%prog %s (%s)' % (program_version, program_build_date)
        #program_usage = '''usage: spam two eggs''' # optional - will be autogenerated by optparse
        program_longdesc = '''''' # optional - give further explanation about what the program does
        program_license = "Copyright 2013 Daniel Warner                                            \
                    Licensed under the Apache License 2.0\nhttp://www.apache.org/licenses/LICENSE-2.0"    

        argv = sys.argv[1:]
        try:
            # setup option parser
            parser = OptionParser(version=program_version_string, epilog=program_longdesc, description=program_license)
            parser.add_option("-i", "--in", dest="infile", help="set input path [default: %default]", metavar="FILE")
            parser.add_option("-o", "--out", dest="outfile", help="set output path [default: %default]", metavar="FILE")
            parser.add_option("-r", "--rxn", dest="rxnfile", help="set reaction path [default: %default]", metavar="FILE")
            parser.add_option("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: %default]")

            # set defaults
            parser.set_defaults(outfile="./mmpa.out", 
                                infile="./mmpa.in",
                                rxnfile="./mmpa.rxn")
            
            # process options
            (opts, args) = parser.parse_args(argv)
            
            if opts.verbose > 0:
                print("verbosity level = %d" % opts.verbose)
            if opts.infile:
                print("infile = %s" % opts.infile)
            if opts.outfile:
                print("outfile = %s" % opts.outfile)
                
            # MAIN BODY #
            
        except Exception, e:
            indent = len(program_name) * " "
            sys.stderr.write(program_name + ": " + repr(e) + "\n")
            sys.stderr.write(indent + "  for help use --help")
        
        # open file handles
        outfile = open(opts.outfile, 'w')
        writer = csv.writer(outfile, delimiter=',')
        writer.writerow(['Context', 'Mol_L', 'Frag_L', 'Mol_R', 'Frag_R'])
        rxnfile = open(opts.rxnfile, 'w')        
        
        infile = open(opts.infile, 'r')
        reader = csv.DictReader(infile, delimiter=',')
        for line in reader:
                
            # define input molecules
            mol1 = Chem.MolFromSmiles(line['Molecule_L'])
            mol2 = Chem.MolFromSmiles(line['Molecule_R'])
                          
            # prepare potential atom-atom mappings and create correspondence graph
            mmp = MMP()
            mmp.setMol1(mol1)
            mmp.setMol2(mol2)
            mmp.createCorrespondence(penalty=3.0)
            
#            networkx.draw(g2)
#            plt.show()            
            
            # calculate maximal cliques and report output
#            print "timer1: ", timeit.timeit('cliques = list(find_cliques(g))', 'from networkx.algorithms.clique import find_cliques; from __main__ import g', number = 1)
             
            # score the cliques and isolate RECS
#            print "timer2: ", timeit.timeit('MMP(mol1,mol2).scoreCliques(list(find_cliques(g)))', 'from networkx.algorithms.clique import find_cliques; from __main__ import mol1, mol2, g, MMP', number = 1)
            cliques = list(find_cliques(mmp))
            mmp.scoreCliques(cliques) 
            mmp.eliminateMCS()
                                    
            # write output
            writer.writerow([line['Context'], line['Molecule_L'], mmp.getFragmentA(), line['Molecule_R'], mmp.getFragmentB()])
                     
#            # create reaction
#            frag1 = Chem.rdmolops.AddHs(mmp.frag1)
#            frag2 = Chem.rdmolops.AddHs(mmp.frag2)
#            reaction = ChemicalReaction()
#            reaction.AddReactantTemplate(frag1)
#            reaction.AddProductTemplate(frag2)
#            reaction.Initialize()
#            
#            # test reaction (this can be slow and run out of memory)
#            mol1 = Chem.rdmolops.AddHs(mol1)
#            for prod in reaction.RunReactants((mol1,))[0]:
#                for atom in prod.GetAtoms():
#                    atom.ClearProp('molAtomMapNumber')
#                prod = Chem.rdmolops.RemoveHs(prod)
#                print Chem.MolToSmiles(prod)
#                
#            # write reaction
#            rxnfile.write(ReactionToRxnBlock(reaction))     
            
#            break       
            
        # close file handles
        infile.close()
        outfile.close()
        rxnfile.close()

            
            