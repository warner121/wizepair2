'''
Created on 24 Nov 2012

@author: Dan
'''

__all__ = []
__version__ = 0.11
__date__ = '2012-11-19'
__updated__ = '2020-09-29'

import csv
import sys
import os
import networkx

from optparse import OptionParser
from rdkit import Chem
from rdkit.Chem.rdChemReactions import ChemicalReaction, ReactionToRxnBlock
from networkx.algorithms.clique import find_cliques
import timeit
    
from mmp import MMP

if __name__ == '__main__':
    
        '''Command line options.'''
        
        program_name = os.path.basename(sys.argv[0])
        program_version = "v%f" %__version__
        program_build_date = "%s" % __updated__
     
        program_version_string = '%%prog %s (%s)' % (program_version, program_build_date)
        program_longdesc = ''''''
        program_license = "Copyright 2013 Daniel Warner                                            \
                    Licensed under the Apache License 2.0\nhttp://www.apache.org/licenses/LICENSE-2.0"    

        argv = sys.argv[1:]
        try:
            # setup option parser
            parser = OptionParser(version=program_version_string, epilog=program_longdesc, description=program_license)
            parser.add_option("-i", "--in", dest="infile", help="set input path [default: %default]", metavar="FILE")
            parser.add_option("-o", "--out", dest="outfile", help="set output path [default: %default]", metavar="FILE")
#            parser.add_option("-r", "--rxn", dest="rxnfile", help="set reaction path [default: %default]", metavar="FILE")
            parser.add_option("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: %default]")

            # set defaults
            parser.set_defaults(outfile="./mmpa.out", 
                                infile="./mmpa.in",
#                                rxnfile="./mmpa.rxn",
                               )
            
            # process options
            (opts, args) = parser.parse_args(argv)
            if opts.verbose is not None and opts.verbose > 0:
                print("verbosity level = %d" % opts.verbose)
            if opts.infile:
                print("infile = %s" % opts.infile)
            if opts.outfile:
                print("outfile = %s" % opts.outfile)
                
        # MAIN BODY #
            
        except Exception as e:
            indent = len(program_name) * " "
            sys.stderr.write(program_name + ": " + repr(e) + "\n")
            sys.stderr.write(indent + "  for help use --help")
        
        # open file handles
        outfile = open(opts.outfile, 'w')
        writer = csv.writer(outfile, delimiter=',')
        writer.writerow(['Context', 'Mol_L', 'Frag_L', 'Mol_R', 'Frag_R'])
#        rxnfile = open(opts.rxnfile, 'w')        
        
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
            
            # score the cliques and isolate RECS
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
#        rxnfile.close()

            
            