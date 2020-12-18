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

from optparse import OptionParser
from rdkit.Chem.rdChemReactions import ReactionFromSmarts, ReactionToRxnBlock
    
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
            parser.add_option("-r", "--rxn", dest="rxnfile", help="set reaction path [default: %default]", metavar="FILE")
            parser.add_option("-v", "--verbose", dest="verbose", action="count", help="set verbosity level [default: %default]")

            # set defaults
            parser.set_defaults(outfile="./mmpa.out", 
                                infile="./mmpa.in",
                                rxnfile="./mmpa.rxn",
                               )
            
            # process options
            (opts, args) = parser.parse_args(argv)
            if opts.verbose is not None and opts.verbose > 0:
                print("verbosity level = %d" % opts.verbose)
            if opts.infile:
                print("infile = %s" % opts.infile)
            if opts.outfile:
                print("outfile = %s" % opts.outfile)
            if opts.rxnfile:
                print("rxnfile = %s" % opts.rxnfile)
                
        # MAIN BODY #
            
        except Exception as e:
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
            
            # prepare potential atom-atom mappings and create correspondence graph
            mmp = MMP(line['Molecule_L'], line['Molecule_R'])
            mmp.createCorrespondence(penalty=3.0)
            mmp.findCliques()
            mmp.eliminateMCS()

            # write output
            writer.writerow([line['Context'], line['Molecule_L'], mmp.getFragment1(), line['Molecule_R'], mmp.getFragment2()])
                 
            # create reaction
            reaction = ReactionFromSmarts(mmp.getSmirks())

            # write reaction
            rxnfile.write(ReactionToRxnBlock(reaction))     
#            break       
            
        # close file handles
        infile.close()
        outfile.close()
        rxnfile.close()

            
            