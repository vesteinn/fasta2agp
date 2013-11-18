#!/usr/bin/env python

# Reads in a fasta file containing scaffolds
# and outputs separated contigs along with agp
# file in accordance with NCBI's AGP Spec. v2.0.
#
# Use at own risk, this script can be found at
# github.com/vesteinn.

import settings

from sys import argv
import re
import textwrap


class Fasta2AGP(object):

    contig_counter = 0
    scaff_name_counter = 1

    def __init__(self, fastafile, outfile):
        self.ifile = open(fastafile)
        self.ffile = open(outfile + '.fasta', 'w')
        self.agpfile = open(outfile + '.agp', 'w')

    def scaffold2parts(self, scaffname, scaffold):
        min_n_string = "n" * settings.MIN_N + '*'
        lower_scaff = scaffold.lower()
        contigs = re.split(min_n_string, lower_scaff)
        gaps = list(re.finditer(min_n_string, lower_scaff))
        contig_start = 1
        outdata = ""
        i = 0
        gapcounter = 0
        contcounter = 0
        totalgaps = len(gaps)
        print str(totalgaps) + " gaps found in scaffold " + scaffname
        components = len(contigs) + len(gaps)
        while i < components:
            if gapcounter <= totalgaps - 1:
                # Gap to the right of the contig
                gap = gaps[gapcounter]
                contig_end = str(gap.start() - 1 + 1)
            else:
                # Contig at the end of a scaffold
                contig_end = str(len(scaffold))
            contig_name = settings.COMPONENT_PRE + str(i + 1)
            print contig_name + " in scaffname."
            # We print the contig info
            outdata += ('\t').join([scaffname,
                                str(contig_start),
                                contig_end,
                                str(i + 1),
                                settings.COMPONENT_TYPE,
                                contig_name,
                                str(1),
                                contig_end,
                                '+\n'])
            i += 1
            self.ffile.writelines('>' + contig_name + '\n')
            self.ffile.writelines("\n".join(textwrap.wrap(contigs[contcounter], 80))+'\n')
            contcounter += 1
            if gapcounter <= totalgaps - 1:
                # We have a gap to print
                outdata += ('\t').join([scaffname,
                                str(gap.start() + 1),
                                str(gap.end()),
                                str(i + 1),
                                settings.GAP_TYPE,
                                str(gap.end() - gap.start()),
                                'scaffold',
                                'yes',
                                settings.LINKAGE_EVIDENCE + '\n'])
                gapcounter += 1
                i += 1
            contig_start = gap.end() + 1
            self.contig_counter +=1
        return outdata


    def fasta2agp(self):
        line = self.ifile.readline()
        while line:
            if line and line[0] == '>':
                # New scaffold
                print "========="
                print "Original scaffold-name:"
                print line[1:-1]
                scaff_name = settings.SCAFFOLD_PRE + str(self.scaff_name_counter)
                scaff = ''
                line = self.ifile.readline().strip()
                print "New name:"
                print scaff_name
                while line and line[0] != '>':
                    scaff += line
                    line = self.ifile.readline().strip()
                self.agpfile.writelines(self.scaffold2parts(scaff_name, scaff))
                self.scaff_name_counter += 1
            else:
                print "Error in input fasta"
                raise


if __name__ == "__main__":
    Fasta2AGP(settings.FASTA_FILE, settings.OUT_PRE).fasta2agp()
