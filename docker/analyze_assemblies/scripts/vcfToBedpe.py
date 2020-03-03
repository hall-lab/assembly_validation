#!/usr/local/bin/python2.6

import sys
from optparse import OptionParser
import pysam

def classifyGenotype(gt):
    if gt[0]==gt[1]:
        return "HOMALT"
    else:
        return "HET"

def parseLine(line, out, minimumSize):
    varLen = len(line.alts[0])-line.rlen
    #out.write(line.alts[0]+"\t"+str(len(line.alts[0]))+"\t"+str(line.rlen)+"\t"+str(varLen)+"\t"+minimumSize+"\n")
    if varLen >= int(minimumSize):
        out.write("\t".join((line.chrom, str(line.pos-1), str(line.pos+1), line.chrom, str(line.pos-1), str(line.pos+1), classifyGenotype(line.samples[0].get("GT"))))+"\n")
    elif varLen <= -1*int(minimumSize):
        varLen = abs(varLen)
        out.write("\t".join((line.chrom, str(line.pos-1), str(line.pos+1), line.chrom, str(line.pos-1+varLen), str(line.pos+varLen+1), classifyGenotype(line.samples[0].get("GT"))))+"\n")

def processVar(opts):
    out = open(opts.outFile, "w")
    out.write("\t".join(("#CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "GENOTYPE"))+"\n")
    with pysam.VariantFile(opts.inFile, 'r') as inVar:
        for line in inVar:
            parseLine(line, out, opts.minimumSize)
    out.close()

def main():

    usage = """%prog -m <minimum size> -i <var.txt file>

varToVcf

Author: Allison Regier	
Description: Converts indels from standard vcf to bedpe
    """
    parser = OptionParser(usage)
    parser.add_option("-i", dest="inFile", help="Variant file as produced by paftools", metavar="FILE")
    parser.add_option("-o", dest="outFile", help="File name to write VCF", metavar="FILE")
    parser.add_option("-m", dest="minimumSize", help="Minimum size of indel to output", metavar="INT", default=1)
    (opts, args) = parser.parse_args()

    try:
        processVar(opts)
    except IOError as err:
        sys.stderr.write("IOError " + str(err) + "\n")
        return

if __name__ == "__main__":
    sys.exit(main())
    (END)
