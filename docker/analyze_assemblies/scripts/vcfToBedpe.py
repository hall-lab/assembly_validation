#!/usr/local/bin/python2.6

import sys
from optparse import OptionParser
import pysam

def classifyType(altLength, refLength):
    if altLength > refLength:
        return "INS"
    elif refLength > altLength:
        return "DEL"
    else:
        return "BALANCED"

def classifyGenotype(gt):
    if gt[0]==gt[1]:
        return "HOMALT"
    else:
        return "HET"

def parseLine(line, out, minimumSize, maximumSize):
    varLen = len(line.alts[0])-line.rlen
    if (varLen >= int(minimumSize) and varLen <= int(maximumSize)) or (varLen >= int(minimumSize) and maximumSize == -1):
        out.write("\t".join((line.chrom, str(line.pos-1), str(line.pos+1), line.chrom, str(line.pos-1), str(line.pos+1), classifyGenotype(line.samples[0].get("GT")), classifyType(len(line.alts[0]), line.rlen)))+"\n")
    elif (varLen <= -1*int(minimumSize) and varLen >= -1*int(maximumSize)) or (varLen <= -1*int(minimumSize) and maximumSize==-1):
        varLen = abs(varLen)
        out.write("\t".join((line.chrom, str(line.pos-1), str(line.pos+1), line.chrom, str(line.pos-1+varLen), str(line.pos+varLen+1), classifyGenotype(line.samples[0].get("GT")), classifyType(len(line.alts[0]), line.rlen)))+"\n")

def processVar(opts):
    out = open(opts.outFile, "w")
    out.write("\t".join(("#CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "GENOTYPE", "SVTYPE"))+"\n")
    with pysam.VariantFile(opts.inFile, 'r') as inVar:
        for line in inVar:
            parseLine(line, out, opts.minimumSize, opts.maximumSize)
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
    parser.add_option("-M", dest="maximumSize", help="Maximum size of indel to output, -1 for no max", metavar="INT", default=-1)
    (opts, args) = parser.parse_args()

    try:
        processVar(opts)
    except IOError as err:
        sys.stderr.write("IOError " + str(err) + "\n")
        return

if __name__ == "__main__":
    sys.exit(main())
    (END)
