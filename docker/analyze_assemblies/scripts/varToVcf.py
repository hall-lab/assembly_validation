#!/usr/local/bin/python2.6

import sys
from optparse import OptionParser
import pysam

counter = 0
def parseLine(line, ref, out, out2, minimumSize, idPrefix):
    varList = line.strip().split('\t')
    if (varList[6] == "-"): #insertion
        if (len(varList[7]) < minimumSize):
            return
        start = int(varList[2])
        padding = ref.fetch(varList[1], start-1, start)
        refBase = padding
        altBase = padding + varList[7]
        start2 = int(varList[9])
        padding2 = ref.fetch(varList[8], start2-1, start2)
        altBase2 = padding2
        refBase2 = padding2 + varList[6]
    elif (varList[7] == "-"): #deletion
        if (len(varList[6]) < minimumSize):
            return
        start = int(varList[2])
        padding = ref.fetch(varList[1], start-1, start)
        altBase = padding
        refBase = padding + varList[6]
        start2 = int(varList[9])
        padding2 = query.fetch(varList[8], start2-1, start2)
        refBase2 = padding2
        altBase2 = padding2 + varList[7]
    else:
        if (minimumSize > 1):
            return
        refBase = varList[6]
        altBase = varList[7]
        start = int(varList[3])
    id = "{}.{}".format(idPrefix, str(counter))
    counter = counter + 1
    out.write("\t".join((varList[1], str(start), id, refBase.upper(), altBase.upper(), varList[5], ".", ";".join(["COV="+varList[4], "QNAME="+varList[8], "QSTART="+varList[9]]), "GT", "1|."))+"\n")
    out2.write("\t".join((varList[8], str(start2, id, refBase2.upper(), altBase2.upper(), varList[5], ".", ";".join(["COV="+varList[4], "QNAME="+varList[1], "QSTART="+varList[2]]), "GT", "1|."))+"\n")

def processVar(opts):
    ref = pysam.Fastafile(opts.refFile)
    out = open(opts.outFile, "w")
    out.write("##fileformat=VCFv4.2\n")
    out.write("\t".join(("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", opts.sample))+"\n")
    out2 = open("{}.2.vcf".format(opts.outFile), "w")
    out2.write("##fileformat=VCFv4.2\n")
    out2.write("\t".join(("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", opts.sample))+"\n")
    with open(opts.inFile) as inVar:
        for line in inVar:
            if (line.startswith("V")):
                parseLine(line, ref, out, out2, opts.minimumSize, opts.idPrefix)
    out.close()
    out2.close()

def main():

    usage = """%prog -r <reference> -i <var.txt file>

varToVcf

Author: Allison Regier	
Description: Converts var format produced by paftools call into VCF
    """
    parser = OptionParser(usage)
    parser.add_option("-r", dest="refFile", help="Reference fasta file indexed by faidx", metavar="FILE")
    parser.add_option("-i", dest="inFile", help="Variant file as produced by paftools", metavar="FILE")
    parser.add_option("-o", dest="outFile", help="File name to write VCF", metavar="FILE")
    parser.add_option("-s", dest="sample", help="Sample name to use in VCF file", metavar="STR")
    parser.add_option("-m", dest="minimumSize", help="Minimum size of indel to output", metavar="INT", default=1)
    parser.add_option("-p", dest="idPrefix", help="Prefix to start variant ids", metavar="STR")
    (opts, args) = parser.parse_args()

    if opts.refFile is None:
        parser.print_help()
        print
    else:
        try:
            processVar(opts)
        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n")
            return

if __name__ == "__main__":
    sys.exit(main())
    (END)
