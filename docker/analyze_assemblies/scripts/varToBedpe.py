#!/usr/local/bin/python2.6

import sys
from optparse import OptionParser
import pysam

def parseLine(line, out, minimumSize):
    varList = line.strip().split('\t')
    if (varList[6] == "-"): #insertion
        ty = "INS"
        if (len(varList[7]) < minimumSize):
            return
        else:
            length=str(len(varList[7]))
    elif (varList[7] == "-"): #deletion
        ty = "DEL"
        if (len(varList[6]) < minimumSize):
            return
        else:
            length="-"+str(len(varList[6]))
    else:
        return
    start = int(varList[2])
    end = int(varList[3])
    filt="PASS"
    if varList[4] != "1":
        filt="contigCount_"+varList[4]
    out.write("\t".join((varList[1], str(start-1), str(start+1), varList[1], str(end-1), str(end+1), ':'.join((varList[1], str(start-1), str(start+1), varList[1], str(end-1), str(end+1), "SVLEN="+length)), varList[5], '+', varList[11], ty, filt, varList[8], '.', '.', '.', '.', '.', "SVTYPE="+ty+";SVLEN="+length, '.', 'GT', "1|."))+"\n")

def processVar(opts):
    out = open(opts.outFile, "w")
    out.write("##fileformat=BEDPE\n")
    out.write("\t".join(("#CHROM_A", "START_A", "END_A", "CHROM_B", "START_B", "END_B", "ID", "QUAL", "STRAND_A", "STRAND_B", "TYPE", "FILTER", "NAME_A", "REF_A", "ALT_A", "NAME_B", "REF_B", "ALT_B", "INFO_A", "INFO_B", "FORMAT", opts.sample))+"\n")
    with open(opts.inFile) as inVar:
        for line in inVar:
            if (line.startswith("V")):
                parseLine(line, out, int(opts.minimumSize))
    out.close()

def main():

    usage = """%prog -r <reference> -i <var.txt file>

varToVcf

Author: Allison Regier	
Description: Converts var format produced by paftools call into VCF
    """
    parser = OptionParser(usage)
    parser.add_option("-i", dest="inFile", help="Variant file as produced by paftools", metavar="FILE")
    parser.add_option("-o", dest="outFile", help="File name to write VCF", metavar="FILE")
    parser.add_option("-s", dest="sample", help="Sample name to use in VCF file", metavar="STR")
    parser.add_option("-m", dest="minimumSize", help="Minimum size of indel to output", metavar="INT", default=1)
    (opts, args) = parser.parse_args()

    if opts.inFile is None:
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
