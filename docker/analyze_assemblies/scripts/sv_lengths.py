import pysam
import sys
from optparse import OptionParser

length_dict = {
    "del": {
        "1 bp": 0,
        "2 bp": 0,
        "3-10 bp": 0,
        "11-50 bp": 0,
        "51-100 bp": 0,
        "101-1000 bp": 0,
        "1-10 kb": 0,
        "10 kb+": 0
    },
    "ins": {
        "1 bp": 0,
        "2 bp": 0,
        "3-10 bp": 0,
        "11-50 bp": 0,
        "51-100 bp": 0,
        "101-1000 bp": 0,
        "1-10 kb": 0,
        "10 kb+": 0
    }
};

def parseLine(line, out):
    varLen = line.rlen - len(line.alts[0])
    if varLen == 1:
        length_dict['del']['1 bp'] += 1
    elif varLen == -1:
        length_dict['ins']['1 bp'] += 1
    elif varLen == 2:
        length_dict['del']['2 bp'] += 1
    elif varLen == -2:
        length_dict['ins']['2 bp'] += 1
    elif varLen > 2 and varLen <= 10:
        length_dict['del']['3-10 bp'] += 1
    elif varLen < -2 and varLen >= -10:
        length_dict['ins']['3-10 bp'] += 1
    elif varLen > 10 and varLen <= 50:
        length_dict['del']['11-50 bp'] += 1
    elif varLen < -10 and varLen >= -50:
        length_dict['ins']['11-50 bp'] += 1
    elif varLen > 50 and varLen <= 100:
        length_dict['del']['51-100 bp'] += 1
    elif varLen < -50 and varLen >= -100:
        length_dict['ins']['51-100 bp'] += 1
    elif varLen > 100 and varLen <= 1000:
        length_dict['del']['101-1000 bp'] += 1
    elif varLen < -100 and varLen >= -1000:
        length_dict['ins']['101-1000 bp'] += 1
    elif varLen > 1000 and varLen <= 10000:
        length_dict['del']['1-10 kb'] += 1
    elif varLen < -1000 and varLen >= -10000:
        length_dict['ins']['1-10 kb'] += 1
    elif varLen > 10000:
        length_dict['del']['10 kb+'] += 1
    elif varLen < -10000:
        length_dict['ins']['10 kb+'] += 1

def processVar(opts):
    out = open(opts.outFile, "w")
    out.write("\t".join(("Type", "Metric", "Value"))+"\n")
    with pysam.VariantFile(opts.inFile, 'r') as inVar:
        for line in inVar:
            parseLine(line, out)
    for key, value in length_dict.items():
        for key2, value2 in value.items():
            out.write("\t".join([key, key2, str(value2)])+"\n")
    out.close()

def main():
    usage = "%prog -v <vcf file> -o <output file>\n"
    parser = OptionParser(usage)
    parser.add_option("-v", dest="inFile", metavar="FILE", help="VCF file")
    parser.add_option("-o", dest="outFile", metavar="FILE", help="Output file")
    (opts, args) = parser.parse_args()
    processVar(opts)

if __name__ == "__main__":
    sys.exit(main())
