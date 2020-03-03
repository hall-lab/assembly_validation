import sys


totals = {}
for line in sys.stdin:
    line = line.strip()
    fields = line.split("\t")
    segment_length = int(fields[2])-int(fields[1])
    if fields[0] in totals:
        if fields[3] in totals[fields[0]]:
            totals[fields[0]][fields[3]]+=segment_length
        else:
            totals[fields[0]][fields[3]] = segment_length
    else:
        totals[fields[0]] = {}
        totals[fields[0]][fields[3]] = segment_length

for chrom, coverage_dict in totals.items():
    for coverage, count in coverage_dict.items():
        print("\t".join([chrom, coverage, str(count)]))
