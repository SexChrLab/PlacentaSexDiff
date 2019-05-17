# What: this script "phases" by assigning Haplotype 1 to contain the genotype where the allele count for that genotype
# over total count is greater than 0.5. We will use Placenta 1 for this.
import sys
def phasing(placenta1_fn):
    pos = []
    hap_1 = ''
    hap_2 = ''
    with open(placenta1_fn, 'r') as f:
        for line in f:
            if not line.startswith("contig"):
                if line.startswith("chrX"):
                    items = line.rstrip("\n").split("\t")
                    if float(items[6])/float(items[7]) > 0.5:
                        pos.append(items[1])
                        hap_1 += items[4]
                        hap_2 += items[3]
                    else:
                        pos.append(items[1])
                        hap_1 += items[3]
                        hap_2 += items[4]

    return pos, hap_1, hap_2

def find_ratio_for_hap(placenta_fn, pos, hap):
    out = {}
    with open(placenta_fn, 'r') as f:
        for line in f:
            if not line.startswith("contig"):
                if line.startswith("chrX"):
                    items = line.rstrip("\n").split("\t")
                    if items[1] in pos:
                        idx = pos.index(items[1])
                        if items[4] == hap[idx]:
                            # print(items[1], idx, hap[idx])
                            out[items[1]] = float(items[6])/float(items[7])
                        else:
                            out[items[1]] = float(items[5]) / float(items[7])
    return out

def main():
    pos, hap_1, hap_2 = phasing(sys.argv[1])
    placenta_1_out = find_ratio_for_hap(sys.argv[1], pos, sys.argv[5]) #sys.argv 5 is hap_1 or hap_2
    placenta_2_out = find_ratio_for_hap(sys.argv[2], pos, sys.argv[5])
    out_1 = open(sys.argv[3], 'w')
    for k, v in placenta_1_out.items():
        out = [str(k), str(v)]
        print ("\t".join(out), file=out_1)

    out_2 = open(sys.argv[4], 'w')
    for k, v in placenta_2_out.items():
        out = [str(k), str(v)]
        print ("\t".join(out), file=out_2)

main()
