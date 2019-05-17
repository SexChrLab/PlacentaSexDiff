# WHAT: This script:
# Phase each gene where gene has more than 2 variants

import sys
from collections import defaultdict

def process_ase(placenta_fn):
    ase_placenta = {}
    with open(placenta_fn, 'r') as f:
        for line in f:
            if not line.startswith('contig'):
                if line.startswith('chrX'):
                    items = line.rstrip("\n").split("\t")
                    ase_placenta[items[1]] = (items[3], items[4], int(items[6]), int(items[7]))  # refAllele, altAllele, altCount, totalCount
    return ase_placenta

def process_genes(genename_vcf_fn, ase_placenta_1, ase_placenta_2):
    chrX_gene_variants = defaultdict(list)  # a dictionary where key is the gene and values is a list of positions on chrX
    with open(genename_vcf_fn, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                items = line.rstrip('\n').split('\t')
                if (items[7].split(';')[0]).split('=')[1] == '1':
                    if items[0] == 'chrX':
                        gene = items[7].split("|")[4]
                        if items[1] in ase_placenta_1 and items[1] in ase_placenta_2:
                            chrX_gene_variants[gene].append(items[1])

    chrX_gene_more_than_2_variants = {k: v for k, v in chrX_gene_variants.items() if len(v) > 1}
    return chrX_gene_more_than_2_variants

def main():
    ase_placenta_1 = process_ase(sys.argv[1]) #OBG0044_placenta_1_hets_totalcountgreater20.csv
    ase_placenta_2 = process_ase(sys.argv[2]) #OBG0044_placenta_2_hets_totalcountgreater20.csv
    chrX_gene_more_than_2_variants = process_genes(sys.argv[3], ase_placenta_1, ase_placenta_2) #OBG0044.gatk.called.raw_vep.vcf

    placenta_1_out = open(sys.argv[4], 'w')
    placenta_2_out = open(sys.argv[5], 'w')
    for gene, variants in chrX_gene_more_than_2_variants.items():
        hap_1 = ''
        hap_2 = ''
        pos = []
        out_placenta_1 = {}
        out_placenta_2 = {}
        for variant in variants:
            placenta_1_ratio = ase_placenta_1[variant][2] / ase_placenta_1[variant][3]
            if placenta_1_ratio > 0.5:
                hap_1 += ase_placenta_1[variant][1]
                hap_2 += ase_placenta_1[variant][0]
                pos.append(variant)
            else:
                hap_1 += ase_placenta_1[variant][0]
                hap_2 += ase_placenta_1[variant][1]
                pos.append(variant)
        for i in pos:
            idx = pos.index(i)
            if ase_placenta_1[i][1] == hap_2[idx]: #hap_1
                out_placenta_1[i] = ase_placenta_1[i][2] / ase_placenta_1[i][3]
            else:
                out_placenta_1[i] = 1 - (ase_placenta_1[i][2] / ase_placenta_1[i][3])

            if ase_placenta_2[i][1] == hap_2[idx]: #hap_1
                out_placenta_2[i] = ase_placenta_2[i][2] / ase_placenta_2[i][3]
            else:
                out_placenta_2[i] = 1 - (ase_placenta_2[i][2] / ase_placenta_2[i][3])

        for k, v in out_placenta_1.items():
            out = [gene, str(k), str(v)]
            print ("\t".join(out), file=placenta_1_out)

        for k, v in out_placenta_2.items():
            out = [gene, str(k), str(v)]
            print("\t".join(out), file=placenta_2_out)

main()

# ase_placenta_1 = {}
# with open("c://Users/tuyen/Documents/postdoc_asu/projects/PlacentaSexDiff/batch_1_exome_processing/analysis/analyze_ase/results/filter_for_hets/OBG0044/OBG0044_placenta_1_hets_totalcountgreater20.csv", "r") as f:
#     for line in f:
#         if not line.startswith("contig"):
#             if line.startswith('chrX'):
#                 items = line.rstrip("\n").split("\t")
#                 ase_placenta_1[items[1]] = (items[3], items[4], int(items[6]), int(items[7])) #refAllele, altAllele, altCount, totalCount
#
# ase_placenta_2 = {}
# with open("c://Users/tuyen/Documents/postdoc_asu/projects/PlacentaSexDiff/batch_1_exome_processing/analysis/analyze_ase/results/filter_for_hets/OBG0044/OBG0044_placenta_2_hets_totalcountgreater20.csv", "r") as f:
#     for line in f:
#         if not line.startswith("contig"):
#             if line.startswith('chrX'):
#                 items = line.rstrip("\n").split("\t")
#                 ase_placenta_2[items[1]] = (items[3], items[4], float(items[6]), float(items[7])) #refAllele, altAllele, altCount, totalCount

# chrX_gene_variants = defaultdict(list) #a dictionary where key is the gene and values is a list of positions on chrX
# with open("c://Users/tuyen/Documents/postdoc_asu/projects/PlacentaSexDiff/batch_1_exome_processing/OBG0044/OBG0044.gatk.called.raw_vep.vcf", "r") as f:
#     for line in f:
#         if not line.startswith("#"):
#             items = line.rstrip("\n").split("\t")
#             if (items[7].split(";")[0]).split("=")[1] == "1":
#                 if items[0] == 'chrX':
#                     gene = items[7].split("|")[4]
#                     if items[1] in ase_placenta_1 and items[1] in ase_placenta_2:
#                         chrX_gene_variants[gene].append(items[1])
#
# chrX_gene_more_than_2_variants = {k: v for k, v in chrX_gene_variants.items() if len(v) > 1}

# for gene, variants in chrX_gene_more_than_2_variants.items():
#     haplotype_1 = ''
#     haplotype_2 = ''
#     pos = []
#     out_placenta_1 = {}
#     out_placenta_2 = {}
#     for variant in variants:
#         placenta_1_ratio = ase_placenta_1[variant][2]/ase_placenta_1[variant][3]
#         if placenta_1_ratio > 0.5:
#             haplotype_1 += ase_placenta_1[variant][1]
#             haplotype_2 += ase_placenta_1[variant][0]
#             pos.append(variant)
#         else:
#             haplotype_1 += ase_placenta_1[variant][0]
#             haplotype_2 += ase_placenta_1[variant][1]
#             pos.append(variant)
#     for i in pos:
#         idx = pos.index(i)
#         if ase_placenta_1[i][1] == haplotype_1[idx]:
#             out_placenta_1[i] = ase_placenta_1[i][2]/ase_placenta_1[i][3]
#         else:
#             out_placenta_1[i] = 1 - (ase_placenta_1[i][2] / ase_placenta_1[i][3])
#
#         if ase_placenta_2[i][1] == haplotype_1[idx]:
#             out_placenta_2[i] = ase_placenta_2[i][2]/ase_placenta_2[i][3]
#         else:
#             out_placenta_2[i] = 1 - (ase_placenta_2[i][2] / ase_placenta_2[i][3])
#
#     for k, v in out_placenta_1.items():
#         print (gene, k, v)