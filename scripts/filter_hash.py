import pysam
import argparse
import functools


parser = argparse.ArgumentParser("Reduce out contigs with redundant non-reference sequences.")
parser.add_argument("--fasta", metavar="fasta", type = str, nargs = 1,
                    help = "Contig fasta file to be reduced")
parser.add_argument("--out", metavar="out", type = str, nargs = 1,
                    help = "Output path for reduced contig fasta file")
parser.add_argument("--blacklist", metavar="blacklist", default = "", type = str, nargs = 1,
                    help = "List of regions to exclude")

args = parser.parse_args()

contig_fasta = pysam.FastaFile(args.fasta[0])
reduced_fasta = open(args.out[0], "w")

blacklist = set()
if len(args.blacklist) > 0:
    with open(args.blacklist[0]) as blacklist_f:
        blacklist.update(blacklist_f.read().splitlines())

hashes = set()

def cmp_contig(x, y):
    (x_chrom, x_start, _, _, _, _, _, _, x_sv_len) = x.split("_")
    (y_chrom, y_start, _, _, _, _, _, _, y_sv_len) = y.split("_")

    x_start, x_sv_len, y_start, y_sv_len = int(x_start), int(x_sv_len), int(y_start), int(y_sv_len)

    if x_chrom == y_chrom:
        if x_start == y_start:
            if x_sv_len < y_sv_len: return 1
            else: return -1
        if x_start < y_start: return -1
        else: return 1
    elif x_chrom < y_chrom: return -1
    else: return 1

for seq_name in sorted(contig_fasta.references, key=functools.cmp_to_key(cmp_contig)):
    print(seq_name)
    (chrom, start, end, ps, hp, num, sample, var_hash, sv_len) = seq_name.split("_")

    # skip if in blacklist
    if "_".join([chrom, start, end]) in blacklist:
        print("Excluded by blacklist")
        continue

    if var_hash not in hashes:
        # output unique variant
        contig_seq = contig_fasta.fetch(seq_name)
        reduced_fasta.writelines([">" + seq_name + "\n", contig_seq + "\n"])
        # never output this hash again
        hashes.add(var_hash)

reduced_fasta.close()
