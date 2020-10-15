import pysam
import argparse


parser = argparse.ArgumentParser("Reduce out contigs with redundant non-reference sequences.")
parser.add_argument("--fasta", metavar="fasta", type = str, nargs = 1,
                    help = "Contig fasta file to be reduced")
parser.add_argument("--out", metavar="out", type = str, nargs = 1,
                    help = "Output path for reduced contig fasta file")

args = parser.parse_args()

contig_fasta = pysam.FastaFile(args.fasta[0])
reduced_fasta = open(args.out[0], "w")

hashes = set()

for seq_name in contig_fasta.references:
    print(seq_name)
    (chrom, start, end, ps, hp, num, sample, var_hash) = seq_name.split("_")
    if var_hash not in hashes:
        # output unique variant
        contig_seq = contig_fasta.fetch(seq_name)
        reduced_fasta.writelines([">" + seq_name + "\n", contig_seq + "\n"])

        # never output this hash again
        hashes.add(var_hash)

reduced_fasta.close()
