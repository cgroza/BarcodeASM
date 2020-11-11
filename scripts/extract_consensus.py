import mappy as mappy
import Bio.Seq
import pysam
import vcfpy
import re
import cigar
import argparse

def make_locus(chrom, st, end):
    return "_".join([chrom, st, end])

def update_genotype(loci, locus, sample, hp, ps):
    # only phased variants at this point
    assert hp != "0"
    # if sample not at this locus yet
    if not sample in loci[locus]:
        loci[locus][sample] = {}

    # no genotype at for this sample at this locus yet
    if len(loci[locus][sample]) == 0:
        # unphased, hp1, hp2, ps1, ps2
        loci[locus][sample] = {"0": "0", "1": "0", "2": "0", "ps": None}

    # update the genotype of locus at sample
    loci[locus][sample][hp] = "1"
    # phase sets are the same for hp1 and hp2
    loci[locus][sample]["ps"] = ps

def collect_genotypes(contig_fasta_path):
    fasta = pysam.FastaFile(contig_fasta_path)
    samples = set()
    loci = {}
    # collect samples from fasta files
    # make first pass to collect all the samples and loci in the contigs
    for contig in fasta.references:
        (chrom, st, end, ps, hp, n, sample, var_hash, var_len) = contig.split("_")
        samples.add(sample)
        loci[make_locus(chrom, st, end)] = {}

    # second pass to fill in sample data
    for contig in fasta.references:
        (chrom, st, end, ps, hp, n, sample, var_hash, var_len) = contig.split("_")
        update_genotype(loci, make_locus(chrom, st, end), sample, hp[2:], ps[2:])

    return (samples, loci)


def extract_consensus_insertions(contig_path, cons_path, ref_fasta_path, vcf_out_path, vcf_template_path, min_insertion_size):
    n_records = 0
    # open input sequences
    cons_fasta = pysam.FastaFile(cons_path)
    ref_fasta = pysam.FastaFile(ref_fasta_path)

    (samples, loci) = collect_genotypes(contig_path)
    print("Found", len(samples), "samples for", len(loci), "phased loci")

    reader = vcfpy.Reader.from_path(vcf_template_path)
    reader.header.samples = vcfpy.SamplesInfos(list(samples))
    writer = vcfpy.Writer.from_path(vcf_out_path, reader.header)

    for contig in cons_fasta.references:
        # parse coordinates
        (chrom, start, end) = contig.split("_")
        (start, end) = int(start), int(end)

        cons_seq  = cons_fasta.fetch(contig)
        ref_seq = ref_fasta.fetch(chrom, start, end)

        aligner = mappy.Aligner(seq = ref_seq, preset = None , k = 10, w = 10, n_threads = 1)
        alignments = list(aligner.map(cons_seq, seq2 = None, cs = True, MD = False))

        if len(alignments) == 0:
            print("No hits in", contig)
            continue

        aln = max(alignments, key = lambda x: x.blen)

        cig = cigar.Cigar(aln.cigar_str)
        ops = list(cig.items())


        cons_pos = aln.q_st
        target_pos = aln.r_st

        strand = "+"
        if aln.strand == -1:
                cons_seq = str(Bio.Seq.Seq(cons_seq).reverse_complement())
                strand = "-"
        # print(contig)
        for op in ops:
            # skip matches
            if op[1] == 'M':
                cons_pos += op[0]
                target_pos += op[0]

            # skip deletions in the query sequence
            elif op[1] == 'D':
                target_pos += op[0]

            # insertions in the query sequence
            elif op[1] == 'I':
                # only interested in large insertions
                if op[0] > min_insertion_size:
                    # Generate pysam.VariantRecord

                    # need to check conversion from 0-based coordinates to 1-based
                    ref_allele = ref_seq[target_pos-1]
                    alt_allele = cons_seq[cons_pos:cons_pos + op[0]]

                    break_point = start + target_pos
                    # output VCF record corresponding to the insertion
                    # print(break_point, (start + end) / 2 )

                    # print(len(loci[contig]), "samples at", contig)

                    # build calls data structure
                    calls = []
                    for sample in samples:
                        sample_gt = "0/0"
                        ps = 0
                        if sample in loci[contig]:
                            sample_gt = loci[contig][sample]["1"] + "|" + loci[contig][sample]["2"]
                            ps = loci[contig][sample]["ps"]
                        sample_call = vcfpy.Call(sample = sample,
                                                 data = vcfpy.OrderedDict(GT = sample_gt, PS = ps))
                        # print(sample_call)
                        calls.append(sample_call)

                    rec = vcfpy.Record(CHROM = chrom, POS = break_point, ID = [contig + "_" + str(cons_pos)],
                                       REF = ref_allele, ALT = [vcfpy.Substitution("INS", ref_allele + alt_allele)],
                                       QUAL = 999, FILTER = ["PASS"],
                                       INFO = vcfpy.OrderedDict(SVLEN = op[0],
                                                                CIGAR = str(cig),
                                                                STRAND = strand,
                                                                CONTIG_START = str(aln.q_st)),
                                       FORMAT = ["GT", "PS"],
                                    calls = calls)

                    # output contig that contains this insertion
                    writer.write_record(rec)

                    # output same contig, but with large flanking sequences
                    # note, the interval is [start, end[
                    n_records += 1

                cons_pos += op[0]
    return n_records

parser = argparse.ArgumentParser("Extract a VCF file from bxlra contig local alignments.")
parser.add_argument("--ref", metavar="ref", type = str, nargs = 1,
                    help = "Reference genome used with bxlra")
parser.add_argument("--cons", metavar="cons", type = str, nargs = 1,
                    help = "Consensus sequences generated with msa.py")
parser.add_argument("--contigs", metavar="contigs", type = str, nargs = 1,
                    help = "Original contigs from which the consensus was generated")
parser.add_argument("--vcf_template", metavar="vcf_template", type = str, nargs = 1,
                    help = "VCF template for you reference genome")
parser.add_argument("--vcf_out", metavar="vcf_out", type = str, nargs = 1,
                    help = "VCF path for the output", default = "out.vcf")
parser.add_argument("--min_insert", metavar="min_insert", type = int, nargs = 1,
                    help = "Minimum insertion size to extract", default = 50)

args = parser.parse_args()

print("Original contigs", args.contigs[0])
print("Consensus contigs", args.cons[0])
print("Reference", args.ref[0])
print("VCF template", args.vcf_template[0])
print("VCF Output path",  args.vcf_out[0])
print("Minimum insertion size ", args.min_insert[0])

records = extract_consensus_insertions(args.contigs[0], args.cons[0], args.ref[0],
                                       args.vcf_out[0], args.vcf_template[0],
                                       args.min_insert[0])

print("Extracted", records, "insertions")
