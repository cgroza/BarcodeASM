import pysam
import pandas
import cigar
import Bio.Seq
import vcfpy
import argparse
from hashlib import sha1

def extract_vcf_records(sample_name,
                        # input paths
                        alignments_path, contigs_path, ref_fasta_path, vcf_template_path,
                        # output paths
                        vcf_out_path, selected_contigs_path, flanked_contigs_path,
                        flank_length, min_insert_size):

    n_records = 0
    ref_fasta = pysam.FastaFile(ref_fasta_path)
    contig_fasta = pysam.FastaFile(contigs_path)

    selected_contig_fasta = open(selected_contigs_path, "w")
    flanked_contig_fasta = open(flanked_contigs_path, "w")

    alns = pandas.read_csv(alignments_path, sep = " ")

    reader = vcfpy.Reader.from_path(vcf_template_path)
    reader.header.samples = vcfpy.SamplesInfos([sample_name])

    writer = vcfpy.Writer.from_path(vcf_out_path, reader.header)

    contig_loci = set()

    # parse each alignment and look for insertions above min_insert_size
    for r in alns.iterrows():
        # skip secondary alignments
        hit =  r[1]["Hit"]
        if hit != "0":
            continue

        query_name = r[1]["QName"]

        # local alignment window in the reference
        ref_chrom, ref_start, ref_end, phase_set, phase, n = query_name.split("_")

        phase_set = phase_set[2:]
        phase = phase[2:]

        # convert to ints
        ref_start, ref_end = (int(ref_start), int(ref_end))

        # alignment start and end for reference sequence
        target_start = r[1]["TStart"]
        target_end = r[1]["TEnd"]

        # alignment start and end for query sequence
        query_start = r[1]["QStart"]
        query_end = r[1]["QEnd"]

        # strand-ness of the query sequence
        strand = r[1]["Strand"]

        # parse cigar for variant extraction
        cig = cigar.Cigar(r[1]["CIGAR"])
        ops = list(cig.items())

        # convert sequences to the positive strand
        query_seq = contig_fasta.fetch(query_name)
        if strand == "-":
            query_seq = str(Bio.Seq.Seq(query_seq).reverse_complement())

        ref_seq =  ref_fasta.fetch(ref_chrom, ref_start, ref_end)

        # initialize iterators for the cigar string
        query_pos = query_start
        target_pos = target_start

        # we are looking to extract insertions larger than 50bp
        for op in ops:
            # skip matches
            if op[1] == 'M':
                query_pos += op[0]
                target_pos += op[0]

            # skip deletions in the query sequence
            elif op[1] == 'D':
                target_pos += op[0]

            # insertions in the query sequence
            elif op[1] == 'I':
                # only interested in large insertions
                if op[0] > min_insert_size:
                    # Generate pysam.VariantRecord

                    # need to check conversion from 0-based coordinates to 1-based
                    ref_allele = ref_seq[target_pos]
                    alt_allele = ref_allele + query_seq[query_pos:query_pos + op[0]]

                    gt = ""
                    if phase == "1":
                        gt = "1|0"
                    elif phase == "2":
                        gt = "0|1"
                    else:
                        gt = "0/1"

                    break_point = ref_start + target_pos
                    # output VCF record corresponding to the insertion
                    rec = vcfpy.Record(CHROM = ref_chrom, POS = break_point + 1, ID = [query_name],
                                    REF = ref_allele, ALT = [vcfpy.Substitution("INS", alt_allele)],
                                    QUAL = 999, FILTER = ["PASS"], INFO = {}, FORMAT =
                                       ["GT", "SVLEN", "PS", "HP", "CIGAR", "STRAND", "CONTIG_START"],
                                    calls = [
                                        vcfpy.Call(sample = sample_name,
                                                   data = vcfpy.OrderedDict(GT = gt, SVLEN = op[0],
                                                                            PS = phase_set, HP = phase,
                                                                            CIGAR = str(cig), STRAND = strand,
                                                                            CONTIG_START = str(query_start)))]
                                    )

                    n_records += 1
                    # output contig that contains this insertion
                    writer.write_record(rec)

                    contig_locus = ">" + query_name + "_" + sample_name
                    contig_hash = sha1("_{chrom}_{pos}_{alt}".format(
                        chrom = ref_chrom, pos = ref_start, alt = alt_allele[1:]).encode()).hexdigest()

                    contig_name = contig_locus + "_" + contig_hash + "_" + str(op[0])

                    if contig_locus not in contig_loci:
                        selected_contig_fasta.writelines([contig_name + "\n",
                                                      query_seq + "\n"])
                        contig_loci.add(contig_locus)


                    # output same insertion, but with flanking sequences
                    # note, the interval is [start, end[
                    if flank_length > 0:
                        left_flank = ref_fasta.fetch(ref_chrom, break_point - flank_length, break_point)
                        right_flank = ref_fasta.fetch(ref_chrom, break_point, break_point + flank_length)
                    else:
                        left_flank = ""
                        right_flank = ""
                    flanked_contig_fasta.writelines([contig_name + "\n",
                                                     left_flank + alt_allele[1:] + right_flank + "\n"])

                query_pos += op[0]
    selected_contig_fasta.close()
    return n_records

parser = argparse.ArgumentParser("Extract a VCF file from bxlra contig local alignments.")
parser.add_argument("--sample", metavar="sample", type = str, nargs = 1,
                    help = "Sample name")
parser.add_argument("--ref", metavar="ref", type = str, nargs = 1,
                    help = "Reference genome used with bxlra")
parser.add_argument("--alns", metavar="alns", type = str, nargs = 1,
                    help = "alignments.tsv from bxlra")
parser.add_argument("--contigs", metavar="contigs", type = str, nargs = 1,
                    help = "contigs.fa from bxlra")
parser.add_argument("--vcf_template", metavar="vcf_template", type = str, nargs = 1,
                    help = "VCF template for you reference genome")
parser.add_argument("--vcf_out", metavar="vcf_out", type = str, nargs = 1,
                    help = "VCF path for the output", default = "out.vcf")
parser.add_argument("--out_contigs", metavar="out_contigs", type = str, nargs = 1,
                    help = "Contains contigs from which variants were selected", default = "selected_contigs.fa")
parser.add_argument("--flanked_inserts", metavar="flanked_inserts", type = str, nargs = 1,
                    help = "Contains insertions large flanks added", default = "flanked_contigs.fa")
parser.add_argument("--flank_length", metavar="flank_length", type = int, nargs = 1,
                    help = "Specify the length of flanks", default = 10000)
parser.add_argument("--min_insert", metavar="min_insert", type = int, nargs = 1,
                    help = "Minimum insertion size to extract", default = 50)

args = parser.parse_args()

print("Sample", args.sample[0])
print("Alignments", args.alns[0])
print("Input Contigs", args.contigs[0])
print("Reference", args.ref[0])
print("VCF template", args.vcf_template[0])
print("VCF Output path",  args.vcf_out[0])
print("Selected contigs fasta output path", args.out_contigs[0])
print("Flanked inserts fasta output path", args.flanked_inserts[0])
print("Flank length", args.flank_length)
print("Minimum insertion size ", args.min_insert[0])

records = extract_vcf_records(args.sample[0], args.alns[0], args.contigs[0], args.ref[0], args.vcf_template[0],  # inputs
                              args.vcf_out[0], args.out_contigs[0], args.flanked_inserts[0],  # outputs
                              args.flank_length[0], args.min_insert[0])  # parameters

print("Extracted", records, "insertions")
