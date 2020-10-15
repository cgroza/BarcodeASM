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
                        vcf_out_path, selected_contigs_path, min_insert_size = 50):
    records = []

    ref_fasta = pysam.FastaFile(ref_fasta_path)
    contig_fasta = pysam.FastaFile(contigs_path)
    selected_contig_fasta = open(selected_contigs_path, "w")
    alns = pandas.read_csv(alignments_path, sep = " ")

    reader = vcfpy.Reader.from_path(vcf_template_path)
    reader.header.samples = vcfpy.SamplesInfos([sample_name])

    writer = vcfpy.Writer.from_path(vcf_out_path, reader.header)

    for r in alns.iterrows():
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

        query_seq = contig_fasta.fetch(query_name)
        if strand == "-":
            query_seq = str(Bio.Seq.Seq(query_seq).reverse_complement())

        # print(strand)

        ref_seq =  ref_fasta.fetch(ref_chrom, ref_start, ref_end)

        query_pos = query_start
        target_pos = target_start

        # we are looking to extract insertions larger than 50bp
        for op in ops:
            # print(op)
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
                    # print(op[0], "insertion from", ops, "on strand", strand)
                    # print(query_start, query_end)
                    # print("query", query_seq)
                    # print("ref", ref_seq)
                    # print("var", query_seq[query_pos:query_pos + op[0]])
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

                    rec = vcfpy.Record(CHROM = ref_chrom, POS = ref_start + target_pos + 1, ID = [query_name],
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
                    records.append(rec)
                    writer.write_record(rec)
                    contig_hash = sha1("_{chrom}_{pos}_{alt}".format(
                        chrom = ref_chrom, pos = ref_start, alt = alt_allele).encode()).hexdigest()
                    selected_contig_fasta.writelines([">" + query_name + contig_hash + "\n", query_seq + "\n"])
                query_pos += op[0]
    selected_contig_fasta.close()
    return records

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
parser.add_argument("--min_insert", metavar="min_insert", type = int, nargs = 1,
                    help = "Minimum insertion size to extract", default = 50)

args = parser.parse_args()
print("Sample", args.sample[0])
print("Alignments", args.alns[0])
print("Input Contigs", args.contigs[0])
print("Reference", args.ref[0])
print("VCF template", args.vcf_template[0])
print("Output path",  args.vcf_out[0])
print("Output Selected Contigs", args.out_contigs[0])
print("Minimum insertion size ", args.min_insert[0])

records = extract_vcf_records(args.sample[0], args.alns[0], args.contigs[0], args.ref[0], args.vcf_template[0],
                              args.vcf_out[0], args.out_contigs[0], args.min_insert[0])

print("Extracted", len(records), "insertions")
