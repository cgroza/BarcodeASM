#module load mugqic/RepeatMasker bcftools
VCF=$1
suffix=$(basename ${VCF})
bcftools norm -N -m- "${VCF}" > norm_${suffix}
bcftools view -H norm_${suffix} | awk '{print(sprintf(">%d_%s %s", NR, $3, $5))}' | tr ' ' '\n' > fasta.fa
mkdir -p out
RepeatMasker -species human -s -pa 40 fasta.fa -dir out

Rscript ~/git/BarcodeAsm/scripts/repmask_vcf.R out/fasta.out norm_${suffix} annot.tsv
bgzip annot.tsv
tabix -s1 -b2 -e2 annot.tsv.gz

echo -e '##INFO=<ID=matching_class,Number=1,Type=String,Description="Repeat family">' > hdr.txt
echo -e '##INFO=<ID=repeat_id,Number=1,Type=String,Description="Repeat name">' >> hdr.txt
echo -e '##INFO=<ID=TE,Number=1,Type=String,Description="If repeat is TE">' >> hdr.txt

bcftools annotate -a annot.tsv.gz -h hdr.txt -c CHROM,POS,INFO/matching_class,INFO/repeat_id,INFO/TE norm_${suffix} | \
    bcftools view -Oz -o annot_${suffix}
tabix annot_${suffix}

