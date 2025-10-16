#!/bin/bash --login

# 1.1 - Remove remnant adapter sequences from raw HiFi reads (hifi.fastq.gz) using HiFiAdapterFilt

export PATH=$PATH:/path/to/HiFiAdapterFilt
export PATH=$PATH:/path/to/HiFiAdapterFilt/DB

bash hifiadapterfilt.sh -p hifi -t 16

# 1.2 - Remove remnant adapter sequences from raw Nanopore reads (ont_runX.fastq) using Porechop_ABI
# Nanopore reads come from 3 seperate sequencing runs. Clean each independently, then concatenate

porechop_abi -abi -i ./ont_run1.fastq -o ./clean_run1.fastq
porechop_abi -abi -i ./ont_run2.fastq -o ./clean_run2.fastq
porechop_abi -abi -i ./ont_run3.fastq -o ./clean_run3.fastq

cat ./clean_run*.fastq > merged_clean_ont.fastq

# 2.1 - Assemble primary and alternate haplotypes using hifiasm

hifiasm --primary -t 48 -o ./ctm --ul ./merged_clean_ont.fastq ./hifi.filt.fastq.gz

# 2.2 - Convert .gfa to .fa, then sort and rename contigs to "contX" by size

awk '/^S/{print ">"$2;print $3}' ./ctm.p_ctg.gfa > ./ctm_p.fa
awk '/^S/{print ">"$2;print $3}' ./ctm.a_ctg.gfa > ./ctm_a.fa

seqkit sort --by-length --reverse ./ctm_p.fa | seqkit replace --pattern '.+' --replacement 'cont{nr}' > ./sorted_ctm_p.fa
seqkit sort --by-length --reverse ./ctm_a.fa | seqkit replace --pattern '.+' --replacement 'cont{nr}' > ./sorted_ctm_a.fa

# 3 - The following analyses and assessments were conducted on the Galaxy web platform (usegalaxy.org)
# Default settings used unless otherwise specified
# 3.1 - QUAST-LG (v5.2.0)
# 3.2 - meryl (v1.3; k=20) followed by GenomeScope2.0 (v2.0.1) and merqury (v1.3)

# 3.3 - BUSCO

busco -i ./sorted_ctm_p.fa -m genome -l sauropsida_odb10 -c 32 --miniprot
busco -i ./sorted_ctm_a.fa -m genome -l sauropsida_odb10 -c 32 --miniprot

# 3.4 - Screen for contamination using FCS-GX

export FCS_DEFAULT_IMAGE=fcs-gx.sif
GXDB_LOC=/path/to/fcs_db

python3 ./fcs.py screen genome --fasta ./sorted_ctm_p.fa --out-dir ./gx_out/ --gx-db "$GXDB_LOC" --tax-id 85608
python3 ./fcs.py screen genome --fasta ./sorted_ctm_a.fa --out-dir ./gx_out/ --gx-db "$GXDB_LOC" --tax-id 85608

# 4.1 - Visualize synteny between current assembly and assembly released by CGEn (Accession GCA_040208395.1) using JupiterPlot

jupiter name=clgu ref=./sorted_ctm_p.fa fa=can_clgu_p.fna t=8 ng=0 m=1000000

# 4.2 - Estimate demographic history of current assembly and assembly released by CGEn (Accession GCA_040208395.1) using PSMC

prefetch SRR27950426 --max-size 35g
fasterq-dump SRR27950426

bwa index ./sorted_ctm_p.fa # Current assembly
bwa index ./can_clgu_p.fna # CGEn/Canadian assembly

bwa mem -t 16 ./sorted_ctm_p.fa ./hifi.filt.fastq.gz | samtools view -Sb | samtools sort -o aligned.bam
bwa mem -t 16 ./can_clgu_p.fna ./SRR27950426.fastq | samtools view -Sb | samtools sort -o can_aligned.bam

bwa index aligned.bam
bwa index can_aligned.bam

bcftools mpileup -C50 -Ou -f ./sorted_ctm_p.fa aligned.bam | bcftools call -c | vcfutils.pl vcf2fq -d 14 -D 80 | gzip > variants.fq.gz
bcftools mpileup -C50 -Ou -f ./can_clgu_p.fna can_aligned.bam | bcftools call -c | vcfutils.pl vcf2fq -d 14 -D 80 | gzip > can_variants.fq.gz

fq2psmcfa -q20 variants.fq.gz > diploid.psmcfa
splitfa diploid.psmcfa > split.psmcfa
psmc -N20 -t20 -r4 -p "50*2+6" -o diploid.psmc diploid.psmcfa
seq 50 | xargs -P 24 -i psmc -N20 -t20 -r4 -b -p "50*2+6" -o round-{}.psmc split.psmcfa | sh
cat diploid.psmc round-*.psmc > combined.psmc
psmc_plot.pl -R -u 1.2e-08 -g 20 combined combined.psmc
psmc_plot.pl -R -u 1.2e-08 -g 20 diploid diploid.psmc
rm round-*.psmc

fq2psmcfa -q20 can_variants.fq.gz > can_diploid.psmcfa
splitfa can_diploid.psmcfa > can_split.psmcfa
psmc -N20 -t20 -r4 -p "50*2+6" -o can_diploid.psmc can_diploid.psmcfa
seq 50 | xargs -P 24 -i psmc -N20 -t20 -r4 -b -p "50*2+6" -o can_round-{}.psmc can_split.psmcfa | sh
cat can_diploid.psmc can_round-*.psmc > can_combined.psmc
psmc_plot.pl -R -u 1.2e-08 -g 20 can_combined can_combined.psmc
psmc_plot.pl -R -u 1.2e-08 -g 20 can_diploid can_diploid.psmc
rm can_round-*.psmc

cat combined*.txt > combined.txt
cat can_combined*.txt > can_combined.txt

rm combined.*.txt
rm can_combined.*.txt
