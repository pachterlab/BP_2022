Downloading data
--
Metadata was fetched with [`ffq`](https://github.com/pachterlab/ffq/)
```bash
ffq -o ffq_GSM5917802.json GSM5917802
```

SRA files for sequencing reads were downloaded and then converted to FASTQ files
```bash
wget --quiet https://sra-download.ncbi.nlm.nih.gov/traces/sra53/SRR/017720/SRR18145553 # illumina
wget --quiet https://sra-download.ncbi.nlm.nih.gov/traces/sra34/SRR/017720/SRR18145555 # ultima

# convert SRA file to FASTQ files
fasterq-dump ./SRR18145553 --split-files --include-technical
fasterq-dump ./SRR18145555 --split-files --include-technical
```

The human pseudoalignment index was downloaded with [`kb-python`](https://github.com/pachterlab/kb_python)
```bash
kb ref -d human -i index -g t2g.txt
```

The human HISAT2 index was downloaded and decompressed
```bash
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xvf grch38_genome.tar.gz
```

The TMSB4X fasta was downloaded with (`gget`)[https://github.com/pachterlab/gget], and indexed with `kallisto` and `hisat2`.
```bash
# download tmsb4x fasta
gget seq -id ENST00000380636.1 -o TMSB4X_ENST00000380636.1.fasta -st gene

# build kallisto index
kallisto index -i index TMSB4X_ENST00000380636.1.fasta
# and t2g map
echo -e "ENST00000380636\tENSG00000205542\tTMSB4X" > t2g.txt

# build hisat2 index
hisat2-build -p 16 TMSB4X_ENST00000380636.1.fasta genome
```

Illumina
--
Illumina reads have the following structure:
- (read 1) 16bp cell barcode
- (read 1) 12bp umi
- (read 2) 55bp cDNA

Reads were pseudoaligned to the human transcriptome and to the TMSB4X gene using (`kb-python`)[https://github.com/pachterlab/kb_python]
```bash
# pseudoalign to transcriptome
kb count --strand unstranded -o out-transcriptome -i human/index -x 10xv3 -g t2g.txt --h5ad -m 16G -t 16 SRR18145553_1.fastq.gz SRR18145553_2.fastq.gz

# pseudoalign to tmsb4x
kb count --strand unstranded -o out-tmsb4x -i tmsb4x/index -x 10xv3 -g t2g.txt --h5ad -m 16G -t 16 SRR18145553_1.fastq.gz SRR18145553_2.fastq.gz
```

HISAT2 alignments and alignment stats were made by 
```bash
# map to tmsb4x
hisat2 -p 32 -x tmsb4x/genome -U SRR18145553_2.fastq.gz -S output.sam

# convert sam to bam and filter unmapped reads
samtools view -bS output.sam | samtools view -b -F 4 > mapped.bam

# indel stats
samtools stats mapped.bam |  grep ^ID | cut -f 2- > indels.txt

# base counting
samtools view mapped.bam | cut -f10 -d$'\t' | tr -d '\n' | wc -c
```

Ultima
--
Ultima reads have the following structure:
- (read1) 22bp insert
- (read1) 26bp cell barcode
- (read1) 12bp umi
- (read1) 12bp polyT homopolymer
- (read1) XXbp cDNA

Reads were pseudoaligned using (`kb-python`)[https://github.com/pachterlab/kb_python], using the entire cDNA
```bash
# entire cDNA
kb count --strand unstranded -o out_full-transcriptome -i human/index -x 10xv3_Ultima -g t2g.txt --h5ad -m 16G -t 16 SRR18145555.fastq.gz

# 55bp maximum length cDNA
kb count --strand unstranded -w 10x_version3_whitelist.txt -o out_trim_55_max-transcriptome -i human/index -x 0,22,38:0,38,50:0,62,117 -g t2g.txt --h5ad -m 16 -t 16 SRR18145555.fastq.gz
```

Note:
```
            bc      umi     min 55bp cDNA
kb count -x 0,22,38:0,38,50:0,62,117
```

The Ultima reads were trimmed so that the minimum cDNA length was 31bp.
```bash
# trim cDNA reads
seqkit subseq -r 63:-1 SRR18145555.fastq.gz | gzip > cdna_63.fastq.gz

# keep cDNA reads longer than or equal to 31 bp
seqkit seq -m 31 cdna_63.fastq.gz | gzip > cdna_31_min.fastq.gz
```

HISAT2 alignments and alignment stats were made by 
```bash
# map to tmsb4x
hisat2 -p 32 -x tmsb4x/genome -U cdna_31_min.fastq.gz -S output.sam

# convert sam to bam and filter unmapped reads
samtools view -bS output.sam | samtools view -b -F 4 > mapped.bam

# indel stats
samtools stats mapped.bam |  grep ^ID | cut -f 2- > indels.txt

# base counting
samtools view mapped.bam | cut -f10 -d$'\t' | tr -d '\n' | wc -c
```
