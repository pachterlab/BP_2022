This Github repository contains the code to reproduce the figures and results from the preprint "Pseudoalignment facilitates assignment of
error-prone Ultima Genomics reads" by Sina Booeshaghi and Lior Pachter.

## Downloading data

Metadata was retrieved with [`ffq`](https://github.com/pachterlab/ffq/)

```bash
ffq -o ffq_GSM6297378.json GSM6297378 # Single-cell Illumina
ffq -o ffq_GSM6297379.json GSM6297379 # Single-cell Ultima
ffq -o ffq_GSM6190598.json GSM6190598 # Perturb-seq Illumina
ffq -o ffq_GSM6190599.json GSM6190599 # Perturb-seq Ultima
```

and downloaded with wget

```bash
# Illumina single-cell count matrices
wget --quiet --show-progress ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6297nnn/GSM6297378/suppl/GSM6297378_cells_counts_Three_Ill.txt.gz
wget --quiet --show-progress ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6297nnn/GSM6297378/suppl/GSM6297378_genes_counts_Three_Ill.txt.gz
wget --quiet --show-progress ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6297nnn/GSM6297378/suppl/GSM6297378_expression_counts_Three_Ill.txt.gz

# Ultima single-cell count matrices
wget --quiet --show-progress ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6297nnn/GSM6297379/suppl/GSM6297379_cells_counts_Three_Ult.txt.gz
wget --quiet --show-progress ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6297nnn/GSM6297379/suppl/GSM6297379_genes_counts_Three_Ult.txt.gz
wget --quiet --show-progress ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6297nnn/GSM6297379/suppl/GSM6297379_expression_counts_Three_Ult.txt.gz

# Perturb-Seq Illumina Cell Line
wget --quiet --show-progress https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6190nnn/GSM6190598/suppl/GSM6190598_illumina_cellranger.tar.gz

# Perturb-Seq Ultima Cell Line
wget --quiet --show-progress https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6190nnn/GSM6190599/suppl/GSM6190599_ultima_cellranger.tar.gz
```

FASTQ files can either be downloaded directly or converted from SRA files

```bash
# download FASTQ files directly
wget --quiet --show-progress ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR200/051/SRR20002551/SRR20002551_1.fastq.gz # Illumina R1
wget --quiet --show-progress ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR200/051/SRR20002551/SRR20002551_2.fastq.gz # Illumina R2
wget --quiet --show-progress ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR200/050/SRR20002550/SRR20002550.fastq.gz   # Ultima R1

# Download SRA files and covert to FASTQ files
wget --quiet https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR20002551/SRR20002551 # illumina
wget --quiet https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR20002550/SRR20002550 # ultima
# convert SRA file to FASTQ files
fasterq-dump ./SRR20002551 --split-files --include-technical
fasterq-dump ./SRR20002550 --split-files --include-technical
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

The TMSB4X fasta was downloaded with [`gget`](https://github.com/pachterlab/gget), and indexed with [`kallisto`](https://github.com/pachterlab/kallisto) and [`hisat2`](https://github.com/DaehwanKimLab/hisat2).

```bash
# download tmsb4x fasta
gget seq -o TMSB4X_ENST00000380636.1.fasta  ENST00000380636.1

# build kallisto index
kallisto index -i index TMSB4X_ENST00000380636.1.fasta
# and t2g map
echo -e "ENST00000380636\tENSG00000205542\tTMSB4X" > t2g.txt

# build hisat2 index
hisat2-build -p 16 TMSB4X_ENST00000380636.1.fasta genome
```

## Illumina

Illumina reads have the following structure:

- (read 1) 16bp cell barcode
- (read 1) 12bp umi
- (read 2) 55bp cDNA

Reads were pseudoaligned to the human transcriptome and to the TMSB4X gene using [`kb-python`](https://github.com/pachterlab/kb_python)

```bash
# pseudoalign to transcriptome
kb count --strand unstranded -o out-transcriptome -i human/index -x 10xv3 -g t2g.txt --h5ad -m 16G -t 16 SRR18145553_1.fastq.gz SRR18145553_2.fastq.gz

# pseudoalign to tmsb4x
kb count --strand unstranded -o out-tmsb4x -i tmsb4x/index -x 10xv3 -g t2g.txt --h5ad -m 16G -t 16 SRR18145553_1.fastq.gz SRR18145553_2.fastq.gz
```

Alignments and alignment stats were made using [`samtools`](https://github.com/samtools/samtools) and [`hisat2`](https://github.com/DaehwanKimLab/hisat2).

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

## Ultima Genomics

Ultima Genomics reads have the following structure:

- (read1) 22bp insert
- (read1) 26bp cell barcode
- (read1) 12bp umi
- (read1) 12bp polyT homopolymer
- (read1) XXbp cDNA

Reads were pseudoaligned using [`kb-python`](https://github.com/pachterlab/kb_python)

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

The Ultima Genomics reads were trimmed so that the minimum cDNA length was 31bp with [`seqkit`](https://github.com/shenwei356/seqkit)

```bash
# trim cDNA reads
seqkit subseq -r 63:-1 SRR18145555.fastq.gz | gzip > cdna_63.fastq.gz

# keep cDNA reads longer than or equal to 31 bp
seqkit seq -m 31 cdna_63.fastq.gz | gzip > cdna_31_min.fastq.gz
```

Alignments and alignment stats were made using [`samtools`](https://github.com/samtools/samtools) and [`hisat2`](https://github.com/DaehwanKimLab/hisat2).

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
