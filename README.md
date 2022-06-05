Illumina
--
(read ) 16bp cell barcode
(read1) 12bp umi
(read2) 55bp cDNA

Ultima
--
(read1) 22bp insert
(read1) 26bp cell barcode
(read1) 12bp umi
(read1) 12bp polyT homopolymer
(read1) XXbp cDNA

Ultima trimming procedure for mapping to TMSB4X generates minimum 31bp cDNA with 
```
seqkit subseq -r 63:-1 SRR18145555.fastq.gz | gzip > cdna_63.fastq.gz
seqkit seq -m 31 cdna_63.fastq.gz | gzip > cdna_31_min.fastq.gz
```

Ultima pseudoalignment discards reads with readlen < 31. For illumina-ultima gene finding with pseudoalignment, set ultima maxreadlen to be 55.
This is done by restricting kb with
```
            bc      umi     min 55bp cDNA
kb count -x 0,22,38:0,38,50:0,62,117
```

Ultima pseudoalignment (max 55bp cDNA)
```
kb count --overwrite --strand unstranded -w ../../reference/10x_version3_whitelist.txt -o out_trim -i ../.
./reference/transcriptome/index -x 0,22,38:0,38,50:0,62,117 -g ../../reference/transcriptome/t2g.txt --h5ad -m 16 -t 16 SRR18145555.fastq.gz
```

Ultima pseudoalignment (no max cDNA)
```
kb count --overwrite --strand unstranded -w ../../reference/10x_version3_whitelist.txt -o out -i ../../reference/transcriptome/index -x 0,22,38:0,38,50:0,62,0 -g ../../reference/transcriptome/t2g.txt --h5ad -m 16 -t 16 SRR18145555.fastq.gz
```

Illumina pseudoalignment
```
kb count --overwrite --strand unstranded -o out -i ../../reference/index -x 10xv3 -g ../../reference/t2g.txt --h5ad -m 16 -t 16 SRR18145553_1.fastq.gz SRR18145553_2.fastq.gz
```
