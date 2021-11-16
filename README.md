# Extract 10xG barcodes
Protocol to detect and extract reads with a specific 10xG barcode

## 0. Required files

A 10xG Chromium library:
```
P8503_1005_S17_L007_I1_001.fastq.gz
P8503_1005_S17_L007_R1_001.fastq.gz
P8503_1005_S17_L007_R2_001.fastq.gz
```

A BAM file of mappings using your favourie mapper:
```
MAPPING.bam
```

## 1. Remove barcodes from 10XG Chromium libraries from the left read and add them to the read ID

Run this longranger command using your preferred name for the --id flag:

```
$ longranger basic --id=sample345 --fastqs=/PATH/TO/10XG/READS
```

It generates a “barcoded.fastq.gz” file. Please, note that the reads in the beginning of the file do not show barcodes! 

## 2. Unshuffle reads

Since we got a single “barcoded.fastq.gz” file with both the left and right reads interleaved:

```
@ST-E00214:160:H277TCCXY:7:2211:28036:31951 BX:Z:AAACACCAGCGATATA-1
ATATCGCTCTTCCGATCTAAACACCAGCGATATATCTAATACTGTTCCTTCGGAATTTGCAAAGACTTACTTTCTGAACAGGCCATAGAAGAACCTCCAGGTGTTTTCCACACTGCTGGAATGTTTGT
+
JJJJJJFJJ<FJJJ7AFJJJJJJJJJ<FJJJJJJ<AF<FJJJFFJJFJJJJJ-7FFJJJJJAJAJJJAFJFJFFAJFJJJJJAFFFJJJJJAJFF7FAFFF-JJ7<-7-<<JJ-AJAJFA-FFF7A<-
@ST-E00214:160:H277TCCXY:7:2211:28036:31951 BX:Z:AAACACCAGCGATATA-1
CCTTCGGAATTTGCAAAGACTATATTCTGTATTGATGGCTAAATGCTATGTGAGCTATATTTGACCTCTAAAAGTAGAAAAAATGTTAATAAACATTCCAGCAGTGTGGAAAACACCTGGAGGTTCTTCTATGGCCTGTTCAGAAAGTAAG
+
AAFFFJJJJJJJJJJJJJFJFFJJJJJJFJJJFFJJFJFJFFFFJ<<7FJJJ<JJAJFJJJJJJJJFF7AFJFJJJJJJJFFF<FFJJJ-<FJJJJJFFJJFFJFF7JJFJJJFJJJJJFJAJ<-<JFFAJJF<JAFJA-A<<<<F<-<<A
@ST-E00214:160:H277TCCXY:7:2118:19471:11945 BX:Z:AAACCCAAGAATCTCC-1
CTTGAGCTACTTTAAGGCTATTGGTTACAATCTCAACAGTAAGCCATATAATAGTTTATATTTGCCTCCTTCTTGAAAATATCTGTGTTTGCTTTTGGGTGTTCTTTTGTATATACTCCAAATCATTA
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJFJFJJJJJJJJJFFJJJJJJJJJJJJJJJF
@ST-E00214:160:H277TCCXY:7:2118:19471:11945 BX:Z:AAACCCAAGAATCTCC-1
CCTTGCAGAACCCCTTGGCCAAGCAGAACAAACAATTACAGAACTCTCTCACGTATTTTCTAATGATTTGGAGTATATACAAAAGAACACCCAAAAGCAAACACAGATATTTTCAAGAAGGAGGCAAATATAAACTATTATATGGCTTACT
+
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJFJFJJJJJJJJJJJJJJJJA<FJ<
```

We need to split it in two files: one for the left reads and other for the right reads. Using, for example this external script https://github.com/ndaniel/fusioncatcher/blob/master/bin/unshuffle.py:

```
$ unshuffle.py -i barcoded.fastq.gz -f barcoded_1.fastq -r barcoded_2.fastq
```

This is the resulting barcoded_1.fastq file:

```
@ST-E00214:160:H277TCCXY:7:2211:28036:31951 BX:Z:AAACACCAGCGATATA-1
ATATCGCTCTTCCGATCTAAACACCAGCGATATATCTAATACTGTTCCTTCGGAATTTGCAAAGACTTACTTTCTGAACAGGCCATAGAAGAACCTCCAGGTGTTTTCCACACTGCTGGAATGTTTGT
+
JJJJJJFJJ<FJJJ7AFJJJJJJJJJ<FJJJJJJ<AF<FJJJFFJJFJJJJJ-7FFJJJJJAJAJJJAFJFJFFAJFJJJJJAFFFJJJJJAJFF7FAFFF-JJ7<-7-<<JJ-AJAJFA-FFF7A<-
@ST-E00214:160:H277TCCXY:7:2118:19471:11945 BX:Z:AAACCCAAGAATCTCC-1
CTTGAGCTACTTTAAGGCTATTGGTTACAATCTCAACAGTAAGCCATATAATAGTTTATATTTGCCTCCTTCTTGAAAATATCTGTGTTTGCTTTTGGGTGTTCTTTTGTATATACTCCAAATCATTA
+
JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJFJFJJJJJJJJJFFJJJJJJJJJJJJJJJF
```

This is the barcoded_2.fastq file:

```
@ST-E00214:160:H277TCCXY:7:2211:28036:31951 BX:Z:AAACACCAGCGATATA-1
CCTTCGGAATTTGCAAAGACTATATTCTGTATTGATGGCTAAATGCTATGTGAGCTATATTTGACCTCTAAAAGTAGAAAAAATGTTAATAAACATTCCAGCAGTGTGGAAAACACCTGGAGGTTCTTCTATGGCCTGTTCAGAAAGTAAG
+
AAFFFJJJJJJJJJJJJJFJFFJJJJJJFJJJFFJJFJFJFFFFJ<<7FJJJ<JJAJFJJJJJJJJFF7AFJFJJJJJJJFFF<FFJJJ-<FJJJJJFFJJFFJFF7JJFJJJFJJJJJFJAJ<-<JFFAJJF<JAFJA-A<<<<F<-<<A
@ST-E00214:160:H277TCCXY:7:2118:19471:11945 BX:Z:AAACCCAAGAATCTCC-1
CCTTGCAGAACCCCTTGGCCAAGCAGAACAAACAATTACAGAACTCTCTCACGTATTTTCTAATGATTTGGAGTATATACAAAAGAACACCCAAAAGCAAACACAGATATTTTCAAGAAGGAGGCAAATATAAACTATTATATGGCTTACT
+
AAFFFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJFJFJJJJJJJJJJJJJJJJA<FJ<
```

## 3. Extract reads

You can do it in two alternative ways.

### 3a. Using a custom script with the pipeline

```
$ extract_barcoded_reads_from_bam.py MAPPING.bam barcoded.fastq.gz @ST-E00
```

### 3b. Running each step of the previous script manually

#### Extract mapped reads in a BAM file

Replace "MAPPING.bam" by your file's name.

```
$ samtools view MAPPING.bam | cut -f1,10,11 | sed 's/^/@/' | sed 's/\134t/\134n/' | sed 's/\134t/\134n+\134n/' >> mapped_reads.fastq
```

#### Get read names from the mapped reads

Change the start of the read ID. Here, it is "@ST-E00".

```
$ grep "@ST-E00" mapped_reads.fastq | sort | uniq | sed 's/@//g' >> mapped_reads_names.txt
```

#### Extract reads from the file to get barcodes

Replace barcoded.fastq.gz by the name of your file after the "longranger basic" command (see above). It also works for gzipped files.

```
$ seqtk subseq barcoded.fastq.gz mapped_reads_names.txt >> mapped_reads_names.fastq
```

#### Get list of uniq barcodes

```
$ grep "BX:Z:" mapped_reads_names.fastq | awk {'print $2'} | sort | uniq >> mapped_reads_barcodes.txt
```

#### Generate a list of read names with a barcode in the previous list

```
$ zgrep -Ff mapped_reads_barcodes.txt barcoded.fastq.gz >> barcodes_names.txt
```

(Alternatively, use grep for non-gzipped barcoded files)

#### Extract read IDs

```
$ awk 'NR%2==1' barcodes_names.txt | awk {'print $1'} | sed 's/@//g' >> barcodes_reads.txt
```

### 4. Assembly with Supernova2

Extract reads from the raw 10xG Chromium reads:

```
$ seqtk subseq P8503_1005_S17_L007_I1_001.fastq.gz barcodes_reads.txt > P8503_1005_S17_L007_I1_001.fastq
$ seqtk subseq P8503_1005_S17_L007_R1_001.fastq.gz barcodes_reads.txt > P8503_1005_S17_L007_R1_001.fastq
$ seqtk subseq P8503_1005_S17_L007_R2_001.fastq.gz barcodes_reads.txt > P8503_1005_S17_L007_R2_001.fastq
```

Move FASTQ files to other location, compress them as GZ and run Supernova2 assembly:

```
$ mkdir reads
$ mv *q reads
$ supernova run --accept-extreme-coverage --maxreads="all" --id=sample345 --fastqs=/PATH/TO/SELECTED/READS
```

### BONUS. Mask positions in the reference with mapped reads in other library

This is a complementary approach to mask position in the reference.

```
$ bedtools bamtobed -i MAPPING.bam > MAPPING.bed
$ bedtools maskfasta -fi reference.fasta -bed MAPPING.bed -fo reference_mask.fasta
```
