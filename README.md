# In silico sequence capture to assembly germline-restricted chromosomes sequences from 10xGenomics Chromium libraries

This protocol was applied for the manuscript by Pei et al. (2021) "Occasional paternal inheritance of the germline-restricted chromosome in songbirds", accepted in PNAS.

It tries to assembly the sequence of few genes in a germline-restricted chromosome (GRC), reducing the interferences of their paralogs in the regular (A) chromosomes.

## 0. Required files

* A 10xG Chromium library for testis (with GRC):
```
P8503_1005_S17_L007_I1_001.fastq.gz
P8503_1005_S17_L007_R1_001.fastq.gz
P8503_1005_S17_L007_R2_001.fastq.gz
```

* Other 10xG Chromium library for a somatic tissue (without GRC) from the same individual:
```
P8503_1004_S14_L006_I1_001.fastq.gz
P8503_1004_S14_L006_R1_001.fastq.gz
P8503_1004_S14_L006_R2_001.fastq.gz
```

* A fasta file with the version of the A chromosomes (Ref) and the same one with GRC-specific alleles (Alt) found using a SNP calling.

```
sequences_ref_alt.fasta
```

## 1. Remove barcodes from both 10XG Chromium libraries from the left read and add them to the read ID

Run this longranger command using your preferred name for the --id flag:

```
$ longranger basic --id=sample345 --fastqs=/PATH/TO/10XG/READS
```

It generates a “barcoded.fastq.gz” file. Please, note that the reads in the beginning of the file do not show barcodes! 

## 2. Unshuffle reads in both libraries

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

After that, we can perform a trimming with Trimmomatic to remove low quality nucleotides and adapters.

## 3. First round of mapping: against the ref and alt sequences

Since the 10xG reads contain introns and the reference sequences do not contain them, we need to use a mapper considering soft clipping for the definitive maps. Here we use SSAHA2, since we can easily control the mapping legth and identity. It can be run in a multithread way with this commands:

```
$ ls barcoded_trim_1.fastq barcoded_trim_2.fastq > list.txt
$ ssaha2_run_multi.py list.txt sequences_ref_alt.fasta 20
```
Where "20" is the number of threads you can choose.

It will generates a sorted and indexed BAM file with the suffix "mapped".

Then, we only keep reads mapping from a BAM file into Alt sequences:

```
$ mapped.bam > bam_list.txt
$ grep "alt" sequences_ref_alt.txt | sed 's/>//g' > alt_list.txt
$ extract_seq_bam.py bam_list.txt alt_list.txt sequences_ref_alt.txt
```
It generates a BAM file with the suffix "mapped.sel.sort"

## 4. Extract reads

You can do it in two alternative ways.

### 4a. Using a custom script with the pipeline

You need to indicate the BAM file, the barcoded fastq file from "longranger basic" and the beginning of the read IDs

```
$ extract_barcoded_reads_from_bam.py mapped.sel.sort.bam barcoded.fastq.gz @ST-E00
```

### 4b. Running each step of the previous script manually

Alternatively, you can run each step from the last script separately.

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

## 5. Extract reads from the raw 10xG Chromium reads:

```
$ seqtk subseq P8503_1005_S17_L007_I1_001.fastq.gz barcodes_reads.txt > P8503_1005_S17_L007_I1_001.fastq
$ seqtk subseq P8503_1005_S17_L007_R1_001.fastq.gz barcodes_reads.txt > P8503_1005_S17_L007_R1_001.fastq
$ seqtk subseq P8503_1005_S17_L007_R2_001.fastq.gz barcodes_reads.txt > P8503_1005_S17_L007_R2_001.fastq
```
Then we prepare the FASTQ files in a separate folder moving them to a new folder and compressin them as GZ:

```
$ mkdir reads
$ cd reads
$ mv ../*q .
$ gzip *q
```

## 6. Assembly with Supernova2

Run Supernova2 assembly with the following options:

```
$ supernova run --accept-extreme-coverage --maxreads="all" --id=sample345 --fastqs=/PATH/TO/SELECTED/READS
```

This will generate a Fasta file with the assembled reads. We can search for contigs matching the genes we used as a reference with BLAST, Exonerate or other local aligner

### 7. (Optional) Second round of mappings, barcode extraction and assembly

This step is only necessary in case we need longer contigs. So we can perform several round of assembly. However, from now on, the mappings are a bit different.

For the second round of mappings, we will use the assembled and selected Supernova2 contigs. In this case, we have introns and intergenic regions identical to the A chromosome paralog. To avoid mappings in such as regions we first map soma reads agaist with more stringent conditions:

```
$ ls barcoded_trim_1.fastq barcoded_trim_2.fastq > list.txt
$ ssaha2_run_multi_id99_len80.py list.txt supernova_sel.fasta 20
```

Then, mask positions in the reference with mapped reads in the soma library with Bedtools

```
$ bedtools bamtobed -i MAPPING.bam > MAPPING.bed
$ bedtools maskfasta -fi supernova_sel.fasta -bed MAPPING.bed -fo supernova_sel_mask.fasta
```

After this, we map the testis reads against the masked reference:

```
$ ls barcoded_trim_1.fastq barcoded_trim_2.fastq > list.txt
$ ssaha2_run_multi_id99_len80.py list.txt supernova_sel_mask.fasta 20
```

Then repeat steps 5 and 6 to extract reads with selected barcodes and assemble them.
