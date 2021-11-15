# Extract 10xG barcodes
Protocol to detect and extract reads with a specific 10xG barcode

## 1. Remove barcodes from 10XG Chromium libraries from the left read and add them to the read ID

Files:<br />
P7359_105_S5_L008_R2_001.fastq.gz<br />
P7359_105_S5_L008_R1_001.fastq.gz<br />
P7359_105_S5_L008_I1_001.fastq.gz<br />

Command:
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

We need to split it in two files: one for the left reads and other for the right reads. Using, for example this commands:

```
$ aunpack barcoded.fastq.gz
$ unshuffle.py list.txt

$ awk '{ if (NR%4==1) { print $1"/1" } else { print } }'  uben_heart_bc_sel_1.fastq > uben_heart_bc_sel_checked_1.fastq
$ awk '{ if (NR%4==1) { print $1"/2" } else { print } }'  uben_heart_bc_sel_2.fastq > uben_heart_bc_sel_checked_2.fastq
```

----

For big files

```
$ split -l 10000000 barcoded.fastq
```

[xaa, xab,...]

```
$ gzip -dc barcoded_1_paired.fastq.gz | split -l 10000000  - split_ --filter='gzip > $FILE.fastq.gz'
```

--------

For compressed files:

awk '{ if (NR%4==1) { print $1"/2" } else { print } }' <(gzip -dc lib_2_unpaired.fastq.gz) | gzip > tgut2_l_s14_unpaired_2.fastq.gz

---------------------------------------------------------------------------------------------------

This is the barcoded_1.fastq file:

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

### 3a. Option1

New script for the next steps!!!!
```
$ extract_barcoded_reads_from_bam.py tgut_alex_testis_barcoded_checked_1_paired_10M_mapped.bam barcoded.fastq
```

---------------------------------------------------------------------------------------------------

### 3b. Option2

#### Extract mapped reads in a BAM file

Replace "MAPPING.bam" by your file's name

```
$ samtools view MAPPING.bam | cut -f1,10,11 | sed 's/^/@/' | sed 's/\134t/\134n/' | sed 's/\134t/\134n+\134n/' >> mapped_reads.fastq
```

#### Get read names from the mapped reads

Change the start of the read id. Here, it is "@ST-E00"

```
$ grep "@ST-E00" mapped_reads.fastq | sort | uniq | sed 's/@//g' >> mapped_reads_names.txt
```

#### Extract reads from the file to get barcodes

Replace barcoded.fastq by the name of your file after the "longranger basic" command (see above). It also works for gzipped files.

```
$ seqtk subseq barcoded.fastq mapped_reads_names.txt >> mapped_reads_names.fastq
```

#### Get list uniq barcodes

```
$ grep "BX:Z:" mapped_reads_names.fastq | awk {'print $2'} | sort | uniq >> mapped_reads_barcodes.txt
```

#### Generate a list of read names with a barcode in the previous list

```
$ grep -Ff mapped_reads_barcodes.txt barcoded.fastq >> barcodes_names.txt
```

(Alternatively, use zgrep for gzipped barcoded files)

#### Extract read IDs

```
$ awk 'NR%2==1' barcodes_names.txt | awk {'print $1'} | sed 's/@//g' >> barcodes_reads.txt
```

#### Extract reads

```
$ seqtk subseq barcoded.fastq barcodes_reads.txt >> barcodes_reads.fastq
```

XXXX

Select reads

```
$ seqtk subseq tgut_alex_liver_barcoded.fastq.gz barcoded_reads.all.txt > barcoded_reads_all.fastq

$ ls barcoded_reads_all.fastq > lili.txt

$ unshuffle.py lili.txt
```

### 3c. Option3
-------
MACROCOMMAND!!!!
```
$ seqtk subseq tgut_alex_testis_barcoded.fastq.gz lista_uniq.txt > tgut_alex_liver_ccnd3.fastq && grep "BX:Z:" med20_mapped_reads_sequence.fastq > sel && awk {'print $2'} sel | sort | uniq > barcodes.txt; zgrep -f barcodes.txt tgut_alex_testis_barcoded.fastq.gz > barcoded_reads.txt; awk {'print $1'} barcoded_reads.txt | sed 's/@//g' | sort | uniq > barcoded_reads.all.txt; seqtk subseq tgut_alex_testis_barcoded.fastq.gz barcoded_reads.all.txt > barcoded_reads_all.fastq; ls barcoded_reads_all.fastq > lala.txt; unshuffle.py lala.txt
```

--------------

### 4. Assembly with Supernova2

```
$ supernova run --accept-extreme-coverage --maxreads="all" --id=sample345 --fastqs=/PATH/TO/SELECTED/READS
```

### 5. Masking positions with mapped reads in soma

```
$ bedtools bamtobed -i heart_id99_len80.bam > heart_id99_len80.bed
$ bedtools maskfasta -fi ugra_round1_supernova_sel.fasta -bed heart_id99_len80.bed -fo ugra_round1_supernova_sel_mask.fasta
```

