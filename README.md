# Extract 10xG barcodes
Protocol to detect and extract reads with a specific 10xG barcode

To remove barcode from 10XG libraries:

Files:
P7359_105_S5_L008_R2_001.fastq.gz
P7359_105_S5_L008_R1_001.fastq.gz
P7359_105_S5_L008_I1_001.fastq.gz

Command:
$ ../../longranger-2.1.6/longranger basic --id=sample345 --fastqs=/mnt/lmigratoria/tgut_alex/ccnd3/test2

It generates a “barcoded.fastq.gz” file. The reads in the beginning of the file do not show barcodes!!!! 

$ aunpack barcoded.fastq.gz

$ unshuffle.py list.txt

$ awk '{ if (NR%4==1) { print $1"/1" } else { print } }'  uben_heart_bc_sel_1.fastq > uben_heart_bc_sel_checked_1.fastq
$ awk '{ if (NR%4==1) { print $1"/2" } else { print } }'  uben_heart_bc_sel_2.fastq > uben_heart_bc_sel_checked_2.fastq

----

For big files

$ split -l 10000000 barcoded.fastq

[xaa, xab,...]

$ gzip -dc barcoded_1_paired.fastq.gz | split -l 10000000  - split_ --filter='gzip > $FILE.fastq.gz'

--------

For compressed files:

awk '{ if (NR%4==1) { print $1"/2" } else { print } }' <(gzip -dc lib_2_unpaired.fastq.gz) | gzip > tgut2_l_s14_unpaired_2.fastq.gz

---------------------------------------------------------------------------------------------------

New script for the next steps!!!!

$ extract_barcoded_reads_from_bam.py tgut_alex_testis_barcoded_checked_1_paired_10M_mapped.bam barcoded.fastq

---------------------------------------------------------------------------------------------------

Extract reads from a BAM file

$ samtools view tgut_alex_liver_checked_1_mapped.bam | cut -f1,10,11 | sed 's/^/@/' | sed 's/\t/\n/' | sed 's/\t/\n+\n/' > mapped_reads_liver.fastq 

Extract names from read pair

$ grep "@ST-E00" mapped_reads_liver.fastq | sort | uniq | sed 's/@//g' > mapped_reads_liver_names_uniq.txt

Select read pairs with barcode 

$ seqtk subseq tgut_alex_liver_barcoded.fastq.gz lista_uniq.txt > tgut_alex_liver_ccnd3.fastq

Select unique barcodes

$ grep "BX:Z:" tgut_alex_liver_ccnd3.fastq | awk {'print $2'} | sort | uniq > barcodes.txt

Select reads with the list of barcodes

$ zgrep -f barcodes.txt tgut_alex_liver_barcoded.fastq.gz > barcoded_reads.txt

Extract read names

$ awk {'print $1'} barcoded_reads.txt | sed 's/@//g' | sort | uniq > barcoded_reads_-all.txt

Select reads

$ seqtk subseq tgut_alex_liver_barcoded.fastq.gz barcoded_reads.all.txt > barcoded_reads_all.fastq

$ ls barcoded_reads_all.fastq > lili.txt

$ unshuffle.py lili.txt

-------
MACROCOMMAND!!!!
$ seqtk subseq tgut_alex_testis_barcoded.fastq.gz lista_uniq.txt > tgut_alex_liver_ccnd3.fastq && grep "BX:Z:" med20_mapped_reads_sequence.fastq > sel && awk {'print $2'} sel | sort | uniq > barcodes.txt; zgrep -f barcodes.txt tgut_alex_testis_barcoded.fastq.gz > barcoded_reads.txt; awk {'print $1'} barcoded_reads.txt | sed 's/@//g' | sort | uniq > barcoded_reads.all.txt; seqtk subseq tgut_alex_testis_barcoded.fastq.gz barcoded_reads.all.txt > barcoded_reads_all.fastq; ls barcoded_reads_all.fastq > lala.txt; unshuffle.py lala.txt

--------------

Masking positions with mapped reads in soma

$ ml BEDTools
$ bedtools bamtobed -i heart_id99_len80.bam > heart_id99_len80.bed
$ bedtools maskfasta -fi ugra_round1_supernova_sel.fasta -bed heart_id99_len80.bed -fo ugra_round1_supernova_sel_mask.fasta

