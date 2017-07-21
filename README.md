# ribo_seq

for mapping of ribosome footprints to a mammalian transcriptome using STAR

## Getting Started
*   Download all of the dependencies and make sure they are in your PATH or PYTHONPATH respectively
*   Make a STAR mapping index for your transcriptome of interest. A full human or mouse genome index require 64Gb of RAM to run.
*   Set --genomeSAsparseD to 2 for a 32Gb system, 3 for 16, and so on. This will increase run time, but it's still faster than bowtie.
*   Make tab-delimited annotation files. Pick a single transcript isoform for each gene (usually the longest coding sequence):
all_coding_Transcripts.tsv : all starts and ends are relative to transcript start site (which is indexed as zero)
tx_id	strand	tx_length	exon_starts	tx_exon_ends	cdsStart	cdsEnd	cdsSeq	geneID	notes
uc007aho.1	-	3526	0,184	183,3525	95	488	ATGAAGTTTTTAGAGAAAGGAGAGCTTGCAAACTTCAGATTCCAAAAGGATTTCTTACGACCTTTTGAACATATAATGAAACGAAACAGGTCTCCAACAATTCGAGATATGGTTGTACGGTGTATAGCACAAATGGTTAATTCTCAGGCTGCAAATATTCGTTCAGGATGGAAGAACATTTTCTCAGTATTCCATCTAGCTGCATCAGACCAAGATGAAAGCATAGTAGAACTTGCATTTCAGACAACAGGGCACATTGTCAGTAAGTATTTTTTAACTATTCAAGTGCAAATAGAAAAGCTGGATGTACTTAGTTGGCAAATGAGGTGCAGTAAGATTATGACTGTAGTATGGCTTTTAGACTTACAAATGTTGTTTTTAAAACTAGAGTGA
all_coding_transcripts.fa: fasta file of coding transcript sequences. names must mathc tx_id in above files.
```
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir mm10_STAR --genomeFastaFiles mm10/mm10.fa --sjdbGTFfile mm10/20170628_known_canonical.gtf --genomeSAsparseD 2
```
*   gzip your fastq files and put them in one folder
*   make a settings file for your experiment. An example file is included, and an annotated version is included, as comments are not allowed in json files.

### Prerequisites
Dependencies
*    skewer read adaptor trimmer downloaded and added to path (https://github.com/relipmoc/skewer) (0.2.2 tested)
*    STAR RNA-seq aligner downlaoded and added to path (tested versionSTAR 020201)
*    samtools (http://www.htslib.org/ version 1.5 tested) compiled and added to path (make install may do this automatically)
*    FASTX-toolkit (http://hannonlab.cshl.edu/fastx_toolkit/) tested with 0.0.14

Python dependencies (for python 2.7) installed with pip install (on mac, a local install with --user may be recommended)
*    simplejson (3.11.1)
*    numpy (1.13.1)
*    scipy (0.19.1)
*    matplotlib (2.0.2)
*    seaborn (0.8)
*    pysam (0.11.2.2)

