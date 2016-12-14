# README

How To Run:
==========
1) Execute make

2) Run 'python stonyccs.py' with the required options. Do '--help' for details

3) To clean, use 'make clean'

Notes:
=====
1) This script must be run only on pacbio bam files sorted by queryname 
   (Use 'samtools sort -n' to sort an unsorted .bam file)

External dependencies:
=====================
1) samtools - http://www.htslib.org/doc/samtools-1.1.html
2) pysam python module - https://github.com/pysam-developers/pysam