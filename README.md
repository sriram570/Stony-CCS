# README

This is a consensus calling tool created as part of the CSE 549 - Computational Biology final project at 
Stony Brook University. It aims to improve the consensus calling methods used to call Pacbio long reads.

Contributors:
1) Swaminathan Sivaraman
2) Sriram Sundar
3) Shyam Sundar Chandrasekaran
4) Prasanth Sridhar

(Stony Brook University, New York, USA)

Details:
=======
The main tool is a python script named stonyccs.py that takes as input a pabio subreads bam(sorted) file and calls
a consensus for each bunch or "well" of reads. It makes use of the poaligner.py and consensus.py modules to do the
consensus alling.

The tool performs ordering heuristics using STAR or a modified "forward-backward" version of STAR before doing a POA.
After doing a POA and generating a graph, the tool can use different scoring and traversal algorithms to generate
the consensus string. Both the ordering heuristic and the scoring and traversal algorithms can be changed using
command-line options.

The script also needs a scoring matrix. Two popular matrices - blosum62 and blosum80 have been included to be used.
(blosum62 is generally preferred for Pacbio data)

How To Run:
==========
1) Execute make

2) Run 'python stonyccs.py' with the required options. Do '--help' for details

3) To clean, use 'make clean'

Notes:
=====
1) This script must be run only on pacbio bam files sorted by queryname 
   (Use 'samtools sort -n' to sort an unsorted .bam file)
   
2) We have included Christopher Lee's POA library in this repository to do the actual Partial-Order Alignment. It
   is included in the 'external' directory (Link - https://sourceforge.net/projects/poamsa/)

External dependencies:
=====================
1) samtools - http://www.htslib.org/doc/samtools-1.1.html

2) pysam python module - https://github.com/pysam-developers/pysam