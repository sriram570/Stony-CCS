# README

This is a consensus calling tool (named stonyccs) created as part of the CSE 549 - Computational Biology final project at 
Stony Brook University. It aims to improve the consensus calling methods used to call Pacbio long reads.

(An included report.pdf file details the methods attempted and results observed)

Contributors:

1) Swaminathan Sivaraman (110951180)

2) Sriram Sundar (110921718)

3) Shyam Sundar Chandrasekaran (110815338)

4) Prasanth Sridhar (110899181)

December 2016, Stony Brook University, New York, USA

Details:
=======
The main tool is a python script named stonyccs.py that takes as input a Pacbio subreads bam(sorted) file and calls
a consensus for each bunch or "well" of reads. It makes use of the included poaligner.py and consensus.py modules 
to do the consensus calling.

The tool performs ordering heuristics using STAR or a modified "forward-backward" version of STAR before doing a POA.
After doing a POA and generating a graph, the tool can use different scoring and traversal algorithms to generate
the consensus string. Both the ordering heuristic and the scoring and traversal algorithms can be changed using
command-line options.

The script also needs a scoring matrix. Two popular matrices - blosum62 and blosum80 have been included to be used.
(blosum62 is generally preferred for Pacbio data)

One can use the sample_bam.bam file provided in the samples_tests directory to
run quick tests of or see the working of the tool.

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

3) blasr and pacbio's ccs tools are from PacbioSciences' GitHub repository

External dependencies:
=====================
1) samtools - http://www.htslib.org/doc/samtools-1.1.html

2) pysam python module - https://github.com/pysam-developers/pysam

Tests:
=====
We have performed tests on the new stonyccs tool and compared it with pacbio's
ccs(baseline) results. There are two kind of tests run:

1) Sample well tests:

   Here, we took a single well of reads (total 6 reads) and ran extensive tests using
   various combinations of stonyccs' algorithms. We mostly used blosum62 and
   sometimes blosum80 as the scoring matrices. All test results are present in
   the sample_tests directory. The bam file used is named as sample_bam.bam. Each
   test case has its own directory named after the test case and inside, there is
   a consensus.fa file, blasr_report.txt and stonyccs_report.txt, which has the
   consensus, the mapping details to the reference and the full details of the
   configuration used respectively. We also ran with two more configurations - 
   no_filters, to not filter the input reads at all, and all_filters, to filter
   the input reads based on read quality, length etc. As no_filters casesperformed
   really poorly, most of the testcases are for the all_filters cases and are in the
   sample_tests/all_filters directory.

NOTE: One can use this sample bam file to perform quick tests on the stonyccs tool.

2) Full data tests:

   Here, we took the pacbio bam file (from https://downloads.pacbcloud.com/public/SequelData/ArabidopsisDemoData/SequenceData/3_C01_customer/)
   and sorted it by queryname and created a sorted bam file. We then chose the
   best-performing configurations from the above sample well tests and ran tests
   for those configurations on this full dataset. All results are filed under
   full_data_tests. The test case details are saved in a manner similar to the
   above sample well cases. Additionally, there is a scores.txt file in each test case
   directory, that contains the list of all match,mismatch,ins,del,sim scores for all
   consensus strings and the average values for each field.

For both cases, the scores and alignment details for Pacbio's original ccs
tool have been included under the baseline directory.

Reference data used for blasr:
ftp://ftp.ensemblgenomes.org/pub/release-32/plants/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.chromosome.1.fa.gz

Testcase organization:
---------------------
1) The sample data test results are in the sample_tests/no_filters/ and sample_tests/all_filters/ directories

2) The full data test results are in the full_data_tests/ directory

3) Each test case is present as a directory. The directory names are shorthand names for the ordering algorithm used, the scoring function
   used and the traversal algorithm used. Finally, it is specified if the blosum62 or blosum80 matrix was used to generate the consensus
   
4) Each test case directory has a consensus.fa file, a stonyccs_report.txt file and a blasr_report.txt file. For full_data test cases,
   there is also a scores.txt file which summarizes multiple scores and gives an average
   
5) The original pbccs (Pacbio's ccs with --noPolish flag) results (baseline) are also saved for both sample and full data cases in the sample_tests/baseline/ and
   the full_data_tests/baseline directories