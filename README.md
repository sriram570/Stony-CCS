# README #

Setup:
-----
1) Do this to use the poaligner python module:
    a. Make the poa program in external/poaV2
    b. `bash --init-file setup_env.sh`

Links and Info:
--------------
Read the following papers in order:
1) http://bioinformatics.oxfordjournals.org/content/18/3/452.full.pdf

2) http://bioinformatics.oxfordjournals.org/content/20/10/1546.full.pdf

3) Details of Consensus builder: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.579.863&rep=rep1&type=pdf
Pysam module docs:
http://pysam.readthedocs.io/en/latest/

Dataset:
https://downloads.pacbcloud.com/public/SequelData/ArabidopsisDemoData/SequenceData/

SAM/BAM file format:
https://samtools.github.io/hts-specs/SAMv1.pdf

Repositories to note:
--------------------
1) https://sourceforge.net/projects/poamsa/
   C++ implementation from the main papers' author - Lee

2) https://github.com/ljdursi/poapy
   Python implementation - Not sure if it has poa-poa alignment

3) https://github.com/mculinovic/cpppoa
   Another C++ implementation - Not sure if it has any new features

SAM/BAM Format:
--------------
https://samtools.github.io/hts-specs/SAMv1.pdf

Optional Fields:
----------------
http://bioinfo-master.ird.fr:8080/smrtanalysis/doc/bioinformatics-tools/pbsamtools/doc/index.html

Details of CCS POA:
---------------
https://github.com/PacificBiosciences/unanimity/blob/master/src/poa/POA-NOTES.txt

Basically they are using the POA-Sequence comparison idea. NOT the POA-POA comparison idea. Do read this. It also explains the scoring function that CCS use - the one that Avi has used in his slide.

The scoring function in code is located inside FindConsensus() in PoaGraphImpl.cpp (https://github.com/PacificBiosciences/unanimity/blob/master/src/poa/PoaGraphImpl.cpp)

Avi's mail:
----------

Copy pasting mail from Avi:

"Hi guys,
Sorry, I took little too long to get back to you but hopefully you have started thinking about the project and directions 
you guys would be moving forward with, in the coming weeks.
As promised you can check Arabidopsis
(https://github.com/PacificBiosciences/DevNet/wiki/Arabidopsis-Sequence-and-Assembly-%28_Arabidopsis-thaliana_-Ler-0%29)
and some other datasets here(https://github.com/PacificBiosciences/DevNet/wiki/Datasets).
This(https://downloads.pacbcloud.com/public/SequelData/ArabidopsisDemoData/SequenceData/1_A01_customer/m54113_160913_184949.subreads.bam)
is the direct link to Arabidopsis sub-reads.bam file.
Genome/Reference for Arabidopsis can be found here(ftp://ftp.ensemblgenomes.org/pub/release-32/plants/fasta/arabidopsis_thaliana/dna/) 

Some libraries which might help:
python: pyfasta, pysam
samtools : A simple tool to convert sam to bam or vice-versa {brew install samtools}
Knowing about tools like awk,sed would be helpful to do fast and dirty testing.

command to see how reads look like:
```samtools view -bS <path to bam file> | head```

I am working on a shell script which will help you install ccs tool from Pac-bio repo, will forward you guys soon. 
Also, like I already told you guys the main repo with all Pac-bio features is unanimity but we can hack our way to the concern files
using the deprecated repo which can be found here (https://github.com/PacificBiosciences/pbccs/tree/28383dc13b1d59659b6ac1d0cce45870bac501a7).
Hope all these details would help.
Let me know if you have any further question."

"There are many ways to install the required software tools and you are more than welcome to explore other methods 
but I have found the following the easiest.
Once you are able to successfully execute the attached script you would be able to run many tools from Pac-Bio but 
most importantly the following two:
1. blasr(https://github.com/PacificBiosciences/blasr) : this tool helps you map reads to reference Genome.
2. ccs(https://github.com/PacificBiosciences/unanimity/blob/master/doc/PBCCS.md#running-ccs): this is our legacy tools which we
	are trying to optimize.

NOTE: remember to use ccs with --noPolish flag when comparing.

To run this tools you need to be in the proper environment, which can be invoked by the command: 
``bash --init-file deployment/setup-env.sh`` after you have executed the attached bash script."

# a small snippet to read a bam file using pysam  
import pysam
bamf = pysam.AlignmentFile('mapt.NA12156.altex.bam', 'rb')
for i, line in enumerate(bamf.fetch(until_eof=True)):
    print line
    if i > 5:
        break
bamf.close()

# Notes:
- Setting up blasr and ccs
1. execute the command from *pitchfork* directory 
*bash --init-file deployment/setup-env.sh*
2. Checking if it is installed:
2.1 *ccs --version* returns *ccs 2.0.5 (commit 1ef34a1)*
2.2 *blasr --version* returns *blasr	5.3.e901e48*
3. then can use *blasr* and *ccs* commands from *any* directory

- the column names for the bam files, printed either from the command-line or through python 
col_names = 'QNAME\tFLAG\tRNAME\tPOS\tMAPQ\tCIGAR\tRNEXT\tPNEXT\tTLEN\tSEQ\tQUAL\tOPTIONAL'

- sam file is used to store sequence data in a csv-ish format. A bam file is a sam file but compressed and stored as binary

- pysam can read both 'bam' and 'sam' files directly without any conversion between the two formats

- Summary of all this: We will write a python script, use the pysam module to read our dataset bam file, reverse reads if necessary and group them based on id. The output can possibly be stored in a simple text file.

Essential Commands:
-------------------

Blasr: https://github.com/PacificBiosciences/blasr/blob/master/README.MANUAL.md

CCS: ccs --noPolish <input subread>.bam <output file>.bam

Blasr output format:
--------------------
https://github.com/PacificBiosciences/blasr/wiki/Blasr-Output-Format

Output of pruned BAM with a random reference genome - Blasr -m 4
---------------------------------------------------
m54113_160913_184949/7864612/ccs/0_3016 2 -6592 87.0339 0 11 1752 3016 1 10752497 10754180 19698289 254
