Problem Description:

Splice junctions are points on a DNA sequence at which `superfluous' DNA is removed during the process of protein creation in higher organisms. The problem posed in this dataset is to recognize, given a sequence of DNA, the boundaries between exons (the parts of the DNA sequence retained after splicing) and introns (the parts of the DNA sequence that are spliced out). This problem consists of two subtasks: recognizing exon/intron boundaries (referred to as EI sites), and recognizing intron/exon boundaries (IE sites). (In the biological community, IE borders are referred to as "acceptors" while EI borders are referred to as "donors".)
This dataset has been developed to help evaluate a "hybrid" learning algorithm (KBANN) that uses examples to inductively refine pre-existing knowledge. Using a "ten-fold cross-validation" methodology on 1000 examples randomly selected from the complete set of 3190, the following error rates were produced by various ML algorithms (all experiments run at the Univ of Wisconsin, sometimes with local implementations of published algorithms).
Attributes predicted: given a position in the middle of a window 60 DNA sequence elements (called "nucleotides" or "base-pairs"), decide if this is a :
a) "intron -> exon" boundary (ie) [These are sometimes called "donors"]
b) "exon -> intron" boundary (ei) [These are sometimes called "acceptors"]
c) neither (n)



To run the project, please use the makefile with the command:
make -f makefile arg1="E" arg2="95" 

Input options: 
arg1="E" to choose computing InfoGain with entropy
arg1="G" to choose computing InfoGain with Gini index
arg2="99" to choose confidence level of 99%
arg2="95" to choose confidence level of 95%

The best accuracy was produced under arg1="E" arg2="95".

When the project terminates, "answer.csv" file will be generated
from testing the "testing.csv" file.


