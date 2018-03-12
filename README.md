<strong>Decision Tree - DNA Classifier</strong><br/>
<br/>
Splice junctions are points on a DNA sequence at which "superfluous" DNA is removed during the process of protein creation in higher organisms. The problem posed in this dataset is to recognize, given a sequence of DNA, the boundaries between exons (the parts of the DNA sequence retained after splicing) and introns (the parts of the DNA sequence that are spliced out). This problem consists of two subtasks: recognizing exon/intron boundaries (referred to as EI sites), and recognizing intron/exon boundaries (IE sites). (In the biological community, IE borders are referred to as "acceptors" while EI borders are referred to as "donors".)<br/>
This dataset has been developed to help evaluate a "hybrid" learning algorithm (KBANN) that uses examples to inductively refine pre-existing knowledge. Using a "ten-fold cross-validation" methodology on 1000 examples randomly selected from the complete set of 3190, the following error rates were produced by various ML algorithms (all experiments run at the Univ of Wisconsin, sometimes with local implementations of published algorithms).<br/>
Attributes predicted: given a position in the middle of a window 60 DNA sequence elements (called "nucleotides" or "base-pairs"), decide if this is a :<br/>
a) "intron -> exon" boundary (ie) [These are sometimes called "donors"]<br/>
b) "exon -> intron" boundary (ei) [These are sometimes called "acceptors"]<br/>
c) neither (n)<br/>
<br/>
<br/>
This project is hard coded without using machine learning libraries.</br>
To run the project, please use the makefile with the command:<br/>
make -f makefile arg1="E" arg2="95" <br/>
<br/>
Input options: <br/>
arg1="E" to choose building a decision tree with InfoGain and Entropy<br/>
arg1="G" to choose building a decision tree with Gini index<br/>
arg2="99" to choose confidence level of 99%<br/>
arg2="95" to choose confidence level of 95%<br/>
<br/>
The best accuracy was produced under arg1="E" arg2="95".<br/>
<br/>
When the project terminates, "answer.csv" file will be generated
from testing the "testing.csv" file.


