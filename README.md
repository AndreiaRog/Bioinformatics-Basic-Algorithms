# BioinformaticsBasicAlgorithms
Practical component of university masters subject.

Viterbi - Receives a sequence X of nucleotids as input. Outputs a path that maximizes P(X|π) over all possible paths π. the paths belong to an HMM model with three states, to identify DNA coding regions. State 1 corresponds to the Start Site signal, state 2 corresponds to an Exon region and State 3 corresponds to an Intro region. Initial probabilities for all the three states are
equal and transitions from all the states to an end state are also equal. Trnasition and Emission probabilities are difined in the code. 

Smith-Waterman - algorithm to compute the local alignment(s) between two amino acid sequences. The inputs are 2 amino acid sequences: S1 and S2 (two strings); the scoring model: gap penalty (integer) and scoring matrix (for isntance BLOSSUM 50).
The gap penalty function is considered to be linear in the number of gaps. It outputs the score of the optimal local alignment(s) between S1 and S2, and prints all possible optimal alignments between S1 and S2. 


