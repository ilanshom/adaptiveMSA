This repository contains the code for adaptiveMSA, a modified version of UPP that uses J-score instead of bitscore, along with a multi-armed bandit adaptive approach for assigning query sequences to HMMs.
The repository for the original UPP code is here: https://github.com/smirarab/sepp

## Installation

In order to install adaptiveMSA, follow the instructions for installing UPP in the file tutorial/upp-tutorial.md with the following modifications.  
First,

`git clone https://github.com/smirarab/sepp.git`

should be replaced with:

`git clone https://github.com/ilanshom/adaptiveMSA.git` 

The installation commands therefore need to be typed in the "adaptiveMSA" folder instead of the "sepp" folder.

Also, after installing this modified version of UPP, move to the adaptiveMSA folder, and type the following command in order to set permissions for executable files.  This resolves a bug in the original UPP installation process on some machines.

`chmod u+x .sepp/bundled-v4.5.2/hmmbuild; chmod u+x .sepp/bundled-v4.5.2/hmmalign; chmod u+x .sepp/bundled-v4.5.2/hmmsearch; chmod u+x .sepp/bundled-v4.5.2/guppy; chmod u+x .sepp/bundled-v4.5.2/pplacer;`

## Running the code

A default run for this modified version of UPP is called using the same command as UPP:

`run_upp.py -A 10 -B 1000 -M -1 -m molecule_type -s input`

We can set the arguments for Algorithm 1 in the paper manually.  The arguments are described briefly when the following command is typed into the terminal:
`run_upp.py -h`

We also describe them in more detail here.

### kmer size
argument: `--kmer_size k` \
This determines the kmer size used by the algorithm for assigning query sequences to HMMs.  In Algorithm 1 in the paper, this is $k$.

### kmers to check
argument: `--kmers_to_check B` \
This determines how many kmers to sample from the query sequence in order to estimate the J-score for each of the candidate HMMs in each round of the adaptive algorithm.  In Algorithm 1 in the paper, this is $B$.

### top sets to check
argument: `--top_sets_to_check T` \
This determines how many of the top scoring HMMs for a query sequence to compute the J-score on after the adaptive estimation rounds are complete.  The HMM from with the best J-score among these top HMMs will be assigned to the query sequence.  In Algorithm 1 in the paper, this is $T$.

### sample rounds
argument: `--sample_rounds R` \
This determines the number of rounds of sampling in the adaptive estimation stage of the algorithm.  In Algorithm 1 in the paper, this is $R$.

### use intersection score
argument: `--use_intersection_score` \
Adding the flag causes the algorithm to use the intersection score instead of the J-score for choosing the best HMM for a query sequence.  The intersection score is given by the numerator in the formula for the J-score.

### exact computation
argument: `--exact_computation` \
Adding this flag causes the algorithm to compute the J-score exactly for every HMM for every query sequence instead of estimating the J-score adaptively.  

### choose K
argument: `--choose_K` \
Adding this flag causes the algorithm to choose K based on how many query sequences share no kmers with any backbone sequences.

### Ks to check
argument: `--Ks_to_check C` \
When K is chosen from a set, this argument determines the set of candidate K values that will be checked.  The K chosen is the largest value of K from the set such that the fraction of query sequences that share no kmers with the any backbone sequence is below a threshold.  The values in the set should be seperated by commas without spaces e.g. "5,9" for {5, 9}.

### sequences unmatched threshold
argument: `--seqs_unmatched_thresh t` \
When K is chosen from a set, this argument determines the threshold for the fraction of query sequences that share no kmers with any backbone sequences.
