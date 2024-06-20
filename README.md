# CodonOptimizer

Simple tool to perform codon optimization for any organism that has a complete, fully-sequenced, single-chromosome genome. Uses the Python package [Biopython](https://biopython.org/).

## Usage
First, create a codon table using the makeCodonTable.py script.
```
makeCodonTable.py <email address> <NCBI nucleotide accession>
```
Next, use the codon table from the previous script to do codon optimization on your protein using the command below.
```
codonOptimizer.py <email address> <output file from makeCodonTable.py>
```
Then, enter the amino acid sequence for your protein of interest when prompted.

Here is a usage example below.
```
makeCodonTable.py foo@bar.edu NC_002570
codonOptimizer.py foo@bar.edu codonFrequencies_NC_002570.tsv
#Prompt // "Enter amino acid sequence: "
PRTEINSEQUENCE
```
