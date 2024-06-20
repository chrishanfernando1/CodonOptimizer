# @author       Chrishan Fernando
# @date         Jun 18, 2024
# @description  Python script to codon optimize a sequence
#               based on a given codon frequency table

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import sys
Entrez.email = sys.argv[1]
codonFrequencyFile = sys.argv[2]
sequence = input("Enter amino acid sequence: ")

def main():

    codonFrequencyDict = {} #dictionary of codon frequencies
    bestCodonDict = {}      #map of amino acid to best codon
    stopCodonDict = {}      #dictionary with stop codon frequencies

    #read in codon frequency dictionary
    with open(codonFrequencyFile, "r") as inFile:
        for line in inFile:
            lineList = line.strip().split("\t")
            codonFrequencyDict[lineList[0]] = float(lineList[1])
    
    #pick best codons
    for codon in codonFrequencyDict:
        aminoAcid = str(Seq(codon).translate())
        if aminoAcid == "*": #if stop codon
            stopCodonDict[codon] = codonFrequencyDict[codon]
        else:                #if not stop codon
            if aminoAcid not in bestCodonDict: #if there isn't a codon assigned for the amino acid in question
                bestCodonDict[aminoAcid] = codon
            elif codonFrequencyDict[codon] > codonFrequencyDict[bestCodonDict[aminoAcid]]: #if the codon being evaluated is better than the current best
                bestCodonDict[aminoAcid] = codon

    #pick best stop codon
    bestStop = ""
    for codon in stopCodonDict:
        if bestStop == "":
            bestStop = codon
        elif codonFrequencyDict[bestStop] < codonFrequencyDict[codon]:
            bestStop = codon

    #perform codon optimization
    optimizedSequence = ""
    for character in sequence:
        optimizedSequence += bestCodonDict[character]
    optimizedSequence += bestStop #add the best stop codon
    print(optimizedSequence)

main()
