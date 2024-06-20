# @author       Chrishan Fernando
# @date         Jun 18, 2024
# @description  Python script to make a codon frequency
#               table for a given NCBI nucleotide accession

from Bio import Entrez, SeqIO
import sys
Entrez.email = sys.argv[1]
accession = sys.argv[2]

def main():

    #get Genbank data via Entrez
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gbwithparts", retmode="text")
    record = SeqIO.read(handle, "gb")
    handle.close()

    totalCodons = 0 #keep track of total codons
    codonDict = {}

    for index in range(len(record.features)): #scan through all features
        if record.features[index].type == "CDS": #stop on protein coding genes
            geneSequence = getGeneSequence(record, index - 1) #get the nucleotide sequence of CDS
            for codonPosition in range(int(len(geneSequence) / 3)): #scan through all codons
                #add codon to codon counts dictionary
                totalCodons += 1
                codon = geneSequence[codonPosition:codonPosition+3]
                if codon in codonDict:
                    codonDict[codon] += 1
                else: 
                    codonDict[codon] = 1
    
    #convert codon dictionary from counts to frequencies
    for element in codonDict:
        codonDict[element] = codonDict[element] / totalCodons

    #print codon dictionary
    with open("codonFrequencies_" + accession + ".tsv", "a") as oFile:
        for element in codonDict:
            oFile.write(element + "\t" + str(codonDict[element]) + "\n")

def getGeneSequence(bioRecord, position):
    # input:   (1) Biopython 'record' object
    #          (2) position of gene 'feature' in record object
    # output:  string containing nucleotide sequence of gene
    if bioRecord.features[position].strand == 1:
        return str(bioRecord.seq[int(bioRecord.features[position].location.start):int(bioRecord.features[position].location.end)])
    else:
        return str(bioRecord.seq[int(bioRecord.features[position].location.start):int(bioRecord.features[position].location.end)].reverse_complement())

main()
