import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
from Bio.Alphabet import DNAAlphabet
import numpy as np
import networkx as nx
from networkx import algorithms

def findCliques():
    #brute force clique finding
    #check each combination of k edges and if there are only k different nodes then its a clique
    fl = open("edgeFile.txt", 'rb')
    g = nx.read_edgelist(fl)
    cliques = []

    for clq in algorithms.find_cliques(g):
        if len(clq) == 3:
            cliques.append(clq)

    return cliques

def consensuses(clique_list):
    conSeqs = []
    instances = []
    for ls in clique_list:
        for node in ls:
            mseq = Seq(node.upper())
            instances.append(mseq) 

        m = motifs.create(instances)
        instances.clear()
        c = m.consensus
        conSeqs.append(c.lower())

    return conSeqs

def exact_appearance(strings, m, L, k): #array of k strings, length of each string, length of motif
    #make sure L<m
    #output motif of length L
    sequence1 = strings[0]
    y = L
    locs = []

    for x in range(m-L+1): #for each substring in s1
        substring = sequence1[x:y]
        for j in range(1,k):
            var = strings[j].find(substring)
            if var >= 0:
                locs.append(var)
            else:
                locs.clear()
                break
    
            if(len(locs) == k):
                return substring
        y+=1

    return ("No exact appearance of a common motif.")

def ldFormulation(strings, m, L, k, d):
    #construct graph from dictionary and lists
    kpDict = {}
    edgeList = []
    eFile = open("edgeFile.txt", "w+")
    for i in range(0,k):
        y = L
        seq = strings[i]
        subs = []
        for x in range(m-L+1):
            substring = seq[x:y]
            subs.append(substring)
            y+=1
        tuple(subs)
        kpDict[i] = subs

    #add edges
    for ls in range(0, k-1):
        x = kpDict.get(ls)
        for i in range(ls+1,k): # for each tuple in kpdict
            comp = kpDict.get(i)
            for xi in x: #for each string in tuple(ls)
                for s in comp: #for each string in tuple(i)
                    mismatch = 0
                    p1 = 0
                    for j in s: #for each letter in tuple(i) string
                        if(xi[p1] != j):
                            mismatch+=1
                        p1+=1
                    if(mismatch <= d):
                        eFile.write(xi)
                        eFile.write(" ")
                        eFile.write(s)
                        eFile.write("\n")
                        edgeList.append((xi, s)) #will be a list of 2 points
    
    eFile.close()

    clique_list = findCliques()
    if clique_list == []:
        return ("No common motif with allowed number of mismatches.")
    else:
        conSeqs = consensuses(clique_list)
        return conSeqs


print("------------------------------------------------------")
print("\nTitle: Bioinformatics Project I - Motif Finding")
print("Author: Rachael Hawthorne")
print("Date: 2/16/2020\n\n")
k = 0
m = 0
input_file = input("Enter the sequence file name: ")
L = int(input("Enter the length of the desired common motif: "))
d = int(input("Enter the allowed number of mismatches: "))
seqs = []

for record in SeqIO.parse(input_file, "fasta"):
    seqs.append(str(record.seq))
    m = len(record)
    k+=1

exact_motif = exact_appearance(seqs, m, L, k)
mismatch_motif = ldFormulation(seqs, m, L, k, d)
print("\n\nThe exact appearance motif of the sequences is: ", exact_motif)
print("\nThe consensus motif(s) of length", L, "with up to", d, "mismatch(es): ", mismatch_motif)
print("\n\n------------------------------------------------------")
