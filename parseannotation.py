__author__ = 'oliverselmoni'

# Given the assembly sequences and the blast hit report of these sequences, this script filters the transcripts keeping only the ones having an annotation

import re
import sys

def rf(fasta): # elaborate a fasta into a dictionary (Key=identifier, value=sequence)

    a = open(fasta, "r").readlines()

    IDpos = []
    counter= -1

    for x in a:

        if x[0] == ">":
            counter+=1
            IDpos+= [counter]
        else:
            counter+=1

    L=range(0, len(IDpos)-1)

    dict = {}

    for x in L:

        cur = IDpos[x]
        nex = IDpos[x+1]

        id = a[cur].split(" ")[0][1:]
        seq= ""

        while cur < (nex-1):

            cur+=1
            seq+=a[cur][:-1]

        dict[id]=seq


    last = len(IDpos)
    lastID = IDpos[last-1]

    id=a[lastID].split(" ")[0][1:]

    L=range(lastID+1, len(a))
    seq=""

    for x in L:

        seq+=a[x][:-1]

    dict[id]=seq

    return dict

def pa(annotation, assembly):

    seqD = rf(assembly) # dictionary of sequences

    ann = open(annotation, "r").readlines() #annotation from Swissprot

    isoD = {}

    outS=[]

    for i in ann: # for each line of the annotation file

        isoID = i.split("\t")[0] # ID of the isoform

        try: isoD[isoID] # isoform already in dictionary? (note: one isoform can appear on more line if the alignment has multiple segments)

        except: # if not, add it...

            isoD[isoID] = ""

            outS+=[">"+isoID+"\n"+seqD[isoID]+"\n"] # add isoform with annotation to new fasta


        else: "" # if yes, do nothing


    o=open("Assembly_filtered.fasta", "w")
    o.writelines(outS)
    o.close()


    print "File Assembly_filtered.fasta contains "+str(len(isoD.keys()))+" transcripts (each of them has one hit vs SP)"



pa(sys.argv[1], sys.argv[2])



