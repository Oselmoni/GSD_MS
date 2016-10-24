# This scripts produces statistics of the Trinity.fasta output of the Trinity de-novo assembler.

__author__ = 'oliverselmoni'

import re
import sys

def ts(fasta):

    a = open(fasta, "r").readlines()

    ids=[]

    for i in a:

         if i[0] == ">":

             ids+=[i] # all transcripts IDs stored in ids

    nbtransc = len(ids) # nb of transcripts

    geneid = {} # will sort transcripts by gene

    del a

    for i in ids:

        r = ">(TR\d+\|c\d+_g\d+)_i\d+"

        gene = re.findall(r, i)[0]

        try: geneid[gene]

        except: geneid[gene] = 0 # if the gene is not in the dictionary yet, it is added to it

        else:   geneid[gene]+=1 # otherwise we increase the count of the nb. of transcripts per gene.


    print "total number of transcripts: "+str(nbtransc)
    print "total number of genes: "+str(len(geneid.keys()))
    print ""

    counts = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

    for i in geneid.keys(): # the count vectors keeps track of the number of isoform per gene.

        if geneid[i] > 14:

            counts[15]+=1

        else:

            counts[geneid[i]]+=1

    for i in range(len(counts)-1):

        print "Number of genes with "+str(i+1)+" isoform: "+str(counts[i])


    print "Number of genes with more than 15 isoforms: "+str(counts[15])


ts(sys.argv[1])
