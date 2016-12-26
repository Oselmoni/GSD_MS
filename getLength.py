__author__ = 'oliverselmoni'

def rf(fasta):

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

import re
import numpy

def GL(fasta):

    seq = rf(fasta)

    lendic = {}

    for i in seq.keys():

        p = "(TR\d+\|c\d+_g\d+)_i\d+"

        geneID= re.findall(p, i)[0]

        try: lendic[geneID]

        except: lendic[geneID]=[len(seq[i])]

        else: lendic[geneID]+=[len(seq[i])]

    out=[]
    for i in lendic.keys():

        leng=numpy.median(lendic[i])

        out+=[i+"\t"+str(leng)+"\n"]


    o = open("geneLenght.txt", "w")
    o.writelines(out)
    o.close()






GL("Assembly_filtered.fasta")