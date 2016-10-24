__author__ = 'oliverselmoni'

# This script converts the output of Kallisto (estimated read count per transcript)
# in a per gene estimated read count. The estimated count per gene is calculated by
# sum of the estimated count of each of the isoforms of the gene.
# -> ***.counts -> table of estimated counts per gene for sample ***
# Input: the abundance.txt output of kallisto.

import sys
import re

def t2g(abundancetxt):

    a = open(abundancetxt, "r").readlines()

    c=0

    gendic = {}

    for i in a:

        if c == 0: c=1 # skip header

        else:

            r= '(TR\d+\|c\d+_g\d+)_i\d+'

            geneID = re.findall(r,i.rstrip())[0]
            counts=i.split("\t")[3]

            try: gendic[geneID] # already in the dictionary?

            except: # if not

                gendic[geneID]=float(counts) # add entry. GENEID as key, values are number of estimated counts
            else: # if yes

                gendic[geneID]+=float(counts) # make sum of effective counts

    outcounts=[]

    for i in gendic.keys():

        outcounts+=[i+"\t"+str(float(gendic[i][0]))+"\n"]

    o=open("sample.count", "w")
    o.writelines(outcounts)
    o.close()



t2g(sys.argv[1])


