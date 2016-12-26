__author__ = 'oliverselmoni'

import re

# maps blast results for isoforms (i.e. more hit per gene are possible) with the custom database corresponding
# to the sequences used for blast. The custom database contains info about protein name and GO terms.
# The output is a table for which for each gene a protein name is associated (taken from the most frequent match within the isoforms)
# and all the GO terms are listed (taken from all the GO terms associated to isoforms). This output can be fed directly into R for DGEA.

def pa(blastres, uniprotDB):

    a = open(blastres, "r").readlines()

    b = open(uniprotDB, "r").readlines()

    c = 0

    uniprotD = {}

    for i in b: #create a dictionnary of uniprot entries (key=ID, value=[protein name, [list of GOs identifier]]

        if c == 0: c=1

        pid=i.split("\t")[0]

        pname=i.split("\t")[2]

        p = "(GO:\d+)"

        pGOstring=i.split("\t")[3]

        pGOs = re.findall(p,pGOstring)

        uniprotD[pid] = [pname, pGOs]


    genedic = {}

    isoformdic = {}

    for i in a: # parse the blast results

        p = "(TR\d+\|c\d+_g\d+_i\d+)"

        isoformid= re.findall(p, i.split("\t")[0])[0] #isoform ID

        p = "(TR\d+\|c\d+_g\d+)_i\d+"

        geneid =  re.findall(p, i.split("\t")[0])[0] #gene ID

        p = "[tr|sp]\|(.+)\|.+"

        pid = re.findall(p, i.split("\t")[1])[0] # entry ID

        try: isoformdic[isoformid] # was this isoform already treated?

        except: # if not, treat it...

            isoformdic[isoformid]=True

            try: genedic[geneid] # was the gene already seen?

            except: genedic[geneid] = [pid] # if not, add entry

            else: genedic[geneid]+= [pid] # if yes, just add the uniprot ID

        else: "" # if yes: do nothing

    out=[]

    f=0

    for i in genedic.keys(): # for each gene:

        iddic = {} # create a dictionarry for count the frequency of associated uniprot entry

        for l in genedic[i]:

                try: iddic[l]
                except: iddic[l]=1
                else: iddic[l]+=1

        golist=[]
        max=0
        maxid=""

        for l in iddic.keys(): # for all of the uniprot entry associated to the gene

            if iddic[l]>max: maxid=uniprotD[l][0] # check which one is the most frequent

            for m in uniprotD[l][1]:

                if (m in golist) == False: golist+=[m] # add all the GO terms observed among the matching uniprot entries

        if len(golist) != 0: f+=1

        out+=[i+"\t"+maxid+"\t"+str(",".join(golist))+"\n"] # write the output, assign the protein name of the most uniprot entry most frequently associated to the gene


    print str(len(genedic.keys()))+" genes, "+str(f)+" annotated with at least 1 GO term"


    o=open("GOs.txt", "w")
    o.writelines(out)
    o.close()

import sys

pa(sys.argv[1], sys.argv[2])


