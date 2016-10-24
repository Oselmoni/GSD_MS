__author__ = 'oliverselmoni'
## this scripts reformat fastq containing reads (sequences must be on one line only!):
## 1) it removes the first ten bases from the read sequences
## 2) removes spaces in the header that and add /1 and /2 to the end of the read (all modifications required by Trinity)

import re
import sys


def ref(fastq):

    a = open(fastq, "r").readlines()

    reg = "(.*).fastq"

    filename = re.findall(reg, str(fastq))[0] #gets filename for output

    out = []

    c=0


    print "converting..." # each read has 4 lines sequeceheader - sequence - qualityheader - quality

    for i in a:

        if (c==0): # sequence header

            reg = "(@.*) (.*)"

            header=re.findall(reg, i)[0][0]

            pair=re.findall(reg, i)[0][1][0]

            out+=[header+"/"+pair+"\n"] # reformat sequence header

            c+=1

        elif c==2: # quality header

            out+=["+"+header[1:]+"/"+pair+"\n"] #reformat the quality header

            c+=1

        else: # sequence or quality

            out+=[i[10:]] # reformat the sequence or quality

            c+=1

        if c==4: c=0 #restart cycle

    print "saving..."

    o=open(filename+"_v2.fastq", "w" )
    o.writelines(out)
    o.close()

    print "done!"

ref(sys.argv[1])
