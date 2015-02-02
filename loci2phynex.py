#!/usr/bin/env python2

import numpy as np

## TODO: Make faster and reduce memory load...
def make(WORK, outname, names, longname, formats):
    " order names "
    names = list(names)
    names.sort()
    
    " read in loci file "
    finalfile = open(WORK+"outfiles/"+outname+".loci").read() 

    " dict for saving the full matrix "
    F = {}
    for name in names:
        F[name] = []

    " remove empty column sites and append edited seqs to dict F "
    for loc in [i for i in finalfile.split("|")[:-1]]:
        anames = [i.split(" ")[0][1:] for i in loc.strip().split("\n")[:-1]]
        array = np.array([tuple(i.split(" ")[-1]) for i in loc.strip().split("\n")][:-1])

        ## which columns are empty
        emps = [i for i in range(len(array.T)) if \
                np.all([j in ['N','-'] for j in array.T[i]])]

        ## delete those columns
        narray = np.delete(array, emps, 1)
        minray = len("".join(narray[0]).replace("x","n").replace("n",""))

        ## append data to dict
        for name in names:
            if name in anames:
                F[name] += "".join(narray[anames.index(name)]).replace("x","n").replace("n","")
            else:
                F[name] += "".join(["N"]*minray)  

    " print out .PHY file "
    superout = open(WORK+"outfiles/"+outname+".phy",'w')
    print >>superout, len(F), len("".join(F[names[0]]).replace("x",""))

    for name in names:
        print >>superout, name+(" "*((longname+3)-len(name)))+"".join(F[name])
    superout.close()


    if 'n' in formats:
        " make nexus output "
        data   = open(WORK+"outfiles/"+outname+".phy").readlines() 
        nexout = open(WORK+"outfiles/"+outname+".nex", 'w')

        ntax,nchar = data[0].strip().split(" ")

        print >>nexout, "#NEXUS"
        print >>nexout, "BEGIN DATA;"
        print >>nexout, "  DIMENSIONS NTAX=%s NCHAR=%s;" % (ntax,nchar)
        print >>nexout, "  FORMAT DATATYPE=DNA MISSING=N GAP=- INTERLEAVE=YES;"
        print >>nexout, "  MATRIX"

        L = {}
        for line in data[1:]:
            a = line.lstrip().rstrip().split(" ")
            L[a[0]] = a[-1]

        n=0
        sz = 100
        while n<len(a[-1]):
            for tax in L:
                print >>nexout, "  "+tax+" "*((longname-len(tax))+3)+L[tax][n:n+sz]
            n += sz
            print >>nexout, ""
        print >>nexout, ';'
        print >>nexout, 'END;'
        nexout.close()
        

if __name__ == "__main__":
    make(WORK, outfile, names, longname, formats)
