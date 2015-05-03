#!/usr/bin/env python2
""" make phylip or nexus formatted data outputs from .loci file """

import numpy as np

def make(params, names, longname, formats):
    """ make outputs """

    ## make sure name is a list and sort it
    names = list(names)
    names.sort()
    
    ## read in the loci file "
    finalfile = open(params["work"]+"outfiles/"+\
                     params["outname"]+".loci").read() 

    ## dict for saving the full matrix
    fdict = {name:[] for name in names}

    ## uncomment and use this if you want to save block 
    ## information for partitioning loci
    ## blocked = 0

    ## remove empty column sites from matrix
    ## and append edited seqs fdict
    for loc in finalfile.split("//")[:-1]:
        anames = [i.split()[0][1:] for i in \
                  loc.strip().split("\n") if ">" in i]

        arrayed = np.array([tuple(i.split()[-1]) for i in \
                               loc.strip().split("\n") if ">" in i])

        ## create mask for columns that are empty or 
        ## that are paired-end separators (compatible w/ pyrad v2 and v3)
        mask = [i for i in range(len(arrayed.T)) if not \
                  np.all([j in list("N-nxX") for j in arrayed.T[i]])]

        ## uncomment to print block info (used to partition by locus)
        #blockend += minray
        #print blockend,

        ## append data to dict
        for name in names:
            if name in anames:
                fdict[name] += "".join(arrayed[anames.index(name), mask])
            else:
                fdict[name] += "".join(["N"]*len(arrayed[0, mask]))

    #############################
    ## print out .PHY file by default
    superout = open(params["work"]+"outfiles/"+\
                    params["outname"]+".phy", 'w')
    print >>superout, len(fdict), len("".join(fdict[names[0]]))
    for name in names:
        print >>superout, name+(" "*((longname+3)-\
                          len(name)))+"".join(fdict[name])
    superout.close()

    #############################
    ## print out interleaved .NEX file
    if 'n' in formats:
        data = open(params["work"]+"outfiles/"+\
                      params["outname"]+".phy").readlines() 
        nexout = open(params["work"]+"outfiles/"+\
                      params["outname"]+".nex", 'w')

        ntax, nchar = data[0].strip().split(" ")

        print >>nexout, "#NEXUS"
        print >>nexout, "BEGIN DATA;"
        print >>nexout, "  DIMENSIONS NTAX=%s NCHAR=%s;" % (ntax, nchar)
        print >>nexout, "  FORMAT DATATYPE=DNA MISSING=N GAP=- INTERLEAVE=YES;"
        print >>nexout, "  MATRIX"

        idict = {}
        for line in data[1:]:
            a = line.lstrip().rstrip().split(" ")
            idict[a[0]] = a[-1]

        nameorder = idict.keys()
        nameorder.sort()

        n = 0
        sz = 100
        while n < len(a[-1]):
            for tax in nameorder:
                print >>nexout, "  "+tax+" "*\
                                ((longname-len(tax))+3)+\
                                idict[tax][n:n+sz]
            n += sz
            print >>nexout, ""
        print >>nexout, ';'
        print >>nexout, 'END;'
        nexout.close()
        

if __name__ == "__main__":
    make(params, names, longname, format)
