#!/usr/bin/env python2
""" make phylip or nexus formatted data outputs from .loci file """

import numpy as np

def make(params, names, longname, formats):
    """ make outputs """

    ## make sure name is a list and sort it
    names = list(names)
    names.sort()
    
    ## read in the loci file " ## use iterwhile to select loci...
    #finalfile = open(params["work"]+"outfiles/"+\
    #                 params["outname"]+".loci").read() 
    locifile = open(params["work"]+"outfiles/"+\
                    params["outname"]+".loci", 'rb') 

    ## dict for saving the full matrix (too much memory!)
    fdict = {name:[] for name in names}

    ## uncomment and use this if you want to save block 
    ## information for partitioning loci
    ## blocked = 0

    ## remove empty column sites from matrix
    ## and append edited seqs fdict
    locus = iter(locifile)
    tnames = []
    seqs = []
    done = 0
    while not done:
        arrayed = np.array([])
        anames = []
        while 1:
            ## get next locus
            try:
                samp = locus.next()
            except StopIteration:
                done = 1
                break
            if "//" in samp:
                #print arrayed
                break
            else:
                name, seq = samp.split()
                anames.append(name[1:])
                seqs.append(seq.strip())
        ## reset
        arrayed = np.array([list(i) for i in seqs])
        seqs = []
        if done:
            break
        ## create mask for columns that are empty or 
        ## that are paired-end separators (compatible w/ pyrad v2 and v3)
        mask = [i for i in range(len(arrayed.T)) if not \
                np.all([j in list("N-nxX") for j in arrayed.T[i]])]
        #masked = arrayed.T[(arrayed.T != "-") & \
        #                   (arrayed.T != "N") & \
        #                   (arrayed.Tq != "X")]

        ## uncomment to print block info (used to partition by locus)
        #blockend += minray
        #print blockend,
        #print loc
        #print arrayed

        ## append data to dict
        for name in names:
            print name, anames
            if name in anames:
                fdict[name].append(arrayed[anames.index(name), mask].tostring())
            else:
                fdict[name].append("N"*len(arrayed[0, mask]))

    #############################
    ## print out .PHY file by default
    superout = open(params["work"]+"outfiles/"+\
                    params["outname"]+".phy", 'wb')
    print fdict[names[0]]
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
