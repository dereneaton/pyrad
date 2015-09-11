#!/usr/bin/env python2

import numpy as np
import sys
import os
import glob


def update(idict, count, WORK, outname):
    """ updates dictionary with the next .5M reads 
    from the super long string phylip file """

    data = iter(open(WORK+"outfiles/"+outname+".phy"))
    ntax, nchar = data.next().strip().split()

    ## read in max N bp at a time                                                                            
    for line in data:
        tax, seq = line.strip().split()
        idict[tax] = idict[tax][100000:]
        idict[tax] += seq[count:count+100000]
    del line

    return idict
    


def makephy(WORK, outname, names, longname, formats):
    " order names "
    names = list(names)
    names.sort()
    
    " read in loci file "
    locifile = open(WORK+"outfiles/"+outname+".loci", 'rb')

    " dict for saving the full matrix "
    fdict = {name:[] for name in names}

    " remove empty column sites and append edited seqs to dict F "
    locus = iter(locifile)
    done = 0
    nloci = 0
    nbases = 0
    while nloci < 50000: #not done:
        seqs = []
        #arrayed = np.array([])
        anames = []
        while 1:
            ## get next locus
            try:
                samp = locus.next()
            except StopIteration:
                done = 1
                break
            if "//" in samp:
                nloci += 1
                #print arrayed
                break
            else:
                name, seq = samp.split()
                anames.append(name[1:])
                seqs.append(seq.strip())
        ## reset
        arrayed = np.array([list(i) for i in seqs])
        if done:
            break
        ## create mask for columns that are empty or 
        ## that are paired-end separators (compatible w/ pyrad v2 and v3)
        #mask = [i for i in range(len(arrayed.T)) if np.any([
        mask = [i for i in arrayed.T if any([j not in list("-Nn") for j in i])]
        masked = np.dstack(mask)[0]
        #not \
        #        np.all([j in list("N-n") for j in arrayed.T[i]])]

        ## uncomment to print block info (used to partition by locus)
        #blockend += minray
        #print blockend,
        #print loc
        #print arrayed

        ## append data to dict
        for name in names:
            if name in anames:
                #fdict[name].append(arrayed[anames.index(name), mask].tostring())
                fdict[name].append(masked[anames.index(name),:].tostring())
            else:
                fdict[name].append("N"*masked.shape[1])
                #fdict[name].append("N"*len(arrayed[0, mask]))
        ## add len to total length
        nbases += len(fdict[name][-1])

        ## after x iterations tmp pickle fdict?
        if not nloci % 1e4:
            ## concat strings
            for name in fdict:
                with open(os.path.join(WORK, "tmp", 
                    "{}_{}.tmp".format(name, nloci)), 'wb') as wout:
                    wout.write("".join(fdict[name]))
            del fdict
            fdict = {name:[] for name in names}

    ## print out .PHY file, if really big, pull form multiple tmp pickle
    superout = open(WORK+"outfiles/"+outname+".phy", 'wb')
    print >>superout, len(names), nbases
    if nloci < 1e4:
        for name in names:
            print >>superout, name+(" "*((longname+3)-\
                              len(name)))+"".join(fdict[name])
    else:
        for name in names:
            superout.write("{}{}{}".format(
                            name,
                            " "*((longname+3)-len(name)),
                            "".join(fdict[name])))
            tmpfiles = glob.glob(os.path.join(WORK, "tmp", name+"*.tmp"))
            tmpfiles.sort()
            for tmpf in tmpfiles:
                with open(tmpf, 'rb') as tmpin:
                    superout.write(tmpin.read())
            superout.write("\n")
    superout.close()



def makenex(WORK, outname, names, longname, formats):
    """ PRINT NEXUS """
    if 1:
        " make nexus output "
        data   = iter(open(WORK+"outfiles/"+outname+".phy"))
        nexout = open(WORK+"outfiles/"+outname+".nex", 'wb')

        ntax, nchar = data.next().strip().split(" ")

        print >>nexout, "#NEXUS"
        print >>nexout, "BEGIN DATA;"
        print >>nexout, "  DIMENSIONS NTAX=%s NCHAR=%s;" % (ntax,nchar)
        print >>nexout, "  FORMAT DATATYPE=DNA MISSING=N GAP=- INTERLEAVE=YES;"
        print >>nexout, "  MATRIX"

        idict = {}

        ## read in max 1M bp at a time
        for line in data:
            tax, seq = line.strip().split()
            idict[tax] = seq[0:100000]
        del line

        nameorder = idict.keys()
        nameorder.sort()

        n=0
        tempn=0
        sz = 100
        while n < len(seq):
            for tax in nameorder:
                print >>nexout, "  "+tax+" "*\
                                 ((longname-len(tax))+3)+\
                                 idict[tax][tempn:tempn+sz]
            n += sz
            tempn += sz
            print >>nexout, ""

            if not n % 100000:
                #print idict[tax][tempn:tempn+sz]
                idict = update(idict, n, WORK, outname)
                tempn -= 100000
            
        print >>nexout, ';'
        print >>nexout, 'END;'
        nexout.close()
        

def make(WORK, outfile, names, longname, formats):
    makephy(WORK, outfile, names, longname, formats)
    makenex(WORK, outfile, names, longname, formats)
    

if __name__ == "__main__":
    make()
