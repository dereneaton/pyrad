#!/usr/bin/env python2

""" script to convert .loci formatted files to .migrate format """

import sys    
from collections import OrderedDict


def make(params, taxadict, minhits, maxnumberloci):
    """ make migrate output file """
    ## outfile
    outfile = open(params["work"]+"/outfiles/"+\
                   params["outname"]+".migrate", 'w')

    ## cleanup taxadict
    taxa = OrderedDict()
    for group in taxadict:
        taxa[group] = []
        for samp in taxadict[group]:
            sampname = samp.split("/")[-1].replace(".consens.gz", "")
            taxa[group].append(sampname)

    ## print to screen the minimums for each group
    print "\t    data set reduced for group coverage minimums"
    for i, j in zip(taxa, minhits):
        print "\t   ", i, taxa[i], "minimum=", j
        if len(taxa[i]) < int(j):
            sys.exit("\t\t**Error: number of samples in group "+i+" is \n"+
                     "\t\tless than the minimum for a locus to be retained\n")

    ## filter data to only the loci that have data
    ## for at least N individuals in each pop
    keep = []
    minzip = zip(taxa.keys(), minhits)

    ## read in data to sample names
    locilines = iter(open(params["work"]+"/outfiles/"+\
                          params["outname"]+".loci", 'rb'))
                     #.\read().strip().split("|")[:-1]

    ## filter loci for mincoverage
    for loc in loci:
        samps = [i.split()[0].replace(">", "") \
                 for i in loc.split("\n") if ">" in i]
        cfilter = []
        for group, mins in minzip:
            cfilter.append(sum([i in samps for i in taxa[group]]) >= int(mins))
        if all(cfilter):
            keep.append(loc)

    if len(keep) == 0:
        sys.exit("\tno loci were found to meet the group assignment criteria")

    ## print data to file
    print >>outfile, str(len(taxa))+" "+str(min(len(keep), maxnumberloci))+\
                     " (npops nloci for data set "+\
                      params["outname"]+".loci)"
    
    ## print all data for each population at a time
    done = 0
    for group in taxadict:
        ## print a list of lengths of each locus
        if not done:
            loclens = [len(loc.split("\n")[1].split()[-1]\
                       .replace("x", "n").replace("n", "")) \
                       for loc in keep]
            print >>outfile, " ".join([str(i) for i in \
                            loclens[:min(len(keep), maxnumberloci)]])
            done += 1

        ## print a list of number of individuals in each locus
        indslist = []
        for loc in keep:
            samps = [i.split()[0].replace(">", "") \
                     for i in loc.split("\n") if ">" in i]
            inds = sum([i in samps for i in taxa[group]])
            indslist.append(inds)
        print >>outfile, " ".join([str(i) for i in \
                         indslist[:min(len(keep), maxnumberloci)]])+" "+group

        ## print sample id, spaces, and sequence data
        #for loc in range(len(keep)):
        for loc in range(min(len(keep), maxnumberloci)):
            seqs = [i.split()[-1] for i in keep[loc].split("\n") if \
                    i.split()[0].replace(">", "") in taxa[group]]
            for i in range(len(seqs)):
                print >>outfile, group[0:7]+"_"+str(i)+\
                      (" "*(10-len(group[0:7]+"_"+str(i))))+\
                      seqs[i].replace("x", "n").replace("n", "")
            
    outfile.close()



if __name__ == "__main__":

    # WORK = "/home/deren/Dropbox/Public/PyRAD_TUTORIALS/tutorial_RAD"
    # outname = "c85m4p3"
    # pops = ['pop1','pop2','pop3']
    # samps = [ ["1A0","1B0","1C0","1D0"],
    #           ["2E0","2F0","2G0","2H0"],
    #           ["3I0","3J0","3K0","3L0"] ]

    # taxadict = OrderedDict(zip(pops,samps))
    # minhits = [4,4,4]
    # seed = 112233
    make(params, taxadict, minhits, seed, maxnumberloci)
