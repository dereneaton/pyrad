#!/usr/bin/env python2

""" create .snp and related outputs """

import numpy as np
from itertools import chain
import alignable

# pylint: disable=E1101


def is_biallelic(arr):
    """ return only columns of snp array that are bi-allelic"""
    cols = []
    for index, col in enumerate(arr.T):
        alls = set(chain(*[alignable.unstruct(i) for i in col if i.strip()]))
        if len(alls) <= 2:
            cols.append(index)
    return arr[:, cols]



def make_dicts(params, names):
    """ make function """
    #WORK, outname, names, formats, seed, ploidy):
    np.random.seed(int(params["seed"]))
    locifile = iter(open(params["work"]+"outfiles/"+\
                         params["outname"]+".loci"))

    ## output .snps and .unlinked_snps"
    snpdict = {name:[] for name in names}   ## snp dict  S
    usnpdict = {name:[] for name in names}  ## unlinked snp dict Si
    counter = {'loci': 0, 
               'snps': 0, 
               'usnps': 0,
               'biusnps': 0}

    ## for each locus select out the SNPs"
    locus = iter(locifile)
    done = start = 0
    while not done:
        seqs = []           
        snpsarray = np.array([])
        biarray = np.array([])        
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
                seq = samp[start:]
                anames.append("SNPS")
                seqs.append(seq.split("|")[0])
                counter["loci"] += 1
                break
            else:
                name, seq = samp.split()
                anames.append(name[1:])
                seqs.append(seq.strip())
                start = samp.rindex(" ")+1                

        ## do this locus
        arrayed = np.array([list(i) for i in seqs])
        if arrayed.size:
            mask = [i for i, j in enumerate(arrayed[-1]) if j.strip()]
            snpsarray = arrayed[:, mask]

        ## store snps and usnps to dicts
        if snpsarray.size:
            ## get bi-allelic array 
            biarray = is_biallelic(snpsarray[:-1,])
            ## get size of snps
            snplen = len(snpsarray[0])
            counter['snps'] += snplen
            counter['usnps'] += 1
            ## fill in dict
            for tax in snpdict:
                if tax in anames:
                    snpstring = snpsarray[anames.index(tax)].tostring()
                    snpdict[tax].append(snpstring)
                else:
                    snpdict[tax].append("N"*snplen)
            ## do biallelic sampling
            if biarray.size:
                rando = np.random.randint(biarray.shape[1])
                counter['biusnps'] += 1
                for tax in snpdict:
                    if tax in anames:
                        usnp = biarray[anames.index(tax), rando]
                        usnpdict[tax].append(usnp)
                    else:
                        usnpdict[tax].append("N")
            else:
                ## do something for no bisnps
                #for tax in snpdict:
                #    usnpdict[tax].append("N")                    
                pass
        else:
            ## do something for no snps
            for tax in snpdict:
                ## don't save last empty string
                if not done:
                    snpdict[tax].append("N")
    ## remove last element of snpdict
    return counter, snpdict, usnpdict



def make_snpfile(params, counter, snpdict, longname):
    """ prints snpdict to file """
    ## names
    snames = list(snpdict.keys())
    snames.sort()

    ## outfile
    snpsout = open(params["work"]+'outfiles/'+\
                   params["outname"]+".snps", 'wb')

    ## print to file
    snpsout.write("## {} taxa, {} loci, {} snps\n".format(
                                            len(snpdict),
                                            counter["loci"],
                                            counter["snps"]))
    for tax in snames:
        snpsout.write("{}{}{}\n".format(
            tax,
            " "*(longname-len(tax)+3),
            " ".join(snpdict[tax])
            ))
    snpsout.close()



def make_usnpfile(params, counter, usnpdict, longname):
    """ prints usnpdict to file """
    ## name order
    snames = list(usnpdict.keys())
    snames.sort()

    ## outfile
    usnpout = open(params["work"]+'outfiles/'+\
                   params["outname"]+".unlinked_snps", 'wb')

    ## print to file
    usnpout.write("{} {}\n".format(len(usnpdict),
                                   counter['usnps']))
    for tax in snames:
        usnpout.write("{}{}{}\n".format(
                tax,
                " "*(longname-len(tax)+3),
                "".join(usnpdict[tax]))
                )
    usnpout.close()



def make_structfile(params, usnpdict, longname):
    """ print structure to file """
    ## name order
    snames = list(usnpdict.keys())
    snames.sort()

    ## outfile
    structout = open(params["work"]+'outfiles/'+\
                  params["outname"]+".str", 'wb')

    stdict = {'A': '0',
              'T': '1',
              'G': '2',
              'C': '3',
              'N': '-9',
              '-': '-9'}

    ## allow haploid output too
    if params["haplos"] > 1:
        for sname in snames:
            #print stdict
            structout.write("{}{}{}\n".format(
                    sname,
                    " "*(longname-len(sname)+3),
                    "\t"*6+"\t".join([stdict[alignable.unstruct(j)[0]] \
                                      for j in usnpdict[sname]])
                    ))
            structout.write("{}{}{}\n".format(
                    sname,
                    " "*(longname-len(sname)+3),
                    "\t"*6+"\t".join([stdict[alignable.unstruct(j)[1]] \
                                      for j in usnpdict[sname]])
                    ))
    else:   
        for sname in snames:
            structout.write("{}{}{}\n".format(
                    sname,
                    " "*(longname-len(sname)+3),
                    "\t"*6+"\t".join([stdict[alignable.unstruct(j)[1]] \
                                      for j in usnpdict[sname]])
                    ))
    structout.close()



def make_genofile(params, snpdict, usnpdict):
    """ make geno output """

    ## name order
    snames = list(usnpdict.keys())
    snames.sort()

    ## outfile for unlinked snps
    genoout = open(params["work"]+'outfiles/'+\
                   params["outname"]+".usnps.geno", 'wb')

    ## print to file
    for snp in range(len(usnpdict.values()[0])):
        getref = 0
        ref = "N"
        while ref == "N":  
            ref = alignable.unstruct(usnpdict[snames[getref]][snp])[0]
            getref += 1

        ## count how many times each taxon has the reference allele
        snprow = []
        for tax in snames:        
            counts = alignable.unstruct(usnpdict[tax][snp]).count(ref)
            if usnpdict[tax][snp] != "N":
                snprow.append(str(2-counts))
            else:
                snprow.append("9")

        ## lines should not be invariable
        if len(set(snprow)) > 1:
            #print >>genoout, snprow 
            genoout.write("{}\n".format("".join(snprow)))
    genoout.close()


    ## outfile for full snps
    genoout = open(params["work"]+'outfiles/'+\
                   params["outname"]+".snps.geno", 'wb')

    for snpstring in range(len(snpdict.values()[0])):
        ## if snpstring isn't empty (Ns)
        if snpdict[snames[0]][snpstring].strip("N").strip():
            for snp in range(len(snpdict[snames[0]][snpstring]\
                                                .strip("N").strip())):
                getref = 0
                ref = "N"
                while ref == "N":
                    base = snpdict[snames[getref]][snpstring][snp]
                    ref = alignable.unstruct(base)[0]
                    getref += 1

                ## count how many times ref is in taxon
                snprow = []
                for tax in snames:
                    counts = alignable.unstruct(
                        snpdict[tax][snpstring][snp]).count(ref)
                    if snpdict[tax][snpstring][snp] != "N":
                        snprow.append(str(2-counts))  
                    else:
                        snprow.append("9")

                if len(set(snprow)) > 1:   
                    genoout.write("{}\n".format("".join(snprow)))
    genoout.close()



def main(params, names):
    """ all functions """
    ## formatting name lengths
    longname = max([len(i) for i in names])

    ## get snps from loci file
    counter, snpdict, usnpdict = make_dicts(params, names)

    ## print output files
    if "s" in params["outform"]:
        make_snpfile(params, counter, snpdict, longname)
    if "u" in params["outform"]:
        make_usnpfile(params, counter, usnpdict, longname)
    if 'k' in params["outform"]:        
        make_structfile(params, usnpdict, longname)        
    if "g" in params["outform"]:
        make_genofile(params, snpdict, usnpdict)

    return counter


if __name__ == "__main__":
    main()


    # statsout = open(params["work"]+"stats/"+\
    #                 params["outname"]+".stats", 'a')
    # print >>statsout, "sampled unlinked SNPs=", len(usnpdict.values()[0])
    # print >>statsout, "sampled unlinked bi-allelic SNPs= "+\
    #                    str(len(usnpdict.values()[0])-bis)
    # statsout.close()
