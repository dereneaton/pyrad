#!/usr/bin/env python2

""" convert .loci to .treemix output """

import numpy as np
import gzip
import sys
try:
    from collections import OrderedDict, Counter
except ImportError:
    from ordereddict import OrderedDict, Counter
import alignable


def make(params, taxadict, minhits, quiet):
    """ make output files """

    ## output files
    outfile = gzip.open(params["work"]+"/outfiles/"+\
                        params["outname"]+".treemix.gz", 'w')

    ## cleanup taxadict to just sample names
    ## that match with the consens file names
    taxa = OrderedDict()
    for group in taxadict:
        taxa[group] = []
        for samp in taxadict[group]:
            consensfile = samp.split("/")[-1].replace(".consens.gz", "")
            taxa[group].append(consensfile)
    if not quiet:
        print "\t    data set reduced for group coverage minimums"        
        for i, j in zip(taxa, minhits):
            print "\t   ", i, taxa[i], 'minimum=', j
            if len(taxa[i]) < int(j):
                sys.exit("\t\t**Error: number of samples in group "+i+" is \n"+
                   "\t\tless than the minimum for a locus to be retained\n")

    ## read in data from unlinked_snps to sample names
    infile = open(params["work"].rstrip("/")+"/outfiles/"+\
                  params["outname"]+".unlinked_snps", 'r')
    dat = infile.readlines()
    nsamp, nsnps = dat[0].strip().split(" ")
    nsamp = int(nsamp)
    nsnps = int(nsnps)
    ndata = np.empty([int(nsamp), int(nsnps)], dtype='object')

    ## read SNP matrix into a numpy.array
    for line in range(len(dat[1:])):
        _, seqs = dat[1:][line].split()
        ndata[line] = list(seqs)
    locs = np.transpose(ndata)

    ## unpack ambiguity bases and find two most common alleles
    ## at every SNP site, save to a list
    alleles = []
    for site in locs:
        ds = []
        for base in site:
            if base in list("RKSYWM"):
                ds.append(alignable.unstruct(base)[0])
                ds.append(alignable.unstruct(base)[1])
            else:
                ds.append(base)
                ds.append(base)
        snp = [base for base in ds if base not in ["N", '-']]
        counts = Counter(snp).most_common(3)
        alleles.append([counts[0][0], counts[1][0]])

    ## create a dictionary mapping sample names to SNPs    
    snps = OrderedDict()
    for line in dat[1:]:
        name, seq = line.split()
        snps[name] = seq

    ## reduce Taxa dict to only samples that are in the unlinkedsnps alignment
    for key in taxa:
        replacement = []
        for val in taxa[key]:
            if val in snps.keys():
                replacement.append(val)
        taxa[key] = replacement

    ## create a dictionary with empty lists for each taxon 
    freq = OrderedDict()
    for tax in taxa:
        freq[tax] = []

    ## fill the FREQ dictionary with SNPs for all 
    ## samples in that taxon
    keeps = []
    for snp in range(int(nsnps)):
        GG = []
        ## if snp meets minhits requirement
        for tax, mins in zip(taxa, minhits):
            GG.append(sum([snps[i][snp] not in ["N", "-"] \
                      for i in taxa[tax]]) >= int(mins))
        if all(GG):
            keeps.append(snp)

    for keep in keeps:
        for tax in freq:
            bunch = []
            for i in taxa[tax]:
                bunch.append(alignable.unstruct(snps[i][keep])[0])
                bunch.append(alignable.unstruct(snps[i][keep])[1])
                #print tax, i, SNPS[i][keep], bunch
            freq[tax].append("".join(bunch))

    ## header
    print >>outfile, " ".join(freq.keys())

    ## data to file
    for i, j in enumerate(keeps):
        a1 = alleles[j][0]
        a2 = alleles[j][1]
        H = [str(freq[tax][i].count(a1))+","+\
             str(freq[tax][i].count(a2)) for tax in freq]
        HH = " ".join(H)

        ## exclude non-biallelic SNPs
        if " 0,0 " not in HH:
            ## exclude invariable sites given this sampling
            if not all([zz.split(",")[1] in '0' for zz in H]):
                print >>outfile, " ".join(H)
        #else:
        #    excludes += 1

    outfile.close()
    

if __name__ == "__main__":
    make(WORK, outname, taxadict)
