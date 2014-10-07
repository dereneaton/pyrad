
import numpy as np
import sys
import gzip
try:
    from collections import Counter
except:
    from counter import Counter
import alignable

def make(WORK, outname, taxadict):

    ## cleanup taxadict
    taxa = {}
    for group in taxadict:
        taxa[group] = []
        for samp in taxadict[group]:
            a = samp.split("/")[-1].replace(".consens.gz","")
            taxa[group].append(a)
    
    ## read in data to sample names
    infile = open(WORK+"/outfiles/"+outname+".unlinked_snps",'r')
    dat = infile.readlines()
    nsamp,nsnps = dat[0].strip().split(" ")
    nsamp = int(nsamp)
    nsnps = int(nsnps)
    NDATA = np.empty([int(nsamp),int(nsnps)],dtype='object')

    ## read SNP matrix into a numpy.array
    for line in range(len(dat[1:])):
        a,b = dat[1:][line].split()
        b = b[0:nsnps]
        NDATA[line] = list(b)
    sites = np.transpose(NDATA)

    ## unpack ambiguity bases and find two most common alleles
    ## at every SNP site, save to a list
    alleles = []
    for site in sites:
        for s in site:
            if s in list("RKSYWM"):
                site = np.append(site,alignable.unstruct(s))
        snp = [s for s in site if s not in list("-NRKSYWM")]
        a = Counter(snp).most_common(2)
        if len(a)>1:
            c,d = a
            alleles.append([c[0][0],d[0][0]])
        else:
            alleles.append(['x','x'])

    ## create a dictionary mapping sample names to SNPs    
    SNPS = {}
    for line in dat[1:]:
        a,b = line.split()
        SNPS[a] = b

    ## create a dictionary with empty lists for each taxon 
    FREQ = {}
    for tax in taxa:
        FREQ[tax] = []

    ## fill the FREQ dictionary with SNPs for all 
    ## samples in that taxon
    for snp in range(int(nsnps)):
        for tax in taxa:
            freq = []
            for samp in taxa[tax]:
                freq.append(SNPS[samp][snp])
            FREQ[tax].append(freq)

    ## create empty dictionary for SNP frequencies
    FREQS = {}
    for tax in taxa:
        FREQS[tax] = []

    ## fill FREQS dictionary with SNP frequencies of 
    ## all samples for each taxon given the two most
    ## common alleles (biallelic SNPs) at each site
    for nloc in range(len(FREQ[FREQ.keys()[0]])):
        for tax in FREQ:
            zz = [m for m in FREQ[tax][nloc] if m not in list("-N")]
            if len(zz)>0:
                pp = []
                for z in zz[0:]:
                    if z in list("RKSYWM"):
                        pp += alignable.unstruct(z)
                    else:
                        pp += [z,z]
                FREQS[tax].append([str(pp).count(alleles[nloc][0]),str(pp).count(alleles[nloc][1])])
            else:
                FREQS[tax].append('xxxx')

    ## write to a file for treemix input file
    writing = [" ".join(FREQS.keys())]
    for posnp in range(int(nsnps)):
        xcheck = [FREQS[tax][posnp] for tax in FREQS]
        if 'xxxx' not in xcheck:
            if [0,0] not in xcheck:
                ## remove sites that are not SNPs given
                ## the subsampling of the data
                if not all([x[1]==0 for x in xcheck]):
                    xcheck = map(str,xcheck)
                    writing.append(" ".join([xc.replace("[","").replace("]","").replace(" ","") for xc in xcheck]))

    outfile = gzip.open(WORK+"/outfiles/"+outname+".treemix.gz",'w')
    outfile.write("\n".join(writing))
    outfile.close()

if __name__ == "__main__":
    make(WORK, outname, taxadict)
