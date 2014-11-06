
import numpy as np
import sys
import gzip
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict
try:
    from collections import Counter
except ImportError:
    from counter import Counter
import alignable


def make(WORK, outname, taxadict, minhits):

    ## cleanup taxadict to just sample names
    taxa = OrderedDict()
    for group in taxadict:
        taxa[group] = []
        for samp in taxadict[group]:
            a = samp.split("/")[-1].replace(".consens.gz","")
            taxa[group].append(a)

    print "\t    data set reduced for group coverage minimums"        
    for i,j in zip(taxa,minhits):
        print "\t   ",i, taxa[i], 'minimum=',j
    
    ## read in data to sample names
    infile = open(WORK.rstrip("/")+"/outfiles/"+outname+".unlinked_snps",'r')
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
    FREQ = OrderedDict()
    for tax in taxa:
        FREQ[tax] = []

    ## fill the FREQ dictionary with SNPs for all 
    ## samples in that taxon
    keep = []
    for snp in range(int(nsnps)):
        GG = []
        for tax,mins in zip(taxa,minhits):
            GG.append( sum([SNPS[i][snp] not in ["N","-"] for i in taxa[tax]]) >= int(mins))
        if all(GG):
            keep.append(snp)
            for tax in FREQ:
                FREQ[tax].append([SNPS[i][snp] for i in taxa[tax]])

    ## output files
    outfile = gzip.open(WORK+"/outfiles/"+outname+".treemix.gz",'w')

    ## print header
    print >>outfile, " ".join(FREQ.keys())

    ## print data
    for i,j in enumerate(keep):
        a1 = alleles[j][0]
        a2 = alleles[j][1]
        H = [str(FREQ[tax][i].count(a1))+","+str(FREQ[tax][i].count(a2)) for tax in FREQ]
        print >>outfile, " ".join(H)

    outfile.close()

if __name__ == "__main__":
    make(WORK, outname, taxadict)
