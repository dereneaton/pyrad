
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
        NDATA[line] = list(b)
    sites = np.transpose(NDATA)

    ## unpack ambiguity bases and find two most common alleles
    ## at every SNP site, save to a list
    alleles = []
    for site in sites:
        ds = []
        for s in site:
            if s in list("RKSYWM"):
                ds.append(alignable.unstruct(s)[0])
                ds.append(alignable.unstruct(s)[1])
            else:
                ds.append(s)
                ds.append(s)
        snp = [s for s in ds if s not in ["N",'-']]
        a = Counter(snp).most_common(3)
        alleles.append([a[0][0],a[1][0]])

    ## create a dictionary mapping sample names to SNPs    
    SNPS = OrderedDict()
    for line in dat[1:]:
        a,b = line.split()
        SNPS[a] = b

    ## create a dictionary with empty lists for each taxon 
    FREQ = OrderedDict()
    for tax in taxa:
        FREQ[tax] = []

    ## fill the FREQ dictionary with SNPs for all 
    ## samples in that taxon
    keeps = []
    for snp in range(int(nsnps)):
        GG = []
        for tax,mins in zip(taxa,minhits):
            GG.append( sum([SNPS[i][snp] not in ["N","-"] for i in taxa[tax]]) >= int(mins))
        if all(GG):
            keeps.append(snp)

    for keep in keeps:
        for tax in FREQ:
            bunch = []
            for i in taxa[tax]:
                #print tax, i, SNPS[i][keep]
                bunch.append(alignable.unstruct(SNPS[i][keep])[0])
                bunch.append(alignable.unstruct(SNPS[i][keep])[1])
            FREQ[tax].append("".join(bunch))

    ## output files
    outfile = gzip.open(WORK+"/outfiles/"+outname+".treemix.gz",'w')

    ## print header
    print >>outfile, " ".join(FREQ.keys())

    ## print data
    for i,j in enumerate(keeps):
        a1 = alleles[j][0]
        a2 = alleles[j][1]
        H = [str(FREQ[tax][i].count(a1))+","+str(FREQ[tax][i].count(a2)) for tax in FREQ]
        #print [FREQ[tax][i] for tax in FREQ], " ".join(H)

        ## exclude non-biallelic SNPs
        if "0,0" not in " ".join(H):
            print >>outfile, " ".join(H)

    outfile.close()

if __name__ == "__main__":
    make(WORK, outname, taxadict)
