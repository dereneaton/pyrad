#!/usr/bin/env python2

""" create .snp and related outputs """

import numpy as np
try:
    from collections import Counter
except ImportError:
    from counter import Counter
import alignable


def make(params, names): 
    """ make function """
    #WORK, outname, names, formats, seed, ploidy):
    np.random.seed(int(params["seed"]))
    finalfile = open(params["work"]+"outfiles/"+\
                     params["outname"]+".loci").read()
    longname = max([len(i) for i in names])

    ## output .snps and .unlinked_snps"
    snpdict = {name:[] for name in names}      ## snp dict  S
    usnpdict = {name:[] for name in names}     ## unlinked snp dict Si
    # for name in list(names):
    #     S[name] = []
    #     Si[name] = []

    ## record bi-allelic snps"
    bis = 0

    ## for each locus select out the SNPs"
    for loc in finalfile.strip().split("|\n"):    #[:-1]:
        pis = ""
        names = []  ## ns
        seqs = []   ## ss
        cov = {}  ## record coverage for each SNP
        for line in loc.split("\n"):
            if ">" in line:
                names.append(line.split()[0].replace(">", ""))
                seqs.append(line.split()[-1])
            else:
                pis = [i[0] for i in enumerate(line) if i[1] in list('*-')]

        
        ## assign snps to S, and record coverage for usnps"
        for tax in snpdict:
            if tax in names:
                if pis:
                    for snpsite in pis:
                        snpsite -= (longname+5)
                        snpdict[tax].append(seqs[names.index(tax)][snpsite])
                        ## weighting of site for "random" selection
                        if snpsite not in cov:
                            cov[snpsite] = 1
                        else:
                            cov[snpsite] += 1
                        ## downweight selection of gap sites "
                        if seqs[names.index(tax)][snpsite] != '-':
                            cov[snpsite] += 1
            else:
                if pis:
                    for snpsite in pis:
                        snpdict[tax].append("N")
                    usnpdict[tax].append("N")

        ## randomly select among snps w/ greatest coverage for unlinked snp "
        maxlist = []
        for j, k in cov.items():
            if k == max(cov.values()):
                maxlist.append(j)

        ## Is bi-allelic ?
        bisnps = []
        for i in maxlist:
            if len(set([seqs[names.index(tax)][i] for tax \
                        in snpdict if tax in names])) < 3:
                bisnps.append(i)

        #rando = pis[np.random.randint(len(pis))]
        #rando -= (longname+5)
        if bisnps:
            rando = bisnps[np.random.randint(len(bisnps))]
        elif maxlist:
            rando = maxlist[np.random.randint(len(maxlist))]
        for tax in snpdict:
            if tax in names:
                if pis:
                    ## if none are bi-allelic "
                    if not bisnps:
                        bis += 1
                    usnpdict[tax].append(seqs[names.index(tax)][rando])
            if pis:
                ## add spacer between loci
                snpdict[tax].append(" ")
            else:
                ## invariable locus
                snpdict[tax].append("_ ")

    ## names
    snames = list(snpdict.keys())
    snames.sort()

    ## print out .SNP file
    if 's' in params["outform"]:
        snpsout = open(params["work"]+'outfiles/'+\
                       params["outname"]+".snps", 'w')
        print >>snpsout, "## %s taxa, %s loci, %s snps" % \
                            (len(snpdict),
                             len("".join(snpdict.values()[0]).split(" "))-1,
                             len("".join(snpdict[snames[0]]).replace(" ", "")))
        for i in snames:
            print >>snpsout, i+(" "*(longname-len(i)+3))+"".join(snpdict[i])
        snpsout.close()


    ## print out .USNP file
    snpout = open(params["work"]+'outfiles/'+\
                  params["outname"]+".unlinked_snps", 'w')
    print >>snpout, len(usnpdict), len("".join(usnpdict.values()[0]))
    for i in snames:
        print >>snpout, i+(" "*(longname-len(i)+3))+"".join(usnpdict[i])
    snpout.close()

    statsout = open(params["work"]+"stats/"+\
                    params["outname"]+".stats", 'a')
    print >>statsout, "sampled unlinked SNPs=", len(usnpdict.values()[0])
    print >>statsout, "sampled unlinked bi-allelic SNPs= "+\
                       str(len(usnpdict.values()[0])-bis)
    statsout.close()

    if 'k' in params["outform"]:
        ## print out .str (structure) file "
        structout = open(params["work"]+'outfiles/'+\
                         params["outname"]+".str", 'w')
        
        stdict = {'A': '0',
                  'T': '1',
                  'G': '2',
                  'C': '3',
                  'N': '-9',
                  '-': '-9'}

        if params["haplos"] > 1:
            for sname in snames:
                #print stdict
                print >>structout, sname+(" "*(longname-len(sname)+3))+\
                      "\t"*6+"\t".join([stdict[alignable.unstruct(j)[0]] \
                                        for j in usnpdict[sname]])
                print >>structout, sname+(" "*(longname-len(sname)+3))+\
                      "\t"*6+"\t".join([stdict[alignable.unstruct(j)[1]] \
                                        for j in usnpdict[sname]])
        else:
            for sname in snames:
                print >>structout, sname+(" "*(longname-len(sname)+3))+\
                      "\t"*6+"\t".join([stdict[alignable.unstruct(j)[1]] \
                                        for j in usnpdict[sname]])
        structout.close()


    if 'g' in params["outform"]:
        ## print out .geno file usnps
        genoout = open(params["work"]+'outfiles/'+\
                       params["outname"]+".usnps.geno", 'w')
        for i in range(len(usnpdict.values()[0])):
            getref = 0
            ref = "N"
            while ref == "N":
                ref = alignable.unstruct(usnpdict[snames[getref]][i])[0]
                getref += 1
            snprow = "".join([str(i) for i in \
                             [alignable.unstruct(usnpdict[j][i]).count(ref) \
                              if usnpdict[j][i] != "N" \
                              else "9" for j in snames]])
            ## print ref,SNProw
            if len(set(snprow)) > 1:
                print >>genoout, snprow 
        genoout.close()

        ## print out .geno file all snps
        genoout = open(params["work"]+'outfiles/'+\
                       params["outname"]+".snps.geno", 'w')
        for i in range(len(snpdict.values()[0])):
            if snpdict[snames[0]][i].strip("_").strip():
                getref = 0
                ref = "N"
                while ref == "N":
                    ref = alignable.unstruct(snpdict[snames[getref]][i])[0]
                    getref += 1
                snprow = "".join([str(i) for i in \
                            [alignable.unstruct(snpdict[j][i]).count(ref) if \
                            snpdict[j][i] != "N" \
                            else "9" for j in snames]])
                ## print ref,SNProw
                if len(set(snprow)) > 1:
                    print >>genoout, snprow 
        genoout.close()


if __name__ == "__main__":
    make(params, names)
    #WORK, outname, names, formats, seed, ploidy)
