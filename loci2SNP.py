
import numpy as np
import sys
import gzip
try:
    from collections import Counter
except ImportError:
    from counter import Counter
from itertools import chain
import alignable


def make(WORK, outname, names, formats, seed):
    np.random.seed(int(seed))
    finalfile = open(WORK+"outfiles/"+outname+".loci").read()
    longname = max(map(len,names))

    " output .snps .snp and .phy files "
    F = {}  ## taxon dict
    S = {}  ## snp dict
    Si = {} ## unlinked snp dict
    for name in list(names):
        F[name] = []
        S[name] = []
        Si[name] = []

    " for each locus select out the SNPs"
    for loc in [i for i in finalfile.strip().split("|\n")[:]]:
        ns = []
        ss = []
        cov = {}  ## record coverage for each SNP
        for line in loc.split("\n"):
            if ">" in line:
                ns.append(line.split(" ")[0].replace(">",""))
                ss.append(line.split(" ")[-1])
            if "//" in line:
                pis = [i[0] for i in enumerate(line) if i[1] in list('*-')]
                
        " assign snps to S, and record coverage for usnps"
        for tax in F:
            if tax in ns:
                F[tax].append(ss[ns.index(tax)])
                if pis:
                    for snpsite in pis:
                        snpsite -= (longname+5)
                        S[tax].append(ss[ns.index(tax)][snpsite])
                        if snpsite not in cov:
                            cov[snpsite] = 1
                        else:
                            cov[snpsite] += 1
                        if ss[ns.index(tax)][snpsite] != '-':
                            cov[snpsite] += 1
            else:
                F[tax].append("N"*len(ss[0]))
                if pis:
                    for snpsite in pis:
                        S[tax].append("N")
                    Si[tax].append("N")
        " randomly select among snps w/ greatest coverage for unlinked snp "
        maxlist = []
        for j,k in cov.items():
            if k == max(cov.values()):
                maxlist.append(j)
        if maxlist:
            rando = pis[np.random.randint(len(pis))]
            rando -= (longname+5)
        for tax in F:
            if tax in ns:
                if pis:
                    Si[tax].append(ss[ns.index(tax)][rando])
            if pis:
                S[tax].append(" ")
            else:
                S[tax].append("_ ")

    " names "
    SF = list(F.keys())
    SF.sort()

    " print out .SNP file "
    if 's' in formats:
        snpsout = open(WORK+'outfiles/'+outname+".snps",'w')
        print >>snpsout, "## %s taxa, %s loci, %s snps" % (len(S),
                                                           len("".join(S.values()[0]).split(" "))-1,
                                                           len("".join(S[SF[0]]).replace(" ","")))
        for i in SF:
            print >>snpsout, i+(" "*(longname-len(i)+3))+"".join(S[i])
        snpsout.close()


    " print out .USNP file "
    snpout = open(WORK+'outfiles/'+outname+".unlinked_snps",'w')
    print >>snpout, len(Si),len("".join(Si.values()[0]))
    for i in SF:
        print >>snpout, i+(" "*(longname-len(i)+3))+"".join(Si[i])
    snpout.close()

    if 'k' in formats:

        "print out .str (structure) file "
        structout = open(WORK+'outfiles/'+outname+".str", 'w')
        
        B = {'A': '0',
             'T': '1',
             'G': '2',
             'C': '3',
             'N': '-9',
             '-': '-9'}

        for line in SF:
            print >>structout, line+(" "*(longname-len(line)+3))+\
                        "\t"*6+"\t".join([B[alignable.unstruct(j)[0]] for j in Si[line]])
            print >>structout, line+(" "*(longname-len(line)+3))+\
                     "\t"*6+"\t".join([B[alignable.unstruct(j)[1]] for j in Si[line]])

        structout.close()
        
        # " array of SNP data "
        # N = np.array([list(Si[i]) for i in SF])

        # " column of names "
        # namescol = list(chain( * [[i+(" "*(longname-len(i)+3)),
        #                            i+(" "*(longname-len(i)+3))] for i in SF] ))
        # " add blank columns "
        # empty = np.array(["" for i in xrange(len(SF)*2)])
        # OUT = np.array([namescol,empty,empty,empty,empty,empty,empty])
        # for col in xrange(len(N[0])):
        #     l = N[:,col]
        #     h = [alignable.unstruct(j) for j in l]
        #     h = list(chain(*h))
        #     bases = list("ATGC")
        #     final = [bases.index(i) if i not in list("-N") else '-9' for i in h]
        #     OUT = np.vstack([OUT, np.array(final)])
        # np.savetxt(WORK+'outfiles/'+outname+".str", OUT.transpose(), fmt="%s", delimiter="\t")


if __name__ == "__main__":
    make(WORK, outname, names, formats, seed)
