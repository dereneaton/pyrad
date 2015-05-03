#!/usr/bin/env python2

""" create VCF format including depth information if available """

import time
import numpy as np
import alignable
try:
    from collections import OrderedDict, Counter
except ImportError:
    from ordereddict import OrderedDict, Counter

def make(params, version, names):
    outfile = open(params["work"]+"/outfiles/"+\
                     params["outname"]+".vcf", 'w')
    inloci = params["work"]+"/outfiles/"+params["outname"]+".loci"
    names = list(names)
    names.sort()

    print >>outfile, "##fileformat=VCFv4.1"
    print >>outfile, "##fileDate="+time.strftime("%Y%m%d")
    print >>outfile, "##source=pyRAD.v."+str(version)
    print >>outfile, "##reference=common_allele_at_each_locus"
    print >>outfile, "##INFO=<ID=NS"+\
                        ",Number=1,Type=Integer,Description="+\
                        "\"Number of Samples With Data\">"
    print >>outfile, "##INFO=<ID=DP"+\
                        ",Number=1,Type=Integer,Description="+\
                        "\"Total Depth\">"
    print >>outfile, "##INFO=<ID=AF"+\
                        ",Number=A,Type=Float,Description="+\
                        "\"Allele Frequency\">"
    print >>outfile, "##INFO=<ID=AA"+\
                        ",Number=1,Type=String,Description="+\
                        "\"Ancestral Allele\">"
    print >>outfile, "##FORMAT=<ID=GT"+\
                        ",Number=1,Type=String,Description="+\
                        "\"Genotype\">"
    print >>outfile, "##FORMAT=<ID=GQ"+\
                        ",Number=1,Type=Integer,Description="+\
                        "\"Genotype Quality\">"
    print >>outfile, "##FORMAT=<ID=DP"+\
                        ",Number=1,Type=Integer,Description="+\
                        "\"Read Depth\">"
    print >>outfile, "\t".join(["#CHROM", "POS", "ID", "REF",
                                "ALT", "QUAL", "FILTER", 
                                "INFO    ", "FORMAT"]+list(names))

    loci = open(inloci).read().split("|")[:-1]
    vcflist = []
    for locusnumber in range(len(loci)):
        samps = [i.split()[0][1:] for i in \
                 loci[locusnumber].strip().split("\n") \
                 if ">" in i]
        loc = np.array([tuple(i.split()[-1]) for i in \
                        loci[locusnumber].strip().split("\n") \
                        if ">" in i])
        ## get number of samples
        vcf_ns = str(len(loc))
        ## get depth info about samples
        vcf_dp = str(params["minsamp"])
        ## 
        ##
        ##
        ##

        ## 
        for site in range(len(loc.T)):
            col = []
            _ = list(loc.T[site])
            esite = list("".join(_).replace("-", "").replace("N", ""))
            if esite:
                for base in esite:
                    if base in list("RKYSWM"):
                        col += alignable.unstruct(base)[0]
                        col += alignable.unstruct(base)[1]
                    else:
                        col += base
                ref = Counter([i for i in col if \
                               i not in list("-RKYSWMN")]).\
                               most_common(1)[0][0]
                alts = set([i for i in col if (i in list("ATGC-N")) \
                           and (i != ref)])
                if alts:
                    geno = [ref]+list(alts)
                    genos = []
                    for samp in names:
                        if samp in samps:
                            idx = samps.index(samp)
                            rbase = alignable.unstruct(loc.T[site][idx])
                            if ('-' in rbase) or ('N' in rbase):
                                genos.append("./.")
                            else:
                                genos.append(str(geno.index(rbase[0]))+\
                                         "|"+str(geno.index(rbase[1])))
                        else:
                            genos.append("./.")
                    vcflist.append("\t".join([str(locusnumber+1),
                                              str(site+1),
                                              '.', ref, ",".join(alts), 
                                              "20", "PASS", 
                                              ";".join(["NS="+vcf_ns,
                                                        "DP="+vcf_dp]),
                                              "GT"]+genos))
        if not locusnumber % 1000:
            outfile.write("\n".join(vcflist)+"\n")
            vcflist = []
                                              
 #print >>outfile, "\t".join([`locusnumber+1`, `base+1`,
 ## '.', REF, ",".join(ALT), "20", "PASS",
 #                            ";".join(["NS="+NS, "DP="+DP]), "GT"]+GENOS)
    

    outfile.write( "\n".join(vcflist) )
    outfile.close()

if __name__ == "__main__":
    make(WORK, version, outname, mindepth, names)
