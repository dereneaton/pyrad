
import time
import numpy as np
from alignable import unstruct
from alignable import most_common


## STILL NEED TO DO SOMETHING ABOUT Ns, and test on real data...

def makeVCF(version, locifile, mindepth, outstring):
    #version = 2.16
    #infile = "/home/deren/Dropbox/Public/PyRAD_TUTORIALS/tutorial_RAD/outfiles/c88d6m4p3.loci"
    #outfile = open("/home/deren/Desktop/test.vcf",'w')

    outfile = open(outstring, 'w')

    samples = list(set([i.strip().split(" ")[0] for i in open(locifile).readlines() if ">" in i]))
    samples.sort()

    print >>outfile, "##fileformat=VCFv4.1"
    print >>outfile, "##fileDate="+time.strftime("%Y%m%d")
    print >>outfile, "##source=pyRAD.v."+str(version)
    print >>outfile, "##reference=common_allele_at_each_locus"
    print >>outfile, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">"
    print >>outfile, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">"
    print >>outfile, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">"
    print >>outfile, "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">"
    print >>outfile, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    print >>outfile, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"
    print >>outfile, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">"
    print >>outfile, "\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO    ","FORMAT"]+list(samples))

    loci = open(locifile).read().strip().split("|")[0:]

    for locusnumber in range(len(loci)):
        samps = [i.split(" ")[0] for i in loci[locusnumber].strip().split("\n")[:-1]]
        loc = np.array([tuple(i.split(" ")[-1]) for i in loci[locusnumber].strip().split("\n")[:-1]])
        NS = str(len(loc))
        DP = str(mindepth)
        for base in range(len(loc.T)):
            col = []
            site = list(loc.T[base])
            site = list("".join(site).replace("-","").replace("N",""))
            if site:
                for bb in site:
                    if bb in list("RKYSWM"):
                        col += unstruct(bb)[0]
                        col += unstruct(bb)[1]
                    else:
                        col += bb
                REF = most_common([i for i in col if i not in list("-RKYSWMN")])
                ALT = set([i for i in col if (i in list("ATGC-N")) and (i!=REF)])
                if ALT:
                    GENO = [REF]+list(ALT)
                    GENOS = []
                    for samp in samples:
                        if samp in samps:
                            idx = samps.index(samp)
                            f = unstruct(loc.T[base][idx])
                            if ('-' in f) or ('N' in f):
                                GENOS.append("./.")
                            else:
                                GENOS.append(str(GENO.index(f[0]))+"|"+str(GENO.index(f[1])))
                        else:
                            GENOS.append("./.")
                    print >>outfile, "\t".join([`locusnumber`, `base+1`, '.', REF, ",".join(ALT), "20", "PASS",
                                                ";".join(["NS="+NS, "DP="+DP]), "GT"]+GENOS)

    outfile.close()


#version = 2.16
#infile = "/home/deren/Documents/Oaks/Virentes/analysis_pyrad/outfiles/virentes_c85d6m20p5noutg.loci"
#mindepth = 6
#outfile = "/home/deren/Documents/Oaks/Virentes/analysis_pyrad/outfiles/virentes_c85d6m20p5noutg.vcf"


if __name__ == "__main__":
    makeVCF(version, infile, mindepth, outfile)
