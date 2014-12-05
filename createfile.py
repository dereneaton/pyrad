#!/usr/bin/env python2

import sys

def main(version):
    output = """
==** parameter inputs for pyRAD version %s  **==========================  affected step ==
./                          ## 1. Working directory                                 (all)
./*.fastq.gz                ## 2. Loc. of non-demultiplexed files (if not line 16)  (s1)
./*.barcodes                ## 3. Loc. of barcode file (if not line 16)             (s1)
TGCAG                       ## 4. Restriction overhang (e.g., C|TGCAG -> TGCAG)     (s1,s2)
2                           ## 5. N processors (parallel)                           (all)
6                           ## 6. Mindepth: min coverage for a cluster              (s4,s5)
4                           ## 7. NQual: max # sites with qual < 20 (line 18)       (s2)
.88                         ## 8. Wclust: clustering threshold as a decimal         (s3,s6)
rad                         ## 9. Datatype: rad,gbs,ddrad,pairgbs,pairddrad,merge   (all)
4                           ## 10. MinCov: min samples in a final locus             (s7)
3                           ## 11. MaxSH: max inds with shared hetero site          (s7)
c88d6m4p3                   ## 12. Prefix name for final output (no spaces)         (s7)
==== optional params below this line ===================================  affected step ==
                       ## 13.opt.: select subset (prefix* only selector)            (s2-s7)
                       ## 14.opt.: add-on (outgroup) taxa (list or prefix*)         (s6,s7)
                       ## 15.opt.: exclude taxa (list or prefix*)                   (s7)
                       ## 16.opt.: loc. of de-multiplexed data                      (s2)
                       ## 17.opt.: maxM: N mismatches in barcodes (def= 1)          (s1)
                       ## 18.opt.: phred Qscore offset (def= 33)                    (s2)
                       ## 19.opt.: filter: def=0=NQual 1=NQual+adapters. 2=strict   (s2)
                       ## 20.opt.: a priori E,H (def= 0.001,0.01, if not estimated) (s5)
                       ## 21.opt.: maxN: max Ns in a cons seq (def=5)               (s5)
                       ## 22.opt.: maxH: max heterozyg. sites in cons seq (def=5)   (s5)
                       ## 23.opt.: ploidy: max alleles in cons seq (def=2;see docs) (s4,s5)
                       ## 24.opt.: maxSNPs: (def=100). Paired (def=100,100)         (s7)
                       ## 25.opt.: maxIndels: within-clust,across-clust (def. 3,99) (s3,s7)
                       ## 26.opt.: random number seed (def. 112233)              (s3,s6,s7)
                       ## 27.opt.: trim overhang left,right on final loci, def(0,0) (s7)
                       ## 28.opt.: output formats: p,n,a,s,v,u,t,m,k,g (see docs)   (s7)
                       ## 29.opt.: call maj. consens if depth < stat. limit (def=0) (s5)
                       ## 30.opt.: keep trimmed reads (def=0). Enter min length.    (s2)
                       ## 31.opt.: max stack size (int), def= max(500,mean+2*SD)    (s3)
                       ## 32.opt.: minDerep: exclude dereps with <= N copies, def=1 (s3)
                       ## 33.opt.: use hierarchical clustering (def.=0, 1=yes)      (s6)
                       ## 34.opt.: repeat masking (def.=1='dust' method, 0=no)      (s3,s6)
                       ## 35.opt.: vsearch threads per job (def.=6; see docs)       (s3,s6)
                       ## 36.opt.: use usearch (enter path to usearch_v.7.0.1090)   (s3,s6)
==== optional: list group/clade assignments below this line (see docs) =================="""  % (version)
    outfile = open("params.txt",'w')
    print >>sys.stderr, "\tnew params.txt file created"
    print >>outfile, "\n".join(output.split("\n")[1:])
    
