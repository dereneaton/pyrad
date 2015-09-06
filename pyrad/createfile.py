#!/usr/bin/env python2

""" create params file for pyrad analyses """

import sys

def main(version):
    """ print params list to file """
    output = """
==** parameter inputs for pyRAD version %s  **======================= affected step ==
./                      ## 1. Working directory                                 (all)
./*.fastq.gz            ## 2. Loc. of non-demultiplexed files (if not line 18)  (s1)
./*_barcodes.txt        ## 3. Loc. of barcode file (if not line 18)             (s1)
vsearch                 ## 4. command (or path) to call vsearch (or usearch)    (s3,s6)
muscle                  ## 5. command (or path) to call muscle                  (s3,s7)
TGCAG                   ## 6. Restriction overhang (e.g., C|TGCAG -> TGCAG)     (s1,s2)
2                       ## 7. N processors (parallel)                           (all)
6                       ## 8. Mindepth: min coverage for a cluster              (s4,s5)
4                       ## 9. NQual: max # sites with qual < 20 (or see line 20)(s2)
.88                     ## 10. Wclust: clustering threshold as a decimal        (s3,s6)
rad                     ## 11. Datatype: rad,gbs,pairgbs,pairddrad,(others:see docs)(all)
4                       ## 12. MinCov: min samples in a final locus             (s7)
3                       ## 13. MaxSH: max inds with shared hetero site          (s7)
c88d6m4p3               ## 14. Prefix name for final output (no spaces)         (s7)
==== optional params below this line =================================  affected step ==
                     ## 15.opt.: select subset (prefix* only selector)            (s2-s7)
                     ## 16.opt.: add-on (outgroup) taxa (list or prefix*)         (s6,s7)
                     ## 17.opt.: exclude taxa (list or prefix*)                   (s7)
                     ## 18.opt.: loc. of de-multiplexed data                      (s2)
                     ## 19.opt.: maxM: N mismatches in barcodes (def= 1)          (s1)
                     ## 20.opt.: phred Qscore offset (def= 33)                    (s2)
                     ## 21.opt.: filter: def=0=NQual 1=NQual+adapters. 2=strict   (s2)
                     ## 22.opt.: a priori E,H (def= 0.001,0.01, if not estimated) (s5)
                     ## 23.opt.: maxN: max Ns in a cons seq (def=5)               (s5)
                     ## 24.opt.: maxH: max heterozyg. sites in cons seq (def=5)   (s5)
                     ## 25.opt.: ploidy: max alleles in cons seq (def=2;see docs) (s4,s5)
                     ## 26.opt.: maxSNPs: (def=100). Paired (def=100,100)         (s7)
                     ## 27.opt.: maxIndels: within-clust,across-clust (def. 3,99) (s3,s7)
                     ## 28.opt.: random number seed (def. 112233)              (s3,s6,s7)
                     ## 29.opt.: trim overhang left,right on final loci, def(0,0) (s7)
                     ## 30.opt.: output formats: p,n,a,s,v,u,t,m,k,g,* (see docs) (s7)
                     ## 31.opt.: maj. base call at depth>x<mindepth (def.x=mindepth) (s5)
                     ## 32.opt.: keep trimmed reads, min length (def=32)          (s2)
                     ## 33.opt.: max stack size (int), def= max(500,mean+2*SD)    (s3)
                     ## 34.opt.: minDerep: exclude dereps with <= N copies, def=1 (s3)
                     ## 35.opt.: use hierarchical clustering (def.=0, 1=yes)      (s6)
                     ## 36.opt.: repeat masking (def.=1='dust' method, 0=no)      (s3,s6)
                     ## 37.opt.: vsearch max threads per job (def.=6; see docs)   (s3,s6)
==== optional: list group/clade assignments below this line (see docs) ================"""  % (version)
    outfile = open("params.txt", 'w')
    print >>sys.stderr, "\tnew params.txt file created"
    print >>outfile, "\n".join(output.split("\n")[1:])
    
