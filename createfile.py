
import sys

def main(version):
    output = """
==== parameter inputs for pyRAD version %s  ============================  affected step ==
./                          ## 1. Working directory                                 (all)
./*.fastq.gz                ## 2. Loc. of non-demultiplexed files (if not line 18)  (s1)
./*.barcodes                ## 3. Loc. of barcode file (if not line 18)             (s1)
usearch7.0.1090_i86linux32  ## 4. command (or path) to call usearch v.7             (s3,s6)
muscle                      ## 5. command (or path) to call muscle                  (s3,s7)
TGCAG                       ## 6. restriction overhang (e.g., C|TGCAG -> TGCAG)     (s1,s2)
2                           ## 7. N processors to use in parallel                   (all)
6                           ## 8. Mindepth: min coverage for a cluster              (s4,s5)
4                           ## 9. NQual: max # sites with qual < 20 (line 20)       (s2)
.88                         ## 10. Wclust: clustering threshold as a decimal        (s3,s6)
rad                         ## 11. Datatype: rad,gbs,ddrad,pairgbs,pairddrad,merge  (all)
4                           ## 12. MinCov: min samples in a final locus             (s7)
3                           ## 13. MaxSH: max inds with shared hetero site          (s7)
c88d6m4p3                   ## 14. prefix name for final output (no spaces)         (s7)
==== optional params below this line ===================================  affected step ==
                       ## 15.opt.: select subset (prefix* selector)                 (s2-s7)
                       ## 16.opt.: add-on (outgroup) taxa (list or prefix*)         (s6,s7)
                       ## 17.opt.: exclude taxa (list or prefix*)                   (s7)
                       ## 18.opt.: Loc. of de-multiplexed data                      (s2)
                       ## 19.opt.: maxM: N mismatches in barcodes (def. 1)          (s1)
                       ## 20.opt.: Phred Qscore offset (def. 33)                    (s2)
                       ## 21.opt.: Filter: 0=NQual 1=NQual+adapters. 2=1+strict     (s2)
                       ## 22.opt.: a priori E,H (def. 0.001,0.01, if not estimated) (s5)
                       ## 23.opt.: maxN: Ns in a consensus seq (def. 5)             (s5)
                       ## 24.opt.: maxH: hetero. sites in consensus seq (def. 5)    (s5)
                       ## 25.opt.: ploidy: max alleles in consens (def. 2) see doc  (s5)
                       ## 26.opt.: maxSNPs: step 7. (def=100). Paired (def=100,100) (s7)
                       ## 27.opt.: maxIndels: within-clust,across-clust (def. 3,99) (s3,s7)
                       ## 28.opt.: random number seed (def. 112233)              (s3,s6,s7)
                       ## 29.opt.: trim overhang left,right on final loci, def(0,0) (s7)
                       ## 30.opt.: add output formats: a,n,s,u (see documentation)  (s7)
                       ## 31.opt.: call maj. consens if dpth < stat. limit (def. 0) (s5)
                       ## 32.opt.: merge/remove paired overlap (def 0), 1=check     (s2)
                       ## 33.opt.: keep trimmed reads (def=0). Enter min length.    (s2)
                       ## 34.opt.: max stack size (int), def= max(500,mean+2*SD)    (s3)
                       ## 35.opt.: minDerep: exclude dereps with <= N copies, def=0 (s3)
                       ## 36.opt.: hierarch. cluster groups (def.=0, 1=yes)         (s6)
==== list hierachical cluster groups below this line ====================================="""  % (version)
    outfile = open("params.txt",'w')
    print >>sys.stderr, "\tnew params.txt file created"
    print >>outfile, "\n".join(output.split("\n")[1:])
    
