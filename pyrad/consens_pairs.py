#!/usr/bin/env python2

import multiprocessing
import glob
import itertools
import sys
import scipy.stats
import scipy.misc
import numpy
import os
import operator
import gzip
from potpour import Worker

from consensdp import binomprobr, simpleconsens, hetero, unhetero, uplow, findalleles,breakalleles, removerepeat_Ns


def stack(D):
    """
    from list of bases at a site D,
    returns an ordered list of counts of bases
    """
    ## TODO: replace with Counter
    L = len(D)
    counts = []
    for i in range(len(D[0])):
        A=C=T=G=N=M=X=n=0
        for nseq in range(L):
            A += D[nseq][i].count("A")
            C += D[nseq][i].count("C")
            T += D[nseq][i].count("T")
            G += D[nseq][i].count("G")
            N += D[nseq][i].count("N")
            M += D[nseq][i].count("-")
            X += D[nseq][i].count("X")
            n += D[nseq][i].count("n")            
        counts.append( [[A,C,T,G],N,M,X,n] )
    return counts


def consensus(infile,E,H,mindepth,maxN,maxH,datatype,
              ploidy,CUT,upperSD,strict,lowcounts):
    """
    from a clust file f,
    reads in all copies at a locus and sorts
    bases at each site, tests for errors at the
    site according to error rate, calls consensus
    """
    f = gzip.open(infile,'r')
    k = itertools.izip(*[iter(f)]*2)
    bases = ['A','C','T','G']
    Dic = {}
    Errors = []
    haplo = []
    Plist = []
    locus = minsamplocus = npoly = P = 0
    while 1:
        try: first = k.next()
        except StopIteration: break
        itera = [first[0],first[1]]
        fname = itera[0].strip().split(";")[0]
        leftjust = rightjust = None

        " lists and variables for this locus"
        S       = []         ## list for sequence data
        S2      = []         ## list for sequence data
        alleles = []         ## for measuring # alleles, detect paralogs
        locus += 1           ## recording n loci
        ploidy = 0           ## for measuring # alleles, detect paralogs
        nHs = 0              ## will record heterozygous sites in this locus
        consensus  = ""       ## empty vector for consensus sequence
        consensus1 = ""       ## empty vector for consensus sequence
        consensus2 = ""       ## empty vector for consensus sequence
        basenumber = 1        ## for recording error locations

        while itera[0] != "//\n":
            nreps = int(itera[0].strip().split(";")[1].replace("size=",""))

            " append sequence * number of dereps "
            for i in xrange(nreps):
                " compatibility from pyrad 2 -> 3 "
                ss = itera[1].strip().replace("X","n")
                S.append(ss)
                S2.append(ss)
            itera = k.next()

        " separate first and second read clusters "
        firsts = [tuple(i.split("n")[0]) for i in S]
        seconds = [tuple(i.split("n")[-1]) for i in S]

        " call first read consensus "
        " Apply depth and paralog filters "
        if (len(firsts) >= min(mindepth,lowcounts)) and (len(firsts) < upperSD):  ## upper limit is meandepth + 2 SD
            minsamplocus += 1
            RAD = stack(firsts)
            for site in RAD:
                nchanged = 0         

                " minimum depth of coverage for base calling at each site"
                depthofcoverage = sum(site[0])                
                if depthofcoverage < min(mindepth,lowcounts):
                    cons = "N"; n1 = depthofcoverage-1; n2=0   ## prevents zero division error.
                else:
                    n1,n2,n3,n4 = sorted(site[0],reverse=True)

                    " speed hack = if diploid exclude if a third base present at > 20% "
                    quickthirdbasetest = 0
                    if ploidy == 2:
                        if float(n3)/(n1+n2+n3+n4) > 0.20:
                            quickthirdbasetest = 1
                    if not quickthirdbasetest:

                        """ if depth > 500 reduce by some factor for base calling """
                        if n1+n2 >= 500:   ## if > 500, random sample 500 
                            firstfivehundred = numpy.array(tuple("A"*n1+"B"*n2))
                            numpy.random.shuffle(firstfivehundred)
                            nchanged = 1
                            oldn1 = n1
                            oldn2 = n2
                            n1 = list(firstfivehundred[:500]).count("A")
                            n2 = list(firstfivehundred[:500]).count("B")

                        """ if lowcounts, make base calls by majority instead of statistics
                        when depth is below mindepth """
                        # if lowcounts:       ## include low count sites or no
                        #     if n1+n2 >= 5:
                        #         P,who = binomprobr(n1,n2,float(E),H)
                        #     else:
                        #         P,who = simpleconsens(n1,n2)
                        # else:
                        #     P,who = binomprobr(n1,n2,float(E),H)
                        """ make base calls using... """
                        if n1+n2 >= mindepth:
                            """ if above stat minimum """
                            P,maf,who = binomprobr(n1,n2,float(E),H)
                        elif n1+n2 >= lowcounts:
                            """ if above maj rule minimum"""
                            P,maf,who = simpleconsens(n1,n2)

                        """ if the base could be called with 95% probability """
                        if float(P) >= 0.95:
                            if who in 'ab':
                                if nchanged:
                                    a = [i for i,l in enumerate(site[0]) if l == oldn1]
                                else:
                                    a = [i for i,l in enumerate(site[0]) if l == n1]
                                if len(a)==2:       ## alleles came up equal freq.
                                    cons = hetero(bases[a[0]],bases[a[1]])
                                    alleles.append(basenumber)
                                else:               ## alleles came up diff freq.
                                    if nchanged:
                                        b= [i for i,l in enumerate(site[0]) if l == oldn2]
                                    else:
                                        b= [i for i,l in enumerate(site[0]) if l == n2]
                                    "if three alleles came up equal, only need if diploid paralog filter off"
                                    if a == b:
                                        cons = hetero(bases[a[0]],bases[a[1]])
                                    else:
                                        cons = hetero(bases[a[0]],bases[b[0]])
                                    alleles.append(basenumber)
                                nHs += 1
                            else:
                                if nchanged:
                                    cons = bases[site[0].index(oldn1)]
                                else:
                                    cons = bases[site[0].index(n1)]
                        else:
                            cons = "N"  ## poor base call 
                    else:
                        cons = "@"  ## third base freq fail
                consensus1 += cons
                basenumber += 1


            if "@" not in consensus1:
                if consensus1.count("N") <= maxN:        ## only allow maxN internal "N"s in a locus
                    if nHs < maxH:                       ## only allow maxH Hs, shortcut if first read fail
                        basenumber += 4  ## separator length
                        
                        " call second read consensus "
                        RAD = stack(seconds)
                        for site in RAD:
                            nchanged = 0         
                            " minimum depth of coverage for base calling at each site"
                            depthofcoverage = sum(site[0])
                            if depthofcoverage < mindepth:
                                cons = "N"; n1 = depthofcoverage-1; n2=0   
                            else:
                                n1,n2,n3,n4 = sorted(site[0],reverse=True)

                                " speed hack = if diploid exclude if a third base present at > 20% "
                                quickthirdbasetest = 0
                                if  ploidy == 2:
                                    if float(n3)/(n1+n2+n3+n4) > 0.20:
                                        quickthirdbasetest = 1
                                if not quickthirdbasetest:

                                    """ if depth > 500 reduce by some factor for base calling """
                                    if n1+n2 >= 500:   ## if > 500, random sample 500 
                                        firstfivehundred = numpy.array(tuple("A"*n1+"B"*n2))
                                        numpy.random.shuffle(firstfivehundred)
                                        nchanged = 1
                                        oldn1 = n1
                                        oldn2 = n2
                                        n1 = list(firstfivehundred[:500]).count("A")
                                        n2 = list(firstfivehundred[:500]).count("B")

                                    """ make base calls using... """
                                    if n1+n2 >= mindepth:
                                        """ if above stat minimum """
                                        P,maf,who = binomprobr(n1,n2,float(E),H)
                                    elif n1+n2 >= lowcounts:
                                        """ if above maj rule minimum"""
                                        P,maf,who = simpleconsens(n1,n2)

                                    """ if the base could be called with 95% probability """
                                    if float(P) >= 0.95:
                                        if who in 'ab':
                                            if nchanged:
                                                a = [i for i,l in enumerate(site[0]) if l == oldn1]
                                            else:
                                                a = [i for i,l in enumerate(site[0]) if l == n1]
                                            if len(a)==2:       ## alleles came up equal freq.
                                                cons = hetero(bases[a[0]],bases[a[1]])
                                                alleles.append(basenumber)
                                            else:               ## alleles came up diff freq.
                                                if nchanged:
                                                    b= [i for i,l in enumerate(site[0]) if l == oldn2]
                                                else:
                                                    b= [i for i,l in enumerate(site[0]) if l == n2]
                                                "if three alleles came up equal, only need if diploid paralog filter off"
                                                if a == b:
                                                    cons = hetero(bases[a[0]],bases[a[1]])
                                                else:
                                                    cons = hetero(bases[a[0]],bases[b[0]])
                                                alleles.append(basenumber)
                                            nHs += 1
                                        else:
                                            if nchanged:
                                                cons = bases[site[0].index(oldn1)]
                                            else:
                                                cons = bases[site[0].index(n1)]
                                    else:
                                        cons = "N"
                                else:
                                    "paralog flag"
                                    cons = "@"
                            consensus2 += cons
                            basenumber += 1
                            

                        "create concatenated consensus sequence from pairs "
                        if "@" not in consensus2:
                            consensus2.replace("-","N")
                            consensus = consensus1 + "nnnn" + consensus2


                        " filter applies to concatenated sequence "
                        if consensus:
                            if "@" not in consensus:
                                " only allow maxH polymorphic sites in a locus "
                                if nHs <= maxH:
                                    " filter for number of 2 alleles - diploids "
                                    if ploidy:
                                        al = []
                                        " only check if more than one hetero site present "
                                        if len(alleles) > 1:
                                            for i in S2:
                                                d = ""
                                                for z in alleles:
                                                    if i[z-1] in unhetero(consensus[z-1]):
                                                        d += i[z-1]+"_"
                                                if "N" not in d:
                                                    if d.count("_") == len(alleles):
                                                        al.append(d.rstrip("_"))

                                            " remove allele if it came up less than one in ten "
                                            " in which case it is likely a true heterozygous site "
                                            " but contains a sequencing error also              "
                                            ## a hack for now. But very conservative.
                                            #if len(al) >= 5:
                                            #    al = [i for i in al if al.count(i) > len(al)/10.]
                                                #TODO allow only 1 bp difference for excludes

                                            AL = sorted(set(al), key=al.count)
                                            diploid = len(AL)

                                            " set correct alleles relative to first polymorphic base"
                                            if AL:
                                                if diploid <= ploidy:
                                                    sss = [zz-1 for zz in alleles]
                                                    consensus = findalleles(consensus,sss,AL)
                                                    ## TODO: incorporate option to output alleles for haplos>2
                                                else:
                                                    consensus += "@E"
                                    else:
                                        None
                                else:
                                    consensus += "@P"
                                    
                            if "@" not in consensus:
                                #print consensus, nHs
                                " strip terminal N's from either end "
                                shortcon1 = consensus1.rstrip("N").replace("-","N")
                                " remove internal - or N, if low count "
                                shortcon1 = removerepeat_Ns(shortcon1)
                                " check for length not trimmed "
                                if (len(shortcon1) >= 32) and (len(consensus2) >= 32):
                                    Dic[fname] = shortcon1 + "nnnn" +consensus2
                                    npoly += nHs
                                        


    if ".gz" in infile[-5:]:
        consens = gzip.open(infile.replace(".clustS",".consens"),'w')
    else:
        consens = open(infile.replace(".clustS",".consens"),'w')
    for i in Dic.items():
        consens.write(str(i[0])+'\n'+str(i[1])+"\n")
    consens.close()
    sys.stderr.write(".")

    if 'pair' in datatype:
        nsites = sum([len(i)-len(CUT)-4 for i in Dic.values()])
    else:
        nsites = sum([len(i)-len(CUT) for i in Dic.values()])
    ldic = len(Dic)
    try: NP = npoly/float(nsites)
    except ZeroDivisionError: NP = 0
    return [infile.split('/')[-1], locus, minsamplocus, ldic, nsites, npoly, round(NP,7)]



def upSD(handle,mindepth):
    " function to calculate mean and SD of clustersize"
    if ".gz" in handle[-5:]:
        infile = gzip.open(handle)
    else:
        infile = open(handle)
    L = itertools.izip(*[iter(infile)]*2)
    a = L.next()[0].strip()
    depth = []
    d = int(a.split(";")[1].replace("size=",""))
    while 1:
        try: a = L.next()[0].strip()
        except StopIteration: break
        if a != "//":
            d += int(a.split(";")[1].replace("size=",""))
        else:
            depth.append(d)
            d = 0
    infile.close()
    keep = [i for i in depth if i>=(mindepth)]
    if keep:
        me = numpy.mean(keep)
        std = numpy.std(keep)
    else:
        me = 0.0
        std = 0.0
    return me, std


def main(Parallel, E, H, ID, mindepth, subset,
         maxN, maxH, ploidy, CUT, datatype,
         lowcounts, strict, WORK, maxstack):

    " find clust.xx directory "
    if not os.path.exists(WORK+'clust'+ID):
        print  "\terror: could not find "+WORK+"clust"+str(ID)+"/ directory,"+ \
                "\n\t\tif you changed the clustering threshold you must transfer *.clustS"+ \
                "\n\t\tfiles to a new directory named clust.xx with xx replaced by new clustering threshold"
        sys.exit()

    " create work queue"
    work_queue = multiprocessing.Queue()

    " iterate over files"
    outfolder = WORK+'clust'+str(ID)
    HH = glob.glob(outfolder+"/"+subset+".clustS*")
    stringout = "\n\tstep 5: created consensus seqs for %i samples, using H=%.5f E=%.5f\n\t" % (len(HH),round(H,5),round(E,5))
    sys.stderr.write(stringout)
    
    if len(HH) > 1:
        " sort files by size"
        for i in range(len(HH)):
            statinfo = os.stat(HH[i])
            HH[i] = HH[i],statinfo.st_size
        HH.sort(key=operator.itemgetter(1))
        FS = [f[0] for f in HH][::-1]
    else: FS = HH
    REMOVE = glob.glob('clust'+ID+"/cat.*")
    FS = [f for f in FS if f not in REMOVE]
    submitted = 0
    for handle in FS:
        if handle.replace('.clustS','.consens').replace('.clust','.consens') not in glob.glob(outfolder+"/*"):
            m,sd = upSD(handle,mindepth)
            if maxstack == "2SD":
                upperSD = max(500,m+(sd*2.5))
            else:
                upperSD = int(maxstack)
            work_queue.put([handle,E,H,mindepth,maxN,maxH,datatype,
                            ploidy,CUT,upperSD,strict,lowcounts])
            submitted += 1
        else:
            print "\tskipping "+handle.replace(".clustS",".consens")+\
                  ', it already exists in '+outfolder+"/"


    " create a queue to pass to workers to store the results"
    result_queue = multiprocessing.Queue()

    " spawn workers"
    jobs = []
    for i in range( min(Parallel,submitted) ):
        worker = Worker(work_queue, result_queue, consensus)
        jobs.append(worker)
        worker.start()
    for j in jobs:
        j.join()

    " get results"
    stats = open(WORK+'stats/s5.consens.txt','a+')
    print >>stats,  "taxon\tnloci\tf1loci\tf2loci\tnsites\tnpoly\tpoly"
    for i in range(submitted):
        a,b,c,d,e,f,g = result_queue.get()
        nn = a.replace(".clustS.gz","")
        print >> stats, "\t".join(map(str,[nn,b,c,d,e,f,g]))
    print >>stats, """
    ## nloci = number of loci
    ## f1loci = number of loci with >N depth coverage
    ## f2loci = number of loci with >N depth and passed paralog filter
    ## nsites = number of sites across f loci
    ## npoly = number of polymorphic sites in nsites
    ## poly = frequency of polymorphic sites"""
    stats.close()




