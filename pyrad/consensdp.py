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


def binomprobr(n1,n2,e,r):
    """
    given two bases are observed at a site
    n1 and n2, and the error rate e, the
    probability the site is truly aa,bb,ab
    is calculated using binomial distribution
    as in Li_et al 2009, 2011, and if
    coverage > 500, 500 reads were randomly
    sampled.
    """
    maf = n1/(n1+n2)
    prior_homo = ((1.-r)/2.)
    prior_het = r
    ab = scipy.misc.comb(n1+n2,n1)/(2.**(n1+n2))
    aa= scipy.stats.binom.pmf(n1,n1+n2,e)
    bb= scipy.stats.binom.pmf(n2,n1+n2,e)
    aa = aa*prior_homo
    bb = bb*prior_homo
    ab = ab*prior_het
    Q = [aa,bb,ab]
    Qn = ['aa','bb','ab']
    P = max(Q)/float(aa+bb+ab)
    return [P,maf,Qn[Q.index(max(Q))]]


def simpleconsens(n1,n2):
    """
    majority consensus calling for sites
    with too low of coverage for
    statistical calling. Only used
    with 'lowcounts' option.
    """
    Qn = ['aa','bb','ab']
    maf = n1/(n1+n2)
    # if not n2:
    #     P = 1.0
    #     aa = 1.0
    #     ab = bb = 0.0
    # else:
    #     P = 0.99
    #     aa = bb = 0.0
    #     ab = 1.0
    ## create an option that saves 
    ## frequencies. Useful for pooled sample data sets.
    #Q = [aa,bb,ab]
    #return [P,Qn[Q.index(max(Q))]]
    return [1.0,maf,'aa']


def hetero(n1,n2):
    """
    returns IUPAC symbol for ambiguity bases,
    used for polymorphic sites.
    """
    D = {('G','A'):"R",
         ('G','T'):"K",
         ('G','C'):"S",
         ('T','C'):"Y",
         ('T','A'):"W",
         ('C','A'):"M"}
    a = D.get((n1,n2))
    b = D.get((n2,n1))
    if a:
        return a
    else:
        return b


def unhetero(amb):
    amb = amb.upper()
    " returns bases from ambiguity code"
    D = {"R":("G","A"),
         "K":("G","T"),
         "S":("G","C"),
         "Y":("T","C"),
         "W":("T","A"),
         "M":("C","A")}
    return D.get(amb)


def uplow(b1):
    " allele precedence "
    " G > T > C > A "
    D = {('G','A'):"G",
         ('A','G'):"G",
         ('G','T'):"G",
         ('T','G'):"G",
         ('G','C'):"G",
         ('C','G'):"G",
         ('T','C'):"T",
         ('C','T'):"T",
         ('T','A'):"T",
         ('A','T'):"T",
         ('C','A'):"C",
         ('A','C'):"C"}
    r = D.get(b1)
    if not r:
        r = b1[0]
    return r
    


def findalleles(consensus,sss,bbb):
    cons = list(consensus)
    " relative to first base "
    bigbase = uplow(tuple([i.split("_")[0] for i in bbb]))
    bigallele = bbb.index([i for i in bbb if i.split("_")[0] == bigbase][0])
    for k in range(1,len(sss)):
        c = uplow(tuple([i.split("_")[k] for i in bbb]))
        which = bbb.index([i for i in bbb if i.split("_")[k] == c][0])
        if bbb[bigallele] != bbb[which]:
            cons[sss[k]] = cons[sss[k]].lower()
    return "".join(cons)


def breakalleles(consensus):
    """ break ambiguity code consensus seqs
    into two alleles """
    a1 = ""
    a2 = ""
    bigbase = ""
    for base in consensus:
        if base in tuple("RKSYWM"):
            a,b = unhetero(base)
            d = set([a,b])
            a1 += uplow((a,b))
            a2 += d.difference(uplow((a,b))).pop()
            if not bigbase:
                bigbase = uplow((a,b))
        elif base in tuple("rksywm"):
            a,b = unhetero(base)
            d = set([a,b])
            a2 += uplow((a,b))
            a1 += d.difference(uplow((a,b))).pop()
        else:
            a1 += base
            a2 += base
    return a1,a2



def stack(D):
    """
    from list of bases at a site D,
    returns an ordered list of counts of bases
    """
    L = len(D)
    counts = []
    for i in range(len(D[0])):
        A=C=T=G=N=M=X=0
        for nseq in range(L):
            A += D[nseq][i].count("A")
            C += D[nseq][i].count("C")
            T += D[nseq][i].count("T")
            G += D[nseq][i].count("G")
            N += D[nseq][i].count("N")
            M += D[nseq][i].count("-")
            X += D[nseq][i].count("X")
        counts.append( [[A,C,T,G],N,M,X] )
    return counts


# def ffmin(x):
#     d = []
#     for i,j in enumerate(x):
#         if j not in ["-","N"]:
#             d.append(i)
#     return min(d)

# def ffmax(x):
#     d = []
#     for i,j in enumerate(x):
#         if j not in ["-","N"]:
#             d.append(i)
#     return max(d)


def removerepeat_Ns(shortcon):
    """
    checks for interior Ns in consensus seqs
    remove those that arise next to *single repeats*
    of at least 3 bases on either side, which may be
    sequencing errors on deep coverage repeats """
    Nlocs = [i for i,j in enumerate(shortcon) if j=="N"]
    repeats = set()
    for n in Nlocs:
        r1 = len(set(list(shortcon)[n-3:n]))
        if r1 < 2:
            repeats.add(n)
        r2  = len(set(list(shortcon)[n+1:n+4]))
        if r2 < 2:
            repeats.add(n)
    return "".join([j for (i,j) in enumerate(shortcon) if i not in repeats])




def consensus(infile,E,H,mindepth,maxN,maxH,datatype,
              haplos,CUT,upperSD,strict,lowcounts):
    """
    from a clust file f,
    reads in all copies at a locus and sorts
    bases at each site, tests for errors at the
    site according to error rate, calls consensus
    """
    f = gzip.open(infile)
    k = itertools.izip(*[iter(f)]*2)
    bases = ['A','C','T','G']
    Dic = {}
    Errors = []
    haplo = []
    #Plist = []
    locus = minsamplocus = npoly = P = 0
    while 1:
        try: first = k.next()
        except StopIteration: break
        itera = [first[0],first[1]]
        fname = itera[0].strip().split(";")[0]
        leftjust = rightjust = None

        " lists and variables for this locus"
        S       = []         ## list for sequence data
        alleles = []         ## for measuring # alleles, detect paralogs
        locus += 1           ## recording n loci
        ploidy = 0          ## for measuring # alleles, detect paralogs
        nHs = 0              ## will record heterozygous sites in this locus
        consensus = ""       ## empty vector for consensus sequence
        basenumber = 1       ## for recording error locations
        lefts = []
        rights = []
        while itera[0] != "//\n":
            " append sequence * number of dereps "
            nreps = int(itera[0].strip().split(";")[1].replace("size=",""))
            for i in xrange(nreps):
                S.append(tuple(itera[1].strip())) 
                #print i, itera[1].strip(), itera[0].strip()[-1], leftjust, rights
                
            " record left and right most index of seed and hits (for GBS) "
            if datatype in ['gbs','merged']:
                " leftjust is seed's left "
                if itera[0].strip()[-1] == ";":
                    leftjust = itera[1].index([i for i in itera[1] if i not in list("-N")][0])

                " rightjust is the shortest reverse hit "
                if itera[0].strip()[-1] == "-":
                    rights.append(max(-1,[itera[1].rindex(i) for i in itera[1] if i in list("ATGC")]))
                    #if rights == -1:
                    #    print itera
                    #lefts.append(itera[1].index([i for i in itera[1] if i not in list("-N")][0]))

            itera = k.next()

        " trim off overhang edges of gbs reads "
        if datatype in ['gbs','merged']:
            if rights:
                " record in name that there was a reverse hit"
                fname = "_".join(fname.split("_")[0:-1])+"_c1"
                try: rightjust = min([min(i) for i in rights])
                except ValueError:
                    S = ""
            
            for s in xrange(len(S)):
                S[s] = S[s][leftjust:]
                if rightjust:
                    #print rights, rightjust, 'right,just'
                    S[s] = S[s][:rightjust+1]
                    
            #for i in S:
            #    print "".join(i)

            #if any([i < leftjust for i in lefts]):
            #    fname = "_".join(fname.split("_")[0:-1])+"_c1"
                #print "".join(list(S[s])), "new"
            
        " Apply depth and paralog filters "
        if (len(S) >= min(lowcounts,mindepth)) and (len(S) < upperSD):  
            minsamplocus += 1
            RAD = stack(S)
            for site in RAD:
                nchanged = 0         

                " minimum depth of coverage for base calling "
                depthofcoverage = sum(site[0])
                if depthofcoverage < min(mindepth,lowcounts):
                    cons = "N"; n1 = depthofcoverage-1; n2=0   ## prevents zero division error.
                else:
                    n1,n2,n3,n4 = sorted(site[0],reverse=True)

                    " speed hack = if diploid exclude if a third base present at > 20% "
                    quickthirdbasetest = 0
                    if haplos == 2:
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
                consensus += cons
                basenumber += 1

            " only allow maxH polymorphic sites in a locus "
            if "@" not in consensus:
                if nHs <= maxH:
                    " filter to limit to N haplotypes (e.g., diploid) "
                    if haplos:
                        al = []
                        if len(alleles) > 1:
                            for i in S:
                                d = ""
                                for z in alleles:
                                    if i[z-1] in unhetero(consensus[z-1]):
                                        d += i[z-1]+"_"
                                if "N" not in d:
                                    if d.count("_") == len(alleles):
                                        al.append(d.rstrip("_"))

                            " remove very rare thirds representing a possible error at a heterozygous site \
                            that changed the base to the alternate allele at that site "
                            #if len(al) >= 50:
                            al = [i for i in al if al.count(i) > len(al)*.25]

                            AL = sorted(set(al), key=al.count)
                            ploidy = len(AL)
                            #Plist.append(ploidy)
                            
                            " set correct alleles relative to first polymorphic base"
                            if AL:
                                if ploidy <= haplos:  
                                    sss = [zz-1 for zz in alleles]
                                    consensus = findalleles(consensus,sss,AL)
                                else:
                                    consensus += "@E"
                                    # print ploidy, haplos
                                    # print alleles
                                    # print "AL", AL
                                    # print "al", al

                        #else: Plist.append(1)
                        
                    if "@" not in consensus:
                        " strip N's from either end "
                        shortcon = consensus.lstrip("N").rstrip("N").replace("-","N")
                        shortcon = removerepeat_Ns(shortcon)
                        if shortcon.count("N") <= maxN:        ## only allow maxN internal "N"s in a locus
                            if len(shortcon) >= 32:            ## minimum length set to 36
                                #print shortcon, 'keep'
                                npoly += nHs
                                Dic[fname] = shortcon


    #with open(infile.replace(".clustS",".ploids"),'w+') as ploidout:
    #    ploidout.write(",".join(map(str,Plist)))
    
    consens = gzip.open(infile.replace(".clustS",".consens"),'w+')
    for i in Dic.items():
        consens.write(str(i[0])+'\n'+str(i[1])+"\n")
    consens.close()
    sys.stderr.write(".")


    if datatype in ['pairgbs','pairddrad']:
        " -4 for the xxxx "
        nsites = sum([len(i)-len(CUT)-4 for i in Dic.values()])
    else:
        nsites = sum([len(i)-len(CUT) for i in Dic.values()])
    ldic = len(Dic)
    try: NP = npoly/float(nsites)
    except ZeroDivisionError: NP = 0
    return [infile.split('/')[-1], locus, minsamplocus, ldic, nsites, npoly, round(NP,7)]



def upSD(handle,mindepth):
    " function to calculate mean and SD of clustersize"
    infile = gzip.open(handle)
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
         maxN, maxH, haplos, CUT, datatype,
         lowcounts, strict, WORK, maxstack):

    " find clust.xx directory "
    if not os.path.exists(WORK+'clust'+ID):
        print  "\terror: could not find "+WORK+"clust"+str(ID)+"/ directory,"+ \
                "\n\t\tif you changed the clustering threshold you must transfer *.clustS"+ \
                "\n\t\tfiles to a new directory named clust.xx with xx replaced by new clustering threshold"
        sys.exit()

    " load up work queue"
    work_queue = multiprocessing.Queue()

    " iterate over files"
    outfolder = WORK+'clust'+str(ID)
    HH = glob.glob(outfolder+"/"+subset+".clustS*")
    stringout = "\n\tstep 5: creating consensus seqs for %i samples, using H=%.5f E=%.5f\n\t" % (len(HH),round(H,5),round(E,5))
    sys.stderr.write(stringout)
    
    if len(HH) > 1:
        " sort files by size"
        for i in xrange(len(HH)):
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
                            haplos,CUT,upperSD,strict,lowcounts])
            submitted += 1
        else:
            print "\tskipping "+handle.replace(".clustS",".consens")+\
                  ', it already exists in '+outfolder+"/"


    " create a queue to pass to workers to store the results"
    result_queue = multiprocessing.Queue()

    " spawn workers"
    jobs = []
    for i in xrange( min(Parallel,submitted) ):
        worker = Worker(work_queue, result_queue, consensus)
        jobs.append(worker)
        worker.start()
    for j in jobs:
        j.join()

    " get results"
    stats = open(WORK+'stats/s5.consens.txt','a+')
    print >>stats,  "taxon          \tnloci\tf1loci\tf2loci\tnsites\tnpoly\tpoly"
    for i in range(submitted):
        a,b,c,d,e,f,g = result_queue.get()
        print >> stats, "\t".join(map(str,[a.replace(".clustS.gz","")+" "*(10-len(a)),b,c,d,e,f,g]))
    print >>stats, """
    ## nloci = number of loci
    ## f1loci = number of loci with >N depth coverage
    ## f2loci = number of loci with >N depth and passed paralog filter
    ## nsites = number of sites across f loci
    ## npoly = number of polymorphic sites in nsites
    ## poly = frequency of polymorphic sites"""
    stats.close()




