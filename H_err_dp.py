import scipy.stats
import scipy.optimize
import numpy
import itertools
import sys
import glob
import multiprocessing
import os
import gzip
from potpour import Worker



def makeP(N):
    """ returns a list of freq. for ATGC"""
    sump = sum([sum(i) for i in N])
    try: p1 = sum([float(i[0]) for i in N])/sump
    except ZeroDivisionError: p1 = 0.0
    try: p2 = sum([float(i[1]) for i in N])/sump
    except ZeroDivisionError: p2 = 0.0
    try: p3 = sum([float(i[2]) for i in N])/sump
    except ZeroDivisionError: p3 = 0.0
    try: p4 = sum([float(i[3]) for i in N])/sump
    except ZeroDivisionError: p4 = 0.0
    return [p1,p2,p3,p4]


def L1(E,P,N):
    """probability homozygous"""
    h = []
    s = sum(N)
    for i,l in enumerate(N):
        p = P[i]
        b = scipy.stats.binom.pmf(s-l,s,E)
        h.append(p*b)
    return sum(h)


def L2(E,P,N):
    """probability of heterozygous"""
    h = []
    s = sum(N)
    for l,i in enumerate(N):
        for j,k in enumerate(N):
            if j>l:
                one = 2.*P[l]*P[j]
                two = scipy.stats.binom.pmf(s-i-k,s,(2.*E)/3.)
                three = scipy.stats.binom.pmf(i,k+i,0.5)
                four = 1.-(sum([q**2. for q in P]))
                h.append(one*two*(three/four))
    return sum(h)


def totlik(E,P,H,N):
    """ total probability """
    lik = ((1-H)*L1(E,P,N)) + (H*L2(E,P,N))
    return lik

def LL(x0,P,Tab):
    """ Log likelihood score given values [H,E] """
    H = x0[0]
    E = x0[1]
    L = []
    if (H <= 0.) or (E <= 0.):
        r = numpy.exp(100)
    else:
        for i in Tab:
            ll = totlik(E,P,H,i[0])
            if ll > 0:
                L.append(i[1] * numpy.log(ll))
        r = -sum(L)
        #print "\t".join(map(str,[r, H, E]))
    return r


def LL_haploid(E,P,Tab):
    """ Log likelihood score given values [H,E] """
    H = 0.
    L = []
    if (E <= 0.):
        r = numpy.exp(100)
    else:
        for i in Tab:
            ll = totlik(E,P,H,i[0])
            if ll > 0:
                L.append(i[1] * numpy.log(ll))
        r = -sum(L)
        #print "\t".join(map(str,[r, H, E]))
    return r



def table_c(N):
    """ makes a dictionary with counts of base counts [x,x,x,x]:x,
    speeds up Likelihood calculation"""
    Tab = {}
    k = iter(N)
    while 1:
        try:
            d = k.next()
        except StopIteration: break
        if tuple(d) in Tab:
            Tab[tuple(d)] += 1
        else:
            Tab[tuple(d)] = 1
    L = []
    for i,j in Tab.items():
        [i,j]
        L.append([i,j])
    return [i for i in L if (0,0,0,0) not in i]


def stack(D):
    """
    from list of bases at a site D,
    returns an ordered list of counts of bases
    """
    L = len(D)
    counts = []
    for i in range(len(D[0])):
        A=C=T=G=N=S=0
        for nseq in range(L):
            A += D[nseq][i].count("A")
            C += D[nseq][i].count("C")
            T += D[nseq][i].count("T")
            G += D[nseq][i].count("G")
            N += D[nseq][i].count("N")
            S += D[nseq][i].count("-")
        counts.append( [[A,C,T,G],N,S] )
    return counts



def consensus(f, minsamp, CUT1, CUT2, datatype):
    """ makes a list of lists of reads at each site """
    f = gzip.open(f)
    k = itertools.izip(*[iter(f)]*2)
    L = []
    locus = 0
    while 1:
        try:
            first = k.next()
        except StopIteration: break
        itera = [first[0],first[1]]
        fname = first[0]
        S = []
        rights = []
        lefts = []
        leftjust = rightjust = None
        while itera[0] != "//\n":
            nreps = int(itera[0].strip().split(";")[1].replace("size=",""))

            " record left and right most for cutting if gbs merge data "
            if datatype in ['mergegbs','gbs']:
                if itera[0].strip().split(";")[-1] == "":
                    leftjust = itera[1].index([i for i in itera[1] if i not in list("-N")][0])
                    rightjust = itera[1].rindex([i for i in itera[1] if i not in list("-N")][0])
                lefts.append(itera[1].index([i for i in itera[1] if i not in list("-N")][0]))
                rights.append(itera[1].rindex([i for i in itera[1] if i not in list("-N")][0]))

            " append sequence * number of dereps "
            for i in range(nreps):
                S.append(tuple(itera[1].strip()))
            itera = k.next()

        " trim off overhang edges of gbs reads "
        if datatype in ['mergegbs','gbs']:
            if any([i < leftjust for i in lefts]):
                rightjust = min(rights)
            if any([i < rightjust for i in rights]):
                leftjust = max(lefts)

            #print "".join(list(S[0]))
            for s in range(len(S)):
                if rightjust:
                    S[s] = S[s][leftjust:rightjust+1]
                if leftjust:
                    S[s] = S[s][leftjust:rightjust+1] ## +1?
            #print "".join(list(S[0]))

        " trim off restriction sites from end/s "
        if datatype in ['mergegbs','mergeddrad','pairddrad','pairgbs','gbs']:
            for s in range(len(S)):
                S[s] = S[s][len(CUT1):-(len(CUT2)+1)]
        else:
            for s in range(len(S)):
                S[s] = S[s][len(CUT1):]
            
        if len(S) >= minsamp:
            " make list for each site in sequences "
            res = stack(S)
            " exclude sites with indels "
            L += [i[0] for i in res if i[2] == 0]      
            locus += 1
    return L





def optim(WORK,handle, minsamp, CUT1, CUT2, datatype, haplos):
    name = handle.split("/")[-1].replace(".clustS.gz","")
    D = consensus(handle, minsamp, CUT1, CUT2, datatype)
    P = makeP(D)
    Tab = table_c(D)
    del D
    #H,E = scipy.optimize.fmin(LL,x0,(P,Tab),maxiter=500,maxfun=200,ftol=0.0001,disp=False,full_output=False)
    if haplos == 1:
        x0 = [0.001]
        H = 0.
        E = scipy.optimize.fmin(LL_haploid,x0,(P,Tab),disp=False,full_output=False)
    else:
        x0 = [0.01,0.001]
        H,E = scipy.optimize.fmin(LL,x0,(P,Tab),disp=False,full_output=False)
    del Tab
    outfile = open(WORK+"stats/."+name+".temp",'w')
    outfile.write("\t".join([name.strip(".gz"),str(round(H,8))[0:10],str(round(E,8))[0:10],"\n"]))
    outfile.close()
    sys.stderr.write(".")




def main(Parallel,ID,minsamp,subset,haplos,WORK,CUT,datatype):
    sys.stderr.write("\n\tstep 4: estimating error rate and heterozygosity\n\t")

    " find clust.xx directory "
    if not os.path.exists(WORK+'clust'+ID):
        print  "\terror: could not find "+WORK+"clust"+str(ID)+"/ directory,"+ \
                "\n\t\tif you changed the clustering threshold you must transfer *.clustS"+ \
                "\n\t\tfiles to a new directory named clust.xx with xx replaced by new clustering threshold"
        sys.exit()


    # warning message for low minsamp
    if minsamp < 5:
        sys.stderr.write("""\n\t warning: Mindepth < 5 is not recommended for this step.\n
                            If you intend to make low coverage base calls use a high mindepth in
                            step 4 to accurately infer H & E parameters, and then use a low mindepth
                            in conjunction with the line 31 params file option to make low coverage
                            base calls""")
        
    # if haploid data
    if haplos == 1:
        sys.stderr.write("""\n\tapplying haploid-based test (infer E while H is fixed to 0)""")

    # if double digest use first cut site
    if "," in CUT:
        CUT1, CUT2 = CUT.strip().split(",")
    else:
        CUT1 = CUT2 = CUT

    # load up work queue
    work_queue = multiprocessing.Queue()

    # iterate over files
    HH = glob.glob(WORK+"clust"+ID+"/"+subset+"*.clustS*")
    submitted = 0
    FS = []
    if len(HH) > 1:
        ## sort files by size
        for i in range(len(HH)):
            statinfo = os.stat(HH[i])
            if statinfo.st_size > 1000:
                FS.append((HH[i],statinfo.st_size))
            else:
                print "excluding ",HH[i],"file is too small\n"
        FS.sort(key=lambda x: x[1])
        FS = [i[0] for i in FS]
    else:
        FS = HH
    REMOVE = glob.glob(WORK+'clust'+ID+"/cat.*")
    FS = [f for f in FS if f not in REMOVE]
    for handle in FS:
        work_queue.put([WORK,handle, minsamp, CUT1, CUT2, datatype, haplos])
        submitted += 1

    " remove temp files if previous run "
    for ff in FS:
        end = ff.split("/")[-1].replace(".clustS.gz","") 
        ff = WORK+"stats/."+end+".temp"
        if os.path.exists(ff):
            os.remove(ff)

    " create a queue to pass to workers to store the results "
    result_queue = multiprocessing.Queue()
    results = []
    
    " spawn workers "
    jobs = []
    for i in range( min(Parallel,submitted) ):
        worker = Worker(work_queue, result_queue, optim)
        worker.start()
        jobs.append(worker)
    for job in jobs:
        job.join()

    " write results to stats file "
    if not os.path.exists(WORK+"stats/Pi_E_estimate.txt"):
        outstats = open(WORK+"stats/Pi_E_estimate.txt",'w')
        outstats.write("taxa\tH\tE\n")
    else:
        outstats = open(WORK+"stats/Pi_E_estimate.txt",'a')
    for ff in FS:
        end = ff.split("/")[-1].replace(".clustS.gz","")
        ft = WORK+"stats/."+end+".temp"
        line = open(ft).readlines()
        outstats.write(line[0])
        os.remove(ft)
        # n,h,e = line[0].strip().split("\t")
        # H.append(float(h))
        # E.append(float(e))
    #outstats.write(" ".join(["mean E =",str(numpy.mean(E))])+"\n")
    #outstats.write(" ".join(["mean H =",str(numpy.mean(H))]))
    outstats.close()
    


