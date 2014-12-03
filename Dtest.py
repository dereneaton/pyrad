#!/usr/bin/env python2

import os
import numpy
import sys
import random
import itertools
import glob
import multiprocessing
import cPickle as pickle
from potpour import Worker


def IUPAC(one):
    """
    returns IUPAC symbol for ambiguity bases,
    used for polymorphic sites.
    """
    D = {"R":['G','A'],
         "K":['G','T'],
         "S":['G','C'],
         "Y":['T','C'],
         "W":['T','A'],
         "M":['C','A']}
    return D[one]


def makefreq(patlist):
    " identify which allele is derived in P3 relative to outgroup "
    " and is the most frequent and use that as the SNP "
    P = {}
    for tax in patlist:
        P[tax] = []

    for tax in patlist:
        for base in patlist[tax]:
            if base in list('ATGC'):
                P[tax].append(base[0])
                P[tax].append(base[0])
            elif base in list("RKSYWM"):
                hh = IUPAC(base[0])
                for i in hh:
                    P[tax].append(i)

    major = [i for i in set(P['p3']) if i not in set(P['o'])]
    " in case of multiple bases "
    if len(major) > 1:
        cc = [P['p3'].count(base) for base in major]
        major = major[cc.index(max(cc))]  ## maybe [0]
    elif not major:
        major = [i for i in set(P['o']) if i in set(P['p3'])]
    elif len(major) == 1:
        major = major[0]

    ret = [float(P[i].count(major))/len(P[i]) for i in ['p1','p2','p3','o']]
    return ret


def Dstat(Loc, pat):
    if pat[0] != pat[1]:
        if pat[0] == pat[3]:
            if pat[1] == pat[2]:
                Loc.abba += 1
        else:
            if pat[0] == pat[2]:
                if pat[1] == pat[3]:
                    Loc.baba += 1
    return Loc


def polyDstat(Loc,patlist):
    ## calculate frequencies
    " look at the P3 taxon first for a derived allele "
    p1,p2,p3,o = makefreq(patlist)   #[a0,a1,a2,a3])
    Loc.abba += ((1.-p1)*p2*p3*(1.-o))
    Loc.baba += (p1*(1.-p2)*p3*(1.-o))
    return Loc
    

def fillin(ll,name,col,ulnames,patlist):
    if len(ll)>1:
        for i in ll:
            patlist[name] = col[ [ulnames.index(i) for i in ll if i in ulnames] ]
    else:
        patlist[name] = col[ [ulnames.index(i) for i in ll ] ]
    return patlist



def IUA(Loc,L):
    Loc.abba = 0.
    Loc.baba = 0.
    for col in Loc.seq.transpose():
        if all(i in list("ATGC") for i in col):
            Loc = Dstat(Loc,col)
    return Loc


def IUAfreq(Loc,L):
    patlist = {}
    Loc.abba = 0.
    Loc.baba = 0.
    for col in Loc.seq.transpose():
        patlist = fillin(L[0],'p1',col,Loc.names,patlist)
        patlist = fillin(L[1],'p2',col,Loc.names,patlist)
        patlist = fillin(L[2],'p3',col,Loc.names,patlist)
        patlist = fillin(L[3],'o', col,Loc.names,patlist)

        if not any([all([i in ["N",'-'] for i in patlist['p1']]),
                    all([i in ["N",'-'] for i in patlist['p2']]),
                    all([i in ["N",'-'] for i in patlist['p3']]),
                    all([i in ["N",'-'] for i in patlist['o']])]):
            if any([i not in patlist['o'] for i in patlist['p3']]):
                Loc = polyDstat(Loc,patlist)
            else:
                None
        else:
            None
    return Loc



def sample_wr(population, k):         
    "used for bootstrap sampling"
    "Chooses k random elements (with replacement) from a population"
    n = len(population)
    _random, _int = random.random, int  # speed hack
    return [_int(_random() * n) for i in itertools.repeat(None, k)]


def bootfreq(Ldict, which):
    Dftop = Dfbot = 0
    while 1:
        try: Lx = Ldict[Ldict.keys()[which.next()]]
        except StopIteration: break
        Dftop += Lx.abba - Lx.baba
        Dfbot += Lx.abba + Lx.baba
    D = 0.
    if Dfbot > 0:
        D = Dftop/float(Dfbot)
    return D
    

def bootfixed(Ldict, which):
    abba = baba = 0
    while 1:
        try: Lx = Ldict[Ldict.keys()[which.next()]]
        except StopIteration: break
        abba += Lx.abba
        baba += Lx.baba
    D = 0.
    if abba+baba > 0:
        D = float(abba-baba)/float(abba+baba)
    return D


def makeSNP(L,snpfreq,loci):
    Ndict = {}
    num = 0
    for loc in loci:
        Loc = Locus()
        Loc.number = num
        " only select loci that have data for all four tiptaxa "
        names = [i.split()[0].replace(">","") for i in loc.lstrip().rstrip().split("\n")[:-1]]
        if snpfreq:
            Loc.names = [i for i in names if i in list(itertools.chain(*L))]
        else:
            Loc.names = L # [i for i in names if i in L]

        " if snpfreq only need one of possibly multiple individuals "
        keep = 0

        if snpfreq:
            for tax in L:
                z = any([tax in names for tax in L[0]])
                y = any([tax in names for tax in L[1]])
                w = any([tax in names for tax in L[2]])
                u = any([tax in names for tax in L[3]])
            if all([z,y,w,u]):
                keep = 1
        else:
            if all(tax in names for tax in Loc.names):
                keep = 1

        if keep:
            N = numpy.array([tuple(i) for i in loc.split("\n")[1:]])
            " only select sites with synapomorphies "
            ## may want to keep autapomorphies in the future, or more
            ## when making a parameterized version of D-statistic
            ## only pyrad 2.1+ finds synapormorphies btwn hetero and fixed sites
            N[-1] = list(N[-1].tostring().replace("-","*"))
            N = N[:, N[-1] == "*"]

            " only select rows with focal taxa"
            if snpfreq:
                Loc.seq = N[[names.index(i) for i in Loc.names],:]
            else:
                Loc.seq = N[[names.index(i) for i in Loc.names],:]
            Ndict[num] = Loc
        num += 1
    return Ndict



class Locus():
    """locus keeps track of position in input file,
    variable sites, and D-statistics"""
    def _init_(self):
        self.number = number
        self.names = names
        self.seq = sequence
        self.abba = abba
        self.baba = baba
    def D(self):
        """ just used to check if abba > baba
        not a global genomic measure of D """
        if self.abba+self.baba > 0:
            return float(self.abba-self.baba)/(self.abba+self.baba)
        else:
            return 0.0
    


def runtest(infile, L, nboots, snpfreq, submitted):
    " print test"
    print L

    " split each locus "
    loci = open(infile).read().strip().split("|")[:-1]
    loci[0] = "\n"+loci[0]

    " returns a {} of Locus objects with data from tiptaxa L"
    Ldict = makeSNP(L,snpfreq,loci)

    " calculate ABBA/BABA for each locus"
    for loc in Ldict:
        if snpfreq:
            Ldict[loc] = IUAfreq(Ldict[loc],L)
        else:
            Ldict[loc] = IUA(Ldict[loc],L)

    " calculate final D "
    dftfinal = sum([Ldict[l].abba-Ldict[l].baba for l in Ldict])
    dbtfinal = sum([Ldict[l].abba+Ldict[l].baba for l in Ldict])
    if dbtfinal > 0:
        Dfinal = float(dftfinal)/dbtfinal
    else:
        Dfinal = 0.

    " proportion of discordant loci "
    try: pdisc = len([i for i in Ldict if Ldict[i].D()]) / float(len(Ldict))
    except ZeroDivisionError:
        pdisc = 0.0

    " do bootstrapping "
    BB = []
    for i in xrange(nboots):
        which = iter(sample_wr(xrange(len(Ldict)),len(Ldict)))
        if snpfreq:
            bb = bootfreq(Ldict, which)
        else:
            bb = bootfixed(Ldict, which)
        BB.append(bb)
    STD = numpy.std(BB)

    " out stats "
    if STD < 0.00001:
        STD = 0.0
    if Dfinal != 0.0:
        if STD != 0.0:
            Z = (abs(Dfinal/STD))
        else:
            Z = 0.0
    else:
        Dfinal = 0.
        Z = 0.

    ABBAloci = [Ldict[l].number for l in Ldict if Ldict[l].D() > 0]
    BABAloci = [Ldict[l].number for l in Ldict if Ldict[l].D() < 0]

    ret = [L,Dfinal,STD,Z,
           len(Ldict),
           sum([Ldict[l].abba for l in Ldict]),
           sum([Ldict[l].baba for l in Ldict]),
           pdisc,submitted,
           ABBAloci,BABAloci, BB ]
    pickle.dump(ret, open(".save.D4temp"+str(submitted),'wb'))




def makesortfiles(outn,locifile,n,loci,outfile,makesort,sub,ps):
    locifile.sort()
    "write to ABBA file all loci indexed in ABBAloci list"
    with open(outfile+"_"+str(sub+1)+"."+outn[0:n]+".txt",'w') as out:
        print >>out, " ".join(ps)
        print >>out, "//"
        print >>out, ",".join(map(str,locifile))
        print >>out, "//"
        if makesort == 2:
            for loc in xrange(len(loci)):
                if loc in locifile:
                    out.write(loci[loc]+"| locus: "+str(loc))





def checktaxa(taxalist,alignfile):
    with open(alignfile) as infile:
        data = infile.readlines()
    taxainfile = set()
    for line in data:
        if ">" in line:
            tax = line.split(" ")[0].replace(">","")
            if tax not in taxainfile:
                taxainfile.add(tax)
    if not set(taxalist).difference(taxainfile):
        return 1





def multiproc_it(tests, alignfile, outfile, nboots, nproc, namelen, makesort, makeboots):

    " submit jobs to processors "
    work_queue = multiprocessing.Queue()
    result_queue = multiprocessing.Queue()
    submitted = 0
    Notes = []
    for rep in tests:
        notes = ""
        if len(rep) == 2:
            rep,notes = rep
        p1,p2,p3,o = rep
        if any(["[" in i for i in rep]):
            p1 = p1[1:-1].split(",")
            p2 = p2[1:-1].split(",")
            p3 = p3[1:-1].split(",")
            o =   o[1:-1].split(",")
            taxalist = list(itertools.chain(*[p1+p2+p3+o]))
            if checktaxa(taxalist,alignfile):
                work_queue.put([alignfile,[p1,p2,p3,o],nboots,1, submitted])
                submitted += 1
            else: 
                print 'a taxon name was found that is not in the sequence file'
        else:
            if checktaxa([p1,p2,p3,o],alignfile):
                work_queue.put([alignfile,[p1,p2,p3,o],nboots,0, submitted])
                submitted += 1
            else: 
                print 'a taxon name was found that is not in the sequence file'

        Notes.append(notes)
    jobs = []
    for i in range(nproc):
        worker = Worker(work_queue, result_queue, runtest)
        jobs.append(worker)
        worker.start()
    for j in jobs:
        j.join()

    " read results back in "
    #Results = [result_queue.get() for i in range(submitted)]
    Results = [pickle.load(open(".save.D4temp"+str(i),'rb')) for i in xrange(submitted)]
    Results.sort(key = lambda x:x[8])

    "setup results file "
    outs = open(outfile+".D4.txt", 'w')
    header = "\t".join([ 'P1'+" "*(namelen[0]-2),
                         'P2'+" "*(namelen[1]-2),
                         'P3'+" "*(namelen[2]-2),
                         'O'+" "*(namelen[3]-1),
                         'D','std(D)','Z',
                         'BABA','ABBA',
                         'nloci','nboot','pdisc', 'notes'])
    print >>outs, header

    for i in range(len(Results)):
        ps,D,STD,Z,nloci,ABBA,BABA,pdisc,sub,ABBAloci,BABAloci,boots = Results[i]
        ps = [str(x).replace("['","[").replace("']","]").replace("', '",",").replace(">","") for x in ps]
        print >>outs, "%s\t%s\t%s\t%s\t%.3f\t%.3f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%.2f\t%s" % (ps[0]+" "*(namelen[0]-len(ps[0])),
                                                                                          ps[1]+" "*(namelen[1]-len(ps[1])),
                                                                                          ps[2]+" "*(namelen[2]-len(ps[2])),
                                                                                          ps[3]+" "*(namelen[3]-len(ps[3])),
                                                                                          D,STD,Z,
                                                                                          BABA,ABBA,
                                                                                          nloci,nboots,
                                                                                          pdisc,Notes[i])



        loci = open(alignfile).read().strip().split("|")[:-1]
        if makesort:
            makesortfiles('ABBA',ABBAloci,4,loci,outfile,makesort,sub,ps)
            makesortfiles('BABA',BABAloci,4,loci,outfile,makesort,sub,ps)            

        if makeboots:
            with open(outfile+"_"+str(sub+1)+".boots",'w') as out:
                out.write(",".join(map(str,boots)))

    for oldpickle in glob.glob(".save.D4temp*"):
        os.remove(oldpickle)
                

def main(tests, alignfile, outfile, nboots, nproc, makesort, makeboots):
    
    P1namelen  = max(map(len,[str(i[0][0]) for i in tests]))
    P2namelen  = max(map(len,[str(i[0][1]) for i in tests]))
    P3namelen  = max(map(len,[str(i[0][2]) for i in tests]))
    Onamelen   = max(map(len,[str(i[0][3]).strip() for i in tests]))
    namelen = [P1namelen,P2namelen,P3namelen,Onamelen]

    multiproc_it(tests,alignfile,outfile,nboots,nproc,namelen,makesort,makeboots)



if __name__ == '__main__':
    main()

