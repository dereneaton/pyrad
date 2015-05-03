#!/usr/bin/env python2

import numpy
import sys
import random
import itertools
import multiprocessing
import cPickle as pickle
from potpour import Worker
from Dtest import IUPAC, sample_wr, fillin, makesortfiles


def most_common(lst):
    return max(set(lst), key=lst.count)


def makefreq(patlist):
    " identify which allele is derived in P3 relative to outgroup "
    " and is the most frequent and use that as the SNP."
    " Also, split up alleles into those that are P3a & vs. or P3b "
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

    " select most common element in outgroup "
    if len(set(P['o'])) > 1:
        minor = most_common(P['o'])
    else:
        minor = P['o'][0]
        
    " select most common element that is not minor "
    bases = list(itertools.chain(*P.values()))
    majors = [i for i in bases if i != minor]
    major = most_common(majors)

    ret = [float(P[i].count(major)) / len(P[i]) for i in ['p1','p2','p3a','p3b','o']]
    ret += [float(P['p3a'].count(major)+P['p3b'].count(major))/(len(P['p3a'])+len(P['p3b']))]
    return ret



def Dstat5(Loc,pat):
    " check site for patterns and add to Locus object if found" 
    if len(set(pat)) < 3:
        if len(set(pat[2:])) > 1:
            minor = pat[-1]
            major = [i for i in [pat[2],pat[3]] if i not in pat[4]][0] 

            o = 0.
            p3ab = 1. if (pat[3] == major) & (pat[2] == major) else 0.
            p3b = 1. if pat[3] == major else 0.
            p3a = 1. if pat[2] == major else 0.
            p2 = 1. if pat[1] == major else 0.
            p1 = 1. if pat[0] == major else 0.

            Loc.abbba += ( (1.-p1)*p2*p3ab*(1.-o) )
            Loc.babba += ( p1*(1.-p2)*p3ab*(1.-o) )

            Loc.abbaa += ( (1.-p1)*p2*p3a*(1.-p3b)*(1.-o) )
            Loc.babaa += ( p1*(1.-p2)*p3a*(1.-p3b)*(1.-o) )

            Loc.ababa += ( (1.-p1)*p2*(1.-p3a)*p3b*(1.-o) )
            Loc.baaba += ( p1*(1.-p2)*(1.-p3a)*p3b*(1.-o) )

    return Loc



def polyDstat5(Loc, pat):
    ## calculate frequencies
    " look at the P3 taxon first for a derived allele "
    p1,p2,p3a,p3b,o,p3ab = makefreq(pat)

    Loc.abbba += ( (1.-p1)*p2*p3ab*(1.-o) )
    Loc.babba += ( p1*(1.-p2)*p3ab*(1.-o) )
    
    Loc.abbaa += ( (1.-p1)*p2*p3a*(1.-p3b)*(1.-o) )
    Loc.babaa += ( p1*(1.-p2)*p3a*(1.-p3b)*(1.-o) )

    Loc.ababa += ( (1.-p1)*p2*(1.-p3a)*p3b*(1.-o) )
    Loc.baaba += ( p1*(1.-p2)*(1.-p3a)*p3b*(1.-o) )

    return Loc



def IUAfreq(Loc, L):
    patlist = {}
    Loc.abbba = 0.
    Loc.babba = 0.
    Loc.abbaa = 0.
    Loc.babaa = 0.
    Loc.ababa = 0.
    Loc.baaba = 0.

    for col in Loc.seq.transpose():
        patlist = fillin(L[0], 'p1', col, Loc.names, patlist)
        patlist = fillin(L[1], 'p2', col, Loc.names, patlist)
        patlist = fillin(L[2], 'p3a', col, Loc.names, patlist)
        patlist = fillin(L[3], 'p3b', col, Loc.names, patlist)
        patlist = fillin(L[4], 'o', col, Loc.names, patlist)

        if not any([ all([i in ["N",'-'] for i in patlist['p1']]),
                     all([i in ["N",'-'] for i in patlist['p2']]),
                     all([i in ["N",'-'] for i in patlist['p3a']]),
                     all([i in ["N",'-'] for i in patlist['p3b']]),
                     all([i in ["N",'-'] for i in patlist['o']]) ]):
            if any([ i not in patlist['o'] for i in numpy.dstack((patlist['p3a'],patlist['p3b']))[0][0] ]):
                Loc = polyDstat5(Loc, patlist)
            else:
                None
        else:
            None
    return Loc



def IUA(Loc,L):
    Loc.abbba = 0.
    Loc.babba = 0.
    Loc.abbaa = 0.
    Loc.babaa = 0.
    Loc.ababa = 0.
    Loc.baaba = 0.
    for col in Loc.seq.transpose():
        if all(i in list("ATGC") for i in col):
            Loc = Dstat5(Loc,col)
    return Loc



def bootfreq(Ldict, which):
    Dft_12 = Dfb_12 = 0
    Dft_1  = Dfb_1  = 0
    Dft_2  = Dfb_2  = 0
    while 1:
        try: Lx = Ldict[Ldict.keys()[which.next()]]
        except StopIteration: break
        Dft_12 += Lx.abbba - Lx.babba
        Dfb_12 += Lx.abbba + Lx.babba
        Dft_1  += Lx.abbaa - Lx.babaa
        Dfb_1  += Lx.abbaa + Lx.babaa
        Dft_2  += Lx.ababa - Lx.baaba
        Dfb_2  += Lx.ababa + Lx.baaba
    D12 = 0.
    if Dfb_12 > 0:
        D12 = Dft_12/float(Dfb_12)
    D1 = 0.
    if Dfb_1 > 0:
        D1 = Dft_1/float(Dfb_1)
    D2 = 0.
    if Dfb_2 > 0:
        D2 = Dft_2/float(Dfb_2)
    return D12, D1, D2

    

def bootfixed(Ldict, which):
    abbba = babba = 0
    abbaa = babaa = 0
    ababa = baaba = 0
    while 1:
        try: Lx = Ldict[Ldict.keys()[which.next()]]
        except StopIteration: break
        abbba += Lx.abbba
        babba += Lx.babba
        abbaa += Lx.abbaa
        babaa += Lx.babaa
        ababa += Lx.ababa
        baaba += Lx.baaba
    D12 = 0.
    if abbba + babba > 0:
        D12 = float(abbba-babba)/(abbba+babba)
    D1 = 0.
    if abbaa + babaa > 0:
        D1 = float(abbaa-babaa)/(abbaa+babaa)
    D2 = 0.
    if ababa + baaba > 0:
        D2 = float(ababa-baaba)/(ababa+baaba)
    return D12, D1, D2



class Locus5():
    """locus keeps track of position in input file,
    variable sites, and D-statistics"""
    def _init_(self):
        self.number = number
        self.taxa = names
        self.seq = sequence
        self.abbba = abbba
        self.babba = abbba
        self.abbaa = abbaa
        self.babaa = babaa
        self.ababa = ababa
        self.baaba = baaba
    " D-stats for an individual locus "
    def D1(self):
        if self.abbaa+self.babaa > 0:
            return float(self.abbaa-self.babaa)/(self.abbaa+self.babaa)
        else:
            return 0.0
    def D2(self):
        if self.ababa+self.baaba > 0:
            return float(self.ababa-self.baaba)/(self.ababa+self.baaba)
        else:
            return 0.0
    def D12(self):
        if self.abbba+self.babba > 0:
            return float(self.abbba-self.babba)/(self.abbba+self.babba)
        else:
            return 0.0






def makeSNP(L, snpfreq, loci):
    Ndict = {}
    num = 0
    for loc in loci:
        Loc = Locus5()
        Loc.number = num

        " only select loci that have data for all five tiptaxa "
        names = [i.split()[0].replace(">","") for i in loc.lstrip().rstrip().split("\n")[:-1]]
        if snpfreq:
            Loc.names = [i for i in names if i in list(itertools.chain(*L))]
        else:
            Loc.names = L #[i for i in names if i in L]

        " if snpfreq only need one of possibly multiple individuals"
        keep = 0

        if snpfreq:
            for tax in L:
                z = any([tax in Loc.names for tax in L[0]])
                y = any([tax in Loc.names for tax in L[1]])
                x = any([tax in Loc.names for tax in L[2]])
                w = any([tax in Loc.names for tax in L[3]])
                u = any([tax in Loc.names for tax in L[4]])
            if all([z,y,x,w,u]):
                keep = 1

        else:
            if all(tax in names for tax in Loc.names):
                keep = 1

        if keep:
            N = numpy.array([tuple(i) for i in loc.split("\n")[1:]])
            " only select sites with synapomorphies "
            # select all variable
            N[-1] = list(N[-1].tostring().replace("-","*"))
            N = N[:, N[-1] == "*"]

            " only select rows with focal taxa "
            Loc.seq = N[[names.index(i) for i in Loc.names],:]
            Ndict[num] = Loc
        num += 1
    return Ndict
                         
        

def runtest(infile, L, nboots, snpfreq, submitted):
    " print test "
    print L

    " split each locus "
    loci = open(infile).read().strip().split("|")[:-1]
    loci[0] = "\n"+loci[0]

    " returns a {} of Locus5 objects with data for tiptaxa L "
    Ldict = makeSNP(L, snpfreq, loci)

    " calculate discordant patterns for each locus "
    for loc in Ldict:
        if snpfreq:
            Ldict[loc] = IUAfreq(Ldict[loc],L)
        else:
            Ldict[loc] = IUA(Ldict[loc],L)
    ################################################

    " final D12 "
    dft_12 = sum([Ldict[l].abbba - Ldict[l].babba for l in Ldict])
    dbt_12 = sum([Ldict[l].abbba + Ldict[l].babba for l in Ldict])
    if dbt_12 > 0:
        D12 = float(dft_12)/dbt_12
    else: D12 = 0.

    " final D1 "
    dft_1 = sum([Ldict[l].abbaa - Ldict[l].babaa for l in Ldict])
    dbt_1 = sum([Ldict[l].abbaa + Ldict[l].babaa for l in Ldict])
    if dbt_1 > 0:
        D1 = float(dft_1)/dbt_1
    else: D1 = 0.

    " final D2 "
    dft_2 = sum([Ldict[l].ababa - Ldict[l].baaba for l in Ldict])
    dbt_2 = sum([Ldict[l].ababa + Ldict[l].baaba for l in Ldict])
    if dbt_2 > 0:
        D2 = float(dft_2)/dbt_2
    else: D2 = 0.

    " proportion of discordant loci "
    try: pdisc = len([i for i in Ldict if any([Ldict[i].D12(),Ldict[i].D1(),Ldict[i].D2()])]) / float(len(Ldict))
    except ValueError:
        pdisc = 0.0
    
    #################################################

    " do bootstrapping "
    BB12 = []
    BB1  = []
    BB2  = []
    for i in xrange(nboots):
        which = iter(sample_wr(xrange(len(Ldict)), len(Ldict)))
        if snpfreq:
            bb12,bb1,bb2 = bootfreq(Ldict, which)
        else:
            #bb12,bb1,bb2 = bootfixed(Ldict, which)
            bb12,bb1,bb2 = bootfreq(Ldict, which)
        BB12.append(bb12)
        BB1.append(bb1)
        BB2.append(bb2)
    STD12 = numpy.std(BB12)
    STD1  = numpy.std(BB1)
    STD2  = numpy.std(BB2)
    ##################################################

    " stats out "
    if STD12 > 0:
        Z12 = (abs(D12/STD12))
    else: Z12 = 0.
    if STD1 > 0:
        Z1 =  (abs(D1/STD1))
    else: Z1 = 0.
    if STD2 > 0:
        Z2 =  (abs(D2/STD2))
    else: Z2 = 0.

    ## make loci files here
    ABBBAloci = [Ldict[l].number for l in Ldict if Ldict[l].D12() > 0]
    BABBAloci = [Ldict[l].number for l in Ldict if Ldict[l].D12() < 0]
    ABBAAloci = [Ldict[l].number for l in Ldict if Ldict[l].D1() > 0]
    BABAAloci = [Ldict[l].number for l in Ldict if Ldict[l].D1() < 0]
    ABABAloci = [Ldict[l].number for l in Ldict if Ldict[l].D2() > 0]
    BAABAloci = [Ldict[l].number for l in Ldict if Ldict[l].D2() < 0]

    " pickle to prevent multiprocessing from freezing on large returns "
    ret = [L,
           D12,Z12,
           D1,Z1,
           D2,Z2,
           len(Ldict),
           sum([Ldict[l].abbba for l in Ldict]),
           sum([Ldict[l].babba for l in Ldict]),
           sum([Ldict[l].abbaa for l in Ldict]),
           sum([Ldict[l].babaa for l in Ldict]),
           sum([Ldict[l].ababa for l in Ldict]),
           sum([Ldict[l].baaba for l in Ldict]),
           pdisc, submitted,
           ABBBAloci, BABBAloci,
           ABBAAloci, BABAAloci,
           ABABAloci, BAABAloci,
           BB12, BB1, BB2]
    pickle.dump(ret, open(".save."+str(submitted),'wb'))



def checktaxa(taxalist,alignfile):
    with open(alignfile) as infile:
        data = infile.readlines()
    taxainfile = set()
    for line in data:
        if ">" in line:
            tax = line.split(" ")[0].replace(">","")
            if tax not in taxainfile:
                taxainfile.add(tax)
    if not set(taxainfile).difference(taxainfile):
        return 1




def multiproc_it(subtests,alignfile,outfile, nboots,nproc,namelen,makesort,makeboots):
    work_queue = multiprocessing.Queue()
    result_queue = multiprocessing.Queue()
    submitted = 0
    Notes = []
    for rep in subtests:
        notes = ""
        if len(rep) == 2:
            rep,notes = rep
        p1,p2,p3a,p3b,o = rep
        if all(["[" in i for i in rep[1:]]):
            p1  = p1[1:-1].split(",")
            p2  = p2[1:-1].split(",")
            p3a = p3a[1:-1].split(",")
            p3b = p3b[1:-1].split(",")
            o   = o[1:-1].split(",")
            if checktaxa([p1,p2,p3a,p3b,o],alignfile):
                work_queue.put([alignfile, [p1,p2,p3a,p3b,o], nboots, 1, submitted])
                submitted += 1
            else: 
                print 'a taxon name was found that is not in the sequence file'
        else:
            if checktaxa([p1,p2,p3a,p3b,o],alignfile):
                work_queue.put([alignfile, [p1,p2,p3a,p3b,o], nboots, 0, submitted])
                submitted += 1
            else: 
                print 'a taxon name was found that is not in the sequence file'
        Notes.append(notes)

    jobs = []
    for i in range(min(submitted,nproc)):
        worker = Worker(work_queue, result_queue, runtest)
        jobs.append(worker)
        worker.start()
    for j in jobs:
        j.join()


    " read results back in "
    #Results = [result_queue.get() for i in range(submitted)]
    Results = [pickle.load(open(".save."+str(i),'rb')) for i in range(submitted)]
    Results.sort(key = lambda x:x[15])


    " setup results file "
    outs = open(outfile+".partD.txt", 'w')
    header = "\t".join([ 'p1'+" "*(namelen[0]-2),
                         'p2'+" "*(namelen[1]-2),
                         'p3_1'+" "*(namelen[2]-4),
                         'p3_2'+" "*(namelen[3]-4),
                         'O'+" "*(namelen[4]-1),
                         'D_12','D_1','D_2',
                         'Z_12','Z_1','Z_2',
                         'BABBA','ABBBA',
                         'BABAA','ABBAA',
                         'BAABA','ABABA',
                         'nloci','pdisc', 'notes'])

    print >>outs, header


    for i in range(len(Results)):
        L,D12,Z12,D1,Z1,D2,Z2,nloc,ABBBA,BABBA,ABBAA,BABAA,ABABA,BAABA,pdisc,sub,ABBBAloci,BABBAloci,ABBAAloci,BABAAloci,ABABAloci,BAABAloci,BB12,BB1,BB2 = Results[i]
        L = [str(x).replace("['","[").replace("']","]").replace("', '",",") for x in L]

        resin = tuple([str(L[0])+" "*(namelen[0]-len(str(L[0]))),
                       str(L[1])+" "*(namelen[1]-len(str(L[1]))),
                       str(L[2])+" "*(namelen[2]-len(str(L[2]))),
                       str(L[3])+" "*(namelen[3]-len(str(L[3]))),
                       str(L[4])+" "*(namelen[4]-len(str(L[4]))),
                       D12, D1, D2, Z12, Z1, Z2, 
                       BABBA, ABBBA, BABAA, ABBAA, BAABA, ABABA,
                       nloc, pdisc, Notes[i]])
        
        print >>outs, "%s\t%s\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%.2f\t%s" % resin 

        loci = open(alignfile).read().strip().split("|")[:-1]
        if makesort:
            makesortfiles("ABBBA",ABBBAloci,5,loci,outfile,makesort,sub,L)
            makesortfiles("BABBA",BABBAloci,5,loci,outfile,makesort,sub,L)            
            makesortfiles("ABBAA",ABBAAloci,5,loci,outfile,makesort,sub,L)            
            makesortfiles("BABAA",BABAAloci,5,loci,outfile,makesort,sub,L)
            makesortfiles("ABABA",ABABAloci,5,loci,outfile,makesort,sub,L)
            makesortfiles("BAABA",BAABAloci,5,loci,outfile,makesort,sub,L)

        if makeboots:
            with open(outfile+"_"+str(sub+1)+".boots_D12",'w') as out:
                out.write(",".join(map(str,BB12)))
            with open(outfile+"_"+str(sub+1)+".boots_D1",'w') as out:
                out.write(",".join(map(str,BB1)))
            with open(outfile+"_"+str(sub+1)+".boots_D2",'w') as out:
                out.write(",".join(map(str,BB2)))
                

def main(tests, alignfile, outfile, nboots, nproc, makesort, makeboots):
    import sys

    P1namelen  = max(map(len,[str(i[0][0]) for i in tests]))
    P2namelen  = max(map(len,[str(i[0][1]) for i in tests]))
    P3anamelen = max(map(len,[str(i[0][2]) for i in tests]))
    P3bnamelen = max(map(len,[str(i[0][3]) for i in tests]))
    Onamelen   = max(map(len,[str(i[0][4]).strip() for i in tests]))
    namelen = [P1namelen,P2namelen,P3anamelen,P3bnamelen,Onamelen]
                
    multiproc_it(tests,alignfile,outfile,nboots,nproc,namelen,makesort,makeboots)


if __name__ == '__main__':
    main()




