import numpy
import sys
import itertools
import multiprocessing
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

    """ select most common element in outgroup if multiple individuals,
    if only one ind but two alleles, select the first one """
    if len(set(P['o'])) > 1:
        minor = most_common(P['o'])
    else:
        minor = P['o'][0]

    " select most common element that is not minor "
    bases = list(itertools.chain(*P.values()))
    majors = [i for i in bases if i != minor]
    major = most_common(majors)

    ret = [float(P[i].count(major)) / len(P[i]) for i in ['p1','p2','p3a','p3b','o']]
    return ret





def Dstatfoil(Loc,pat):
    " check site for patterns and add to Locus object if found" 
    if len(set(pat)) < 3:
        " only allow biallelic "
        minor = pat[-1]
        " select only alternative to the outgroup allele "
        major = [i for i in pat if i!= pat[-1]][0]

        o = 0.
        p3b = 1. if pat[3] == major else 0.
        p3a = 1. if pat[2] == major else 0.
        p2 = 1. if pat[1] == major else 0.
        p1 = 1. if pat[0] == major else 0.

        ## from partitioned D-stat
        Loc.abbba += ( (1.-p1)*p2*p3a*p3b*(1.-o) )                #                       DFI[5]   DOL[5]
        Loc.babba += ( p1*(1.-p2)*p3a*p3b*(1.-o) )                #                       DFI[1]   DOL[1]

        Loc.abbaa += ( (1.-p1)*p2*p3a*(1.-p3b)*(1.-o) )           #     DFO[6]   DIL[0]   DFI[4]   DOL[2]
        Loc.babaa += ( p1*(1.-p2)*p3a*(1.-p3b)*(1.-o) )           #     DFO[0]   DIL[6]   DFI[0]   DOL[6]

        Loc.ababa += ( (1.-p1)*p2*(1.-p3a)*p3b*(1.-o) )           #     DFO[2]   DIL[4]   DFI[2]   DOL[4]
        Loc.baaba += ( p1*(1.-p2)*(1.-p3a)*p3b*(1.-o) )           #     DFO[4]   DIL[2]   DFI[6]   DOL[0]

        ## new to foil, contrast of xxbba
        Loc.bbbaa += ( p1*p2*p3a*(1.-p3b)*(1.-o) )                #     DFO[1]   DIL[1]
        Loc.bbaba += ( p1*p2*(1.-p3a)*p3b*(1.-o) )                #     DFO[5]   DIL[5]

        ## terminal branch patterns
        if not Loc.noterminals:
            Loc.aaaba += ( (1.-p1)*(1.-p2)*(1.-p3a)*p3b*(1.-o) )      #     DFO[3]   DIL[3]
            Loc.aabaa += ( (1.-p1)*(1.-p2)*p3a*(1.-p3b)*(1.-o) )      #     DFO[7]   DIL[7]  
            Loc.abaaa += ( (1.-p1)*p2*(1.-p3a)*(1.-p3b)*(1.-o) )      #                       DFI[3]   DOL[3]
            Loc.baaaa += ( p1*(1.-p2)*(1.-p3a)*(1.-p3b)*(1.-o) )      #                       DFI[7]   DOL[7]
    return Loc



def polyDstatfoil(Loc, pat):
    ## calculate frequencies
    " look at the P3 taxon first for a derived allele "
    p1,p2,p3a,p3b,o = makefreq(pat)
    # else:
    #     pat = [1. if base!=pat[-1] else 0. for base in pat]
    #     p1,p2,p3a,p3b,o = pat

    ## from partitioned D-stat
    Loc.abbba += ( (1.-p1)*p2*p3a*p3b*(1.-o) )                #                       DFI[5]   DOL[5]
    Loc.babba += ( p1*(1.-p2)*p3a*p3b*(1.-o) )                #                       DFI[1]   DOL[1]

    Loc.abbaa += ( (1.-p1)*p2*p3a*(1.-p3b)*(1.-o) )           #     DFO[6]   DIL[0]   DFI[4]   DOL[2]
    Loc.babaa += ( p1*(1.-p2)*p3a*(1.-p3b)*(1.-o) )           #     DFO[0]   DIL[6]   DFI[0]   DOL[6]

    Loc.ababa += ( (1.-p1)*p2*(1.-p3a)*p3b*(1.-o) )           #     DFO[2]   DIL[4]   DFI[2]   DOL[4]
    Loc.baaba += ( p1*(1.-p2)*(1.-p3a)*p3b*(1.-o) )           #     DFO[4]   DIL[2]   DFI[6]   DOL[0]

    ## new to foil, contrast of xxbba
    Loc.bbbaa += ( p1*p2*p3a*(1.-p3b)*(1.-o) )                #     DFO[1]   DIL[1]
    Loc.bbaba += ( p1*p2*(1-p3a)*p3b*(1.-o) )                 #     DFO[5]   DIL[5]
    
    ## terminal branch patterns
    if not Loc.noterminals:
        Loc.aaaba += ( (1.-p1)*(1.-p2)*(1.-p3a)*p3b*(1.-o) )      #     DFO[3]   DIL[3]
        Loc.aabaa += ( (1.-p1)*(1-p2)*p3a*(1.-p3b)*(1.-o) )       #     DFO[7]   DIL[7]  
        Loc.abaaa += ( (1.-p1)*p2*(1.-p3a)*(1.-p3b)*(1.-o) )      #                       DFI[3]   DOL[3]
        Loc.baaaa += ( p1*(1.-p2)*(1.-p3a)*(1.-p3b)*(1.-o) )      #                       DFI[7]   DOL[7]
    
    return Loc



def IUAfreq(Loc, L):
    patlist = {}
    Loc.abbba = 0.
    Loc.babba = 0.
    Loc.abbaa = 0.
    Loc.babaa = 0.
    Loc.ababa = 0.
    Loc.baaba = 0.

    Loc.bbbaa = 0.
    Loc.bbaba = 0.
    Loc.aaaba = 0.
    Loc.aabaa = 0.
    Loc.abaaa = 0.
    Loc.baaaa = 0.

    for col in Loc.seq.transpose():
        patlist = fillin(L[0], 'p1', col, Loc.names, patlist)
        patlist = fillin(L[1], 'p2', col, Loc.names, patlist)
        patlist = fillin(L[2], 'p3a', col, Loc.names, patlist)
        patlist = fillin(L[3], 'p3b', col, Loc.names, patlist)
        patlist = fillin(L[4], 'o', col, Loc.names, patlist)

        " exclude sites with missing data "
        if not any([ all([i in ["N",'-'] for i in patlist['p1']]),
                     all([i in ["N",'-'] for i in patlist['p2']]),
                     all([i in ["N",'-'] for i in patlist['p3a']]),
                     all([i in ["N",'-'] for i in patlist['p3b']]),
                     all([i in ["N",'-'] for i in patlist['o']]) ]):
            " if site in not invariable "
            if len(set(col)) > 1:
                " look for patterns in site "
                Loc = polyDstatfoil(Loc, patlist)
    return Loc



def IUA(Loc,L):
    Loc.abbba = 0.
    Loc.babba = 0.
    Loc.abbaa = 0.
    Loc.babaa = 0.
    Loc.ababa = 0.
    Loc.baaba = 0.

    Loc.bbbaa = 0.
    Loc.bbaba = 0.
    Loc.aaaba = 0.
    Loc.aabaa = 0.
    Loc.abaaa = 0.
    Loc.baaaa = 0.

    for col in Loc.seq.transpose():
        " exclude heterozygous sites "
        if all(i in list("ATGC") for i in col):
            " if site is not invariable "
            if len(set(col)) > 1:
                " look for patterns in site "
                Loc = Dstatfoil(Loc,col)
    return Loc



def bootfreq(Ldict, which):
    Dfo_t = Dfo_b = 0.
    Dil_t = Dil_b = 0.
    Dfi_t = Dfi_b = 0.
    Dol_t = Dol_b = 0.
    while 1:
        try: Lx = Ldict[Ldict.keys()[which.next()]]
        except StopIteration: break
        " iterate over loci summing top and bottom values of Ds"
        Dfo_t  += Lx.DFO_t()
        Dfo_b  += Lx.DFO_b()
        Dil_t  += Lx.DIL_t()
        Dil_b  += Lx.DIL_b()
        Dfi_t  += Lx.DFI_t()
        Dfi_b  += Lx.DFI_b()
        Dol_t  += Lx.DOL_t()
        Dol_b  += Lx.DOL_b()
    " take top over bottom values to calc Ds "
    DFO = 0.
    if Dfo_b > 0:
        DFO = Dfo_t/float(Dfo_b)
    DIL = 0.
    if Dil_b > 0:
        DIL = Dil_t/float(Dil_b)
    DFI = 0.
    if Dfi_b > 0:
        DFI = Dfi_t/float(Dfi_b)
    DOL = 0.
    if Dol_b > 0:
        DOL = Dol_t/float(Dol_b)

    return DFO,DIL,DFI,DOL



class Locusfoil():
    """locus keeps track of position in input file,
    variable sites, and D-statistics"""
    def _init_(self):
        self.number = number
        self.names = names
        self.seq = seq
        self.noterminals = noterminals
        
        self.abbba = abbba
        self.babba = abbba
        self.abbaa = abbaa
        self.babaa = babaa
        self.ababa = ababa
        self.baaba = baaba

        self.bbbaa = bbbaa
        self.bbaba = bbaba

        self.aaaba = aaaba
        self.aabaa = aabaa
        self.abaaa = abaaa
        self.baaaa = baaaa

    """ per-locus top or bottom values of Dstats """
    def DFO_t(self):
        part1 = [self.babaa,self.bbbaa,self.ababa,self.aaaba]
        part2 = [self.baaba,self.bbaba,self.abbaa,self.aabaa]
        if self.noterminals:
            part1 = part1[:-1]
            part2 = part2[:-1]
        return float(sum(part1)-sum(part2))

    def DFO_b(self):
        part1 = [self.babaa,self.bbbaa,self.ababa,self.aaaba]
        part2 = [self.baaba,self.bbaba,self.abbaa,self.aabaa]
        if self.noterminals:
            part1 = part1[:-1]
            part2 = part2[:-1]
        return float(sum(part1)+sum(part2))


    def DIL_t(self):
        part1 = [self.abbaa,self.bbbaa,self.baaba,self.aaaba]
        part2 = [self.ababa,self.bbaba,self.babaa,self.aabaa]
        if self.noterminals:
            part1 = part1[:-1]
            part2 = part2[:-1]
        return float(sum(part1)-sum(part2))


    def DIL_b(self):
        part1 = [self.abbaa,self.bbbaa,self.baaba,self.aaaba]
        part2 = [self.ababa,self.bbaba,self.babaa,self.aabaa]
        if self.noterminals:
            part1 = part1[:-1]
            part2 = part2[:-1]
        return float(sum(part1)+sum(part2))


    def DFI_t(self):
        part1 = [self.babaa,self.babba,self.ababa,self.abaaa]
        part2 = [self.abbaa,self.abbba,self.baaba,self.baaaa]
        if self.noterminals:
            part1 = part1[:-1]
            part2 = part2[:-1]
        return float(sum(part1)-sum(part2))


    def DFI_b(self):
        part1 = [self.babaa,self.babba,self.ababa,self.abaaa]
        part2 = [self.abbaa,self.abbba,self.baaba,self.baaaa]
        if self.noterminals:
            part1 = part1[:-1]
            part2 = part2[:-1]
        return float(sum(part1)+sum(part2))


    def DOL_t(self):
        part1 = [self.baaba,self.babba,self.abbaa,self.abaaa]
        part2 = [self.ababa,self.abbba,self.babaa,self.baaaa]
        if self.noterminals:
            part1 = part1[:-1]
            part2 = part2[:-1]
        return float(sum(part1)-sum(part2))


    def DOL_b(self):
        part1 = [self.baaba,self.babba,self.abbaa,self.abaaa]
        part2 = [self.ababa,self.abbba,self.babaa,self.baaaa]
        if self.noterminals:
            part1 = part1[:-1]
            part2 = part2[:-1]
        return float(sum(part1)+sum(part2))
        




def makeSNP(L, snpfreq, loci, noterminals):
    Ndict = {}
    num = 0
    for loc in loci:
        Loc = Locusfoil()
        Loc.noterminals = noterminals
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
            # select all variable sites
            N[-1] = list(N[-1].tostring().replace("-","*"))
            N = N[:, N[-1] == "*"]

            " only select rows with focal taxa "
            Loc.seq = N[[names.index(i) for i in Loc.names],:]
            Ndict[num] = Loc
        num += 1
    return Ndict
                         
        

def runtest(infile, L, nboots, snpfreq, submitted, noterminals):
    " print test "
    print L

    " split each locus "
    loci = open(infile).read().strip().split("|")[:-1]
    loci[0] = "\n"+loci[0]

    " returns a {} of Locusfoil objects with data for tiptaxa L "
    Ldict = makeSNP(L, snpfreq, loci, noterminals)

    " calculate discordant patterns for each locus "
    for loc in Ldict:
        if snpfreq:
            Ldict[loc] = IUAfreq(Ldict[loc],L)
        else:
            Ldict[loc] = IUA(Ldict[loc],L)
    ################################################

    " final DFO "
    DFO_t = sum([(Ldict[l].babaa + Ldict[l].bbbaa + Ldict[l].ababa + Ldict[l].aaaba) -\
                 (Ldict[l].baaba + Ldict[l].bbaba + Ldict[l].abbaa + Ldict[l].aabaa) for l in Ldict])
    DFO_b = sum([(Ldict[l].babaa + Ldict[l].bbbaa + Ldict[l].ababa + Ldict[l].aaaba) + \
                 (Ldict[l].baaba + Ldict[l].bbaba + Ldict[l].abbaa + Ldict[l].aabaa) for l in Ldict])
    if DFO_b > 0:
        DFO = float(DFO_t)/DFO_b
    else: DFO = 0.
    
    " final DIL "
    DIL_t = sum([(Ldict[l].abbaa + Ldict[l].bbbaa + Ldict[l].baaba + Ldict[l].aaaba) - \
                 (Ldict[l].ababa + Ldict[l].bbaba + Ldict[l].babaa + Ldict[l].aabaa) for l in Ldict])
    DIL_b = sum([(Ldict[l].abbaa + Ldict[l].bbbaa + Ldict[l].baaba + Ldict[l].aaaba) + \
                 (Ldict[l].ababa + Ldict[l].bbaba + Ldict[l].babaa + Ldict[l].aabaa) for l in Ldict])
    if DIL_b > 0:
        DIL = float(DIL_t)/DIL_b
    else: DIL = 0.

    " final DFI "
    DFI_t = sum([(Ldict[l].babaa + Ldict[l].babba + Ldict[l].ababa + Ldict[l].abaaa) - \
                 (Ldict[l].abbaa + Ldict[l].abbba + Ldict[l].baaba + Ldict[l].baaaa) for l in Ldict])
    DFI_b = sum([(Ldict[l].babaa + Ldict[l].babba + Ldict[l].ababa + Ldict[l].abaaa) + \
                 (Ldict[l].abbaa + Ldict[l].abbba + Ldict[l].baaba + Ldict[l].baaaa) for l in Ldict])
    if DFI_b > 0:
        DFI = float(DFI_t)/DFI_b
    else: DFI = 0.

    " final DOL "
    DOL_t = sum([(Ldict[l].baaba + Ldict[l].babba + Ldict[l].abbaa + Ldict[l].abaaa) - \
                 (Ldict[l].ababa + Ldict[l].abbba + Ldict[l].babaa + Ldict[l].baaaa) for l in Ldict])
    DOL_b = sum([(Ldict[l].baaba + Ldict[l].babba + Ldict[l].abbaa + Ldict[l].abaaa) + \
                 (Ldict[l].ababa + Ldict[l].abbba + Ldict[l].babaa + Ldict[l].baaaa) for l in Ldict])
    if DOL_b > 0:
        DOL = float(DOL_t)/DOL_b
    else: DOL = 0.

    " proportion of discordant loci "
    #try: pdisc = len([i for i in Ldict if any([Ldict[i].D12(),Ldict[i].D1(),Ldict[i].D2()])]) / float(len(Ldict))
    #except ValueError:
    #    pdisc = 0.0

    " TODO "
    pdisc = 0.0
    
    #################################################

    " do bootstrapping "
    BBFO = []
    BBIL = []
    BBFI = []
    BBOL = []
    for i in xrange(nboots):
        which = iter(sample_wr(xrange(len(Ldict)), len(Ldict)))
        bbfo,bbil,bbfi,bbol = bootfreq(Ldict, which)
        BBFO.append(bbfo)
        BBIL.append(bbil)
        BBFI.append(bbfi)
        BBOL.append(bbol)
    STDfo  = numpy.std(BBFO)
    STDil  = numpy.std(BBIL)
    STDfi  = numpy.std(BBFI)
    STDol  = numpy.std(BBOL)
    ##################################################

    " stats out "
    if STDfo > 0:
        ZFO = (abs(DFO/STDfo))
    else: ZFO = 0.
    if STDil > 0:
        ZIL =  (abs(DIL/STDil))
    else: ZIL = 0.
    if STDfi > 0:
        ZFI =  (abs(DFI/STDfi))
    else: ZFI = 0.
    if STDol > 0:
        ZOL =  (abs(DOL/STDol))
    else: ZOL = 0.

    ## make loci files here
    #ABBBAloci = [Ldict[l].number for l in Ldict if Ldict[l].D12() > 0]
    #BABBAloci = [Ldict[l].number for l in Ldict if Ldict[l].D12() < 0]
    #ABBAAloci = [Ldict[l].number for l in Ldict if Ldict[l].D1() > 0]
    #BABAAloci = [Ldict[l].number for l in Ldict if Ldict[l].D1() < 0]
    #ABABAloci = [Ldict[l].number for l in Ldict if Ldict[l].D2() > 0]
    #BAABAloci = [Ldict[l].number for l in Ldict if Ldict[l].D2() < 0]

    return [L,
            DFO,ZFO,
            DIL,ZIL,
            DFI,ZFI,
            DOL,ZOL,
            len(Ldict),
            sum([Ldict[l].babba for l in Ldict]),
            sum([Ldict[l].abbba for l in Ldict]),
            sum([Ldict[l].babaa for l in Ldict]),
            sum([Ldict[l].abbaa for l in Ldict]),
            sum([Ldict[l].baaba for l in Ldict]),
            sum([Ldict[l].ababa for l in Ldict]),
            sum([Ldict[l].bbbaa for l in Ldict]),
            sum([Ldict[l].bbaba for l in Ldict]),
            sum([Ldict[l].aabaa for l in Ldict]),
            sum([Ldict[l].aaaba for l in Ldict]),
            sum([Ldict[l].baaaa for l in Ldict]),
            sum([Ldict[l].abaaa for l in Ldict]),
            pdisc, submitted,
            BBFO, BBIL, BBFI, BBOL]


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




def multiproc_it(subtests,alignfile,outfile, nboots,nproc,namelen,makesort,makeboots,noterminals):
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
                work_queue.put([alignfile, [p1,p2,p3a,p3b,o], nboots, 1, submitted, noterminals])
                submitted += 1
            else: 
                print 'a taxon name was found that is not in the sequence file'
        else:
            if checktaxa([p1,p2,p3a,p3b,o],alignfile):
                work_queue.put([alignfile, [p1,p2,p3a,p3b,o], nboots, 0, submitted, noterminals])
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
    Results = [result_queue.get() for i in range(submitted)]
    Results.sort(key = lambda x:x[15])



    " setup results file "
    if noterminals: 
        outs = open(outfile+".Dfoilalt.txt", 'w')
    else:
        outs = open(outfile+".Dfoil.txt", 'w')
    header = "\t".join([ 'p1'+" "*(namelen[0]-2),
                         'p2'+" "*(namelen[1]-2),
                         'p3'+" "*(namelen[2]-2),
                         'p4'+" "*(namelen[3]-2),
                         'O'+" "*(namelen[4]-1),
                         'Dfo','Dil','Dfi','Dol',
                         'Z_fo','Z_il','Z_fi','Z_ol',
                         'BABBA','ABBBA',
                         'BABAA','ABBAA',
                         'BAABA','ABABA',
                         'BBBAA','BBABA',
                         'AABAA','AAABA',
                         'BAAAA','ABAAA',
                         'nloci','sign', 'notes'])
    print >>outs, header

    for i in range(len(Results)):
        L,DFO,ZFO,DIL,ZIL,DFI,ZFI,DOL,ZOL,nloc,BABBA,ABBBA,BABAA,ABBAA,BAABA,ABABA,BBBAA,BBABA,AABAA,AAABA,BAAAA,ABAAA,pdisc,sub,BBFO,BBIL,BBFI,BBOL = Results[i]
        L = [str(x).replace("['","[").replace("']","]").replace("', '",",") for x in L]

        sign = []
        for s,d in zip([ZFO,ZIL,ZFI,ZOL],[DFO,DIL,DFI,DOL]):
            if s>3.0:
                if d>0:
                    sign.append("+")
                else:
                    sign.append("-")
            else:
                sign.append("0")
        #print sign

        resin = tuple([str(L[0])+" "*(namelen[0]-len(str(L[0]))),
                       str(L[1])+" "*(namelen[1]-len(str(L[1]))),
                       str(L[2])+" "*(namelen[2]-len(str(L[2]))),
                       str(L[3])+" "*(namelen[3]-len(str(L[3]))),
                       str(L[4])+" "*(namelen[4]-len(str(L[4]))),
                       DFO,DIL,DFI,DOL,
                       ZFO,ZIL,ZFI,ZOL,
                       BABBA,ABBBA,BABAA,ABBAA,BAABA,ABABA,BBBAA,BBABA,AABAA,AAABA,BAAAA,ABAAA,
                       nloc, "".join(sign), Notes[i]])
        
        print >>outs, "%s\t%s\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%s\t%s" % resin 

        loci = open(alignfile).read().strip().split("|")[:-1]
        if makesort:
            None
            # makesortfiles("ABBBA",ABBBAloci,5,loci,outfile,makesort,sub,L)
            # makesortfiles("BABBA",BABBAloci,5,loci,outfile,makesort,sub,L)            
            # makesortfiles("ABBAA",ABBAAloci,5,loci,outfile,makesort,sub,L)            
            # makesortfiles("BABAA",BABAAloci,5,loci,outfile,makesort,sub,L)
            # makesortfiles("ABABA",ABABAloci,5,loci,outfile,makesort,sub,L)
            # makesortfiles("BAABA",BAABAloci,5,loci,outfile,makesort,sub,L)

        if makeboots:
            None
            # with open(outfile+"_"+str(sub+1)+".boots_D12",'w') as out:
            #     out.write(",".join(map(str,BB12)))
            # with open(outfile+"_"+str(sub+1)+".boots_D1",'w') as out:
            #     out.write(",".join(map(str,BB1)))
            # with open(outfile+"_"+str(sub+1)+".boots_D2",'w') as out:
            #     out.write(",".join(map(str,BB2)))
                

def main(tests, alignfile, outfile, nboots, nproc, makesort, makeboots,noterminals):
    import sys

    P1namelen  = max(map(len,[str(i[0][0]) for i in tests]))
    P2namelen  = max(map(len,[str(i[0][1]) for i in tests]))
    P3anamelen = max(map(len,[str(i[0][2]) for i in tests]))
    P3bnamelen = max(map(len,[str(i[0][3]) for i in tests]))
    Onamelen   = max(map(len,[str(i[0][4]).strip() for i in tests]))
    namelen = [P1namelen,P2namelen,P3anamelen,P3bnamelen,Onamelen]
                
    multiproc_it(tests,alignfile,outfile,nboots,nproc,namelen,makesort,makeboots,noterminals)


if __name__ == '__main__':
    main()




