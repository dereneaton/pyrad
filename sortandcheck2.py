#!/usr/bin/env python2

import gzip
import itertools
import sys
import glob
import os
import cPickle as pickle
import multiprocessing
from potpour import Worker


def combinefiles(GLOB):
    "combines first and second reads file names"
    if len(glob.glob(GLOB)) > 1:
        FS = glob.glob(GLOB)
    else:
        FS = glob.glob(GLOB)
    firsts = [i for i in FS if "_R1_" in i]
    seconds = [ff.replace("_R1_","_R2_") for ff in firsts]
    if len(firsts) != len(seconds):
        raise Exception("different numbers of first and second read files.\
              Check that the names of files are correct")
    return zip(firsts,seconds)


def revcomp(s):
    "returns reverse complement of a string"
    ss = s[::-1].strip().replace("A","t").replace("T","a").\
         replace("C","g").replace("G","c").upper()
    return ss


def matching(a,b, maxmismatch):
    "allows for N base difference between barcodes"
    if len(a) == len(b):
        t = [a[i]==b[i] for i in range(len(a))]
        if t.count(False) <= maxmismatch:
            return 1
        else:
            return 0
    else:
        return 0


def unambig(seq):
    """ returns both resolutions of a cut site
    that has an ambiguous base in it """
    resos = []
    D = {"R":("G","A"),
         "K":("G","T"),
         "S":("G","C"),
         "Y":("T","C"),
         "W":("T","A"),
         "M":("C","A")}
    for base in list("RKSYWM"):
        if base in seq:
            resos.append(seq.replace(base,D[base][0]))
            resos.append(seq.replace(base,D[base][1]))
    return resos


def findbcode(CUT,longB,l):
    barcode = 'N'*20
    " in case ambiguous base in CUT "
    if any([i in CUT for i in list("RKYSWM")]):
        CUT = unambig(CUT)
        Bs = []
        for cut in CUT:
            if l[1][0:longB+len(cut)].count(cut) == 1:
                barcode = l[1].split(cut)[0].strip()
            elif l[1][0:longB+len(cut)].count(cut) == 2:
                barcode = cut.join(l[1].split(cut)[0:2]).strip()
            else:
                barcode = ""
            Bs.append(barcode)
        longestbar = Bs[map(len,Bs).index(max(map(len,Bs)))]
        return longestbar
    else:
        if l[1][0:longB+len(CUT)].count(CUT) == 1:
            barcode = l[1].split(CUT)[0].strip()
        elif l[1][0:longB+len(CUT)].count(CUT) == 2:
            barcode = CUT.join(l[1].split(CUT)[0:2]).strip()
        else:
            barcode = ""
    return barcode



def barmatch(C, Raws, CUT, datatype, num, maxmismatch, WORK, longB):
    """matches reads to barcodes in barcode file
    and writes to individual temp files, after all
    read files have been split, temp files are collated
    into .fq files"""
    
    #CUT1 = CUT = unambig(CUT)[0]
    locus = 0
    match = 0
    match2 = 0
    barcodehits = set()

    " dictionary to record barcode misses" 
    M = {}
    M['_'] = 0
    
    "read in paired end read files"
    if 'pair' in datatype:
        if '.gz' in Raws[0][-3:]:
            fr1 = gzip.open(Raws[0])
        else:
            fr1 = open(Raws[0])
        if '.gz' in Raws[1][-3:]:
            fr2 = gzip.open(Raws[1])
        else:
            fr2 = open(Raws[1])
        R1 = itertools.izip(*[iter(fr1)]*4)
        R2 = itertools.izip(*[iter(fr2)]*4)
    else:
        "read in single end read file"
        if '.gz' in Raws[-3:]:
            fr1 = gzip.open(Raws)
        else:
            fr1 = open(Raws)
        R1 = itertools.izip(*[iter(fr1)]*4)

    D = {}
    DD = {}
    while 1:
        try: r1 = R1.next()
        except StopIteration: break

        "match paired end reads together, for now"
        if 'pair' in datatype:
            r2 = R2.next()
            l = [r.strip() for r in r1]
            l = [l[0],l[1],l[2],l[3]]
            ll = [r.strip() for r in r2]
            ll = [ll[0],ll[1],ll[2],ll[3]]
        else:
            "make list of four fastq line elements"
            l = [r.strip() for r in r1]
            l = [l[0],l[1],l[2],l[3]]

        locus += 1

        if 'pair' in datatype:
            if longB[1] == 'same':
                " bars are all same length"
                barcode = l[1][:longB[0]]
            else:
                " find barcodes"
                barcode = findbcode(CUT,longB[0],l)
            if barcode:
                if barcode in M:
                    M[barcode] += 1
                else:
                    M[barcode] = 1

                "exclude the read if no cutsite/barcode found"
                if barcode in D:
                    l[1] = l[1][len(barcode):]   #barcode.join(l[1].split(barcode)[1:])
                    l[3] = l[3][len(barcode):]
                    D[barcode].append("\n".join(l).strip())
                    DD[barcode].append("\n".join(ll).strip())
                    match += 1
                else:
                    l[1] = l[1][len(barcode):]   #barcode.join(l[1].split(barcode)[1:])
                    l[3] = l[3][len(barcode):]                        
                    D[barcode] = l
                    DD[barcode] = ll
                    match += 1

            else:
                M["_"] += 1

        else:
            if longB[1] == 'same':
                if datatype=='2brad':
                    barcode = l[1][-longB[0]:]
                else:
                    barcode = l[1][:longB[0]]
            else:
                barcode = findbcode(CUT,longB[0],l)
            if barcode:
                " tracker of number of occurrences of each barcode"
                if barcode in M:
                    M[barcode] += 1
                else:
                    M[barcode] = 1

                "exclude the read if no cutsite/barcode found"
                "saves reads from barcodes to a dictionary D"
                if barcode in D:
                    #l[1] = CUT1+l[1][len(barcode)+len(CUT):]
                    if datatype=='2brad':
                        l[1] = l[1][:-len(barcode)]
                        l[3] = l[3][:-len(barcode)]
                    else:
                        l[1] = l[1][len(barcode):]
                        l[3] = l[3][len(barcode):]
                    D[barcode].append("\n".join(l).strip())
                    match += 1
                else:
                    l[1] = l[1][len(barcode):]
                    l[3] = l[3][len(barcode):]                        
                    D[barcode] = l
                    match += 1
            else:
                M["_"] += 1


        "write to file every 50Kth read"
        " only writes reads that match to a barcode in C by less than some N differences "
        if not locus % 50000:
            for bar in C:
                outF1 = gzip.open(WORK+"fastq/."+C[bar]+'.temp_R1_'+str(num)+'.gz','ab')
                if 'pair' in datatype:
                    outF2 = gzip.open(WORK+"fastq/."+C[bar]+'.temp_R2_'+str(num)+'.gz','ab')
                for barcode in D:
                    if matching(bar,barcode,maxmismatch):
                        barcodehits.add(barcode)
                        if D[barcode]:
                            match2 += len(D[barcode])  ## -3
                            outF1.write("\n".join(D[barcode])+'\n')
                        if 'pair' in datatype:
                            if DD[barcode]:
                                outF2.write("\n".join(DD[barcode])+'\n')
                        D[barcode] = []
                        DD[barcode] = []
                outF1.close()
                if 'pair' in datatype:
                    outF2.close()
                D[bar] = []
                DD[bar] = []

                
    "write the remaining reads to file"
    for bar in C:
        outF1 = gzip.open(WORK+"fastq/."+C[bar]+'.temp_R1_'+str(num)+'.gz','ab')
        if 'pair' in datatype:
            outF2 = gzip.open(WORK+"fastq/."+C[bar]+'.temp_R2_'+str(num)+'.gz','ab')
        for barcode in D:
            if matching(bar,barcode,maxmismatch):
                barcodehits.add(barcode)
                if D[barcode]:
                    match2 += len(D[barcode])  ## -3
                    outF1.write("\n".join(D[barcode])+'\n')
                if 'pair' in datatype:
                    if DD[barcode]:
                        outF2.write("\n".join(DD[barcode])+'\n')
                D[barcode] = []
                DD[barcode] = []
        outF1.close()
        if 'pair' in datatype:
            outF2.close()
        D[bar] = []
        DD[bar] = []


    sys.stderr.write(".")
    fr1.close()
    if 'pair' in datatype:
        fr2.close()
                        
    "writes statistics out" 
    statout = open(WORK+"stats/s1.sorting.txt",'a')
    if 'pair' in datatype:
        name = Raws[0].split("/")[-1].replace("_R1_","_")
    else:
        name = Raws.split("/")[-1].replace("_R1_","_")

    match2 = sum([M[i] for i in M if i in barcodehits])
    writeit = "%s\t%i\t%i\t%i\n" % (name, locus, match, match2)
    statout.write(writeit)
    statout.close()
    pickout = open(WORK+"fastq/."+name+".pickle","wb")
    pickle.dump( M, pickout)
    pickout.close()


def writefunc(GLOB,Parallel,Bcode,CUT,datatype,maxmismatch,WORK):
    "create barcode dictionary"
    codetable = open(Bcode, 'r')
    codes = [line.strip().split() for line in codetable.readlines()]
    C = {}
    for line in codes:
        if line[0]:
            C[line[1].strip().upper()] = line[0]

    " find longest barcode "
    keylens = map(len,C.keys())
    if len(set(keylens)) == 1:
        longB = (keylens[0],'same')
    else:
        longB = (max(keylens),'diff')

    " check for CUT in barcodes "
    CCC = unambig(CUT)
    if len(CCC)>1:
        for cut in CCC:
            if any([cut in i for i in C.keys()]):
                print "\n\twarning: CUT site matches within one of the barcodes, "+\
                "I suggest double \n\tchecking the file to make sure it properly demultiplexes"
    else:
        if any([CUT in i for i in C.keys()]):
            print "\n\twarning: CUT site matches within one of the barcodes, "+\
            "I suggest double \n\tchecking the file to make sure it properly demultiplexes"

    " read in sequence files "
    if len(glob.glob(GLOB)) > 1:
        FS = [f for f in glob.glob(GLOB)]
    else:
        FS = glob.glob(GLOB)
    if 'pair' in datatype:
        Raws = combinefiles(GLOB)
    else:
        Raws = FS

    "send jobs to multiprocess queue"
    num = 0
    work_queue = multiprocessing.Queue()
    submitted = 0
    for fs in Raws:
        if 'pair' in datatype:
            work_queue.put([C, [fs[0],fs[1]], CUT, datatype, num, maxmismatch, WORK, longB])
            submitted += 1
        else:
            work_queue.put([C, fs, CUT, datatype, num, maxmismatch, WORK, longB])
            submitted += 1
        num += 1

    result_queue = multiprocessing.Queue()

    "spawn workers, give function"
    jobs = []
    for i in range( min(Parallel,submitted) ):
        worker = Worker(work_queue, result_queue, barmatch)
        worker.start()
        jobs.append(worker)
    for job in jobs:
        job.join()

    Ms = {}

    if len(glob.glob(WORK+"fastq/.*.pickle")) > 1:
        for pick in glob.glob(WORK+"fastq/.*.pickle"):
            pickin = open(pick, "rb")
            M = pickle.load( pickin )
            pickin.close()
            for key in M:
                if key not in Ms:
                    Ms[key] = M[key]
                else:
                    Ms[key] += M[key]
            os.remove(pick)
    elif len(glob.glob(WORK+"fastq/.*.pickle")) == 1:
        pick = glob.glob(WORK+"fastq/.*.pickle")[0]
        pickin = open(pick, 'rb')
        Ms = pickle.load( pickin )
        pickin.close()
        os.remove(pick)
    else:
        print "\nno stats file generated"

    Mkeys = Ms.keys()
    Mkeys.sort(key=lambda x: Ms[x], reverse=True)

    statout = open(WORK+"stats/s1.sorting.txt",'a')
    statout.write("\n\n")
    statout.write("sample\ttrue_bar\tobs_bars\tN_obs\n")

    Cnames = C.keys()
    Cnames.sort()
    print Ms.values()
    try: maxl = max(map(len,map(str,Ms.values())))
    except ValueError: maxl = 5

    hits = []
    for bar in Cnames:
        for barcode in Mkeys:
            if matching(bar, barcode, maxmismatch):
                print >>statout, "%s    \t%s    \t%s\t%s" % (C[bar], bar, barcode,
                                                             str(Ms[barcode])+" "*(maxl+3-len(str(Ms[barcode]))))
                hits.append(barcode)

    statout.write("\n")
    maxl = max(map(len,Mkeys))
    for barcode in Mkeys:
        if barcode not in hits:
            print >>statout, "nomatch  \t%s    \t%i" % (barcode+" "*(maxl+3-len(barcode)), Ms[barcode])
    statout.close()



def main(Bcode, GLOB, CUT, datatype, Parallel, maxmismatch, WORK):

    if not len(glob.glob(GLOB)) > 0:
        sys.stderr.write("\tno data found in "+GLOB+" fix path to the data files\n")
        sys.exit()

    "check for previous output"
    if not os.path.exists(WORK+'stats'):
        os.makedirs(WORK+'stats')
    if os.path.exists(WORK+'fastq'):
        if os.listdir(WORK+'fastq'):
            print ("\n\tfastq/ directory in working directory contains data, move/remove it before running step 1\n")
            sys.exit()
    else:
        os.makedirs(WORK+'fastq')

    if "*" in Bcode:
        if len(glob.glob(Bcode)) == 1:
            Bcode = glob.glob(Bcode)[0]

    sys.stderr.write("\n\tstep 1: sorting reads by barcode\n\t ")

    " seperate double digest cut sites, only need first read one for now "
    if "," in CUT:
        CUT,CUT2 = CUT.split(",")

    statout = open(WORK+"stats/s1.sorting.txt",'w')
    statout.write("\t".join(["file    ","Nreads","cut_found","bar_matched"])+"\n")
    statout.close()

    " DO THE BARCODE SORTING "
    writefunc(GLOB, Parallel, Bcode, CUT, datatype, maxmismatch, WORK)
    names = [line.split()[0] for line in open(Bcode).readlines()]

    " remove tiny sorted temp files "
    if len(glob.glob(GLOB)) > 1:
        for name in names:
            if len(glob.glob(WORK+"fastq/."+name+"*")) > 0:
                "remove very small files, probably errors"
                for ff in glob.glob(WORK+'fastq/.'+name+"*"):
                    statinfo = os.stat(ff)
                    s = statinfo.st_size
                    if s < 1000:
                        os.remove(ff)
                        

    " concatenate temp files "
    for name in names:
        if len(glob.glob(WORK+"fastq/."+name+"*")) > 0:
            os.system("/bin/cat "+WORK+"fastq/."+name+".temp_R1_*.gz > "+WORK+"fastq/"+name+"_R1.fq.gz")
            if datatype in ['pairgbs','pairddrad']:
                os.system("/bin/cat "+WORK+"fastq/."+name+".temp_R2_*.gz > "+WORK+"fastq/"+name+"_R2.fq.gz")

    if len(glob.glob(WORK+"fastq/*")) > 0:
        os.system("/bin/ls "+WORK+"fastq/.*temp_* | xargs /bin/rm" )
    if len(glob.glob(WORK+"fastq/*.pickle")) > 0:
        os.system("/bin/ls "+WORK+"fastq/.*pickle | xargs /bin/rm" )
