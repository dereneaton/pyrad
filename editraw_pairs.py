#!/usr/bin/env python2

import multiprocessing
import itertools
import sys
import os
import glob
import operator
import gzip
from potpour import Worker
from sortandcheck2 import unambig


def revcomp(s):
    ss = s[::-1].strip().replace("A","t").replace("T","a").\
         replace("C","g").replace("G","c").replace("n","Z").upper().replace("Z","n")
    return ss


def unambar(CUT):
    if any([i in CUT for i in list("RKYSWM")]):
        CUTa, CUTb = unambig(CUT)
        return [CUTa,CUTb]
    else:
        return False


def Afilter(CUT,s,strict,read):
    read1 = read==1
    a = b = lookfor = wheretocut = None
    " lookfor cut site "

    " if ambiguity in cutter "
    if unambar(CUT):
        CUTa,CUTb = unambar(CUT)
        if strict == 2:
            if read1:
                lookfor1 = CUTa+"A"
                lookfor2 = CUTb+"A"
            else:
                lookfor1 = CUTa
                lookfor2 = CUTb
        else:
            if read1:
                lookfor1 = CUTa+"AGAT"
                lookfor2 = CUTb+"AGAT"
            else:
                lookfor1 = "A"*50
                lookfor2 = "A"*50
        if lookfor1 in s:
            a = s.rindex(lookfor1)
        if lookfor2 in s:
            b = s.rindex(lookfor2)
        if (a or b):
            wheretocut = min([i for i in [a,b] if i])
        else:
            wheretocut = None
    else:
        "no ambiguity in cutter "
        if strict == 2:
            if read1:
                lookfor1 = CUT+"A"
            else:
                lookfor1 = CUT
        else:
            if read1:
                lookfor1 = CUT+"AGA"
            else:
                lookfor1 = "A"*50
        if lookfor1 in s:
            wheretocut = s.rindex(lookfor1)
        else:
            wheretocut = None
                
    " look for adapter sequence "
    if not wheretocut:
        if strict == 2:
            lookfor1 = "AGATCG"
        else:
            lookfor1 = "AGATCGGA"
        if lookfor1 in s:
            if read1:
                wheretocut = s.rindex(lookfor1)-(len(CUT)+1)
            else:
                wheretocut = s.rindex(lookfor1)-(len(CUT)+6)
        else:
            wheretocut = None

    " look for CUT and end of seq "
    if not wheretocut:
        if CUT in s[-(len(CUT)+5):]:
            wheretocut = s.rindex(CUT)
    return wheretocut
            



def rawedit(WORK, infile, CUT, pN, trimkeep, strict, Q, datatype):
    """ three functions:
    (1) replaces low quality base calls with Ns,
    (2) checks for adapter sequence and xxbarcodesxx if strict set to 1 or 2 
    (3) concatenate paired reads with a separator and write to file """

    if CUT:
        if "," in CUT:
            CUT1,CUT2 = CUT.split(",")
            CUT2=revcomp(CUT2)
            x = 0   
        else:
            CUT1=CUT
            CUT2=revcomp(CUT1)
            x = 1   ## trims garbage base off gbs

    " create iterators for R1 and R2 files "
    if ".gz" in infile:
        f1 = gzip.open(infile, 'rb')
        if ".forward." in infile:
            f2 = gzip.open(infile.replace(".forward.",".reverse."), 'r')
        else:
            f2 = gzip.open(infile.replace("_R1.","_R2."), 'r')
    else:
        f1 = open(infile,'r')
        if ".forward." in infile:
            f2 = open(infile.replace(".forward.",".reverse."), 'r')
        else:
            f2 = open(infile.replace("_R1.","_R2."), 'r')
    n = str(infile.split('/')[-1])
    while n.split(".")[-1] in ["fastq","fastQ","gz","fq","FastQ","nomerge"]:
        n = n.replace('.'+n.split(".")[-1], "")
    if '.forward' in n:
        n = n.split(".unassembled")[0]
    else:
        n = n.replace("_R1","")

    k1 = itertools.izip(*[iter(f1)]*4)
    k2 = itertools.izip(*[iter(f2)]*4)
    writing_r = []
    writing_c = []

    orig = keep = keepcut = 0
    handle = WORK+'edits/'+str(n)+".edit"

    "iterate over paired reads, edits 1st, if OK, append both to .edit file"
    while 1:
        try: d = k1.next()
        except StopIteration: break
        dd = k2.next()

        orig += 1
        SS = d[1].strip()
        ph = map(ord,d[3].strip())
        offset = int(Q)
        phred = map(lambda x:x-offset,ph)
        seq = ["N"]*len(phred)
        for base in range(len(phred)):
            if base >= len(CUT1):              ## don't quality check cut site
                if phred[base] >= 20:          ## quality threshold
                    try: seq[base] = SS[base]
                    except IndexError:
                        None
                else:
                    seq[base] = "N"
            else:
                if unambar(CUT1):
                    seq[base] = unambar(CUT1)[0][base]
                else:
                    seq[base] = CUT1[base]

        s = "".join(seq)
        wheretocut1 = None

        " apply filters for adapter sequences "
        " if GBS CUT2 = revcomp(CUT1)   ex: CTGCA"
        " if ddRAD CUT2 = revcomp(CUT2) ex: AATT "
        if strict:
            wheretocut1 = Afilter(CUT2,s,strict,1)

        if s.count("N") <= pN:              ## max allowed Ns
            if len(s) >= max(36,trimkeep):  ## if trimmed read1 length atleast t

                " first read is (maybe) good, now filter second reads "
                SS = dd[1].strip()
                ph = map(ord,dd[3].strip())
                " if PEAR filtered then seqs are reversed "
                if '.forward' in infile:
                    SS = revcomp(SS)
                    ph = ph[::-1]
                
                offset = int(Q)
                phred = map(lambda x:x-offset,ph)
                seq = ["N"]*len(phred)
                for base in range(len(phred)):
                    if base > len(CUT2):               ## don't quality check cut site
                        if phred[base] >= 20:          ## quality threshold
                            try: seq[base] = SS[base]
                            except IndexError: None
                        else:
                            seq[base] = "N"
                    else:
                        try: seq[base] = SS[base]
                        except IndexError: None
                s2 = "".join(seq)

                " filter for gbs read2s, b/c they will be clustered"
                badread = 0
                if datatype == "pairgbs":
                    s2 = s2[:len(s)]
                    if s2.count("N")>pN:
                        badread = 1

                " apply adapter filter to read 2 "
                wheretocut2 = None
                if strict:
                    wheretocut2 = Afilter(revcomp(CUT1),s2,strict,2)

                if (wheretocut1 and wheretocut2):
                    cutter = min(wheretocut1,wheretocut2)
                else:
                    cutter = max(wheretocut1,wheretocut2)
                if strict:
                    if not cutter:
                        if (revcomp(CUT1) in s2[-16:]) or (CUT2 in s[-10:]):
                            cutter = len(s)-16

                if not badread:
                    if cutter:
                        if cutter > max(36,trimkeep):
                            sout = ">"+n+"_"+str(keepcut)+"_trim1"+"\n"+s[:cutter]+\
                                   "\n"+d[2]+d[3][:cutter]+"\n"
                            writing_c.append(sout)
                            keepcut += 1
                            sout = ">"+n+"_"+str(keepcut)+"_trim2"+"\n"+revcomp(s2[x:cutter+5])+\
                                   "\n"+d[2]+d[3][x:cutter+5]+"\n"
                            writing_c.append(sout)
                            keepcut += 1
                    else:
                        sout = ">"+n+"_"+str(keep)+"_pair"+"\n"+s[:-1]+"nnnn"+revcomp(s2[x:])+"\n"
                        writing_r.append(sout)
                        keep += 1

        if not orig % 5000:
            if trimkeep:
                with open(WORK+'mergedreads/'+str(n)+"M.fq",'a') as outfile:
                    outfile.write("".join([z for z in writing_c]))
            " writes only full length reads "
            with open(WORK+'edits/'+str(n)+".edit",'a') as outfile:
                outfile.write("".join([z for z in writing_r]))
            writing_r = []
            writing_c = []

    if trimkeep:
        with open(WORK+'mergedreads/'+str(n)+"M.fq",'a') as outfile:
            outfile.write("".join([z for z in writing_c]))
    " writes only full length reads "
    with open(WORK+'edits/'+str(n)+".edit",'a') as outfile:
        outfile.write("".join([z for z in writing_r]))
    writing_r = []
    writing_c = []

    f1.close()
    f2.close()
    sys.stderr.write(".")
    if not trimkeep:
        keepcut = 0
    return [handle.split("/")[-1].replace(".edit",""),str(orig),str(keepcut/2),str(keep)]



def main(Parallel, WORK, FQs, CUT, pN, Q, strict, trimkeep, datatype):

    print >>sys.stderr, "\n\tstep 2: quality filtering \n\t",

    " create output directories "
    if not os.path.exists(WORK+'stats'):
        os.makedirs(WORK+'stats')
    if not os.path.exists(WORK+'edits'):
        os.makedirs(WORK+'edits')

    " load up work queue "
    submitted = 0
    work_queue = multiprocessing.Queue()

    " do not select merged or discarded reads if PEAR was used on data"
    FQs = glob.glob(FQs)
    fqs = [i for i in FQs if not any([j in i for j in ["discarded",".assembled."]])]
    
    if len(fqs) > 1:
        " subselect only the first reads "
        if any([".unassembled.forward." in i for i in fqs]):
            FS = [i for i in fqs if '.forward.' in i]
        else:
            FS = [i for i in fqs if '_R1.' in i]
        
        " order files by size "
        for i in range(len(FS)):
            statinfo = os.stat(FS[i])
            FS[i] = FS[i],statinfo.st_size
        FS.sort(key=operator.itemgetter(1))
        FS = [i[0] for i in FS][::-1]

        " submit jobs to queue "
        for handle in FS:
            n = handle.split('/')[-1]
            while n.split(".")[-1] in ["fastq","fastQ","gz","fq","FastQ","nomerge"]:
                n = n.replace('.'+n.split(".")[-1], "")
            if '.forward.' in n:
                n = n.split(".unassembled")[0]
            else:
                "_".join(n.split('_R')[:-1])
            if WORK+"edits/"+n+".edit" not in glob.glob(WORK+"edits/*"):
                if os.stat(handle).st_size > 0:     ## exclude empty files
                    args = [WORK, handle, CUT, float(pN), trimkeep, strict, Q, datatype]
                    work_queue.put(args)
                    submitted += 1
            else:
                print "\t"+n+'.edit'+" already in edits/"
    else:
        print "no de-multiplexed files found."

    " create a queue to pass to workers to store the results "
    result_queue = multiprocessing.Queue()

    " spawn workers, give function "
    jobs = []
    for i in range( min(Parallel,submitted) ):
        worker = Worker(work_queue, result_queue, rawedit)
        worker.start()
        jobs.append(worker)
    for job in jobs:
        job.join()


    " collect the results off the queue "
    outstats = open(WORK+"stats/s2.rawedit.txt",'w')
    print >> outstats, "\t".join(["sample","Nreads","exclude","trimmed","passed"])
    for i in range(submitted):
        a,b,c,d = result_queue.get()
        print >> outstats, "\t".join([a,b, str(int(b)-int(d)), c, d])

    print >>outstats, """
    Nreads = total number of reads for a sample
    exclude = reads that were excluded
    trimmed = reads that had adapter trimmed but were kept
    passed = total kept reads
    "
    #note: you can set minimum length of kept trimmed reads on line 32 the params file 
    #kept trimmed (fragment) reads written to""" ,WORK+"mergedreads/\n"
    outstats.close()

