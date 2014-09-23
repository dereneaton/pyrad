import multiprocessing
import itertools
import sys
import os
import glob
import operator
import gzip
from potpour import Worker
from sortandcheck2 import unambig
from cluster_cons7_shuf import comp



def unambar(CUT):
    if any([i in CUT for i in list("RKYSWM")]):
        CUTa, CUTb = unambig(CUT)
        return [CUTa,CUTb]
    else:
        return False


def Afilter(CUT,s,strict):
    a = b = wheretocut = None
    " lookfor cut site "
    if unambar(CUT):
        " if ambiguity in cutter "
        CUTa,CUTb = unambar(CUT)
        if strict == 2:
            lookfor1 = CUTa+"A"
            lookfor2 = CUTb+"A"
        else:
            lookfor1 = CUTa+"AGA"
            lookfor2 = CUTb+"AGA"
        if lookfor1 in s:
            a = s.rindex(lookfor1)
        if lookfor2 in s:
            b = s.rindex(lookfor2)
        if (a or b):
            wheretocut = min([i for i in [a,b] if i])
        else:
            wheretocut = None
    else:
        if strict == 2:
            lookfor = CUT+"A"
        else:
            lookfor = CUT+"AGA"
        if lookfor in s:
            wheretocut = s.rindex(lookfor)
        else:
            wheretocut = None
                
    if not wheretocut:
        " look for adapter sequence "
        if strict == 2:
            lookfor1 = "AGATCG"
        else:
            lookfor1 = "AGATCGGA"
        if lookfor1 in s:
            wheretocut = s.rindex(lookfor1)-(len(CUT)+1)
        else:
            wheretocut = None

    " look for CUT at end of seq "
    if not wheretocut:
        if CUT in s[-len(CUT)-5:]:
            wheretocut = s.rindex(CUT)
    return wheretocut
            
                



def rawedit(WORK, infile, CUT, pN, trimkeep, strict, Q):
    """ three functions:
    (1) replaces low quality base calls with Ns,
    (2) checks for adapter sequence if strict set to 1 or 2 """

    if "," in CUT:
        CUT1,CUT2 = CUT.split(',')
    else:
        CUT1=CUT2=CUT
        
    if ".gz" in infile:
        f = gzip.open(infile, 'r')
    else:
        f = open(infile,'r')
    n = str(infile.split('/')[-1]).replace("_R1.",".")
    while n.split(".")[-1] in ["fastq","fastQ","gz","fq","FastQ"]:
        n = n.replace('.'+n.split(".")[-1], "")
    k = itertools.izip(*[iter(f)]*4)
    writing_r = []
    writing_c = []

    orig = keep = keepcut = 0
    handle = WORK+'edits/'+str(n)+".edit"

    while 1:
        try: d = k.next()
        except StopIteration: break
        orig += 1 
        SS = d[1].strip()
    
        ph = map(ord,d[3].strip('\n'))      
        offset = int(Q) 
        phred = map(lambda x:x-offset,ph)
        seq = ["N"]*len(phred)
        for base in range(len(phred)):
            if base >= len(CUT1):              ## don't quality check cut site
                if phred[base] >= 20:         ## quality threshold
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
                #try: seq[base] = SS[base]
                #except IndexError:
                #    None
            
        if not orig % 5000:
            if trimkeep:
                " write full length and fragment reads "
                with open(WORK+'edits/'+str(n)+".edit",'a') as outfile:
                    outfile.write("".join([z for z in writing_r]))
                    outfile.write("".join([z for z in writing_c]))
            else:
                " write only full length reads "
                with open(WORK+'edits/'+str(n)+".edit",'a') as outfile:
                    outfile.write("".join([z for z in writing_r]))
            writing_r = []
            writing_c = []

        s = "".join(seq)
        wheretocut1 = None
        if strict:
            wheretocut1 = Afilter(comp(CUT2)[::-1],s,strict)
            s = s[:wheretocut1]

        if s.count("N") <= pN:             ## max allowed Ns
            if len(s) >= max(36,trimkeep): ## if read is trimmed, must be minlen long
                if wheretocut1:            ## if it was trimmed...
                    writing_c.append(">"+n+"_"+str(keepcut)+"_c1"+"\n"+s+"\n")
                    keepcut += 1
                else:
                    writing_r.append(">"+n+"_"+str(keep)+"_r1"+"\n"+s+"\n")
                    keep += 1

    if trimkeep:
        with open(WORK+'edits/'+str(n)+".edit",'a') as outfile:
            outfile.write("".join([z for z in writing_r]))
            outfile.write("".join([z for z in writing_c]))
    else:
        with open(WORK+'edits/'+str(n)+".edit",'a') as outfile:
            outfile.write("".join([z for z in writing_r]))
    writing_r = []
    writing_c = []

    f.close()
    sys.stderr.write(".")
    if not trimkeep:
        keepcut = 0
    return [handle.split("/")[-1].replace(".edit",""),str(orig),str(keep),str(keepcut)]



def main(Parallel, WORK, FQs, CUT, pN, Q, strict, trimkeep):
    print >>sys.stderr, "\tstep 2: editing raw reads \n\t",

    " create output directories "
    if not os.path.exists(WORK+'stats'):
        os.makedirs(WORK+'stats')
    if not os.path.exists(WORK+'edits'):
        os.makedirs(WORK+'edits')

    " used to find if files already exist "
    lookfor = ".edit"

    " load up work queue "
    submitted = 0
    work_queue = multiprocessing.Queue()
    if len(glob.glob(FQs)) > 1:
        FS = [f for f in glob.glob(FQs)]

        " order files by size "
        for i in range(len(FS)):
            statinfo = os.stat(FS[i])
            FS[i] = FS[i],statinfo.st_size
        FS.sort(key=operator.itemgetter(1))
        FS = [i[0] for i in FS][::-1]

        " submit jobs to queue "
        for handle in FS:
            finder = WORK+'edits/'+handle.split("/")[-1]
            while finder.split(".")[-1] in ["fastq","fastQ","gz","fq","FastQ"]:
                finder = finder.replace('.'+finder.split(".")[-1], "").replace("_R1","")
            if finder+".edit" not in glob.glob(WORK+"edits/*"):
                if os.stat(handle).st_size > 0:   ## exclude empty files
                    args = [WORK, handle, CUT, float(pN), trimkeep, strict, Q]
                    work_queue.put(args)
                    submitted += 1
            else:
                print "\t"+finder+" already in edits/"

    elif len(glob.glob(FQs)) == 1:
        " if only one file "
        work_queue.put([WORK, glob.glob(FQs)[0], CUT, float(pN), trimkeep, strict, Q])
        submitted += 1

    else:
        print "\tNo demultiplexed files found. Check path."

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
    outstats = open(WORK+"stats/s2.rawedit.txt",'a')
    print >> outstats, "\t".join(["sample ","Nreads","passed","passed.w.trim","passed.total"])
    STATS = []
    for i in range(submitted):
        STATS.append(result_queue.get())

    STATS.sort(key = lambda x: x[0])
    for i in range(submitted):
        a,b,c,d = STATS[i]
        print >> outstats, "\t".join([a,b,c,d,str(int(c)+int(d))])

    print >>outstats, """
    Nreads = total number of reads for a sample
    passed = retained reads that passed quality filtering at full length
    passed.w.trim= retained reads that were trimmed due to detection of adapters
    passed.total  = total kept reads of sufficient length
    note: you can set the option in params file to include trimmed reads of xx length. """
    outstats.close()

    #" zip files to save size "
    #for ff in glob.glob(WORK+"edits/*"):
    #    os.system("gzip "+ff)
