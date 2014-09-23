
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
        return CUT


def most_common(L):
  return max(itertools.groupby(sorted(L)), key=lambda(x, v):(len(list(v)),-L.index(x)))[0]


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

    " remove name suffix"
    n = str(infile.split('/')[-1]).replace("_R1.",".")
    while n.split(".")[-1] in ["fastq","fastQ","gz","fq","FastQ"]:
        n = n.replace('.'+n.split(".")[-1], "")

    " read infile 4 lines at a time, setup counters and lists"
    k = itertools.izip(*[iter(f)]*4)
    writing_r = []
    writing_c = []
    orig = keep = keepcut = 0
    handle = WORK+'edits/'+str(n)+".edit"


    " do a test run on first 1000 reads to find if extra bases on right end of reads"
    rightend = []
    while len(rightend) < 1000:
        try: d = k.next()
        except StopIteration: break
        s = "".join(d[1].strip())
        
        " cutters "
        find1 = CUT1
        find2 = comp(CUT2)[::-1]

        " are cutters found on both ends? A correct merge"
        a = s[:len(find1)]
        b = s[-len(find2)-2:]  ## w/ wiggle room
        if (find1 in a) and (find2 in b) :
            xtra = s.rindex(find2)+len(find2)
            rightend.append(len(s)-xtra)
            
    " find most common element in rightend "
    if rightend:
        a = most_common(rightend)
        if a>3:
            Roffset = 0
        else:
            Roffset = a
    else:
        Roffset = 0

    " reset iterable "
    if ".gz" in infile:
        f = gzip.open(infile, 'r')
    else:
        f = open(infile,'r')
    k = itertools.izip(*[iter(f)]*4)

    " iterate over each read "
    while 1:
        try: d = k.next()
        except StopIteration: break
        orig += 1 
        SS = d[1].strip()
    
        " apply Phred Q filter "
        ph = map(ord,d[3].strip('\n'))      
        offset = int(Q) 
        phred = map(lambda x:x-offset,ph)
        seq = ["N"]*len(phred)
        for base in range(len(phred)):
            "don't quality check cut sites "
            if (base >= len(CUT1)) and (base < len(phred)-len(CUT2)):
                if phred[base] >= 20:       
                    try: seq[base] = SS[base]
                    except IndexError:
                        None
                else:
                    seq[base] = "N"
            else:
                try: seq[base] = SS[base]
                except IndexError:
                    None

        " write to file "    
        if not orig % 5000:
            with open(WORK+'edits/'+str(n)+".edit",'a') as outfile:
                outfile.write("".join([z for z in writing_r]))
            writing_r = []

        s = "".join(seq)

        wheretocut = [None,None,None]
        " filter for N"
        if s.count("N") <= pN:

            " apply filter for Adapters "
            find1 = CUT1
            find2 = comp(CUT2)[::-1]

            if "trim" in d[0]:
                " filters for non-merged, trimmed reads from s2 "
                if (find1 in s[:len(find1)]) or (find2 in s[len(find2)-2:]):
                    None
                else:
                    " CUT1 rarely missing, CUT2 sometimes missing"
                    s = s[:-len(CUT2)-Roffset]

            else:
                " merged reads. Are cutters found on both ends? A correct merge"
                a = s[:len(find1)]
                b = s[-len(find2)-2:]  ## w/ wiggle room
                if (find1 in a) and (find2 in b) :
                    " find end of read2 "
                    xtra = s.rindex(find2)+len(find2)
                    wheretocut = [None, len(s)-Roffset, 'complete']
                else:
                    " look for CUT2 from right side "
                    if find2 in s[len(s)/2:]:   ## check that this is a good general number...
                        a = s.rindex(find2)+len(find2)
                        wheretocut = [None, a, 'find2 in s']
                    else:
                        "couldn't find cut2, maybe has error, look for adapter"
                        if 'AGATCG' in s:
                            a = s.rindex('AGATCG')-len(CUT2)
                            wheretocut = [None, a, 'AGATCG in s']
                        else:
                            if "TCGGAAG" in s:
                                a = s.rindex('TCGGAAG')-len(CUT2)-3
                                wheretocut = [None, a, 'TCGGAAG in s']
                            else:
                                "no sign of overshoot to right --->"
                                " look for overshoot on left <---- "
                                wheretocut = [None, len(s)-Roffset, "None"]

                    " look for CUT1 from left side "
                    if CUT1 in s:
                        a = s.index(CUT1)
                        wheretocut[0] = a
                    else:
                        "exclude read"
                        wheretocut[0] = wheretocut[1]

            w1,w2,reason = wheretocut
            if len(s[w1:w2]) > trimkeep:
                #print s[w1:w2], reason, len(s[w1:w2]), trimkeep
                s = s[w1:w2]
            else:
                s = ""
                
            if len(s) >= max(36,trimkeep): ## if read is trimmed, must be minlen long
                writing_r.append(">"+n+"_"+str(keep)+"_r1"+"\n"+s+"\n")
                keep += 1

    with open(WORK+'edits/'+str(n)+".edit",'a') as outfile:
        outfile.write("".join([z for z in writing_r]))
    writing_r = []

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
