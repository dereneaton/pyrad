#!/usr/bin/env python2

""" edits reads based on quality scores. Can be used
to check for adapters and primers, but is not optimized 
for all types of cutters """

import multiprocessing
import itertools
import sys
import os
import glob
import operator
import gzip
from potpour import Worker
from sortandcheck2 import unambig
from cluster_cons7_shuf import fullcomp


def unambar(cut):
    """ returns both resolutions of an ambiguous cutter """
    if any([i in cut for i in list("RKYSWM")]):
        cuta, cutb = unambig(cut)
        return [cuta, cutb]
    else:
        return False


def afilter(cut, seq, strict):
    """ applies filter for primers & adapters """
    check1 = check2 = wheretocut = None
    ## lookfor cut site
    if unambar(cut):
        ## if ambiguity in cutter
        cuta, cutb = unambar(cut)
        if strict == 2:
            lookfor1 = cuta+"A"
            lookfor2 = cutb+"A"
        else:
            lookfor1 = cuta+"AGA"
            lookfor2 = cutb+"AGA"
        if lookfor1 in seq:
            check1 = seq.rindex(lookfor1)
        if lookfor2 in seq:
            check2 = seq.rindex(lookfor2)
        if check1 or check2:
            wheretocut = min([i for i in [check1, check2] if i])
        else:
            wheretocut = None
    else:
        if strict == 2:
            lookfor = cut+"A"
        else:
            lookfor = cut+"AGA"
        if lookfor in seq:
            wheretocut = seq.rindex(lookfor)
        else:
            wheretocut = None
                
    if not wheretocut:
        ## look for adapter sequence "
        if strict == 2:
            lookfor1 = "AGATCG"
        else:
            lookfor1 = "AGATCGGA"
        if lookfor1 in seq:
            wheretocut = seq.rindex(lookfor1)-(len(cut)+1)
        else:
            wheretocut = None

    ## look for CUT at end of seq
    if not wheretocut:
        if cut in seq[-len(cut)-5:]:
            wheretocut = seq.rindex(cut)
    return wheretocut


def rawedit(params, infile, quiet):
    """ three functions:
    (1) replaces low quality base calls with Ns,
    (2) checks for adapter sequence if strict set to 1 or 2 """
    ## WORK, infile, CUT, pN, trimkeep, strict, Q, datatype):
    if "," in params["cut"]:
        cut1, cut2 = params["cut"].split(',')
    else:
        cut1 = cut2 = params["cut"]
        
    ## read in the demultiplexed reads file
    if infile.endswith(".gz"):
        dems = gzip.open(infile, 'r')
    else:
        dems = open(infile, 'r')

    ## get sample name from the file name
    fname = str(infile.split('/')[-1]).replace("_R1.", ".")
    while fname.split(".")[-1] in ["fastq", "fastQ", 
                                   "gz", "fq", "FastQ"]:
        fname = fname.replace('.'+fname.split(".")[-1], "")

    ## iterator that grabs 4 lines at a time
    k = itertools.izip(*[iter(dems)]*4)

    ## lists for stroring full (r) and trimmed (c) reads
    writing_r = []
    writing_c = []

    ## counters
    orig = keep = keepcut = 0
    handle = params["work"]+'edits/'+str(fname)+".edit"

    while 1:
        try: 
            quart = k.next()
        except StopIteration: 
            break
        orig += 1 
        iseq = quart[1].strip()   ## SS
        qscore = [ord(i) for i in quart[3].strip('\n')]  ## ph
        offset = int(params["Q"])
        phred = [x-offset for x in qscore]
        seq = ["N"]*len(phred)
        for base in range(len(phred)):
            if base >= len(cut1):              ## don't quality check cut site
                if phred[base] >= 20:          ## quality threshold
                    try: 
                        seq[base] = iseq[base]
                    except IndexError:
                        pass
                else:
                    seq[base] = "N"
            else:
                if unambar(cut1):
                    seq[base] = unambar(cut1)[0][base]
                else:
                    seq[base] = cut1[base]

        if not orig % 5000:
            if params["trimkeep"]:
                ## write full length and fragment reads
                with open(params["work"]+'edits/'+str(fname)+\
                          ".edit", 'a') as outfile:
                    outfile.write("".join([z for z in writing_r]))
                    outfile.write("".join([z for z in writing_c]))
            else:
                ## write only full length reads "
                with open(params["work"]+'edits/'+str(fname)+\
                          ".edit", 'a') as outfile:
                    outfile.write("".join([z for z in writing_r]))
            writing_r = []
            writing_c = []

        sseq = "".join(seq)  ## s
        wheretocut1 = None
        if params["strict"]:
            wheretocut1 = afilter(fullcomp(cut2)[::-1],
                                  sseq, params["strict"])
            sseq = sseq[:wheretocut1]

        if params["datatype"] == 'merged':
            ## remove extra forward base so forwards match reverse length
            sseq = sseq[:-1]

        ## max allowed Ns
        if sseq.count("N") <= params["maxN"]:  
            ## if read is trimmed, must be minlen long
            if len(sseq) >= max(32, params["trimkeep"]): 
                ## if it was trimmed
                if wheretocut1:     
                    writing_c.append(">"+fname+"_"+str(keepcut)+\
                                     "_c1"+"\n"+sseq+"\n")
                    keepcut += 1
                else:
                    writing_r.append(">"+fname+"_"+str(keep)+\
                                     "_r1"+"\n"+sseq+"\n")
                    keep += 1

    if params["trimkeep"]:
        with open(params["work"]+'edits/'+str(fname)+".edit", 'a') as outfile:
            outfile.write("".join([z for z in writing_r]))
            outfile.write("".join([z for z in writing_c]))
    else:
        with open(params["work"]+'edits/'+str(fname)+".edit", 'a') as outfile:
            outfile.write("".join([z for z in writing_r]))
    writing_r = []
    writing_c = []

    dems.close()
    if not quiet:
        sys.stderr.write(".")
    if not params["trimkeep"]:
        keepcut = 0

    return [handle.split("/")[-1].replace(".edit", ""), 
            str(orig), str(keep), str(keepcut)]


def main(params, fastqs, quiet):
    """ runs the main script """
    ##Parallel, WORK, FQs, CUT, pN, Q, strict, trimkeep, datatype):

    if not quiet:
        print >>sys.stderr, "\n\tstep 2: editing raw reads \n\t",

    ## create output directories "
    if not os.path.exists(params["work"]+'stats'):
        os.makedirs(params["work"]+'stats')
    if not os.path.exists(params["work"]+'edits'):
        os.makedirs(params["work"]+'edits')

    ## load up work queue "
    submitted = 0
    work_queue = multiprocessing.Queue()
    if len(glob.glob(fastqs)) > 1:
        ffiles = glob.glob(fastqs)  ## FS

        ## order files by size 
        for ffile in range(len(ffiles)):
            statinfo = os.stat(ffiles[ffile])
            ffiles[ffile] = ffiles[ffile], statinfo.st_size
        ffiles.sort(key=operator.itemgetter(1))
        ffiles = [ffile[0] for ffile in ffiles][::-1]

        ## submit jobs to queue 
        for ffile in ffiles:
            finder = params["work"]+'edits/'+ffile.split("/")[-1]
            while any([finder.endswith(i) for i in ["fastq", "fastQ",
                                                    "gz", "fq", "FastQ"]]):
                finder = finder.replace('.'+finder.split(".")[-1], "").\
                         replace("_R1", "")
            if finder+".edit" not in glob.glob(params["work"]+"edits/*"):
                ## exclude empty files
                if os.stat(ffile).st_size > 0:  
                    work_queue.put([params, ffile, quiet])
                    submitted += 1
                else:
                    print "skipping", ffile, ", file is empty"
            else:
                print "\t"+finder+" already in edits/"

    elif len(glob.glob(fastqs)) == 1:
        ## if only one file "
        work_queue.put([params, glob.glob(fastqs)[0], quiet])
        submitted += 1

    else:
        print "\tNo demultiplexed files found. Check path."
        sys.exit()

    ## create a queue to pass to workers to store the results "
    result_queue = multiprocessing.Queue()

    ## spawn workers, give function "
    jobs = []
    for _ in range(min(params["parallel"], submitted)):
        worker = Worker(work_queue, result_queue, rawedit)
        worker.start()
        jobs.append(worker)
    for job in jobs:
        job.join()

    ## collect the results off the queue "
    outstats = open(params["work"]+"stats/s2.rawedit.txt", 'a')
    print >> outstats, "\t".join(["sample", "Nreads", "passed",
                                  "passed.w.trim", "passed.total"])
    stats = []
    for _ in range(submitted):
        stats.append(result_queue.get())

    stats.sort(key=lambda x: x[0])
    for job in range(submitted):
        a, b, c, d = stats[job]
        print >> outstats, "\t".join([a, b, c, d, str(int(c)+int(d))])

    print >>outstats, """
    Nreads = total number of reads for a sample
    passed = retained reads that passed quality filtering at full length
    passed.w.trim= retained reads that were trimmed due to detection of adapters
    passed.total  = total kept reads of sufficient length
    note: you can set the option in params file to include 
    trimmed reads of xx length.\n"""
    outstats.close()


if __name__ == "__main__":
    PARAMS = {}
    FASTQS = []
    QUIET = 0
    main(PARAMS, FASTQS, QUIET)
    