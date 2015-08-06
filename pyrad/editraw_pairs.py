#!/usr/bin/env python2

""" filter and edit reads for paired data """

import multiprocessing
import itertools
import sys
import os
import glob
import operator
import gzip
from editraw_rads import unambar
from potpour import Worker
from cluster_cons7_shuf import fullcomp


def afilter(cut, seq, strict, read):
    read1 = (read == 1)
    check1 = check2 = wheretocut = None
    ## lookfor cut site "

    if unambar(cut):
        ## if ambiguity in cutter "        
        cuta, cutb = unambar(cut)
        if strict == 2:
            if read1:
                lookfor1 = cuta+"A"
                lookfor2 = cutb+"A"
            else:
                lookfor1 = cuta
                lookfor2 = cutb
        else:
            if read1:
                lookfor1 = cuta+"AGAT"
                lookfor2 = cutb+"AGAT"
            else:
                lookfor1 = "A"*50
                lookfor2 = "A"*50
        if lookfor1 in seq:
            check1 = seq.rindex(lookfor1)
        if lookfor2 in seq:
            check2 = seq.rindex(lookfor2)
        if check1 or check2:
            wheretocut = min([i for i in [check1, check2] if i])
        else:
            wheretocut = None
    else:
        ## no ambiguity in cutter "
        if strict == 2:
            if read1:
                lookfor1 = cut+"A"
            else:
                lookfor1 = cut
        else:
            if read1:
                lookfor1 = cut+"AGA"
            else:
                lookfor1 = "A"*50
        if lookfor1 in seq:
            wheretocut = seq.rindex(lookfor1)
        else:
            wheretocut = None
                
    ## look for adapter sequence "
    if not wheretocut:
        if strict == 2:
            lookfor1 = "AGATCG"
        else:
            lookfor1 = "AGATCGGA"
        if lookfor1 in seq:
            if read1:
                wheretocut = seq.rindex(lookfor1)-(len(cut)+1)
            else:
                wheretocut = seq.rindex(lookfor1)-(len(cut)+6)
        else:
            wheretocut = None

    ## look for CUT and end of seq 
    if not wheretocut:
        if cut in seq[-(len(cut)+5):]:
            wheretocut = seq.rindex(cut)
    return wheretocut


def rawedit(params, infile, quiet):
    """ three functions:
    (1) replaces low quality base calls with Ns,
    (2) checks for adapter sequence and xxbarcodesxx if strict set to 1 or 2 
    (3) concatenate paired reads with a separator and write to file """

    ## get cutters
    if "," in params["cut"]:
        cut1, cut2 = params["cut"].strip().split(",")
        cut2 = fullcomp(cut2)
    else:
        cut1 = params["cut"]
        cut2 = fullcomp(cut1)

    ## create iterators for R1 and R2 files "
    if infile.endswith(".gz"):
        dems1 = gzip.open(infile, 'rb')
        if ".forward." in infile:
            dems2 = gzip.open(infile.replace(".forward.", ".reverse."), 'r')
        else:
            dems2 = gzip.open(infile.replace("_R1.", "_R2."), 'r')
    else:
        dems1 = open(infile, 'r')
        if ".forward." in infile:
            dems2 = open(infile.replace(".forward.", ".reverse."), 'r')
        else:
            dems2 = open(infile.replace("_R1.", "_R2."), 'r')
    name = str(infile.split('/')[-1])
    while name.split(".")[-1] in ["fastq", "fastQ", "gz", "fq",
                                  "FastQ", "nomerge"]:
        name = name.replace('.'+name.split(".")[-1], "")
    if '.forward' in name:
        name = name.split(".forward")[0]
    else:
        name = name.replace("_R1", "")

    k1 = itertools.izip(*[iter(dems1)]*4)
    k2 = itertools.izip(*[iter(dems2)]*4)
    writing_r = []
    writing_c = []

    orig = keep = keepcut = 0
    handle = params["work"]+'edits/'+str(name)+".edit"

    ## iterate over paired reads, edits 1st, if OK, append both to .edit file"
    while 1:
        try: 
            itera1 = k1.next()
        except StopIteration: 
            break
        itera2 = k2.next()

        orig += 1
        iseq = itera1[1].strip()
        qscore = [ord(i) for i in itera1[3].strip()]
        offset = int(params["Q"])
        phred = [x-offset for x in qscore]
        seq = ["N"]*len(phred)

        ## fix cut sites to be error free
        liseq = list(iseq)
        liseq[:len(cut1)] = list(cut1)
        iseq = "".join(liseq)
        if "merge" in params["datatype"]:
            iseq[-len(cut2):] = fullcomp(cut2)

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
                if unambar(cut1)[1]:
                    seq[base] = unambar(cut1)[0][base]
                else:
                    seq[base] = cut1[base]

        sseq1 = "".join(seq)
        wheretocut1 = None

        ## apply filters for adapter sequences "
        ## if GBS CUT2 = revcomp(CUT1)   ex: CTGCA"
        ## if ddRAD CUT2 = revcomp(CUT2) ex: AATT "
        if params["strict"]:
            wheretocut1 = afilter(cut2, sseq1, params["strict"], 1)
        
        ## max allowed Ns
        if sseq1.count("N") <= int(params["maxN"]):
            ## if trimmed read1 length atleast t
            if len(sseq1) >= max(32, params["trimkeep"]):     
                ## first read is (maybe) good, now filter second reads "
                iseq = itera2[1].strip()
                qscore = [ord(i) for i in itera2[3].strip()]
                ## if PEAR filtered then seqs are revcomps "
                if '.forward' in infile:
                    iseq = revcomp(iseq)
                    qscore = qscore[::-1]
                offset = int(params["Q"])
                phred = [x-offset for x in qscore]
                seq = ["N"]*len(phred)
                for base in range(len(phred)):
                    ## don't quality check cut site                    
                    if base > len(cut2):
                        ## quality threshold                        
                        if phred[base] >= 20:          
                            try: 
                                seq[base] = iseq[base]
                            except IndexError: 
                                pass
                        else:
                            seq[base] = "N"
                    else:
                        try:
                            seq[base] = iseq[base]
                        except IndexError: 
                            pass
                sseq2 = "".join(seq)

                ## filter for gbs read2s, b/c they will be clustered"
                badread = 0
                if params["datatype"] == "pairgbs":
                    sseq2 = sseq2[:len(sseq1)]
                    if sseq2.count("N") > params["maxN"]:
                        badread = 1

                ## apply adapter filter to read 2
                wheretocut2 = None
                if params["strict"]:
                    wheretocut2 = afilter(revcomp(cut1), sseq2, 
                                          params["strict"], 2)

                if wheretocut1 and wheretocut2:
                    cutter = min(wheretocut1, wheretocut2)
                else:
                    cutter = max(wheretocut1, wheretocut2)
                ## extra strict check for undigested cutsites
                if params["strict"]:
                    if not cutter:
                        if (cut1 in sseq2[-16:]) or \
                           (cut2 in sseq1[-10:]):
                            cutter = len(sseq1)-16

                if not badread:
                    if cutter:
                        ## second read was trimmed
                        if cutter > max(32, params["trimkeep"]):
                            ## include only the first read, with an 
                            ## N placeholder for read2
                            ## since it was trimmed off
                            sout = ">"+name+"_"+str(keepcut)+\
                                   "_trim1"+"\n"+sseq1[:cutter]+\
                                   "nnnnN\n"
                            writing_c.append(sout)
                            keepcut += 1
                            ## cannot keep trimmed second read in 
                            ## pairddrad method but can in pairgbs
                            if params["datatype"] == 'pairgbs':
                                sout = ">"+name+"_"+str(keepcut)+\
                                       "_trim2\nNnnnn"+\
                                       fullcomp(sseq2[x:cutter+5])+"\n"
                                writing_c.append(sout)
                                keepcut += 1
                    else:
                        ## second read is good, not trimmed
                        sout = ">"+name+"_"+str(keep)+\
                               "_pair"+"\n"+sseq1[:-1]+"nnnn"+\
                               fullcomp(sseq2[x:])+"\n"
                        writing_r.append(sout)
                        keep += 1
                else:
                    ##...
                    pass
            else:
                ## ...
                pass

        if not orig % 50000:
            ## writes only full length reads "
            with open(params["work"]+'edits/'+\
                      str(name)+".edit", 'a') as outfile:
                outfile.write("".join([z for z in writing_r]))
            ## writes only full length reads "
            with open(params["work"]+'edits/'+\
                      str(name)+".edit", 'a') as outfile:
                outfile.write("".join([z for z in writing_c]))
            writing_r = []
            writing_c = []

    ## writes only full length reads
    with open(params["work"]+'edits/'+str(name)+".edit", 'a') as outfile:
        outfile.write("".join([z for z in writing_r]))
    ## writes only full length reads "
    with open(params["work"]+'edits/'+str(name)+".edit", 'a') as outfile:
        outfile.write("".join([z for z in writing_c]))
    writing_r = []
    writing_c = []

    dems1.close()
    dems2.close()
    if not quiet:
        sys.stderr.write(".")

    if params["datatype"] == 'pairgbs':
        keepcut = keepcut*2

    return [handle.split("/")[-1].replace(".edit", ""), 
            str(orig), str(keepcut), str(keep)]



def main(params, fastqs, quiet): 
    """ call the main script """
    ##Parallel, WORK, FQs, CUT, pN, Q, strict, trimkeep, datatype):

    if not quiet: 
        sys.stderr.write("\n\tstep 2: quality filtering \n\t")

    ## create output directories
    if not os.path.exists(params["work"]+'stats'):
        os.makedirs(params["work"]+'stats')
    if not os.path.exists(params["work"]+'edits'):
        os.makedirs(params["work"]+'edits')

    ## load up work queue
    submitted = 0
    work_queue = multiprocessing.Queue()

    ## do not select merged or discarded reads if PEAR was used on data
    fastqs = glob.glob(fastqs)
    fastqs = [i for i in fastqs if not \
              any([j in i for j in ["discarded", ".assembled."]])]
    
    ## if not just one file
    if len(fastqs) > 1:
        ## subselect only the first reads
        if any([".unassembled.forward." in i for i in fastqs]):
            infiles = [i for i in fastqs if '.forward.' in i] ## FS
        else:
            infiles = [i for i in fastqs if '_R1.' in i]
        
        ## order files by size
        for i in range(len(infiles)):
            statinfo = os.stat(infiles[i])
            infiles[i] = infiles[i], statinfo.st_size
        infiles.sort(key=operator.itemgetter(1))
        infiles = [i[0] for i in infiles][::-1]

        ## submit jobs to queue
        for handle in infiles:
            ## trim name endings
            name = handle.split('/')[-1]
            while name.split(".")[-1] in ["fastq", "fastQ", "gz", 
                                          "fq", "FastQ", "nomerge"]:
                name = name.replace('.'+name.split(".")[-1], "")
            if '.forward.' in name:
                name = name.split(".forward")[0]
            else:
                "_".join(name.split('_R')[:-1])

            if params["work"]+"edits/"+name+".edit" \
               not in glob.glob(params["work"]+"edits/*"):
                if os.stat(handle).st_size > 0:     ## exclude empty files
                    args = [params, handle, quiet]
                    work_queue.put(args)
                    submitted += 1
                else:
                    print 'skipping', handle, ", file is empty"
            else:
                print "\t"+name+'.edit'+" already in edits/"
    elif len(fastqs) == 1:
        ## if only one file "
        work_queue.put([params, glob.glob(infiles)[0], quiet])
        submitted += 1

    else:
        print "no _paired_ de-multiplexed files found in this location."
        sys.exit()

    ## create a queue to pass to workers to store the results "
    result_queue = multiprocessing.Queue()

    ## spawn workers, give function "
    jobs = []
    for i in range(min(params["parallel"], submitted)):
        worker = Worker(work_queue, result_queue, rawedit)
        worker.start()
        jobs.append(worker)
    for job in jobs:
        job.join()

    ## collect the results off the queue
    outstats = open(params["work"]+"stats/s2.rawedit.txt", 'a')
    print >> outstats, "\t".join(["sample", "Nreads", "exclude",
                                  "trimmed", "passed"])
    for i in range(submitted):
        a, b, c, d = result_queue.get()
        print >> outstats, "\t".join([a,b, str(int(b)-int(d)), c, d])

    print >>outstats, """
    Nreads = total number of reads for a sample
    exclude = reads that were excluded
    trimmed = reads that had adapter trimmed but were kept
    passed = total kept reads
    """
    outstats.close()

if __name__ == "__main__":
    main()
