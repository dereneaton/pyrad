
import glob
import multiprocessing
import gzip
import subprocess
import sys
import os
import itertools
from potpour import Worker


def mergepairs(WORK, UCLUST, handle, match, Q):

    handle1 = handle
    handle2 = handle.replace("_R1.","_R2.")
    outfile = handle.replace("_R1.","M.")

    while outfile.split(".")[-1] in ["fastq","fastQ","gz","fq","FastQ"]:
        outfile = outfile.replace('.'+outfile.split(".")[-1], "")
    outfile = outfile.split("/")[-1]
    outfile = WORK+"mergedreads/"+outfile+".fq"

    if [handle1 and handle2]:
        if ".gz" in handle1[-4:]:
            k1 = itertools.izip(*[iter(gzip.open(handle1))]*4)
            k2 = itertools.izip(*[iter(gzip.open(handle2))]*4)
            thandle1 = WORK+"mergedreads/"+handle1.split("/")[-1].replace(".gz",".temp2")
            thandle2 = WORK+"mergedreads/"+handle2.split("/")[-1].replace(".gz",".temp2")
            numout1 = open(thandle1, 'w')
            numout2 = open(thandle2, 'w')
        else:
            k1 = itertools.izip(*[iter(open(handle1))]*4)
            k2 = itertools.izip(*[iter(open(handle2))]*4)
            thandle1 = WORK+"mergedreads/"+handle1.split("/")[-1]+".temp2"
            thandle2 = WORK+"mergedreads/"+handle2.split("/")[-1]+".temp2"
            numout1 = open(thandle1, 'w')
            numout2 = open(thandle2, 'w')
    else:
        print "pair missing"
        sys.exit()

    N1 = []
    N2 = []
    cnt = 0
    
    while 1:
        try: d = k1.next()
        except StopIteration: break
        e = k2.next()
        N1.append("".join([d[0].strip()+"_"+str(cnt)+"\n",d[1],d[2],d[3]]))
        N2.append("".join([e[0].strip()+"_"+str(cnt)+"\n",e[1],e[2],e[3]]))
        cnt+=1
        if not cnt % 50000:
            numout1.write("".join(N1))
            numout2.write("".join(N2))
            N1 = []
            N2 = []
    numout1.write("".join(N1))
    numout2.write("".join(N2))
    numout1.close()
    numout2.close()

    cmd = UCLUST+\
          " -fastq_mergepairs "+thandle1 +\
          " -reverse "+thandle2 +\
          " -fastq_maxdiffs 6 " +\
          " -fastq_truncqual 2 " +\
          " -fastq_minlen 36 " +\
          " -fastq_minmergelen 50 "+\
          " -fastqout "+outfile +\
          " -fastq_allowmergestagger" +\
          " -quiet "
    subprocess.call(cmd, shell=True)


    stats = statsout(thandle1, thandle2, outfile, WORK)
    sys.stderr.write(".")
    return stats


def statsout(h1,h2,m,WORK):
    " remove merged reads from 1st & 2nd read files "

    " stat counters "
    cnt = 0
    mcnt = 0

    " create list of merged IDs "
    MIDS = []
    if os.path.exists(m): 
        merged = open(m, 'r') 
        for line in itertools.izip(*[iter(merged)]*4):
            MIDS.append(int(line[0].strip().split("_")[-1]))
        merged.close()
    ## if not...

    if ".gz" in h1[-5:]:
        hand1 = gzip.open(h1, 'rb')
        hand2 = gzip.open(h2, 'rb')
    else:
        hand1 = open(h1, 'r')
        hand2 = open(h2, 'r')
        
    r1 = itertools.izip(*[iter(hand1)]*4)
    r2 = itertools.izip(*[iter(hand2)]*4)

    " lists to write output "
    ONE = []
    TWO = []

    " outfile names for mergeless reads "
    outname = WORK+"fastq/"+h1.split("/")[-1].replace(".temp2",".nomerge")+".gz"

    if os.path.exists(outname):
        os.remove(outname)
    outname2 = outname.replace("_R1.","_R2.")
    if os.path.exists(outname2):
        os.remove(outname2)

    while 1:
        try: one = r1.next()
        except StopIteration: break
        two = r2.next()
        cnt += 1
        find = int(one[0].strip().split("_")[-1])
        if MIDS:
            if find == MIDS[0]:
                "reads were merged, don't write to file"
                mcnt += 1
                MIDS.pop(0)
            else:
                ONE.append(one) #[i.strip() for i in one])
                TWO.append(two) #[i.strip() for i in two])
        else:
            ONE.append(one) #[i.strip() for i in one])
            TWO.append(two) #[i.strip() for i in two])

        if not cnt % 10000:
            outfile = gzip.open(outname, 'ab')
            outfile.write("".join(["".join(i) for i in ONE]))
            outfile.close()
            outfile2 = gzip.open(outname2, 'ab')
            outfile2.write("".join(["".join(i) for i in TWO]))
            outfile2.close()    
            ONE = []
            TWO = []


    if os.path.exists(h1):
        cmd1 = "/bin/rm "+h1
        cmd2 = "/bin/rm "+h2
        subprocess.call(cmd1, shell=True)
        subprocess.call(cmd2, shell=True)

    outfile = gzip.open(outname, 'ab')
    outfile.write("".join(["".join(i) for i in ONE]))
    outfile.close()
    outfile2 = gzip.open(outname2, 'ab')
    outfile2.write("".join(["".join(i) for i in TWO]))
    outfile2.close()
    sys.stderr.write(".")
    return [outname,mcnt]

        

def main(WORK, UCLUST, FQs, match, Q, Parallel):

    " create output directories " 
    if not os.path.exists(WORK+'fastq/'):
        os.makedirs(WORK+'fastq')
    if not os.path.exists(WORK+'mergedreads'):
        os.makedirs(WORK+'mergedreads')
    if not os.path.exists(WORK+'stats'):
        os.makedirs(WORK+'stats')


    submitted = 0
    work_queue = multiprocessing.Queue()

    names = [i for i in glob.glob(FQs) if "_R1.fq" in i]

    " submit jobs to queue "
    if len(names) > 1:
        for handle in names:
            if "nomerge." not in handle:
                n = str(handle.split('/')[-1]).replace("_R1.",".")
                while n.split(".")[-1] in ["fastq","fastQ","gz","fq","FastQ"]:
                    n = n.replace('.'+n.split(".")[-1], "")
                finder = WORK+'edits/'+n+".edit"
                if finder not in glob.glob(WORK+"edits/*"):
                    if os.stat(handle).st_size > 0:   ## exclude empty files
                        if os.path.exists(handle.replace("_R1.","_R2.")):
                            if not os.path.exists(handle.replace(".fq",".nomerge.fq")):
                                args = [WORK, UCLUST, handle, match, Q]
                                work_queue.put(args)
                                submitted += 1
                            else:
                                print "merge file already created for", handle.split("/")[-1]
                        else:
                            print "cannot find 2nd read file for", handle.split("/")[-1]
                    else:
                        print "\t"+finder+" already in edits/"
    else:
        if not names:
            if [i for i in glob.glob(FQs) if "_R1_." in i]:
                print "\n\tfile names should have _R1. not _R1_."
            print "\n\tcannot find input files"
            sys.exit()
        else:
            work_queue.put([WORK, UCLUST, names[0], match, Q])
            submitted += 1

    " create a queue to pass to workers to store the results "
    result_queue = multiprocessing.Queue()


    " spawn workers, give function "
    jobs = []
    for i in range( min(Parallel,submitted) ):
        worker = Worker(work_queue, result_queue, mergepairs)
        worker.start()
        jobs.append(worker)
    for job in jobs:
        job.join()

    if submitted > 0:
        statout = open(WORK+"stats/s2.mergedreads.txt",'w')
        print >>statout, "\t".join(["taxon","mergedreads"])

        for i in range(submitted):
            stat = result_queue.get()
            a,b = stat
            n = a.strip().split("/")[-1].replace(".nomerge.gz","")
            print >>statout, "\t".join([n,str(b)])
        print >>statout, "\nmerged reads written to", WORK+"mergedreads/ "
        statout.close()
