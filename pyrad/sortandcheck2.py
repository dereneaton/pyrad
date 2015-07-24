#!/usr/bin/env python2

""" demultiplex raw sequence data given a barcode map """

import gzip
import itertools
import sys
import glob
import os
import subprocess
import cPickle as pickle
import multiprocessing
from potpour import Worker


def combinefiles(filepath):
    """combines first and second reads file names """
    ## unpack seq files in filepath
    if len(glob.glob(filepath)) > 1:
        fastqs = glob.glob(filepath)
    else:
        fastqs = glob.glob(filepath)
    firsts = [i for i in fastqs if "_R1_" in i]

    ## check names
    if len(firsts) < 1:
        print "\n\tFirst read files names must contain '_R1_'."
        sys.exit()
    seconds = [ff.replace("_R1_", "_R2_") for ff in firsts]

    ## check pair matching
    if len(firsts) != len(seconds):
        print "Different numbers of first and second read files."+\
              " Check that the names of files are correct"
        sys.exit()
    return zip(firsts, seconds)


def revcomp(sequence):
    "returns reverse complement of a string"
    sequence = sequence[::-1].strip().replace("A", "t").replace("T", "a").\
                                 replace("C", "g").replace("G", "c").upper()
    return sequence


def matching(bcd, hit, maxmismatch):
    "allows for N base difference between barcodes"
    if len(bcd) == len(hit):
        diffs = [bcd[i] == hit[i] for i in range(len(bcd))]
        if diffs.count(False) <= maxmismatch:
            return 1
        else:
            return 0
    else:
        return 0


def unambig(seq):
    """ returns both resolutions of a cut site
    that has an ambiguous base in it """
    resos = []
    ambigs = {"R":("G", "A"),
              "K":("G", "T"),
              "S":("G", "C"),
              "Y":("T", "C"),
              "W":("T", "A"),
              "M":("C", "A")}
    for base in list("RKSYWM"):
        if base in seq:
            resos.append(seq.replace(base, ambigs[base][0]))
            resos.append(seq.replace(base, ambigs[base][1]))
    return resos


def findbcode(cut, longbar, l):
    """ find barcode sequence in the beginning of read """

    ## default barcode string
    barcode = 'N'*20

    ## in case ambiguous base in CUT "
    if any([i in cut for i in list("RKYSWM")]):
        cuts = unambig(cut)
        barcodes = [] 
        for icut in cuts:
            if l[1][0:longbar+len(icut)].count(icut) == 1:
                barcode = l[1].split(icut)[0].strip()
            elif l[1][0:longbar+len(icut)].count(icut) == 2:
                barcode = cut.join(l[1].split(icut)[0:2]).strip()
            else:
                barcode = ""
            barcodes.append(barcode)
        longestbar = barcodes[[len(i) for i in barcodes].\
                    index(max([len(j) for j in barcodes]))]
        return longestbar
    else:
        if l[1][0:longbar+len(cut)].count(cut) == 1:
            barcode = l[1].split(cut)[0].strip()
        elif l[1][0:longbar+len(cut)].count(cut) == 2:
            barcode = cut.join(l[1].split(cut)[0:2]).strip()
        else:
            barcode = ""
    return barcode


def barmatch(params, raws, bcdmap, longbar, num, localcut, quiet):
    """matches reads to barcodes in barcode file
    and writes to individual temp files, after all
    read files have been split, temp files are collated
    into .fq files"""
    ## C, Raws, CUT, datatype, num, maxmismatch, WORK, longB):

    ## counters for stats output
    locus = 0
    match = 0     ## ...?cut site found
    match2 = 0    ## ...?bar matches 
    barcodehits = set()

    ## dictionary to record barcode misses
    misses = {}
    misses['_'] = 0
    
    ## read in paired end read files"
    if 'pair' in params["datatype"]:
        if raws[0].endswith(".gz"):
            fr1 = gzip.open(raws[0])
        else:
            fr1 = open(raws[0])
        if raws[1].endswith(".gz"):
            fr2 = gzip.open(raws[1])
        else:
            fr2 = open(raws[1])
        ## create iterators to sample 4 lines at a time
        quart1 = itertools.izip(*[iter(fr1)]*4)
        quart2 = itertools.izip(*[iter(fr2)]*4)
    else:
        ## read in single end read file"
        if raws.endswith('.gz'):
            fr1 = gzip.open(raws)
        else:
            fr1 = open(raws)
        quart1 = itertools.izip(*[iter(fr1)]*4)

    ## dictionaries to store first and second reads
    dsort1 = {}  ## M
    dsort2 = {}

    ## go until end of the file
    while 1:
        try: 
            read1 = quart1.next()
        except StopIteration: 
            break

        ## match paired end reads together, for now
        if 'pair' in params["datatype"]:
            read2 = quart2.next()
            l = [r.strip() for r in read1]
            l = [l[0], l[1], l[2], l[3]]
            ll = [r.strip() for r in read2]
            ll = [ll[0], ll[1], ll[2], ll[3]]
        else:
            ## make list of four fastq line elements
            l = [r.strip() for r in read1]
            l = [l[0], l[1], l[2], l[3]]
        locus += 1

        if 'pair' in params["datatype"]:
            ## if all bars are the same length, 
            ## demultiplexing is really easy
            if longbar[1] == 'same':
                barcode = l[1][:longbar[0]]
            else:
                ## find barcodes in the reads"
                barcode = findbcode(localcut, longbar[0], l)

            ## tracker of number of occurrences of each barcode
            if barcode:
                if barcode in misses:
                    misses[barcode] += 1
                else:
                    misses[barcode] = 1

                ## exclude the read if no cutsite/barcode found"
                if barcode in dsort1:
                    l[1] = l[1][len(barcode):]  
                    l[3] = l[3][len(barcode):]
                    ## 
                    dsort1[barcode].append("\n".join(l).strip())
                    dsort2[barcode].append("\n".join(ll).strip())
                    match += 1
                else:
                    l[1] = l[1][len(barcode):] 
                    l[3] = l[3][len(barcode):]                        
                    dsort1[barcode] = l
                    dsort2[barcode] = ll
                    match += 1

            else:
                misses["_"] += 1

        else:
            ## NOT paired data
            if longbar[1] == 'same':
                if params["datatype"] == '2brad':
                    barcode = l[1][-longbar[0]:]
                else:
                    barcode = l[1][:longbar[0]]
            else:
                barcode = findbcode(localcut, longbar[0], l)

            ## tracker of number of occurrences of each barcode"
            if barcode:
                if barcode in misses:
                    misses[barcode] += 1
                else:
                    misses[barcode] = 1

                ## exclude the read if no cutsite/barcode found"
                if barcode in dsort1:
                    if params["datatype"] == '2brad':
                        l[1] = l[1][:-len(barcode)]
                        l[3] = l[3][:-len(barcode)]
                    else:
                        l[1] = l[1][len(barcode):]
                        l[3] = l[3][len(barcode):]
                    dsort1[barcode].append("\n".join(l).strip())
                    match += 1
                else:
                    l[1] = l[1][len(barcode):]
                    l[3] = l[3][len(barcode):]                        
                    dsort1[barcode] = l
                    match += 1
            else:
                misses["_"] += 1


        ## write to file every 500Kth read"
        ## only writes reads that match to a barcode in 
        ## C by less than some N differences "
        if not locus % 500000:
            barcodehits, match2, dsort1, dsort2 = writetofile(params,
                                                      bcdmap, barcodehits,
                                                      dsort1, dsort2, 
                                                      match2, num)

    ## write the remaining reads to file"
    barcodehits, match2, dsort1, dsort2 = writetofile(params,
                                                      bcdmap, barcodehits,
                                                      dsort1, dsort2, 
                                                      match2, num)

    if not quiet:
        sys.stderr.write(".")
    fr1.close()
    if 'pair' in params["datatype"]:
        fr2.close()
                        
    ## writes statistics out
    statout = open(params["work"]+"stats/s1.sorting.txt", 'a')
    if 'pair' in params["datatype"]:
        outname = raws[0].split("/")[-1].replace("_R1_", "_")
    else:
        outname = raws.split("/")[-1].replace("_R1_", "_")

    match2 = sum([misses[i] for i in misses if i in barcodehits])
    writeit = "%s\t%i\t%i\t%i\n" % (outname, locus, match, match2)
    statout.write(writeit)
    statout.close()
    pickout = open(params["work"]+"fastq/."+outname+".pickle", "wb")
    pickle.dump(misses, pickout)
    pickout.close()



def writetofile(params, bcdmap, barcodehits, dsort1,
                dsort2, match2, num):
    """ write to file duh """

    for bcd in bcdmap:
        outf1 = gzip.open(params["work"]+"fastq/."+bcdmap[bcd]+\
                          '.temp_R1_'+str(num)+'.gz', 'ab')

        if 'pair' in params["datatype"]:
            outf2 = gzip.open(params["work"]+"fastq/."+bcdmap[bcd]+\
                             '.temp_R2_'+str(num)+'.gz', 'ab')

        for barcode in dsort1:
            if matching(bcd, barcode, params["maxmismatch"]):
                barcodehits.add(barcode)

                if dsort1[barcode]:
                    match2 += len(dsort1[barcode])
                    outf1.write("\n".join(dsort1[barcode])+'\n')

                if 'pair' in params["datatype"]:
                    if dsort2[barcode]:
                        outf2.write("\n".join(dsort2[barcode])+'\n')
                ## remove seqs that are written to file                        
                dsort1[barcode] = []
                dsort2[barcode] = []
        ## close the read1 file
        outf1.close()
        if 'pair' in params["datatype"]:
            ## close the read2 file
            outf2.close()
        ## remove seqs that are not written to file
        dsort1[bcd] = []
        dsort2[bcd] = []
    return barcodehits, match2, dsort1, dsort2



def writefunc(params, quiet, localcut):
    "create barcode dictionary"

    ## get barcode map
    codetable = open(params["bcode"], 'r')
    codes = [line.strip().split() for line in codetable.readlines()]
    bcdmap = {}
    for line in codes:
        if line:
            bcdmap[line[1].strip().upper()] = line[0]

    ## find longest barcode "
    keylens = [len(i) for i in bcdmap.keys()]
    if len(set(keylens)) == 1:
        longbar = (keylens[0], 'same')
    else:
        longbar = (max(keylens), 'diff')

    ## check for CUT in barcodes "
    cuts = unambig(localcut)

    if len(cuts) > 1:
        for cut in cuts:
            if any([cut in i for i in bcdmap.keys()]):
                print "\n\twarning: CUT site matches within one"+\
                      " of the barcodes, I suggest double \n\tchecking"+\
                      " the file to make sure it properly demultiplexes"
    else:
        if any([localcut in i for i in bcdmap.keys()]):
            print "\n\twarning: CUT site matches within one of "+\
                  "the barcodes, I suggest double \n\tchecking the "+\
                  "file to make sure it properly demultiplexes"

    ## read in sequence files "
    if len(glob.glob(params["glob"])) > 1:
        fastqs = glob.glob(params["glob"]) ## [f for f in glob.glob(GLOB)]
    else:
        fastqs = glob.glob(params["glob"])
    if 'pair' in params["datatype"]:
        raws = combinefiles(params["glob"])
    else:
        raws = fastqs

    ## send jobs to multiprocess queue"
    num = 0
    work_queue = multiprocessing.Queue()
    submitted = 0
    for fastq in raws:
        if 'pair' in params["datatype"]:
            work_queue.put([params, fastq, bcdmap,
                            longbar, num, localcut, quiet])
            submitted += 1
        else:
            work_queue.put([params, fastq, bcdmap,
                            longbar, num, localcut, quiet])
            submitted += 1
        num += 1

    result_queue = multiprocessing.Queue()

    ## spawn workers, give function"
    jobs = []
    for _ in range(min(params["parallel"], submitted)):
        worker = Worker(work_queue, result_queue, barmatch)
        worker.start()
        jobs.append(worker)
    for job in jobs:
        job.join()

    ## a dictionary for ...?
    allpics = {}    ## Ms

    ## combine pickled results from all files
    if len(glob.glob(params["work"]+"fastq/.*.pickle")) > 1:
        for pick in glob.glob(params["work"]+"fastq/.*.pickle"):
            pickin = open(pick, "rb")
            picdic = pickle.load(pickin)  ## M
            pickin.close()
            for key in picdic:
                if key not in allpics:
                    allpics[key] = picdic[key]
                else:
                    allpics[key] += picdic[key]
            os.remove(pick)
    elif len(glob.glob(params["work"]+"fastq/.*.pickle")) == 1:
        pick = glob.glob(params["work"]+"fastq/.*.pickle")[0]
        pickin = open(pick, 'rb')
        allpics = pickle.load(pickin)
        pickin.close()
        os.remove(pick)
    else:
        print "\nno stats file generated"

    allkeys = allpics.keys()
    allkeys.sort(key=lambda x: allpics[x], reverse=True)

    statout = open(params["work"]+"stats/s1.sorting.txt", 'a')
    statout.write("\n\n")
    statout.write("sample\ttrue_bar\tobs_bars\tN_obs\n")

    bcdnames = bcdmap.keys()
    bcdnames.sort()
    try: maxl = max(map(len,map(str, allpics.values())))
    except ValueError: maxl = 5

    hits = []
    for bcd in bcdnames:
        for barcode in allkeys:
            if matching(bcd, barcode, params["maxmismatch"]):
                print >>statout, "%s    \t%s    \t%s\t%s" % (
                                       bcdmap[bcd],
                                       bcd, barcode,
                                       str(allpics[barcode])+\
                                       " "*(maxl+3-len(str(allpics[barcode]))))
                hits.append(barcode)

    ## look back at old one to find the diff btwn allpics and allkeys (M, Ms?)

    statout.write("\n")
    maxl = max([len(i) for i in allkeys])
    for barcode in allkeys:
        if barcode not in hits:
            print >>statout, "nomatch  \t%s    \t%i" % (barcode+\
                                               ""*(maxl+3-len(barcode)), 
                                               allpics[barcode])
    statout.close()



def main(params, quiet):
    """ the main script """

    ## seperate double digest cut sites, only need first read one for now "
    if "," in params["cut"]:
        localcut = params["cut"].split(",")[0]
    else:
        localcut = params["cut"]

    if not glob.glob(params["glob"]):
        sys.exit("\tNo data found in "+params["glob"]+\
                 ". \n\tFix path to the data files\n")

    ## check for previous output"
    if not os.path.exists(params["work"]+'stats'):
        os.makedirs(params['work']+'stats')
    if os.path.exists(params["work"]+'fastq'):
        if os.listdir(params["work"]+'fastq'):
            sys.exit("\n\tfastq/ directory in working directory contains"+\
                      " data, move/remove it before running step 1\n")
    else:
        os.makedirs(params["work"]+'fastq')

    if "*" in params["bcode"]:
        if len(glob.glob(params["bcode"])) == 1:
            params["bcode"] = glob.glob(params["bcode"])[0]

    if not quiet:
        sys.stderr.write("\n\tstep 1: sorting reads by barcode\n\t ")

    ## write stats output
    with open(params["work"]+"stats/s1.sorting.txt", 'w') as statout:
        statout.write("\t".join(["file  ", "Nreads",
                      "cut_found", "bar_matched"])+"\n")

    ## do barcode sorting
    writefunc(params, quiet, localcut)
    names = [line.split()[0] for line \
             in open(params["bcode"]).readlines() if line.strip()]

    ## concatenate temp files "
    for name in names:
        cmd = "cat "+params["work"]+"fastq/."+name+".temp_R1_* > "+\
                     params["work"]+"fastq/"+name+"_R1.fq.gz"
        subprocess.call(cmd, shell=True)
        if "pair" in params["datatype"]:
            cmd = "cat "+params["work"]+"fastq/."+name+".temp_R2_* > "+\
                         params["work"]+"fastq/"+name+"_R2.fq.gz"
            subprocess.call(cmd, shell=True)

    ## remove tempfiles
    for removefile in glob.glob(params["work"]+"fastq/.*"):
        os.remove(removefile)

if __name__ == "__main__":
    main()
    