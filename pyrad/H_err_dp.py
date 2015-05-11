#!/usr/bin/env python2

""" jointly infers heterozygosity and error rate from stacked sequences """

import scipy.stats
import scipy.optimize
import numpy as np
import itertools
import sys
import glob
import multiprocessing
import os
import gzip
from potpour import Worker

try:
    from collections import OrderedDict, Counter
except ImportError:
    from ordereddict import OrderedDict, Counter


def get_freqs(stack):
    """ returns a list as frequencies for ATGC"""
    sump = sum([sum(cc.values()) for cc in stack])
    #sump = sum([sum(i) for i in site])
    totalcount = Counter()
    for stackcount in stack:
        totalcount += stackcount
    basefreqs = np.array([totalcount["A"],
                          totalcount["T"],
                          totalcount["G"],
                          totalcount["C"]])/float(sump)
    return basefreqs


def likelihood1(errors, base_frequencies, stacks):
    """probability homozygous"""
    ## make sure base_frequencies are in the right order
    #print uniqstackl.sum()-uniqstack, uniqstackl.sum(), 0.001
    tot = np.array([stacks.sum(axis=1)]*4).T
    b = scipy.stats.binom.pmf(tot-stacks, tot, errors)
    return np.sum(base_frequencies*b, axis=1)


def likelihood2(errors, base_frequencies, stacks):
    """probability of heterozygous"""

    returnlist = []
    for stackl in stacks:
        spair = list(itertools.combinations(stackl, 2))
        bpair = list(itertools.combinations(base_frequencies, 2))
        one = 2.*np.product(bpair, axis=1)
        tot = stackl.sum() #np.sum(spair, axis=1)
        atwo = tot - np.array([i[0] for i in spair]) -\
                     np.array([i[1] for i in spair])
        two = scipy.stats.binom.pmf(atwo, tot, (2.*errors)/3.)
        three = scipy.stats.binom.pmf(np.array([i[0] for i in spair]),
                                      np.array([i[0]+i[1] for i in spair]),
                                      0.5)
        four = 1.-np.sum(base_frequencies**2)
        returnlist.append(np.sum(one*two*(three/four)))
    return np.array(returnlist)

    # total = sum(uniqstackl)
    # for num, thisbase in enumerate(uniqstackl):
    #     for j, k in enumerate(uniqstackl):
    #         if j > num:
    #             one = 2.*base_frequencies[num]*base_frequencies[j]
    #             two = scipy.stats.binom.pmf(total-thisbase-k, 
    #                                         total,
    #                                         (2.*errors)/3.)
    #             three = scipy.stats.binom.pmf(thisbase, k+thisbase, 0.5)
    #             four = 1.-(sum([q**2. for q in base_frequencies]))
    #             hetero.append(one*two*(three/four))
    # return sum(hetero)



def get_diploid_lik(starting_params, base_frequencies, stacks, stackcounts):
    """ Log likelihood score given values [H,E] """
    hetero = starting_params[0]
    errors = starting_params[1]
    if (hetero <= 0.) or (errors <= 0.):
        score = np.exp(100)
    else:
        ## get likelihood for all sites
        lik1 = (1.-hetero)*likelihood1(errors, base_frequencies, stacks)
        lik2 = (hetero)*likelihood2(errors, base_frequencies, stacks)
        liks = lik1+lik2
        logliks = np.log(liks[liks>0])*stackcounts[liks>0]
        #logliks = np.log(liks)*stackcounts
        score = -logliks.sum()
    return score


def get_haploid_lik(errors, base_frequencies, tabled_stacks):
    """ Log likelihood score given values [E] """
    hetero = 0.
    listofliks = []
    if errors <= 0.:
        score = np.exp(100)
    else:
        for uniqstack in tabled_stacks:
            loglik = ((1.-hetero)*\
                     likelihood1(errors, base_frequencies, uniqstack))+\
                     (hetero*likelihood2(errors, base_frequencies, uniqstack))
            if loglik > 0:
                listofliks.append(tabled_stacks[uniqstack]*np.log(loglik))
        score = -sum(listofliks)
    return score


def table_c(stack):
    """ makes a dictionary with counts of base counts [x,x,x,x]:x,
        greatly speeds up Likelihood calculation"""
    countedstacks = []
    for stackcount in stack:
        ## convert to string for use with counter
        countedstacks.append(str([stackcount["A"],
                                  stackcount["T"],
                                  stackcount["G"],
                                  stackcount["C"]]))
    return Counter(countedstacks)


def consensus(params, handle, cut1, cut2):
    """ makes a list of lists of reads at each site """
    #f, minsamp, CUT1, CUT2, datatype):
    infile = gzip.open(handle)
    duo = itertools.izip(*[iter(infile)]*2)
    stacked = []
    while 1:
        try:
            first = duo.next()
        except StopIteration:
            break
        itera = [first[0], first[1]]
        thisstack = []
        rights = []
        lefts = []
        leftjust = rightjust = None
        while itera[0] != "//\n":
            nreps = int(itera[0].strip().split(";")[1].replace("size=", ""))

            ## record left and right most for cutting if gbs merge data "
            if params["datatype"] in ['merged', 'gbs', 'pairgbs']:
                if itera[0].strip().split(";")[-1] == "":
                    leftjust = itera[1].index([i for i in itera[1] \
                                        if i not in list("-N")][0])
                    rightjust = itera[1].rindex([i for i in itera[1] \
                                        if i not in list("-N")][0])
                lefts.append(itera[1].index([i for i in itera[1] \
                                      if i not in list("-N")][0]))
                rights.append(itera[1].rindex([i for i in itera[1] \
                                        if i not in list("-N")][0]))

            ## append sequence * number of dereps
            for _ in range(nreps):
                thisstack.append(tuple(itera[1].strip()))
            itera = duo.next()

        ## trim off overhang edges of gbs reads "
        if params["datatype"] in ['merged', 'gbs', 'pairgbs']:
            if any([i < leftjust for i in lefts]):
                rightjust = min(rights)
            if any([i < rightjust for i in rights]):
                leftjust = max(lefts)

            for seq in range(len(thisstack)):
                if rightjust:
                    thisstack[seq] = thisstack[seq][leftjust:rightjust+1]
                if leftjust:
                    thisstack[seq] = thisstack[seq][leftjust:rightjust+1]

        ## trim off restriction sites from end/s
        if params["datatype"] in ['merged', 'pairddrad', 'pairgbs', 'gbs']:
            for seq in range(len(thisstack)):
                thisstack[seq] = thisstack[seq][len(cut1):-(len(cut2)+1)]
        else:
            for seq in range(len(thisstack)):
                thisstack[seq] = thisstack[seq][len(cut1):]
        
        if len(thisstack) >= params["minsamp"]:
            arrayed = np.array(thisstack)
            ## make list for each site in sequences
            res = [Counter(seq) for seq in arrayed.T]
            ## exclude sites with indels
            stacked += [i for i in res if "-" not in i]
    return stacked



def optim(params, handle, cut1, cut2, quiet):
    """ fun scipy optimize to find best parameters"""
    ##WORK,handle, minsamp, CUT1, CUT2, datatype, haplos):
    name = handle.split("/")[-1].replace(".clustS.gz", "")

    ## make a list of Counter objects for each site in each stack
    stacked = consensus(params, handle, cut1, cut2)

    ## get base frequencies
    base_frequencies = get_freqs(stacked)

    ## get tabled counts of base patterns
    tabled_stacks = table_c(stacked)

    def toarray(uniqstack):
        """ converts string lists to arrays"""
        return np.array([int(i) for i in uniqstack[1:-1].strip().split(',')])

    stacks = np.array([toarray(i) for i in tabled_stacks.keys()])
    stackcounts = np.array(tabled_stacks.values())

    del stacked

    ## if data are haploid fix H to 0
    if params["haplos"] == 1:
        starting_params = [0.001]
        hetero = 0.
        errors = scipy.optimize.fmin(LL_haploid, 
                                starting_params,
                                (base_frequencies, 
                                    stacks,
                                    stackcounts),
                                disp=False,
                                full_output=False)
    else:
        starting_params = [0.01, 0.001]
        hetero, errors = scipy.optimize.fmin(get_diploid_lik,
                                             starting_params,
                                             (base_frequencies, 
                                                stacks,
                                                stackcounts),
                                             disp=False,
                                             full_output=False)
    outfile = open(params["work"]+"stats/."+name+".temp", 'w')
    outfile.write("\t".join([name.strip(".gz"),
                             str(round(hetero, 8))[0:10],
                             str(round(errors, 8))[0:10],
                             "\n"]))
    outfile.close()
    if not quiet:
        sys.stderr.write(".")


def main(params, quiet, mindepth):
    """ calls the main functions """

    ## assign tempmindepth to params
    params["mindepth"] = mindepth

    ##Parallel,ID,minsamp,subset,haplos,WORK,CUT,datatype):
    if not quiet:
        sys.stderr.write("\n\tstep 4: estimating error rate "+\
                         "and heterozygosity\n\t")

    ## find clust.xx directory
    if not os.path.exists(params["work"]+'clust'+params["wclust"]):
        sys.exit("\n\terror: could not find "+params["work"]+"clust"+\
                            str(params["wclust"])+"/ directory,"+ \
                            "\n\t\tif you changed the clustering threshold"+\
                            " you must transfer *.clustS"+\
                            "\n\t\tfiles to a new directory named clust.xx "+\
                            "with xx replaced by new clustering threshold")

    # warning message for low minsamp
    if params["mindepth"] < 5:
        sys.stderr.write("\n\twarning: Mindepth < 5 is not recommended\n"+\
               "\t\tfor this step. If you intend to make low\n"+\
               "\t\tcoverage base calls use a high mindepth in\n"+\
               "\t\tstep 4 to accurately infer H & E parameters, \n"+\
               "\t\tand then use a low mindepth in conjunction \n"+\
               "\t\twith the line 31 params file option to make\n"+\
               "\t\tlow coverage base calls")
        
    # if haploid data
    if params["haplos"] == 1:
        sys.stderr.write("\n\tapplying haploid-based test (infer E"+\
                         "while H is fixed to 0)\n\t")

    # if double digest use first cut site
    if "," in params["cut"]:
        cut1, cut2 = params["cut"].strip().split(",")
    else:
        cut1 = cut2 = params["cut"]

    # load up work queue
    work_queue = multiprocessing.Queue()

    # iterate over files
    clustsfiles = glob.glob(params["work"]+"clust"+params["wclust"]+\
                          "/"+params["subset"]+"*.clustS*")
    submitted = 0
    funcfiles = []
    if len(clustsfiles) > 1:
        ## sort files by size
        for clustfile in range(len(clustsfiles)):
            statinfo = os.stat(clustsfiles[clustfile])
            if statinfo.st_size > 1000:
                funcfiles.append((clustsfiles[clustfile], statinfo.st_size))
            else:
                sys.stderr.write("\n\texcluding ", clustsfiles[clustfile]+\
                                 "file is too small\n")
        funcfiles.sort(key=lambda x: x[1])
        funcfiles = [i[0] for i in funcfiles]
    else:
        funcfiles = clustsfiles
    removes = glob.glob(params["work"]+'clust'+params["wclust"]+"/cat.*")
    funcfiles = [f for f in funcfiles if f not in removes]
    for handle in funcfiles:
        work_queue.put([params, handle, cut1, cut2, quiet])
        submitted += 1

    ## remove temp files if previous run 
    for handle in funcfiles:
        end = handle.split("/")[-1].replace(".clustS.gz", "") 
        tempfile = params["work"]+"stats/."+end+".temp"
        if os.path.exists(tempfile):
            os.remove(tempfile)

    ## create a queue to pass to workers to store the results
    result_queue = multiprocessing.Queue()
    
    ##  spawn workers 
    jobs = []
    for _ in range(min(params["parallel"], submitted)):
        worker = Worker(work_queue, result_queue, optim)
        worker.start()
        jobs.append(worker)
    for job in jobs:
        job.join()

    ## write results to stats file
    if not os.path.exists(params["work"]+"stats/Pi_E_estimate.txt"):
        outstats = open(params["work"]+"stats/Pi_E_estimate.txt", 'w')
        outstats.write("taxa\tH\tE\n")
    else:
        outstats = open(params["work"]+"stats/Pi_E_estimate.txt", 'a')

    ## remove stats temp files
    for handle in funcfiles:
        end = handle.split("/")[-1].replace(".clustS.gz", "")
        tempfile = params["work"]+"stats/."+end+".temp"
        line = open(tempfile).readlines()
        outstats.write(line[0])
        os.remove(tempfile)
    outstats.close()


if __name__ == "__main__":
    main()
