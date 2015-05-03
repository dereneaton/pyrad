#!/usr/bin/env python2

"""
#######################################################
##    pyRAD for filtering and aligning RAD tags      ##
##    for phylogenetic and introgression analyses    ##
#######################################################
"""

import os
import sys
import glob
import multiprocessing
import subprocess
from optparse import OptionParser
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict
#-------------
import sortandcheck2
import editraw_rads
import editraw_pairs
import H_err_dp
import cluster7dp
import consensdp
import consens_pairs
import cluster_cons7_shuf
import tier2clust
import alignable
#------------
import Dtest
import Dtest_5
import Dtest_foil
import createfile
from potpour import Worker
#------------


def parseparams(paramfile):
    """ parse the params file """

    ## read in the file by line number to a list
    readin = [line.strip().split('##')[0].strip() \
              for line in open(paramfile).readlines()]

    ## check that params file is correct version format
    if "==** " not in str(readin[0]):
        sys.exit("\n\twarning: update params input file"+\
                 "format to latest version\n")

    ## get required params
    params = {
        "work" : str(readin[1]),
        "glob" : str(readin[2]),
        "bcode" : str(readin[3]),
        "vsearch" : str(readin[4]),
        "muscle" : str(readin[5]),
        "cut" : str(readin[6]),
        "parallel" : int(readin[7]),
        "mindepth" : int(readin[8]),
        "maxN" : str(readin[9]),
        "wclust" : str(readin[10]),
        "datatype" : str(readin[11]),
        "minsamp" : int(readin[12]), 
        "maxpoly" : str(readin[13]), 
        "outname" : str(readin[14]),
        ###########################
        ## 15 is separator line
        ###########################
        "subset" : str(readin[16]),
        "outgroup" : str(readin[17]),
        "exclude" : str(readin[18]),
        "floc" : str(readin[19])
        }

    ## get optional params
    if readin[20]:
        params["maxmismatch"] = int(readin[20])
    else:
        params["maxmismatch"] = 1
    if readin[21]:
        params["Q"] = int(readin[21])
    else:
        params["Q"] = 33        
    if readin[22]:
        params["strict"] = int(readin[22])
    else:
        params["strict"] = 0
    if readin[23]:
        params["E"], params["H"] = [float(i) for i in \
                            str(readin[23]).strip().split(",")]
    else:
        params["E"] = ""
        params["H"] = ""
    if readin[24]:
        params["maxN"] = int(readin[24])
    else:
        params["maxN"] = 5
    if readin[25]:
        params["maxH"] = int(readin[25])
    else:
        params["maxH"] = 5
    if readin[26]:
        params["haplos"] = int(readin[26])
    else:
        params["haplos"] = 2
    if readin[27]:
        params["maxSNP"] = str(readin[27])
    else:
        params["maxSNP"] = "99"
    if readin[28]:
        params["maxinserts"] = str(readin[28])
    else:
        params["maxinserts"] = "3"
    if readin[29]:
        params["seed"] = int(readin[29])
    else:
        params["seed"] = 112233
    if readin[30]:
        params["overhang"] = [int(i) for i in str(readin[30]).\
                              strip().split(",")]
    else:
        params["overhang"] = [0, 0]
    if readin[31]:
        params["outform"] = str(readin[31])
    else:
        params["outform"] = ""
    if readin[32]:
        params["lowcounts"] = int(readin[32])
    else:
        params["lowcounts"] = params["mindepth"]
    if readin[33]:
        params["trimkeep"] = int(readin[33])
    else:
        params["trimkeep"] = 0
    if readin[34]:
        params["maxstack"] = int(readin[34])
    else:
        params["maxstack"] = "2SD"
    if readin[35]:
        params["minuniq"] = int(readin[35])
    else:
        params["minuniq"] = 0
    if readin[36]:
        params["hierarch"] = int(readin[36])
    else:
        params["hierarch"] = 0
    if readin[37]:
        params["mask"] = int(readin[37])
        if params["mask"] == 1:
            params["mask"] = 'dust'
        else:
            params["mask"] = 'none'            
    else:
        params["mask"] = 'dust'
    if readin[38]:
        params["threads"] = int(readin[38])
    else:
        params["threads"] = 6
    try: 
        params["clustprefix"] = readin[40:]
        params["clustprefix"] = [i for i in params["clustprefix"] if i]
    except IndexError:
        params["clustprefix"] = ""        

    return params


def expander(namepath):
    """ expand ./ ~ and ../ designators in location names """        
    if "~" in namepath:
        namepath = namepath.replace("~", os.path.expanduser("~"))
    if "../" in namepath:
        a, b = namepath.split("../")
        namepath = os.path.abspath(os.path.join(os.path.dirname(""), '..', b))
    elif "./" in namepath:
        a, b = namepath.split("./")
        namepath = os.path.abspath("")+"/"+b
    return namepath


def expand_params(params):
    """ expand paths in params args """
    if not params["work"]:
        params["work"] = os.path.abspath("")+"/"
    else:
        params["work"] = expander(params["work"])

    if params["work"][-1] != "/":
        params["work"] = params["work"]+"/"

    stripped = 0
    if params["floc"]:
        if params["floc"][0] == "@":
            stripped = 1
            params["floc"] = expander(params["floc"][1:])
        else:
            params["floc"] = expander(params["floc"])

    for path in ["glob", "bcode", "vsearch"]:
        if params[path]:
            params[path] = expander(params[path])
    return params, stripped


def cmd_exists(cmd):
    """ check if dependency program is there """
    return subprocess.call("type " + cmd,
                           shell=True, 
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE) == 0


def check_params(params):
    """ check for dependency programs params before starting """

    ## check for numpy and scipy
    try: 
        import numpy
    except ImportError:
        sys.exit("\n\tError: numpy is not installed/loaded")
    try: 
        import scipy.stats
    except ImportError:
        sys.exit("\n\tError: scipy is not installed/loaded")
        sys.exit()

    ## find location of vsearch (or usearch) and muscle
    if not cmd_exists(params["vsearch"]):
        sys.exit("\tcannot find vsearch (or usearch), edit path in param file")
    if not cmd_exists(params["muscle"]):
        sys.exit("\tcannot find muscle, edit path in input file")

    ## set threads = 1 for usearch"
    if 'vsearch' not in params["vsearch"]:
        params["threads"] = 1

    ## currently allowed data types
    datopts = ['rad', 'gbs', 'ddrad',
               'pairgbs', 'pairddrad',
               'merged', '2brad']

    ## check data types
    if params["datatype"] not in datopts:
        sys.exit("\t datatype argument not recognized")

    ## parse max_inserts argument "
    ## defaults
    params["w1"] = 3
    params["w2"] = 6
    params["a1"] = params["a2"] = 99
    if 'pair' in params["datatype"]:
        if "," in params["maxinserts"]:
            wargs = params["maxinserts"].strip().split(",")
            if len(wargs) == 2:
                params["w1"] = params["w2"] = wargs[0]
                params["a1"] = params["a2"] = wargs[1]
            elif len(wargs) == 4:
                params["w1"], params["w2"], params["a1"], params["a2"] = wargs
            else:
                sys.exit("\n\tmax_inserts parameter not recognized."+\
                         " see documentation")
    else:
        if "," in params["maxinserts"]:
            params["w1"], params["a1"] = [int(i) for i in \
                                          params["maxinserts"].split(",")]

    ## parse maxSNP argument "
    if 'pair' in params["datatype"]:
        if "," in params["maxSNP"]:
            params["s1"], params["s2"] = [int(i) for i in \
                                          params["maxSNP"].split(",")]
        else:
            params["s1"] = params["s2"] = int(params["maxSNP"])
    else:
        if "," in params["maxSNP"]:
            params["s1"] = int(params["maxSNP"][0])
        else:
            params['s1'] = params['s2'] = params["maxSNP"]

    return params


def get_hierarchical_groups(params):
    """ expand clustprefix cluster groups if present """
    gids = []
    groups = []
    minhits = []

    ## hierarchical clustering "
    for line in params["clustprefix"]:
        gid, hits, inds = line.strip().split()
        gids.append(gid)
        minhits.append(hits)
        if "," in inds:
            thisgroup = []
            ii = inds.split(",")
            for i in ii:
                if "*" in i:
                    expanded = glob.glob(params["work"]+"clust"+\
                                         params["wclust"]+"/"+i+".consens*")
                    for _ in expanded:
                        thisgroup.append(_)
                else:
                    thisgroup.append(params["work"]+"clust"+\
                                     params["wclust"]+"/"+i+".consens.gz")
            groups.append(thisgroup)
        else:
            if "*" in inds:
                expanded = glob.glob(params["work"]+"clust"+\
                                     params["wclust"]+"/"+inds+".consens*")
                groups.append(expanded)
            else:
                inds = inds.split(",")
                groups.append([params["work"]+"clust"+\
                               params["wclust"]+"/"+i+".consens.gz" \
                               for i in inds])
    ## TODO check for size=1 "
    if not gids:
        gids = ""
    return gids, groups, minhits


def step1(params, quiet):
    """ demultiplex reads using barcodes file """
    ## expand Barcode file name if necessary "
    if "*" in params["bcode"]:
        try: 
            params["bcode"] = glob.glob(params["bcode"])[0]
        except IndexError:
            sys.exit("\tcould not find barcodes file "+ params["bcode"]+\
                     "\n\tcomment out line 3 of params file"+\
                     " or edit path to barcodes file")
    if params["floc"]:
        if not quiet:
            sys.stderr.write("\tskipping step 1: line 18 of input file"+\
                             " shows seqs already sorted")
    else:
        ## if directory as input select all inside
        if params["glob"]:
            if params["glob"].endswith("/"):
                params["glob"] = params["glob"]+"*"
    sortandcheck2.main(params, quiet)


def step2(params, stripped, quiet):
    """ filter based on quality scores (and check for short fragments) """

    ## check alternative fastq location
    if params["floc"]:
        print >>sys.stderr, "\tsorted .fastq from "+\
                            params["floc"]+" being used"
        if len(glob.glob(params["floc"])) < 1:
            sys.stderr.write("\t... no files found in line 18 location"+\
                             ", check required file name formatting\n")
            sys.exit()
        else:
            fastqs = params["floc"]

        ## if barocodes are already stripped off
        if stripped:
            print "\tbarcode & restriction site are already "+\
                 "stripped off of sequences"
            params['cut'] = ""
            if params["strict"]:
                print "\tApplying step 2 filter is not recommended "+\
                      "for data that is stripped (w/ @) \n"
    else:
        ## default location
        fastqs = params["work"]+"fastq/"+params["subset"]+"*.fq.gz"

    ## if directory as input select all inside
    if fastqs.endswith("/"):
        fastqs += "*"

    ## if not paired filter only read 1 "
    if 'pair' not in params["datatype"]: 
        ## extra check on cut sites for ddrad data
        if params["datatype"] == 'ddrad':
            if "," not in params["cut"]:
                sys.exit("\n\tyou must enter two restriction sites "+\
                        "for pair ddRAD data")
        editraw_rads.main(params, fastqs, quiet)

    else: 
        if params["datatype"] == 'pairddrad':
            if "," not in params["cut"]:
                sys.exit("\n\tyou must enter two restriction "+\
                         "sites for pair ddRAD data")
        editraw_pairs.main(params, fastqs, quiet)


def step3(params, quiet, remake):
    """ dereplicate, cluster, and align reads """
    cluster7dp.main(params, quiet, remake)


def step4(params, quiet):
    """ get parameters for making consensus base calls """
    ## even if using low depth option limit mindepth for statistical calls
    if params["mindepth"] < 5:
        tempmindepth = 5
    else:
        tempmindepth = params["mindepth"]
    H_err_dp.main(params, quiet, tempmindepth) 


def step5(params, quiet):
    """ make consensus base calls, keep alleles for haploid and diploids, 
    keep base counts for all sites """
    ## get H and E from the stats output of step 4, or else use defaults
    if not params["E"]:
        try: 
            inferp = open(params["work"]+"stats/Pi_E_estimate.txt").readlines()
        except IOError: 
            inferp = ""
        if inferp:
            all_e = []
            all_h = []
            for line in inferp[1:]:
                try: 
                    _, het, err = line.strip().split("\t")
                    all_h.append(float(het))
                    all_e.append(float(err))
                except IndexError:
                    pass

            if len(all_h) == 0:
                print "\n\terror in step 4, "+\
                      "no estimates in file stats/Pi_E_estimate.txt"
                sys.exit()
            params["H"] = sum(all_h)/len(all_h)
            params["E"] = sum(all_e)/len(all_e)
        else:
            ## use default values
            params["E"] = 0.001
            params["H"] = 0.01
            print "\n\tstep 4 values not detected, using E=0.001, H=0.01"

    if 'pair' in params["datatype"]:
        ## call consensus on each pair separately
        consens_pairs.main(params, quiet)
    else:
        ## call consensus on single end clusters "
        consensdp.main(params, quiet)


def step6(params, gids, groups, minhits, quiet):
    """ perform across sample clustering. Do step-wise (heirarchical)
        approach if hierarch groups are provided """

    if not params["hierarch"]:
        if "," in params["subset"]:
            inlist = [params["work"]+"clust"+params["wclust"]+\
                      "/"+i+".consens*" for i in params["subset"]\
                      .strip().split(",")]
        else:
            inlist = glob.glob(params["work"]+"clust"+params["wclust"]+\
                               "/"+params["subset"]+"*.consens*")

        ## cluster consens files in inlist
        cluster_cons7_shuf.main(params, inlist, "", "", "", quiet)
        ## vsearch, wclust, datatype, 
        ## outgroup, seed, gids, minhits, 
        ## inlist, WORK, MASK, 0)
        if not quiet:
            sys.stderr.write("\n\tfinished clustering\n")

    else:
        print gids
        print groups
        print minhits

        if not os.path.exists(params["work"]+"prefix/"):
            os.makedirs(params["work"]+"prefix/")

        ## queue up jobs
        work_queue = multiprocessing.Queue()
        result_queue = multiprocessing.Queue()

        ## submit jobs
        for (gid, minhit) in zip(gids, minhits):
            inlist = groups[gid]
            work_queue.put([params, inlist, gid, minhit, "", quiet])
    #     inlist = Hgroups[Hgid]
    #     work_queue.put([vsearch, wclust, datatype, 
    #         outgroup, seed,
    #         Hgid, Hminhit, inlist,
    #         WORK, MASK, 1 ])
                        
        ## execute first tier jobs "    
        jobs = []
        for i in range(params["parallel"]):
            worker = Worker(work_queue,
                            result_queue,
                            cluster_cons7_shuf.main)
            jobs.append(worker)
            worker.start()
            for j in jobs:
                j.join()

        ## cluster second tier
        tier2clust.main(params)
        #  Hgids, seed, WORK, MASK)
        print "\n\tfinished clustering\n"



def step7(params, gids, groups, minhits,
          version, quiet):
    """ align final cluster file, unless it already exists, and create 
    additional output formats """

    ## do not ouput singletons
    if params["minsamp"] < 2:
        print "\n\tminimum minCov setting is <2: changing to 2"
        params["minsamp"] = 2

    ## if group ids are provided check for clust file in prefix/
    ## directory first, then look in normal clust dir/
    if params["hierarch"]:
        inclustfile = params["work"]+"prefix/cat.clust_.gz"
    else:
        inclustfile = params["work"]+'clust'+params["wclust"]+"/cat.clust_.gz"

    ## if clust file was found, print to screen 
    if os.path.exists(inclustfile):
        if not quiet:
            sys.stderr.write("\n\tUsing across-sample cluster input file: "+\
                             "\n\t"+inclustfile+"\n\n")
    else:
        sys.exit("\tclustfile not found")

    ## group assignments for certain outputs
    taxadict = OrderedDict(zip(gids, groups))

    ## final function
    alignable.main(params, inclustfile, taxadict, 
                   minhits, version, quiet)


def step8hidden(params):
    """ recreate step 3 cluster files from unfinished vsearch run"""
    cluster7dp.main(params, remake=1) 



def parseD(dtestfile):
    """ parses dtest input file and returns parameters """

    ## read in the input file
    readin = [line.strip() for line in open(dtestfile).readlines()]
    
    ## parse params
    dparams = {}
    dparams[nboots] = int(readin[0].split("##")[0].strip())
    alignfile = str(readin[1].split("##")[0].strip())
    outfile = str(readin[2].split("##")[0].strip())
    ntax = str(readin[3].split("##")[0].strip())
    nproc = int(readin[4].split("##")[0].strip())
    makesort =  int(readin[5].split("##")[0].strip())
    makeboots = int(readin[6].split("##")[0].strip())
    
    ## read in the tests    
    tests = []
    for line in readin[8:]:
        if line:
            notes = ""
            if "##" in line:
                tax, notes = [line.strip().split("##")[0],
                             line.strip().split("##")[-1]]
                if tax:
                    tests.append([tax.strip().split(),
                                  notes.strip()])
            else:
                tests.append(line.strip().split()) 

    ## parse complex taxon assignments
    ## 
    ## 

    ## check that params are correct
    ## tests should only have four taxa for d4 test, etc.
    ## 
    ##

    ## return paramslist


############################################################


if __name__ == "__main__":

    ## parse the command line arguments
    PARSER = OptionParser(prog="pyRAD", 
                          usage="%prog [options]",
                          version="%prog 3.1.0a0")

    PARSER.add_option('-p', action="store", 
                            type="string",
                            dest="params",
                 help="input parameter file and start analysis\n")

    PARSER.add_option('-s', action="store", 
                            dest="steps",
                 help="""with -p performs step-wise parts of analysis\n
                      1 = barcode sorting                        \
                      2 = filter/edit raw sequences              \
                      3 = within-sample clustering               \
                      4 = estimate pi and e                      \
                      5 = consensus calling                      \
                      6 = cluster consensus                      \
                      7 = align & create output files """)

    PARSER.add_option('-d', action="store", 
                            type="string", 
                            dest="dtest",
                      help="""input file for D-test of introgression,
                              can iterate over multiple samples """)

    PARSER.add_option('-n', action="store_true",
                            dest="newparamsfile",
                      help="""creates a new empty input params.txt file """)

    PARSER.add_option('-D', action="store_true",
                            dest="newDtestfile",
                      help="""creates a new empty Dtest input file """)

    PARSER.add_option('-q', action="store_true",
                            dest="quiet",
                      help="do not print to screen")

    (OPTIONS, ARGS) = PARSER.parse_args()

    ## check entered opts
    if not any([OPTIONS.params, OPTIONS.dtest, 
                OPTIONS.newparamsfile, OPTIONS.newDtestfile]):
        print "\n\tmust include option of -p, -d, -D or -n\n"
        sys.exit()

    ## create new params file if asked
    if OPTIONS.newparamsfile:
        if os.path.exists("./params.txt"):
            print "\tfile params.txt already exists"
            sys.exit()
        else:
            createfile.main(PARSER.version.split(" ")[1])

    ## create D-stat input file if asked
    if OPTIONS.newDtestfile:
        OUTSTRING = """200                          ## N bootstrap replicates
        test.loci                    ## loc/path to input .loci file
        dstats/test1_res             ## output file path/name (no suffix)
        4                            ## which test: (4,part,foil,foilalt)
        2                            ## N cores (execute jobs [lines below] in parallel)
        0                            ## output ABBA/BABA loci to files (0=no,1,2=verbose)
        0                            ## output bootstrap Ds to files (0=no,1=yes)
        -----------------------------------------------------------\n"""
        sys.stdout.write(OUTSTRING)

    ## do D-stat measurements if asked
    if OPTIONS.dtest:
        ## parse the dstat input file
        DPARAMS = parseD(OPTIONS.dtest)
        ## call the appropriate test 
        if DPARAMS["test"] == '4':
            Dtest()
        elif DPARAMS["test"] == 'part':
            Dtest()
        elif 'foil' in DPARAMS["test"]:
            Dtest()
        else:
            sys.exit("error in input file")

        # if ntax == '4':
        #     Dtest.main(tests,alignfile,outfile,nboots,nproc,makesort,makeboots)
        # elif ntax == 'part':
        #     Dtest_5.main(tests,alignfile,outfile,nboots,nproc,makesort,makeboots)
        # elif ntax == 'foil':
        #     Dtest_foil.main(tests,alignfile,outfile,nboots,nproc,makesort,makeboots,0)
        # elif ntax == 'foilalt':
        #     Dtest_foil.main(tests,alignfile,outfile,nboots,nproc,makesort,makeboots,1)
        # else:
        #     print "error in input file"

        ## run the appropriate step(s) 

    ## print splash page if not quiet
    if OPTIONS.params:
        if not OPTIONS.quiet:
            sys.stderr.write('\n\n'+' '*5+'---'*20+'\n'+' '*6+\
               'pyRAD : RADseq for phylogenetics & introgression analyses\n'+\
               ' '*5+'---'*20+'\n\n')

        ## get params from the params file
        PARAMS = parseparams(OPTIONS.params)

        ## check params 
        PARAMS, STRIPPED = expand_params(PARAMS)
        PARAMS = check_params(PARAMS)
        GIDS, GROUPS, MINHITS = get_hierarchical_groups(PARAMS)

        ## steps of the analysis
        K = tuple('1234567')
        if OPTIONS.steps:
            K = tuple(str(OPTIONS.steps))

        if '1' in K:
            step1(PARAMS, OPTIONS.quiet)
        if '2' in K:
            step2(PARAMS, STRIPPED, OPTIONS.quiet)
        if '3' in K:
            step3(PARAMS, OPTIONS.quiet, remake=0)
        if '4' in K:
            step4(PARAMS, OPTIONS.quiet)
        if '5' in K:
            step5(PARAMS, OPTIONS.quiet)
        if '6' in K:
            step6(PARAMS, GIDS, GROUPS, MINHITS, OPTIONS.quiet)
        if '7' in K:
            step7(PARAMS, GIDS, GROUPS, MINHITS, 
                  PARSER.version.split(" ")[1], OPTIONS.quiet)

##############################################





