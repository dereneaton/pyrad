#!/usr/bin/env python2

#######################################################
##    pyRAD for filtering and aligning RAD tags      ##
##    for phylogenetic and introgression analyses    ##
#######################################################

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
#-----------   check for numpy and scipy
try: import numpy
except ImportError:
    print "\n\tError: numpy is not installed/loaded"
    sys.exit()
try: import scipy.stats
except ImportError:
    print "\n\tError: scipy is not installed/loaded"
    sys.exit()
#-----------
import sortandcheck2
import editraw_rads, editraw_pairs
import H_err_dp
import cluster7dp
import consensdp, consens_pairs
import cluster_cons7_shuf
import tier2clust
import alignable
import loci2vcf, loci2phynex, loci2treemix
import loci2SNP, loci2mig, loci2gphocs
#------------
import Dtest, Dtest_5, Dtest_foil
import createfile
import potpour
from potpour import Worker
#------------


def main():
    parser = OptionParser(prog="pyRAD", usage="%prog [options]", version="%prog 3.0.64")
    parser.add_option('-p', action="store", type="string", dest="params",
                      help="input file for within sample filtering and clustering\n")
    parser.add_option('-s', action="store", dest="steps",
                      help="""perform step-wise parts of within analysis\n
                      1 = barcode sorting                        \
                      2 = filter/edit raw sequences              \
                      3 = within-sample clustering               \
                      4 = estimate pi and e                      \
                      5 = consensus calling                      \
                      6 = cluster consensus                      \
                      7 = align & create output files """ )
    parser.add_option('-d', action="store", type="string", dest="dtest",
                      help="""input file for D-test of introgression,
                              can iterate over multiple samples """ )
    parser.add_option('-n', action="store_true", dest="newparamsfile",
                      help="""creates a new empty input params.txt file """ )
    parser.add_option('-D', action="store_true", dest="newDtestfile",
                      help="""creates a new empty Dtest input file """ )


    (options, args) = parser.parse_args()

    if not any([options.params,options.dtest,options.newparamsfile,options.newDtestfile]):
        print "\n\tmust include option of -p, -d, -D or -n\n"
        sys.exit()

    if options.params:
        sys.stderr.write('\n\n'+' '*5+'---'*20+'\n'+\
                         ' '*6+'pyRAD : RADseq for phylogenetics & introgression analyses\n'+\
                         ' '*5+'---'*20+'\n\n')
        
        readin = [line.strip().split('##')[0].strip() for line in open(options.params).readlines()]
        if "==** " not in str(readin[0]):
            print "\n\twarning: update params input file format to latest version\n"; sys.exit()

        WORK     = str(readin[1])
        GLOB     = str(readin[2])
        Bcode    = str(readin[3])
        vsearch  = str(readin[4])
        muscle   = str(readin[5])
        CUT      = str(readin[6])  
        parallel = int(readin[7])
        mindepth = int(readin[8])
        pN       = str(readin[9])    
        wclust   = str(readin[10])   
        datatype = str(readin[11])   
        minsamp  = int(readin[12])
        maxpoly  = str(readin[13])
        outname  = str(readin[14])
        ###########################
        ## 15 is separator line
        ###########################
        subset   = str(readin[16])
        outgroup = str(readin[17])
        exclude  = str(readin[18])
        Floc     = str(readin[19])
        try: maxmismatch = int(readin[20])
        except (ValueError,IndexError): maxmismatch = 1
        try: Q = int(readin[21])
        except (ValueError,IndexError): Q = 33
        try: strict     = int(readin[22])
        except (ValueError, IndexError): strict = 0
        try: E,H      = str(readin[23]).strip().split(",")
        except ValueError: E = ""; H = ""
        try: maxN     = int(readin[24])
        except ValueError: maxN = 5
        try: maxH     = int(readin[25])
        except ValueError: maxH = 5
        try: haplos   = int(readin[26])
        except ValueError: haplos = 2
        maxSNP   = str(readin[27])
        if maxSNP == "": maxSNP = "99"
        max_inserts = str(readin[28])
        if max_inserts == "": max_inserts = "3"
        try: seed     = int(readin[29])
        except ValueError: seed = 112233
        try: overhang    = [int(i) for i in str(readin[30]).strip().split(',')]
        except (ValueError,IndexError): overhang = [0,0]
        try: outform   = str(readin[31])
        except (ValueError,IndexError): outform = ""
        try: lowcounts   = int(readin[32])
        except (ValueError, IndexError): lowcounts = mindepth
        ##mergepairs = str(readin[31])
        ##if mergepairs in [0,""]: mergepairs = 0
        try: trimkeep = int(readin[33])
        except ValueError: trimkeep = 0
        try: maxstack = int(readin[34])  
        except ValueError: maxstack = "2SD"
        try: minuniq = int(readin[35])  
        except ValueError: minuniq = 0
        try: hierarch = int(readin[36])  
        except ValueError: hierarch = 0
        try: MASK = int(readin[37])
        except ValueError: MASK = 'dust'
        if MASK == 1: MASK='dust'
        else: MASK='none'
        try: threads = int(readin[38])
        except ValueError: threads = 6
        ###############################
        ## 39 is separator line
        ###############################
        try: clustprefix = readin[40:]
        except IndexError: clustprefix = ""
        clustprefix = [i for i in clustprefix if i]
        

        """ expand ./ ~ and ../ designators in location names """
        def expander(namepath):
            if "~" in namepath:
                namepath = namepath.replace("~",os.path.expanduser("~"))
            if "../" in namepath:
                a,b = namepath.split("../")
                namepath = os.path.abspath(os.path.join(os.path.dirname( "" ), '..', b))
            elif "./" in namepath:
                a,b = namepath.split("./")
                namepath = os.path.abspath("")+"/"+b
            return namepath
            
        if WORK == "":
            WORK = os.path.abspath("")+"/"
        else:
            WORK = expander(WORK) 
        if WORK[-1] != "/":
            WORK = WORK+"/"
        stripped = 0
        if Floc:
            if Floc[0] == "@":
                stripped = 1
                Floc = expander(Floc[1:])
            else:
                Floc = expander(Floc)
        if GLOB:   GLOB = expander(GLOB)
        if Bcode:  Bcode = expander(Bcode)
        if vsearch: vsearch = expander(vsearch)
        if options.dtest: options.dtest = expander(options.dtest)

        """ find location of vsearch (or usearch) and muscle """
        def cmd_exists(cmd):
            return subprocess.call("type " + cmd, shell=True, 
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

        # " check platform: mac v linux "
        # if 'linux' in sys.platform:
        #     vsearch = "vsearch-1.0.3-linux-x86_64"
        # else:
        #     vsearch = "vsearch-1.0.3-mac-x86_64"

        # " find vsearch and muscle in user's lib/"
        # PYRADPATH = os.path.dirname(os.path.realpath(__file__))
        # vsearch = PYRADPATH+"/lib/"+vsearch
        # muscle = PYRADPATH+"/lib/muscle"

        " threads = 1 for usearch"
        if 'vsearch' not in vsearch:
            threads = 1
    
        if not cmd_exists(vsearch):
            print "\tcannot find vsearch (or usearch), edit path in param file"
            sys.exit()
        if not cmd_exists(muscle):
            print "\tcannot find muscle, edit path in input file"
            sys.exit()

        """ expand clustprefix cluster groups """
        gids = []
        groups = []
        minhits = []
        "hierarchical clustering "
        for line in clustprefix:
            gid, hits, inds = line.strip().split()
            gids.append(gid)
            minhits.append(hits)
            if "," in inds:
                thisgroup = []
                ii = inds.split(",")
                for i in ii:
                    if "*" in i:
                        expanded = glob.glob(WORK+"clust"+wclust+"/"+i+".consens*")
                        [thisgroup.append(i) for i in expanded]
                    else:
                        thisgroup.append(WORK+"clust"+wclust+"/"+i+".consens.gz")
                groups.append(thisgroup)
            else:
                if "*" in inds:
                    expanded = glob.glob(WORK+"clust"+wclust+"/"+inds+".consens*")
                    groups.append(expanded)
                else:
                    inds = inds.split(",")
                    groups.append([WORK+"clust"+wclust+"/"+i+".consens.gz" for i in inds])
            "TODO check for size=1 "
        if not gids:
            gids = ""


        " step of the analysis "
        k = tuple('1234567')
        if options.steps:
            k = tuple(str(options.steps))

        " check that the data type was entered correctly "
        datopts = ['rad','gbs','ddrad','pairgbs','pairddrad','merged','2brad']
        if datatype not in datopts:
            print "\t datatype argument (line 11) not recognized "
            sys.exit()
        # if datatype == 'merged':
        #     print "specify mergetype in params file, ex: mergeddrad or mergegbs "
        #     sys.exit()

        " parse max_inserts argument "
        w1=3
        w2=6
        a1=a2=99
        if 'pair' in datatype:
            if "," in max_inserts:
                wargs = max_inserts.strip().split(",")
                if len(wargs) == 2:
                    w1 = w2 = wargs[0]
                    a1 = a2 = wargs[1]
                elif len(wargs) == 4:
                    w1,w2,a1,a2 = wargs
                else:
                    print "\n\tmax_inserts parameter not recognized. see documentation"
                    sys.exit()
        else:
            if "," in max_inserts:
                w1,a1 = map(int,max_inserts.split(","))


        #########  Begin analysis  ###################################################
        if '1' in k:
            " expand Barcode file name if necessary "
            if "*" in Bcode:
                try: Bcode = glob.glob(Bcode)[0]
                except IndexError:
                    print "\tcould not find barcodes file ",Bcode,
                    "\n\tcomment out line 3 of params file or edit path to barcodes file"
                    sys.exit()
            if Floc:
                print "\tskipping step 1: line 18 of input file shows seqs already sorted"
            else:
                " if directory as input select all inside"
                if GLOB:
                    if GLOB[-1] == "/":
                        GLOB = GLOB+"*"
                sortandcheck2.main(Bcode,GLOB,CUT,datatype,parallel,maxmismatch,WORK)


        ### step 2 ###################
        if '2' in k:
            if Floc:
                print >>sys.stderr, "\tsorted .fastq from %s being used" % Floc
                if len(glob.glob(Floc))<1:
                    sys.stderr.write("\t... no files found in line 18 location, check required file name formatting\n")
                    sys.exit()
                FQs = Floc
                if stripped:
                    print "\tbarcode & restriction site are already stripped off of sequences"
                    CUT = ""
                    if strict:
                        print "\tApplying step 2 filter (param 19) is not recommended for data that is stripped (w/ @) \n"
            else:
                " default location "
                FQs = WORK+"fastq/"+subset+"*.fq.gz"

            " if directory as input select all inside"
            if FQs[-1] == "/":
                FQs = FQs+"*"

            " if not paired filter only read 1 "
            if 'pair' not in datatype:  # in ['rad','ddrad','gbs','merged','2brad']:
                editraw_rads.main(parallel, WORK, FQs, CUT,
                                  pN, Q, strict, trimkeep, datatype)

            else:   #elif datatype in ['pairddrad','pairgbs']:
                " check for both CUT sites in pairddrad"
                if datatype == 'pairddrad':
                    if "," not in CUT:
                        print "\n\tyou must enter two restriction sites for pair ddRAD data"
                        sys.exit()
                editraw_pairs.main(parallel, WORK, FQs, CUT, 
                                   pN, Q, strict, trimkeep, datatype)

            #elif "merge" in datatype:
            #    editraw_merges.main(parallel, WORK, FQs, CUT,
            #                       pN, Q, strict, trimkeep)



        ### step 3  ####################
        if '3' in k:
            cluster7dp.main(WORK, parallel, wclust, mindepth,
                            subset, datatype, w1, w2, minuniq,
                            MASK, muscle, vsearch, threads, remake=0)


        ### step 4  ####################
        if '4' in k:
            " if using low depth option still use a reasonable limit for parameter estimates"
            if mindepth < 5:
                tempmindepth = 5
            else:
                tempmindepth = mindepth
            H_err_dp.main(parallel, wclust, tempmindepth, subset,
                          haplos, WORK, CUT, datatype)


        ### step 5  ####################
        if '5' in k:
            if not E:
                try: Pi = open(WORK+"stats/Pi_E_estimate.txt").readlines()
                except IOError: Pi = ""
                if Pi:
                    El = []
                    Hl = []
                    for line in Pi[1:]:
                        try: _,h,e = line.strip().split("\t")
                        except IndexError:
                            None
                        Hl.append(float(h))
                        El.append(float(e))
                    if len(Hl) == 0:
                        print "\n\terror in step 4, no estimates in file stats/Pi_E_estimate.txt"
                        sys.exit()
                    H = sum(Hl)/len(Hl)
                    E = sum(El)/len(El)
                else:
                    E = 0.001
                    H = 0.01
                    print "\n\tstep 4 values not detected, using E=0.001, H=0.01"
            if 'pair' in datatype:
                " call consensus on each pair separately "
                consens_pairs.main(parallel, float(E), float(H), wclust, mindepth, subset+"*",
                                   maxN, maxH, haplos, CUT, datatype,
                                   lowcounts, strict, WORK, maxstack)
            else:
                " call consensus on single end clusters "
                consensdp.main(parallel, float(E), float(H), wclust, mindepth, subset+"*",
                               maxN, maxH, haplos, CUT, datatype,
                               lowcounts, strict, WORK, maxstack)


        ### step 6  ####################
        if '6' in k:
            if not hierarch:
                gids = ""
                if "," in subset:
                    inlist = [WORK+"clust"+wclust+"/"+i+".consens*" for i in subset.strip().split(",")]
                else:
                    inlist = glob.glob(WORK+"clust"+wclust+"/"+subset+"*.consens*")
                cluster_cons7_shuf.main(vsearch, wclust, datatype, 
                                        outgroup, seed, gids, minhits, 
                                        inlist, WORK, MASK, 0)
                print "\n\tfinished clustering"
            else:
                """ re-expand clustprefix cluster groups in case no -s """
                Hgids = []
                Hgroups = {}
                Hminhits = []
                "hierarchical clustering "
                for line in clustprefix:
                    Hgid, Hhits, Hinds = line.strip().split()
                    Hgids.append(Hgid)
                    Hminhits.append(Hhits)
                    Hgroups[Hgid] = []
                    if "," in Hinds:
                        Hinds = Hinds.split(",")
                        for Hind in Hinds:
                            if "*" in Hind:
                                expanded = glob.glob(WORK+"clust"+wclust+"/"+Hind+".consens*")
                                Hgroups[Hgid] += expanded #.append(expanded)
                            else:
                                Hgroups[Hgid].append(WORK+"clust"+wclust+"/"+Hind+".consens.gz")
                    else:
                        if "*" in Hinds:
                            expanded = glob.glob(WORK+"clust"+wclust+"/"+Hinds+".consens*")
                            Hgroups[Hgid] += expanded #.append(expanded)
                        else:
                            Hgroups[Hgid].append(WORK+"clust"+wclust+"/"+Hinds+".consens.gz")

                for i,j in zip(Hgids,Hminhits):
                    for cons in Hgroups[i]:
                        if cons not in glob.glob(WORK+"clust"+wclust+"/*.consens.gz"):
                            print "\n\tsample name",cons,"in group",i,"does not match any filenames"
                            sys.exit()

                preclusts = []
                for i in Hgroups.values():
                    preclusts += i

                for cons in glob.glob(WORK+"clust"+wclust+"/*.consens.gz"):
                    if cons not in preclusts:
                        print "\n\twarning: sample",cons,"not assigned to a cluster group"

                #if not gids:
                #    gids = ""
                    
                " make prefix directory "
                if not os.path.exists(WORK+'prefix/'):
                    os.makedirs(WORK+'prefix')


                ########### TODO ####################################
                # if os.path.exists(WORK+"prefix/cat.clust_.gz"):
                #     print "\tRemaking clusters from existing clustprefix files "+\
                #           "using minmatches: ",minmatch
                #     print "\t(To completely re-start hierarchical clustering delete the prefix/ directory)\n"
                #    
                #     for (gid,minhit,inlist) in zip(gids,minhits,groups):
                #         handle = WORK+"clust"+wclust+"/cat.haplos_"+gid
                #         #cluster_cons7_shuf.makeclust(handle, datatype, pre, pre, minm, WORK, 1)
                #     #tier2clust.makeclust(wclust, datatype, WORK)
                #######################################################

                " queue up jobs "
                work_queue = multiprocessing.Queue()
                result_queue = multiprocessing.Queue()

                " submit jobs "
                for (Hgid,Hminhit) in zip(Hgids,Hminhits):
                    inlist = Hgroups[Hgid]
                    work_queue.put([vsearch, wclust, datatype, 
                                    outgroup, seed,
                                    Hgid, Hminhit, inlist,
                                    WORK, MASK, 1 ])
                        
                " execute first tier jobs "    
                jobs = []
                for i in range(parallel):
                    worker = Worker(work_queue, result_queue, cluster_cons7_shuf.main)
                    jobs.append(worker)
                    worker.start()
                for j in jobs:
                    j.join()

                " cluster second tier "
                tier2clust.main(vsearch, wclust, datatype,
                                Hgids, seed, WORK, MASK)

                print "\n\tfinished clustering\n"

            " cleanup "
            #for ff in glob.glob(WORK+"clust"+wclust+"/cat.consens_*.gz"):
            #    os.remove(ff)
            #for ff in glob.glob(WORK+"clust"+wclust+"/cat.u*"):
            #    os.remove(ff)


        if '7' in k:
            if minsamp < 2:
                print "\n\tminimum minCov setting is <2: changing to 2"
                minsamp = 2
                
            if gids:
                inclustfile = WORK+"prefix/cat.clust_.gz"
            else:
                inclustfile = WORK+'clust'+wclust+"/cat.clust_.gz"

            if not os.path.exists(inclustfile):
                #sys.stderr.write("\n\t didn't find hierarchically clustered subset: \n\t"+inclustfile)
                #sys.stderr.write("\n\t looking for default full cluster file")
                if os.path.exists(WORK+'clust'+wclust+"/cat.clust_.gz"):
                    inclustfile = WORK+'clust'+wclust+"/cat.clust_.gz"
                    sys.stderr.write("\n\tCluster input file: using \n\t"+inclustfile+"\n\n")
                else:
                    print "\t{} not found".format(inclustfile)
                    #print "\tcat.clust_ file is selected based on line 15 subset argument "
                    #print "\n\t if you wish to exclude samples from an existing cat.clust file "+\
                    #      "\n\t in your output alignments list exclude names on line 17 of the params file.\n "
                    sys.exit()
            #if any([i in outform for i in ['t','m']]):
            #    if gids:
            #        print "\tgroups for 't' or 'm' outputs:", gids
            taxadict = OrderedDict(zip(gids,groups))
            alignable.main(outgroup, minsamp, outname,
                           inclustfile, maxpoly, parallel,
                           maxSNP, muscle, exclude, overhang,
                           outform, WORK, gids, CUT,
                           a1, a2, datatype, subset,
                           parser.version.split(" ")[1],
                           mindepth, taxadict, minhits, seed, haplos)

        if '8' in k:
            cluster7dp.main(WORK, parallel, wclust, mindepth,
                            subset, datatype, w1, w2, minuniq,
                            MASK, muscle, vsearch, threads, remake=1)

    if options.dtest:
        readin = [line.strip() for line in open(options.dtest).readlines()]

        nboots =    int(readin[0].split("##")[0].strip())
        alignfile = str(readin[1].split("##")[0].strip())
        outfile   = str(readin[2].split("##")[0].strip())
        ntax =      str(readin[3].split("##")[0].strip())
        nproc =     int(readin[4].split("##")[0].strip())
        makesort =  int(readin[5].split("##")[0].strip())
        makeboots = int(readin[6].split("##")[0].strip())
        
        tests = []
        for line in readin[8:]:
            if line:
                notes = ""
                if "##" in line:
                    tax,notes = line.strip().split("##")[0], line.strip().split("##")[-1], 
                    if tax:
                        tests.append([tax.strip().split(), notes.strip()])   #.split("\t"),notes.strip()])
                else:
                    tests.append(line.strip().split()) # "\t"))
        if ntax == '4':
            Dtest.main(tests,alignfile,outfile,nboots,nproc,makesort,makeboots)
        elif ntax == 'part':
            Dtest_5.main(tests,alignfile,outfile,nboots,nproc,makesort,makeboots)
        elif ntax == 'foil':
            Dtest_foil.main(tests,alignfile,outfile,nboots,nproc,makesort,makeboots,0)
        elif ntax == 'foilalt':
            Dtest_foil.main(tests,alignfile,outfile,nboots,nproc,makesort,makeboots,1)
        else:
            print "error in input file"

    if options.newparamsfile:
        if os.path.exists("./params.txt"):
            print "\tfile params.txt already exists"
            sys.exit()
        else:
            createfile.main(parser.version.split(" ")[1])

    if options.newDtestfile:
        outstring = """200                          ## N bootstrap replicates
test.loci                    ## loc/path to input .loci file
dstats/test1_res             ## output file path/name (no suffix)
4                            ## which test: 4,part,foil,foilalt
2                            ## N cores (execute jobs [lines below] in parallel
0                            ## output ABBA/BABA loci to files (0=no,1,2=verbose)
0                            ## output bootstrap Ds to files (0=no,1=yes)
-----------------------------------------------------------\n"""
        sys.stdout.write(outstring)

if __name__ == "__main__":
    main()







