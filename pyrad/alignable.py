#!/usr/bin/env python2

""" final alignment and output formatting """

import os
import sys
import glob
import subprocess
import multiprocessing
import gzip
import time
import cPickle as pickle
from numpy import array
from itertools import izip
from copy import copy
from potpour import Worker
from cluster_cons7_shuf import breakalleles

try:
    from collections import Counter
except ImportError:
    from ordereddict import Counter

import loci2phynex
import loci2vcf
import loci2treemix
import loci2SNP
import loci2mig
import loci2gphocs
import loci2cat


def unstruct(amb):
    " returns bases from ambiguity code"
    amb = amb.upper()
    trans = {"R":["G", "A"],
             "K":["G", "T"],
             "S":["G", "C"],
             "Y":["T", "C"],
             "W":["T", "A"],
             "M":["C", "A"],
             "A":["A", "A"],
             "T":["T", "T"],
             "G":["G", "G"],
             "C":["C", "C"],
             "N":["N", "N"],
             "-":["-", "-"]}
    return trans.get(amb)


def alignfast(params, pronum, names, seqs):
    """ inputs data to muscle and returns aligned data
        as a string that needs to be parsed """
    inputstring = "\n".join('>'+i+'\n'+j[0] for i, j in zip(names, seqs))
    
    ##  if inputstring is very large it needs to be written to 
    ##  file, otherwise the process can just be piped
    if len(inputstring) > 100000:
        fstring = params["work"]+".tempalign_"+pronum
        with open(fstring, 'w') as instring:
            print >>instring, inputstring
        cmd = params["muscle"]+" -quiet -in "+fstring #+" -out "+ostring
    else:
        cmd = "/bin/echo '"+inputstring+"' | "+\
                            params["muscle"]+" -quiet -in -"

    ## RUN muscle... is this faster than .call? need to test speed.
    fout = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    returnstring = fout.stdout.read()
    #NB: check if fout.communicate()[0] might be faster
    return returnstring


def polyfilter(seqs, maxpoly):
    """ filter for maximum number of polymorphic sites
        in an aligned locus with user supplied max """

    ## how many columns (sites) are polymorphic """
    arrayed = array([tuple(seq) for seq in seqs])
    counted = [Counter(site) for site in arrayed.T]

    ## parse maxpoly if it is a proportion
    if 'p' in str(maxpoly):
        maxpoly = len(seqs)*float(maxpoly.replace("p", ""))

    ## are there more H sites than allowed?
    notpass = []
    for count in counted:
        if any([count[ambig] > maxpoly for ambig in list("RKSYWM")]):
            notpass.append(count)
    if notpass:
        return 1
    else:
        return 0


def sortalign(stringnames):
    """ gets (names, seqs) from muscle output """
    ## splits string into separate names\nseqs
    splitstring = stringnames.split("\n>")  
    ## removes a bunch of newline characters within seqs and > chars
    elements = [i.split("\n")[0].replace(">", "")+"\n"+\
                "".join(i.split('\n')[1:]) for i in splitstring]
    aligned = [i.split("\n") for i in elements]
    names = [">"+i[0] for i in aligned]
    seqs = [i[1] for i in aligned]
    return names, seqs


def screensamples(ingroup, outgroup, exclude):
    """ print sampling strategy to stderr """    
    ## print includes and excludes to screen
    toprint = [i for i in list(ingroup) if i not in exclude]
    toprint.sort()
    i = 0
    sys.stderr.write('\tingroup '+", ".join(toprint[i:i+4])+"\n")
    i += 4    
    while i < len(toprint):
        sys.stderr.write("\t        "+", ".join(toprint[i:i+4])+"\n")
        i += 4

    toprint = [i for i in outgroup if i not in exclude]
    toprint.sort()
    i = 0
    sys.stderr.write('\taddon   '+", ".join(toprint[i:i+4])+"\n")
    i += 4    
    while i < len(toprint):
        sys.stderr.write("\t        "+", ".join(toprint[i:i+4])+"\n")
        i += 4

    toprint = exclude
    toprint.sort()
    i = 0
    sys.stderr.write('\texclude '+", ".join(toprint[i:i+4])+"\n")
    i += 4    
    while i < len(toprint):
        sys.stderr.write("\t        "+", ".join(toprint[i:i+4])+"\n")
        i += 4
    sys.stderr.write("\t")


def trimmer(params, names, tseqs):
    """ overhang trim or keep """
    fm1 = sm1 = None   ## FM, SM
    fm2 = sm2 = None

    ## only trim if more than three samples with
    ## setting 1, unless the minsamp allows fewer 
    ## than 4, in which case make the min for trimming
    ## equal to minsamp
    mintcov = 4
    if params["minsamp"] < 4:
        mintcov = params["minsamp"]
    
    ## exclude singletons
    if params["minsamp"] > 1:
        if 'pair' in params["datatype"]:
            firsts = [i.split("n")[0] for i in tseqs]
            seconds = [i.split("n")[-1] for i in tseqs]

            ## treat each read separately, or do the total read?
            if len(params["overhang"]) == 4:
                t1left, t1right, t2left, t2right = params["overhang"]
            elif len(params["overhang"]) == 2:
                t1left, t2left = params["overhang"]
                t1right = t1left
                t2right = t2left

            ## trim 1st read
            leftlimit = [edger(i, 'min') for i in firsts]
            rightlimit = [edger(i, 'max') for i in firsts]

            ## trim1 means at least four samples have data at that site
            ## trim2 means that all samples with data have data at that site
            if t1left == 1:
                ## trim1 1st left overhang
                fm1 = min([i for i in xrange(len(firsts[0])) if \
                          [j <= i for j in leftlimit].count(True) \
                             >= mintcov])
            elif t1left == 2:
                try: 
                    fm1 = min([i for i in xrange(len(firsts[0])) if \
                              [j <= i for j in leftlimit].count(True) \
                                 == len(names)])
                except ValueError: 
                    fm1 = 1
                        
            if t1right == 1:
                ## trim 1st right overhang
                sm1 = max([i for i in xrange(len(firsts[0])) if \
                          [j >= i for j in rightlimit].count(True) \
                             >= mintcov])
            elif t1right == 2:
                try: 
                    sm1 = max([i for i in xrange(len(firsts[0])) if \
                              [j >= i for j in rightlimit].count(True) \
                                 == len(names)])
                except ValueError: 
                    sm1 = 1

            ## trim 2nd read
            leftlimit = [edger(i, 'min') for i in seconds]
            rightlimit = [edger(i, 'max') for i in seconds]

            if t2left == 1:
                ## trim1 2nd left overhang
                try: 
                    fm2 = min([i for i in xrange(len(seconds[0])) if \
                              [j <= i for j in leftlimit].count(True) \
                                 >= mintcov])
                except ValueError:
                    ## no sites where 4 samples have data
                    fm2 = 1
                    ## empty2 = 1
                        
            elif t2left == 2:
                ## trim1 2nd left overhang
                try: 
                    fm2 = min([i for i in xrange(len(seconds[0])) if \
                              [j <= i for j in leftlimit].count(True) \
                                 == len(names)])
                except ValueError:
                    ## no sites where all samples have data
                    fm2 = 1
                    #empty2 = 1

            if t2right == 1:
                ## trim 2nd right overhang
                sm2 = max([i for i in xrange(len(seconds[0])) if \
                          [j >= i for j in rightlimit].count(True) \
                             >= mintcov])
            elif t2right == 2:
                sm2 = max([i for i in xrange(len(seconds[0])) if \
                          [j >= i for j in rightlimit].count(True) \
                             == len(names)])+1

            ## put pair back together
            ## TODO check this...
            tseqs = [i+"nnnn"+j for i, j in zip(firsts, seconds)]
                        
        else:
            leftlimit = [edger(i, 'min') for i in tseqs]
            rightlimit = [edger(i, 'max') for i in tseqs]

            ## trim left overhang
            if params["overhang"][0] == 1:
                fm1 = min([i for i in xrange(len(tseqs[0])) if \
                          [j <= i for j in leftlimit].count(True) \
                             >= mintcov])
            elif params["overhang"][0] == 2:
                fm1 = max(leftlimit)
            else:
                fm1 = min(leftlimit)
                        
            ## trim right overhang
            if params["overhang"][1] == 1:
                sm1 = max([i for i in xrange(len(tseqs[0])) if \
                          [j >= i for j in rightlimit].count(True) \
                             >= mintcov])+1
            elif params["overhang"][1] == 2:
                sm1 = min(rightlimit)+1
            else:
                sm1 = max(rightlimit)+1
    return [fm1, fm2, sm1, sm2]


def getvariable(tseqs):
    """ returns a snpstring that records synapomophies
        with * and autapomorphies with - and is placed
        under each locus in the .loci output file """
    ## create list of bases at each site
    arrayed = array([tuple(seq) for seq in tseqs])
    basenumber = 0
    snpsite = [" "]*len(tseqs[0])

    ## put in split for pairs
    for i in range(len(tseqs[0])):
        if tseqs[0][i] == "n":
            snpsite[i] = "n"

    ## record a string for variable sites in snpsite
    for site in arrayed.T:
        ## if site is variable
        reals = [i for i in site if i not in list("N-")]
        if len(set(reals)) > 1:
            ## convert ambiguity bases to reals
            for i in xrange(len(reals)):
                if reals[i] in list("RWMSYK"):
                    for j in unstruct(reals[i]):
                        reals.append(j)
            reals = [i for i in site if i not in list("N-")]                    
            ## if not an autapomorphy
            if sorted([reals.count(i) for i in set(reals)], 
                       reverse=True)[1] > 1:
                ## mark PIS for outfile
                snpsite[basenumber] = "*"
            else:
                snpsite[basenumber] = "-"
        basenumber += 1
    return snpsite


def alignfunc(params, infile, ingroup, exclude, longname, quiet):
    """ align each cluster in infile using muscle """

    ## split cutters
    if "," in params["cut"]:
        cut1, cut2 = params["cut"].split(",")
    else:
        cut1 = cut2 = params["cut"]

    ## assign number to files for this process
    pronum = str("".join(infile.split("_")[-1]))

    ## create temp out files for aligned clusters
    aout = open(params["work"]+".align_"+pronum, 'w')
    nout = open(params["work"]+".not_"+pronum, 'w')

    ## read in clust file 2 lines at a time
    clusts = open(infile)   ## f
    duo = izip(*[iter(clusts)]*2)

    while 1:
        try: 
            itera = duo.next()  ## d
        except StopIteration:
            break
        ## local lists
        names = []
        cnames = []
        onames = []
        seqs = []
        nameiter = 0
        while "//\n" not in itera:
            ## record names and seqs, remove # at end"
            ## record the name into locus name. "
            nam = "_".join(itera[0].split(">")[1].split("_")[:-2])
            if nam not in exclude:
                cnames.append("_".join(itera[0].split(">")\
                                       [1].split("_")[:-2]))
                names.append("_".join(itera[0].split(">")[1].split("_")\
                                       [:-2])+"_"+str(nameiter))
                onames.append(itera[0].strip().split("_")[-1])
                seqs.append(itera[1].strip())
            itera = duo.next()
            nameiter += 1
        ## get old locus id
        olocus = itera[0].strip()
        ## apply duplicate filter "
        ## no grouping un-clustered copies from same taxon
        if not len(cnames) != len(set(cnames)):      

            ## apply minsamp filter checking if
            ## too few ingroup samples in locus
            ingroupnames = [i for i in cnames if i in ingroup]
            if len(ingroupnames) >= (params["minsamp"]):
            
                ## align read1 separate from read2
                if 'pair' in params["datatype"]:
                    ## compatibility from pyrad 2 -> 3
                    seqs = [i.replace("X", 'n') for i in seqs]
                    firsts = [[i.split("nnnn")[0]] for i in seqs]
                    seconds = [[i.split("nnnn")[-1]] for i in seqs]

                    ## align first reads
                    stringnames = alignfast(params, pronum, names, firsts)
                    names, seqs1 = sortalign(stringnames)
                    read1dic = {}
                    for i in range(len(names)):
                        read1dic[names[i]] = seqs1[i]

                    ## reorder keys by name
                    keys = read1dic.keys()
                    keys.sort(key=lambda x: int(x.split("_")[-1]), reverse=True)

                    ## align second reads
                    stringnames = alignfast(params, pronum, names, seconds)
                    names, seqs2 = sortalign(stringnames)
                    read2dic = {}
                    for i in range(len(names)):
                        read2dic[names[i]] = seqs2[i]
                    names = keys 
                    seqs = [read1dic[key]+"nnnn"+read2dic[key] for key in keys]

                else:
                    ## align reads
                    seqs = [[i] for i in seqs]
                    stringnames = alignfast(params, pronum, names, seqs)
                    names, seqs = sortalign(stringnames)  ## nn, sss
                
                ## now strip off cut sites
                if params["datatype"] == "merged":
                    tseqs = [i[len(cut1):-len(cut2)] for i in seqs]
                ## TODO: double check that pair in names is right here
                elif ("c1" in onames) or ("pair" in onames):
                    tseqs = [i[len(cut1):-len(cut2)] for i in seqs]
                else:
                    tseqs = [i[len(cut1):] for i in seqs]

                ## trim off numbers that were added to names
                names = ["_".join(i.split("_")[:-1]) for i in names]

                ## apply paralog filter
                if not polyfilter(tseqs, params["maxpoly"]):

                    ## tupled
                    zz = zip(names, tseqs)

                    ## record variable sites
                    snpsite = getvariable(tseqs)

                    ## get trimmed edges
                    fm1, fm2, sm1, sm2 = trimmer(params, names, tseqs)

                    ## alphabetize names
                    zz.sort()

                    ## filter for duplicates or paralogs, then SNPs and Indels
                    ffilter = ""
                    ffilter = sandi_filter(params, zz, snpsite, fm1, sm1, sm2)

                    if not ffilter:
                        ## write aligned loci to temp files for later 
                        ## concatenation into the .loci file
                        writetokeep(params, snpsite, zz, longname,
                                    fm1, fm2, sm1, sm2, aout, olocus)
                    else:
                        #filterlist.append(ffilter)
                        writetoexclude(params, snpsite, zz,
                                       longname, fm1, fm2, 
                                       sm1, sm2, ffilter, nout)
                else:
                    #filterlist.append("P")
                    writetoexclude(params, snpsite, zz,
                                   longname, fm1, fm2, 
                                   sm1, sm2, "P", nout)                   
            else:
                #filterlist.append("L") ## g4 += 1
                zz = zip(cnames, seqs)                
                writetoexclude(params, "---unaligned---", zz,
                               longname, None, None, None, None,
                               "L", nout)
        else:
            #filterlist.append("D") 
            zz = zip(cnames, seqs)
            writetoexclude(params, "---unaligned---", zz,
                           longname, None, None, None, None,
                           "D", nout)
    nout.close()
    aout.close()
    if not quiet:
        sys.stderr.write('.')


def writetokeep(params, snpsite, zz, longname,
                fm1, fm2, sm1, sm2, aout, olocus):
    """ writes aligned locus to output file to be further
        formatted later """

    if 'pair' in params["datatype"]:
        snp1, snp2 = "".join(snpsite).split("nnnn")
        for name, seq in zz:
            first, second = seq.split("nnnn")
            space = ((longname+5)-len(name))
            print >>aout, name+" "*space+first[fm1:sm1].upper()+\
                              'nnnn'+second[fm2:sm2].upper()
        print >>aout, '//'+' '*(longname+3)+snp1[fm1:sm1]+\
                      "    "+snp2[fm2:sm2]+"|"+olocus+"|"
    else:
        for name, seq in zz:
            space = ((longname+5)-len(name))
            print >>aout, name+" "*space + seq[fm1:sm1].upper()
        print >>aout, '//'+' '*(longname+3)+\
                      "".join(snpsite[fm1:sm1])+"|"+olocus+"|"



def writetoexclude(params, snpsite, zz, longname,
                   fm1, fm2, sm1, sm2, thisfilter, nout):
    """ write to exclude file with letter to designate filter """

    if 'pair' in params["datatype"]:
        snp1, snp2 = "".join(snpsite).split("nnnn")
        for name, seq in zz:
            first, second = seq.split("nnnn")
            space = ((longname+5)-len(name))
            print >>nout, name+" "*space+first[fm1:sm1].upper()+\
                          'nnnn'+second[fm2:sm2].upper()
        print >>nout, '//'+thisfilter+' '*(longname+3-len(thisfilter))+\
                           snp1[fm1:sm1]+"    "+snp2[fm2:sm2]+"|"#+notes

    else:
        for name, seq in zz:
            space = ((longname+5)-len(name))
            print >>nout, name+" "*space+seq[fm1:sm1].upper()
        print >>nout, '//'+thisfilter+' '*(longname+3-len(thisfilter))+\
                      "".join(snpsite[fm1:sm1])+"|"#+notes
                


def sandi_filter(params, zz, snpsite, fm1, sm1, sm2):
    """ filter for the max allowed number of snps in a locus"""

    ffilter = ""
    ## apply SNP filter
    if 'pair' in params["datatype"]:
        ## if paired, apply separate filter to first and second reads
        snp1, snp2 = "".join(snpsite).split("nnnn")
        snp1 = snp1.replace("*", "-")
        if snp1.count("-") > int(params["s1"]):
            ffilter = "S"
        else:
            snp2 = snp2.replace("*", "-")
            if snp2.count("-") > int(params["s2"]):
                ffilter = "S"
    else:
        if "".join(snpsite[fm1:sm1]).replace("*", "-")\
                       .count('-') > int(params["s1"]):
            ffilter = "S"

    ## Indel filter
    if not ffilter:
        if "pair" in params["datatype"]:
            spacer = zz[0][0].index("n")
            if any([seq[fm1:spacer].count("-") > int(params["a1"]) \
                                            for _, seq in zz]):
                ffilter = "I"
            elif any([seq[spacer:sm2].count("-") > int(params["a2"]) \
                                            for _, seq in zz]):
                ffilter = "I"
        else:
            if any([seq[fm1:sm1].count("-") > int(params["a1"]) \
                                        for _, seq in zz]):
                ffilter = "I"
    return ffilter


def edger(seq, minmax):
    """ finds leftmost or rightmost base in an alignment """
    if minmax == 'max':
        try:
            fmost = max([i for i, j in enumerate(seq) if j != "-"])
        except ValueError:
            fmost = 1
    elif minmax == 'min':
        try: 
            fmost = min([i for i, j in enumerate(seq) if j != "-"])
        except ValueError:
            ## only Ns and -s"
            fmost = len(seq)
    return fmost


def splitandalign(params, infile, ingroup, exclude, longname, quiet):
    """ split cluster file into smaller files depending on the number
    of processors and align each file separately using alignfunc function."""

    ## break input file into chunks for n threads
    for i in glob.glob(params["work"]+".align*"):
        os.remove(i)
    for i in glob.glob(params["work"]+".chunk*"):
        os.remove(i)

    ## read infile, split into chunks for aligning, nchuncks
    ## depends on number of available processors
    data = gzip.open(infile, 'rb').read().strip().split("//\n")
    minpar = max(3, params["parallel"])  ## pp
    chunks = [0+(len(data)/minpar)*i for i in range(minpar)]
    for i in range(len(chunks)-1):
        dat = open(params["work"]+".chunk_"+str(i), 'w')
        #print "//\n\n".join(data[chunks[i]])+"//\n\n"
        dat.write("//\n".join(data[chunks[i]:chunks[i+1]])+"//\n")
        dat.close()
    dat = open(params["work"]+".chunk_"+str(i+1), 'w')
    dat.write("//\n".join(data[chunks[i+1]:])+"\n")
    dat.close()

    ## set up parallel
    work_queue = multiprocessing.Queue()
    result_queue = multiprocessing.Queue()
    for handle in glob.glob(params["work"]+".chunk*"):
        work_queue.put([params, handle, ingroup, 
                        exclude, longname, quiet])

    ## spawn workers
    jobs = []
    for i in range(minpar):
        worker = Worker(work_queue, result_queue, alignfunc)
        jobs.append(worker)
        worker.start()
    for j in jobs:
        j.join()


def makelocifile(params):
    """ concatenate split aligned files into .loci file """
    ## keep track of how many loci passed filters
    locicounter = 1

    ## dictionary to map old locus names to new ones
    new2olddict = {}

    ## location of aligned loci
    aligns = glob.glob(params["work"]+".align*")

    ## filename to write to
    locifile = open(params["work"]+"outfiles/"+\
                    params["outname"]+".loci", "w")
    
    ## iterate over aligned loci files and write
    ## aligned loci to the .loci output file
    ## and delete temp files
    for chunkfile in aligns:
        chunkdata = open(chunkfile, "r")
        for lines in chunkdata:
            if lines.startswith("//"):
                ## compatibility between pyrad versions 3.0--3.1
                if "|\n" in lines:
                    lines = lines.replace("|\n", "|"+str(locicounter)+"\n", 1)
                else:
                    ## map filtered locus numbers to unfiltered locus numbers
                    snpstring, olocus = lines.strip().split("|")
                    #print olocus, locicounter
                    lines = snpstring+"|"+str(locicounter)+"\n"
                    #dataobj = Locobj(olocus, indeldict)
                    new2olddict[locicounter] = olocus
                locicounter += 1
            locifile.write(lines)
        chunkdata.close()
        os.remove(chunkfile)
    locifile.close()

    ## concatenate all of the bad loci into excluded_loci file    
    ## and delete temp files
    unaligns = glob.glob(params["work"]+".not*")
    excluded_loci_file = open(params["work"]+"outfiles/"+\
                              params["outname"]+".excluded_loci", "w")
    filterlist = []
    for excludechunk in unaligns:
        excludedata = open(excludechunk, "r")
        for lines in excludedata:
            if "//" in lines:
                filterlist.append(lines[2])
            excluded_loci_file.write(lines)
        excludedata.close()
        os.remove(excludechunk)
    excluded_loci_file.close()

    ## remove chunk files
    for chunkfile in glob.glob(params["work"]+".chunk*"):
        os.remove(chunkfile)

    return locicounter-1, filterlist, new2olddict

    

def dostats(params, ingroup, outgroup, longname,
            locus, filterlist, version, quiet):
    """ get final stats """

    ## print to screen the location of final stats files
    if not quiet:
        sys.stderr.write("\n\tfinal stats written to:\n\t "+\
              params["work"]+"stats/"+params["outname"]+".stats")
        sys.stderr.write("\n\toutput files being written to:\n\t "+\
               params["work"]+"outfiles/ directory\n")

    ## open stats file for writing
    statsout = open(params["work"]+"stats/"+\
                     params["outname"]+".stats", 'w')   
    finalfile = open(params["work"]+"outfiles/"+\
                     params["outname"]+".loci").read() 


    #print Counter(filterlist), "counter filterlist in dostats"
    ## print header for stats output
    stat = params["outname"]
    print >>statsout, stat+" "*(20-len(stat))+"## Named outputs"
    stat = time.strftime("%Y/%m/%d %Hh:%Mm")
    print >>statsout, stat+" "*(20-len(stat))+"## Date"
    stat = "pyRAD.v."+str(version)
    print >>statsout, stat+" "*(20-len(stat))+"## Source\n"

    ## print how many loci passed filtering
    stat = str(locus+len(filterlist))
    print >>statsout, stat+" "*(12-len(stat))+\
                      "## total clusters excluding singletons"

    stat = str(locus+len([i for i in filterlist if i not in list("D")]))
    print >>statsout, stat+" "*(12-len(stat))+\
                      "## loci [filters: noDups]"

    stat = str(locus+len([i for i in filterlist if i not in list("DL")]))
    print >>statsout, stat+" "*(12-len(stat))+\
                      "## loci [filters: noDups, >minCov]"

    stat = str(locus+len([i for i in filterlist if i not in list("DLP")]))
    print >>statsout, stat+" "*(12-len(stat))+\
                      "## loci [filters: noDups, >minCov, <maxSH]"

    stat = str(locus+len([i for i in filterlist if i not in list("DLPS")]))
    print >>statsout, stat+" "*(12-len(stat))+\
                      "## loci [filters: noDups, >minCov, <maxSH, <maxSNP]"

    stat = str(locus+len([i for i in filterlist if i not in list("DLPSI")]))
    print >>statsout, stat+" "*(12-len(stat))+\
                      "## loci [filters: noDups, >minCov, <maxSH, "+\
                                              "<maxSNP, <maxIndel]"

    ## print columns for how many loci were found in each sample "
    print >>statsout, "\n## number of loci recovered in final "+\
                      "data set for each taxon."
    names = list(ingroup)+outgroup
    names.sort()
    print >>statsout, '\t'.join(['taxon'+" "*(longname-5), 'nloci'])
    for name in names:
        print >>statsout, name+" "*(longname-len(name))+\
                          "\t"+str(finalfile.count(">"+name+" "))

    ## print distribution of coverage
    print >>statsout, "\n## nloci = number of loci with data for exactly ntaxa"
    print >>statsout, "## ntotal = number of loci for which at least "+\
                         "ntaxa have data"
    print >>statsout, '\t'.join(['ntaxa', 'nloci', 'saved', 'ntotal'])
    coverage = [i.count(">") for i in finalfile.strip().split("//")[:-1]]
    if not coverage:
        sys.exit("\twarning: no loci meet 'minCov' setting"+\
                 " (line 11)\n\tno results written")
    coverage.sort()
    tot = locus
    print >>statsout, str(1)+"\t-"  
    for i in range(2, max(set(coverage))+1):
        if i >= params["minsamp"]:
            tot -= coverage.count(i-1)
            print >>statsout, str(i)+"\t"+\
                              str(coverage.count(i))+\
                              "\t*\t"+str(tot)
        else:
            print >>statsout, str(i)+"\t-\t\t-"

    ## print distribution of snps and pis
    print >>statsout, "\n## nvar = number of loci containing n "+\
                      "variable sites (pis+autapomorphies)."
    print >>statsout, "## sumvar = sum of variable sites (SNPs)."
    print >>statsout, "## pis = number of loci containing n "+\
                      "parsimony informative sites."
    print >>statsout, "## sumpis = sum of parsimony informative sites."    
    print >>statsout, '\t'.join(['n', 'nvar', 'sumvar', 'PIS', 'sumPIS'])
    snps = [line.count("-")+line.count("*") for line in \
             finalfile.split("\n") if "|" in line]
    pis = [line.count("*") for line in finalfile.split("\n") if "|" in line]
    zero = sum([line.count("*")+line.count("-") == 0 for line in \
                 finalfile.split("\n") if "|" in line])

    print >>statsout, str(0)+"\t"+str(zero)+"\t"+str(0)+\
                      "\t"+str(pis.count(0))+"\t"+str(0)

    for i in range(1, max(snps)+1):
        sumvar = sum([(j)*snps.count(j) for j in range(1, i+1)])
        sumpis = sum([(j)*pis.count(j) for j in range(1, i+1)])
        print >>statsout, str(i)+"\t"+str(snps.count(i))+"\t"+\
                          str(sumvar)+"\t"+str(pis.count(i))+"\t"+str(sumpis)
    totalvar = sum(snps) #+sum(pis)
    print >>statsout, "total var=", totalvar
    print >>statsout, "total pis=", sum(pis)



def makehaplos(params, longname):
    """ split loci with heterozygous sites into the 
        two phased alleles for diploids """

    ## outfile to write to
    outfile = open(params["work"]+"outfiles/"+\
                   params["outname"]+".alleles", 'w')
    ## data from .loci file
    lines = open(params["work"]+"outfiles/"+\
                 params["outname"]+".loci").readlines()
    ## store data until ready to write
    writing = []
    ## counter to know when to write
    loc = 0
    for line in lines:
        if ">" in line:
            a, b = line.split(" ")[0], line.split(" ")[-1]
            a1, a2 = breakalleles(b.strip())
            writing.append(a+"_0"+" "*(longname-len(a)+3)+a1)
            writing.append(a+"_1"+" "*(longname-len(a)+3)+a2)
        else:
            writing.append(line.strip())
        loc += 1

        ## print every 10K loci "
        if not loc % 10000:
            outfile.write("\n".join(writing)+"\n")
            writing = []
            
    outfile.write("\n".join(writing))
    outfile.close()



def main(params, infile, taxadict, minhits, version, quiet):
    """ setup to call functions """

    ## remove old temp files
    for i in glob.glob(params["work"]+".chunk_*"):
        os.remove(i)
    for i in glob.glob(params["work"]+".align_*"):
        os.remove(i)
    for i in glob.glob(params["work"]+".not_*"):
        os.remove(i)

    ## create output directory "
    if not os.path.exists(params["work"]+'outfiles'):
        os.makedirs(params["work"]+'outfiles')
    if not os.path.exists(params["work"]+'stats'):
        os.makedirs(params["work"]+'stats')

    ## read names from file
    temp = gzip.open(infile, 'r').readlines()
    names = set(["_".join(i.split(">")[1].split("_")[:-2]) \
                for i in temp if ">" in i])

    ## find subset names "
    params["subset"] = set([i for i in names if params["subset"] in i])

    ## remove excludes and outgroups from list "
    if params["exclude"]:
        exclude = params["exclude"].strip().split(",")
    else:
        exclude = []
    exclude += list(names.difference(params["subset"]))
    if params["outgroup"]:
        outgroup = params["outgroup"].strip().split(",")
    else:
        outgroup = []
    for i in exclude:
        names.discard(i)
    ingroup = copy(names)
    for i in outgroup:
        if i in ingroup:
            ingroup.remove(i)

    if not quiet:
        screensamples(ingroup, outgroup, exclude)

    if len(ingroup) < 2:
        sys.exit("\n\twarning: must have at least two samples "+\
                 "selected for inclusion in the data set ")

    ## dont allow more processors than available on machine
    if params["parallel"] > multiprocessing.cpu_count():
        params["parallel"] = multiprocessing.cpu_count()

    ## find longest name for prettier output files
    longname = max([len(i) for i in list(ingroup)+list(outgroup)])

    ## check if output files already exist with this outname prefix
    if os.path.exists(params["work"]+"outfiles/"+params["outname"]+".loci"):
        sys.stderr.write("\n\tWarning: data set "+params["outname"]+\
                 ".loci already exists"+\
                 "\n\t  Skipping re-alignment. Creating extra data formats "+\
                 "\n\t  from the existing .loci file. To create a new .loci"+\
                 "\n\t  file and stats output move/delete "+params["outname"]+\
                 ".loci\n\t  or change the outname prefix in the params file\n")

    else:
        ## split up clusters and align on different nodes
        splitandalign(params, infile, ingroup,
                      exclude, longname, quiet)

        ## make .loci file and get filtered data
        locus, filterlist, new2olddict = makelocifile(params)
        ## pickle the new2olddict
        pickleloc = gzip.open(params["work"]+"clust"+\
                              params["wclust"]+"/"+
                              params["outname"]+".new2olddict", 'wb')
        pickle.dump(new2olddict, pickleloc)
        pickleloc.close()

        ## make stats output
        dostats(params, ingroup, outgroup, longname,
                locus, filterlist, version, quiet)

    ## make other formatted files "
    if "*" in params["outform"]:
        params["outform"] = ",".join(list("pnasvutmkgfc"))
    formats = params["outform"].split(",")

    ## make alleles output
    if "a" in formats:
        if not quiet:
            sys.stderr.write("\twriting alleles file\n")
        makehaplos(params, longname)

    ## make phy and/or nex
    if any([i in formats for i in ['n', 'p']]):
        if not quiet:
            if 'n' in formats:
                sys.stderr.write("\twriting nexus file\n")
            if 'p' in formats:
                sys.stderr.write("\twriting phylip file\n")
        loci2phynex.make(params, names, longname, formats)

    ## make gphocs format
    if 'f' in formats:
        if not quiet:        
            sys.stderr.write("\twriting gphocs file\n")
        loci2gphocs.make(params["work"], params["outname"])

    ## formats that depend on snp files
    if any([i in formats for i in list("usktg")]):
        if not quiet:
            if 's' in formats:
                sys.stderr.write("\twriting full SNPs file\n")
            if 'u' in formats:
                sys.stderr.write("\twriting unlinked SNPs file\n")
            if 'k' in formats:
                sys.stderr.write("\twriting STRUCTURE file\n")
            if 'g' in formats:
                sys.stderr.write("\twriting geno file\n")
        loci2SNP.make(params, names)

    ## make treemix output (uses the usnp files above)
    if "t" in formats:
        if taxadict.keys():
            if not quiet:
                sys.stderr.write("\t  + writing treemix file\n")
            loci2treemix.make(params, taxadict, minhits, quiet)
        else:
            print "\t  ** must enter group/clade assignments"+\
                  "for treemix output "

    ## make migrate output 
    if 'm' in formats: #[i[0] for i in formats]:
        migstring = [i for i in formats if i[0] == 'm'][0]
        if len(migstring) > 1:
            maxnumberloci = int(migstring.strip()[1:])
        else:
            maxnumberloci = int(9e6)
        if taxadict.keys():
            if not quiet:
                sys.stderr.write("\twriting migrate-n file\n")
            loci2mig.make(params, taxadict, minhits, maxnumberloci)
        else:
            print "\t  ** must enter group/clade assignments for migrate-n output "

    ## make vcf
    if 'v' in formats:
        if not quiet:
            sys.stderr.write("\twriting vcf file\n")
        loci2vcf.make(params, version, names)
    
    ## make cat output 
    if 'c' in formats:
        if not quiet:
            sys.stderr.write("\twriting cat file\n")
        ## check for dependencies...
        loci2cat.make(params, names, quiet)



if __name__ == "__main__":
    main()