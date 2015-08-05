#! /usr/bin/env python2

""" 
de-replicates edit files and clusters de-replciated reads 
by sequence similarity using vsearch
"""

import multiprocessing
import sys
import os
import numpy
import itertools
import glob
import subprocess
import operator
import gzip
import re
import fileinput
from potpour import Worker


def makederepclust(params, outfolder, handle):
    """ combines information from .u and ._temp files 
    to create .clust files, which contain un-aligned clusters """

    ## if .u files are present this read them in as userout
    if os.path.exists(outfolder+"/"+\
                      handle.split("/")[-1].\
                      replace(".edit", ".u")):
        userout = open(outfolder+"/"+handle.split("/")[-1].\
                       replace(".edit", ".u"), 'r').readlines()
    else:
        userout = []
        print "\n\tSkipping: no '.u' file available for sample %s" % (
               handle.split("/")[-1])

    ## create an output file to write cluster to        
    outfile = gzip.open(outfolder+"/"+\
                        handle.split("/")[-1].\
                        replace(".edit", ".clust.gz"), 'w')

    ## load reads into a dictionary"
    hits = {}  ## D = {}
    dereps = open(handle.replace(".edit", ".derep")).read()  ## f
    for line in dereps.split(">")[1:]:
        # Workaround to comply with both vsearch and usearch
        line = line.replace(";\n", "\n", 1) 
        a, b, c = line.replace("\n", ";", 1).replace("\n", "").split(";")
        hits[">"+a+";"+b+";"] = [int(b.replace("size=", "")), c.strip()]

    ## create dictionary of .u file cluster hits info "
    udic = {}  ## U
    for uline in [line.split("\t") for line in userout]:
        # Workaround to comply with both vsearch and usearch
        if uline[1].endswith(";") == False: 
            uline[1] += ";"
            uline[0] += ";"
            uline[5] = re.sub("\..*\n", "\n", uline[5])
        if ">"+uline[1] in udic:
            udic[">"+uline[1]].append([">"+uline[0], uline[4],
                                           uline[5].strip(),
                                           uline[3]])
        else:
            udic[">"+uline[1]] = [[">"+uline[0], uline[4],
                                       uline[5].strip(),
                                       uline[3]]]

    ## map sequences to clust file in order
    seq = ""
    seqslist = []  ## SEQS
    for key, values in udic.items():
        seq = key+"\n"+hits[key][1]+'\n'  
        matchnames = [i[0] for i in values]       ## names of matches (S)
        forwrev = [i[1] for i in values]       ## + or - for strands (R)
        cov = [int(float(i[2])) for i in values]  ## query coverage (overlap)
        ins = [int(i[3]) for i in values]  ## indels

        ## allow only 'w1' indels in hits to seed
        if not any([int(i) > int(params["w1"]) for i in ins]):
            for i in range(len(matchnames)):
                if forwrev[i] == "+":

                    ## only match forward reads if high Cov
                    if cov[i] >= 90:
                        seq += matchnames[i]+'+\n'+\
                               hits[matchnames[i]][1]+"\n"
                else:
                    ## name change for reverse hits
                    ## allow low Cov for reverse hits
                    seq += matchnames[i].replace("_r1;", "_c1;")+\
                               '-\n'+comp(hits[matchnames[i]][1][::-1])+"\n"
        seqslist.append(seq)
    outfile.write("//\n//\n".join(seqslist))

    ## make Dict. from seeds (_temp files) 
    seedsdic = {}     ##  I
    invar = open(outfolder+"/"+handle.split("/")[-1].\
                               replace(".edit", "._temp"), 'r')
    invarlist = [(">"+i.replace("\n", "")).split(';') for i \
                      in invar.read().split(">")[1:]]
    invar.close()

    ## fill the seedsdic dictionary with seeds from invarlist
    for i in invarlist:
        try: 
            seedsdic[i[0]+';'+i[1]+';'] = i[2]
        except IndexError:
            pass  ## skip binary errors 
    del invarlist

    ## create a set for keys in I not in seedsdic
    set1 = set(seedsdic.keys())       ## temp file (no hits) seeds
    set2 = set(udic.keys())       ## u file (with hits) seeds
    diff = set1.difference(set2)  ## seeds in 'temp not matched to in 'u
    if len(diff) > 1:
        for i in diff:
            hits[i][1].replace("n", "Z").upper().replace("Z", "n")
            outfile.write("//\n//\n"+i+"\n"+hits[i][1]+'\n')
    else:
        if diff:
            popped = diff.pop()
            hits[hits.keys()[0]][1].replace("n", "Z").upper().replace("Z", "n")            
            outfile.write(popped+"\n"+hits[hits.keys()[0]][1]+'\n')
    outfile.write("//\n//\n")
    outfile.close()
    del dereps
    del userout


def derep(params, handle):
    """ dereplicates reads and write to .step file """
    if 'vsearch' in params["vsearch"]:
        threads = " "
    else:
        threads = " -threads 1"
    if params["datatype"] in ['pairgbs', 'gbs', 'merged']:
        reverse = " -strand both "
    else:
        reverse = " "
    if params["minuniq"]:
        mins = " -minuniquesize "+str(params["minuniq"])
    else:
        mins = " "
    cmd = params["vsearch"]+\
        " -derep_fulllength "+handle+\
        reverse+\
        mins+\
        " -output "+handle.replace(".edit", ".step")+\
        " -sizeout "+\
        threads
    subprocess.call(cmd, shell=True, 
                         stderr=subprocess.STDOUT, 
                         stdout=subprocess.PIPE)


def sortbysize(params, handle):
    """ sorts dereplicated file (.step) so reads that were highly
    replicated are at the top, and singletons at bottom, writes
    output to .derep file """
    cmd = params["vsearch"]+\
          " -sortbysize "+handle.replace(".edit", ".step")+\
          " -output "+handle.replace(".edit", ".derep")
    subprocess.call(cmd, shell=True, 
                         stderr=subprocess.STDOUT,
                         stdout=subprocess.PIPE)
    if os.path.exists(handle.replace(".edit", ".step")):
        os.remove(handle.replace(".edit", ".step"))


## SKIP FOR VSEARCH
def splitbigfilesforderep(dfiles, params):
    """ work around 4GB limit of 32-bit usearch
    by splitting files for derep then rejoining """
    for handle in dfiles:
        if not os.path.exists(handle.replace(".edit", ".derep")):
            ## break every 7.5M reads
            bsize = 15000000
            statinfo = os.stat(handle)
            size = statinfo.st_size
            ## if file size is over 3G
            if size > 300000000:
                with open(handle, 'r') as infile:
                    alllines = infile.readlines()
                breaks = len(alllines)/bsize
                lines = 0
                if breaks > 1:
                    for brake in range(breaks+1):
                        out = open(handle+"_piece_"+str(brake), 'w')
                        out.write("".join(alllines[lines:lines+bsize]))
                        out.close()
                        derep(params, handle+"_piece_"+str(brake))
                        lines += bsize
                else:
                    brakes = 0
                    with open(handle+"_piece_"+str(brakes), 'w') as out:
                        out.write("".join(alllines[lines:lines+bsize]))
                    derep(params, handle+"_piece_"+str(brakes))
                    lines += bsize

                with open(handle+"_piece_"+str(breaks+1), 'w') as out:
                    out.write("".join(alllines[lines:]))
                del alllines
                derep(params, handle+"_piece_"+str(brakes+1))

                sublist = glob.glob(handle.replace(".edit", ".step")+"_piece*")
                if len(sublist) > 0:
                    fout = open(handle.replace(".edit", ".derep"))
                    for line in fileinput.input(sublist):
                        fout.write(line)        
                for rmfile in glob.glob(handle.\
                              replace(".edit", ".step")+"_piece*"):
                    os.remove(rmfile)
        else:
            print 'skipping derep of '+handle.replace(".edit", ".derep")+\
                                                        ', aleady exists'


def fullcluster(params, outfolder, handle):
    """ calls vsearch for clustering """
    if params["datatype"] == 'pairddrad':
        comm = " -cluster_smallmem "+handle.replace(".edit", ".firsts")
    else:
        comm = " -cluster_smallmem "+handle.replace(".edit", ".derep")
    if params["datatype"] in ['gbs', 'merged']:
        reverse = " -strand both "
        cov = " -query_cov .35 " 
    elif params["datatype"] == 'pairgbs':
        reverse = " -strand both "
        cov = " -query_cov .60 " 
    else:     ## rad, ddrad, ddradmerge
        reverse = " -leftjust "
        cov = " -query_cov .60"
    ## if vsearch and not usearch
    if 'vsearch' not in params["vsearch"]:
        masker = " "
    else:
        masker = " -qmask "+params["mask"]
    cmd = params["vsearch"]+\
        comm+\
        reverse+\
        cov+\
        masker+\
        " -id "+params["wclust"]+\
        " -userout "+outfolder+"/"+\
                     handle.split("/")[-1].replace(".edit", ".u")+\
        " -userfields query+target+id+gaps+qstrand+qcov"+\
        " -maxaccepts 1"+\
        " -maxrejects 0"+\
        " -minsl 0.5"+\
        " -fulldp"+\
        " -threads "+str(params["threads"])+\
        " -usersort "+\
        " -notmatched "+outfolder+"/"+handle.split("/")[-1].\
                                     replace(".edit", "._temp")
    subprocess.call(cmd, shell=True,
                         stderr=subprocess.STDOUT,
                         stdout=subprocess.PIPE)



def stats(params, outfolder, handle, multihits, quiet):
    """ return stats output after clustering is finished """
    temphandle = outfolder+"/"+handle.split("/")[-1].\
                               replace(".edit", ".clustS.gz")
    infile = gzip.open(temphandle)
    duo = itertools.izip(*[iter(infile)]*2)
    try: 
        ifas = duo.next()[0]   ## a
    except StopIteration: 
        print "no clusters found in ", temphandle+"\n\n"
        sys.exit()
    depth = []
    thisdepth = int(ifas.split(";")[1].replace("size=", ""))
    while 1:
        try: 
            ifas = duo.next()[0]
        except StopIteration: 
            break
        if ifas != "//\n":
            thisdepth += int(ifas.split(";")[1].replace("size=", ""))
        else:
            depth.append(thisdepth)
            thisdepth = 0
    infile.close()

    keep = [i for i in depth if i >= params["mindepth"]]
    namecheck = temphandle.split("/")[-1].replace(".clustS.gz", "")
    if depth:
        me = round(numpy.mean(depth), 3)
        std = round(numpy.std(depth), 3)
    else:
        me = std = 0.0
    if keep:
        mek = round(numpy.mean(keep), 3)
        stdk = round(numpy.std(keep), 3)
    else:
        mek = stdk = 0.0
    out = [namecheck, len(depth),
           me, std, len(keep), 
           mek, stdk, multihits]

    bins = [0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 100, 250, 500, 99999]
    ohist, edges = numpy.histogram(depth, bins)
    hist = [float(i)/sum(ohist) for i in ohist]
    hist = [int(round(i*30)) for i in hist]

    if not quiet:
        sys.stderr.write("\tsample "+handle.split("/")[-1].\
                         split(".")[0]+" finished, "+str(len(depth))+" loci\n")
    del depth, keep
    return out, edges, hist, ohist


def comp(seq):
    """ returns a seq with small complement"""
    return seq.replace("A", 't')\
           .replace('T', 'a')\
           .replace('C', 'g')\
           .replace('G', 'c')\
           .replace('n', 'Z')\
           .upper().replace("Z", "n")


def sortalign(stringnames):
    """ parses muscle output from a string to two list """
    objs = stringnames.split("\n>")
    seqs = [i.split("\n")[0].replace(">", "")+"\n"+\
              "".join(i.split('\n')[1:]) for i in objs]
    aligned = [i.split("\n") for i in seqs]
    newnames = [">"+i[0] for i in aligned]
    seqs = [i[1] for i in aligned]     
    return newnames, seqs


def alignfast(names, seqs, muscle):
    """ performs muscle alignments on cluster and returns output as string"""
    inputstring = "\n".join('>'+i+'\n'+j for i, j in zip(names, seqs))
    cmd = "/bin/echo '"+inputstring+"' | "+muscle+" -quiet -in -"
    piped = subprocess.Popen(cmd, shell=True, 
                                  stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE, 
                                  stderr=subprocess.STDOUT,
                                  close_fds=True)
    (_, fout) = (piped.stdin, piped.stdout)
    return fout.read()


def alignwrappair(params, handle):
    """ same as alignwrap but for pairddrads,
        feeds in first and second reads separately """
    ## iterator for 2 lines at a time
    infile = gzip.open(handle)
    duo = itertools.izip(*[iter(infile)]*2)

    ## lists for storing first and second aligned loci
    out = []
    outp = []
    cnts = 0
    multihits = 0

    while 1:
        try:
            first = duo.next()
        except StopIteration:
            break
        itera = [first[0], first[1]]
        names = []   
        seqs = []
        stack = []
        badpair = []
        nameiter = 0

        ##read in all data for this stack "
        while itera[0] != "//\n":
            names.append(itera[0].strip()+"_"+str(nameiter))
            seqs.append(itera[1].strip())  ## .upper())
            itera = duo.next()
            nameiter += 1

        ## if longer than 1 it needs aligning "
        if len(names) > 1:
            firsts = [i.split("n")[0] for i in seqs]
            seconds = [i.split("n")[-1] for i in seqs]

            ## align first reads "
            stringnames = alignfast(names[0:200], 
                                    firsts[0:200],
                                    params["muscle"])
            anames1, aseqs1 = sortalign(stringnames)
            somedic1 = {}
            for i in range(len(anames1)):
                somedic1[anames1[i]] = aseqs1[i]
            
            ## reorder keys by nameiter order
            keys = somedic1.keys()
            keys.sort(key=lambda x: int(x.split(";")[1].\
                                        replace("size=", "")),
                                        reverse=True)

            ## align second reads "
            stringnames = alignfast(names[0:200],
                                    seconds[0:200],
                                    params["muscle"])
            anames2, aseqs2 = sortalign(stringnames)
            somedic2 = {}
            for i in range(len(anames2)):
                somedic2[anames2[i]] = aseqs2[i]
            
            ## check that second reads do not align poorly 
            badpair = any([somedic2[i].count("-") > int(params["w2"]) \
                             for i in somedic2 if '_trim' not in i])

            if not badpair:
                for key in keys:
                    stack.append("_".join(key.split("_")[:-1])+'\n'+\
                                     somedic1[key]+"nnnn"+somedic2[key])
            else:
                for key in keys:
                    badpair.append("_".join(key.split("_")[:-1])+'\n'+\
                                     somedic1[key]+"nnnn"+somedic2[key])
                multihits += 1
                
        else:
            if seqs:  ## sequence could have been trimmed
                stack.append("_".join(names[0].split("_")[:-1])+'\n'+seqs[0])

        cnts += 1
        if stack:
            out.append("\n".join(stack))
        if badpair:
            outp.append("\n".join(badpair))

        if not cnts % 5000:
            if out:
                outfile = gzip.open(infile.replace(".clust", ".clustS"), 'a')
                outfile.write("\n//\n//\n".join(out)+"\n//\n//\n")
                outfile.close()
            out = []

    outfile = gzip.open(infile.replace(".clust", ".clustS"), 'a')
    if out:
        outfile.write("\n//\n//\n".join(out)+"\n//\n//\n")
    outfile.close()
    if outp:
        outbads = gzip.open(infile.replace(".clust", ".badpairs"), 'a')
        outbads.write("\n//\n//\n".join(outp)+"\n//\n//\n")
        outbads.close()

    return multihits
    
    
def alignwrap(params, handle):
    """ splits clusters and feeds them into alignfast function """
    ## iterator for 2 lines at a time
    infile = gzip.open(handle)
    duo = itertools.izip(*[iter(infile)]*2)

    ## list for storing until writing
    out = []
    cnts = 0
    while 1:
        try: 
            first = duo.next()
        except StopIteration:
            break
        itera = [first[0], first[1]]
        stack = []
        names = []   
        seqs = []
        while itera[0] != "//\n":
            names.append(itera[0].strip())
            seqs.append(itera[1].strip().replace("nnnn", ""))
            itera = duo.next()
        if len(names) > 1:
            ## keep only the 200 most common dereps, 
            ## aligning more is surely junk
            stringnames = alignfast(names[0:200], seqs[0:200], params["muscle"])
            anames, aseqs = sortalign(stringnames)
            ## a dictionary for names2seqs post alignment and indel check
            somedic = {}
            leftlimit = 0
            for i in range(len(anames)):
                ## apply filter for number of indels again, post-alignment,
                ## this affects indels of sequences relative to each other, not
                ## just relative to the seed sequence """
                if aseqs[i].rstrip("-").lstrip("-").count("-") <= params["w1"]:
                    somedic[anames[i]] = aseqs[i]

                ## do not allow seqeuence to the left of the 
                ## seed (may include adapter/barcodes)"
                if not anames[i].split(";")[-1]:
                    leftlimit = min([aseqs[i].index(j) for j in\
                                              aseqs[i] if j != "-"])
                    
            ## reorder keys by derep number "
            keys = somedic.keys()
            keys.sort(key=lambda x: int(x.split(";")[1].\
                          replace("size=", "")), reverse=True)
            for key in keys:
                stack.append(key+'\n'+somedic[key][leftlimit:])
        else:
            if names:
                stack = [names[0]+"\n"+seqs[0]]

        if stack:
            out.append("\n".join(stack))

        cnts += 1
        ## only write to file after 5000 aligned loci
        if not cnts % 5000:
            if out:
                outfile = gzip.open(handle.replace(".clust", ".clustS"), 'a')
                outfile.write("\n//\n//\n".join(out)+"\n//\n//\n")
                outfile.close()
            out = []

    outfile = gzip.open(handle.replace(".clust", ".clustS"), 'a')
    if out:
        outfile.write("\n//\n//\n".join(out)+"\n//\n//\n")
    outfile.close()


def orderseqs(names, seqs):
    """ reorders cluster by derep number because muscle output does
        not retain the order """
    try: 
        dereps = [int(i.split(";")[1].replace("size=", "")) for i in names]
    except IndexError:
        print names
    ordered = sorted(range(len(dereps)), key=lambda a: dereps[a], reverse=True)
    nnames = [names[i] for i in ordered]
    sseqs = [seqs[i] for i in ordered]
    return nnames, sseqs



def splitter(handle):
    """splits paired reads and writes firsts to a file """
    ## read in the derep file
    lines = iter(open(handle.replace(".edit", ".derep")).\
                             read().strip().split(">")[1:])
    firsts = []
    cnts = 0
    for line in lines:
        obj = line.split('\n')
        name = obj[0]
        seq = "".join(obj[1:])
        ## legacy fix for pyrad2 -> pyrad 3
        seq = seq.replace("XXXX", "nnnn")
        ## split on nn separator
        seq = seq.split("nn")[0]
        firsts.append(">"+name+"\n"+seq)
        cnts += 1
        if not cnts % 100000:
            orderfirsts = open(handle.replace(".edit", ".firsts"), 'a')
            orderfirsts.write("\n".join(firsts))
            orderfirsts.close()
            firsts = []
    orderfirsts = open(handle.replace(".edit", ".firsts"), 'a')
    orderfirsts.write("\n".join(firsts))
    orderfirsts.close()
    

def final(params, outfolder, handle, fileno, remake, quiet):
    """ run the full script """

    multihits = 0
    if not remake:
        ## de-replicate the reads if not done by big file method"
        if handle.replace(".edit", ".derep") not in \
                          glob.glob(params["work"]+"edits/*"):
            derep(params, handle)
            sortbysize(params, handle)

        if params["datatype"] == 'pairddrad':
            if handle.replace(".edit", ".firsts") not in \
                              glob.glob(params["work"]+"edits/*"):
                splitter(handle)

        ## cluster the reads "
        fullcluster(params, outfolder, handle)

    ## build cluster files from .u & .temp files
    makederepclust(params, outfolder, handle)

    # " thread each align job x2 to reach ~100% "
    # " split file in half"
    # f = gzip.open(outfolder+"/"+handle.split("/")[-1].replace(".edit",".clust.gz"), 'rb').read().strip().split("//\n")
    # chunk1 = f/2
    # ff1 = gzip.open(outfolder+"/"+handle.split("/")[-1].replace(".edit",".clust.gz.1"))
    # ff1.write("//\n\n".join(f[:chunk1])+"//\n\n")
    # ff1.close()
    # ff2 = gzip.open(outfolder+"/"+handle.split("/")[-1].replace(".edit",".clust.gz.2"))
    # ff2.write("//\n\n".join(f[chunk1:])+"//\n\n")
    # ff2.close()

    # threads = []
    # for ff in range(1,3):
    #     if 'pair' in datatype:
    #         margs = (outfolder+"/"+handle.split("/")[-1].replace(".edit",".clust.gz."+str(ff)),
    #                               mindepth, muscle, w2,)
    #         t = threading.Thread(target=alignwrapPAIR, args=margs)
    #     else:
    #         margs = (outfolder+"/"+handle.split("/")[-1].replace(".edit",".clust.gz."+str(ff)),
    #                  mindepth, muscle, w1,)
    #         t = threading.Thread(target=alignwrap, args=margs)
    #     threads.append(t)
    #     t.start()
        
    # " combine split alignment files"
    # ff1 = outfolder+"/"+handle.split("/")[-1].replace(".edit",".clust.gz.1")
    # ff2 = outfolder+"/"+handle.split("/")[-1].replace(".edit",".clust.gz.2")
    # ff3 = outfolder+"/"+handle.split("/")[-1].replace(".edit",".clust.gz")
                                                      
    # cmd = "/bin/cat "+ff1+" "+ff2+" > "+ff3
    # os.system(cmd)
    # os.remove(ff1)
    # os.remove(ff2)    

    ## align clusters w/ muscle "
    if 'pair' in params["datatype"]:
        multihits = alignwrappair(params, 
                           outfolder+"/"+handle.split("/")[-1].\
                           replace(".edit", ".clust.gz"))

    else:
        alignwrap(params, outfolder+"/"+handle.split("/")[-1].\
                          replace(".edit", ".clust.gz"))

    ## get stats 
    out, edges, hist, ohist = stats(params,
                                    outfolder, 
                                    handle, 
                                    multihits, 
                                    quiet)

    end = handle.split("/")[-1].replace(".edit", "")
    outwrite = open(outfolder+"/.temp."+end, 'w')
    outwrite.write("\t".join([str(i) for i in out]))

    ## print histograms to file
    print >>outwrite, "\nbins\tdepth_histogram\tcnts"
    print >>outwrite, "   :\t0------------50-------------100%"

    for i, j, k in zip(edges, hist, ohist):
        firststar = " "
        if k > 0:
            firststar = "*"
        print >>outwrite, i, '\t', firststar+"*"*j + " "*(34-j), k   ## HERE
    outwrite.close()



def main(params, quiet, remake):
    """ calls the main script """

    #WORK, parallel, wclust, mindepth,
    #     subset, datatype, w1, w2, minuniq,
    #     MASK, muscle, vsearch, threads, remake):

    ## find .edit files in edits/ directory
    if not os.path.exists(params["work"]+'edits/'):
        sys.exit("\terror: could not find edits/ folder in working directory")

    ## make output folder for clusters
    if not os.path.exists(params["work"]+'clust'+params["wclust"]):
        os.makedirs(params["work"]+'clust'+params["wclust"])
    outfolder = params["work"]+'clust'+str(params["wclust"])
    if not os.path.exists(params["work"]+'stats'):
        os.makedirs(params["work"]+'stats')

    ## remake option... in development"
    if remake:
        for ufile in glob.glob(outfolder+"/*.u"):
            #infile = open(ufile).readlines()
            ## delete the last line in the file that was probably
            ## not completely written
            cmd = "/bin/sed '$d' < " + ufile + " > tempfile"
            subprocess.call(cmd, shell=True) ##os.system(cmd)
            ## make a backup b/c this isn't tested enough yet
            cmd = "/bin/mv "+ufile+" "+ufile+".backup"
            os.system(cmd)
            ## replace original that is backed up with the new file
            ## that has the last line removed.
            cmd = "/bin/mv tempfile "+ufile
            os.system(cmd)

    dfiles = []   ## FS

    ## if not only 1 sample "
    if len(glob.glob(params["work"]+"edits/"+params["subset"]+"*.edit*")) > 1:  
        for efile in glob.glob(params["work"]+"edits/"+\
                               params["subset"]+"*.edit*"):
            ## append files to list if not already clustered or empty"
            if not os.path.exists(outfolder+"/"+efile.\
                   replace(".edit", ".clustS.gz")):
                size = os.stat(efile)
                if size.st_size > 0:
                    dfiles.append(efile)
                else:
                    print "excluding "+str(efile)+" file is empty"
            else:
                print efile.replace(".edit", ".clustS")+" already exists"
        ## " arranges files by decreasing size for fast clustering order"
        for efile in range(len(dfiles)):
            statinfo = os.stat(dfiles[efile])
            dfiles[efile] = dfiles[efile], statinfo.st_size
        dfiles.sort(key=operator.itemgetter(1), reverse=True)
        dfiles = [i[0] for i in dfiles]

    ## if only one files
    elif len(glob.glob(params["work"]+"edits/"+params["subset"]+\
                                               "*.edit*")) == 1:
        dfiles = glob.glob(params["work"]+"edits/"+params["subset"]+"*.edit*")
        size = os.stat(dfiles[0])
        ## check that the file is not empty
        if size.st_size > 0:
            pass
        else:
            print "excluding "+dfiles[0]+" file is empty"
    else:
        print "\tNo .edit files found in edits/ dir."

    if not quiet:
        sys.stderr.write("\n\tde-replicating files for clustering...\n")

    ## do not split big files if using 64-bit Usearch,
    ## or if using Vsearch, else do it to avoid 4GB limit of 32-bit usearch"""

    if "vsearch" not in params["vsearch"]:
        print '\n\tsplitting big files'
        splitbigfilesforderep(dfiles, params)

    ## load work queue"
    work_queue = multiprocessing.Queue()
    result_queue = multiprocessing.Queue()

    ## perform function 'final' on files in dfiles list "
    submitted = {}
    fileno = 1

    ## if not reconstructing unfinished clusters
    if not remake:
        if params["threads"] == 0:
            nthreads = 'all'
        else:
            nthreads = params["threads"]
        nproc = min(params["parallel"], len(dfiles))
        if not quiet:
            sys.stderr.write("\n\tstep 3: within-sample clustering of "+\
                         str(len(dfiles))+" samples at \n\t        "+\
                         str(params["wclust"])+\
                         " similarity. Running "+str(nproc)+\
                         " parallel jobs\n\t"+\
                         " \twith up to "+str(nthreads)+" threads per job."+\
                         " If needed, \n\t\tadjust to avoid CPU and"+\
                         " MEM limits\n\n")
    else:
        sys.stderr.write("\n\tstep 3: rebuilding clusters "+\
                          "from unfinished step 3 files\n")

    for handle in dfiles:
        if outfolder+"/"+handle.split("/")[-1].replace(".edit", ".clustS.gz")\
                    not in glob.glob(outfolder+"/*"):
            work_queue.put([params, outfolder, handle, fileno, remake, quiet])
            # vsearch,outfolder,handle,wclust,mindepth,
            #                 parallel,muscle,datatype,fileno, w1, w2, 
            #                 WORK, minuniq, MASK, threads, remake])
            submitted[handle] = 1
            fileno += 1
        else:
            print "\tskipping "+handle.split("/")[-1].\
                                replace(".edit", ".clustS.gz")+\
                  ' already exists in '+params["work"]+outfolder.split("/")[-1]

    ## create a queue to pass to workers to store the results"
    jobs = []
    for _ in range(min(submitted, params["parallel"])):
        worker = Worker(work_queue, result_queue, final)
        jobs.append(worker)
        worker.start()
    for j in jobs:
        j.join()

    ## output statistics on depth of coverage"
    outstats = open(params["work"]+"stats/s3.clusters.txt", 'a')
    print >>outstats, '\n'+'\t'.join(['taxa', 'total', 'dpt.me',
                                      'dpt.sd', 'd>'+\
                                      str(params["mindepth"]-1)+'.tot',
                                      'd>'+str(params["mindepth"]-1)+'.me',
                                      'd>'+str(params["mindepth"]-1)+'.sd',
                                      'badpairs'])

    res = []
    histo = []
    for ffile in dfiles:
        end = ffile.split("/")[-1].replace(".edit", "")
        ffile = outfolder+"/.temp."+end
        if os.path.exists(ffile):
            line = open(ffile).readlines()
            res.append(line[0].strip().split("\t"))
            histo.append([line[0].split("\t")[0], "".join(line[1:])])
            os.remove(ffile)
    res.sort(key=lambda x: x[0])
    histo.sort(key=lambda x: x[0])
    
    for i in res:
        print >>outstats, "\t".join(i)
    
    print >>outstats, """
    ## total = total number of clusters, including singletons
    ## dpt.me = mean depth of clusters
    ## dpt.sd = standard deviation of cluster depth
    ## >N.tot = number of clusters with depth greater than N
    ## >N.me = mean depth of clusters with depth greater than N
    ## >N.sd = standard deviation of cluster depth for clusters with depth greater than N
    ## badpairs = mismatched 1st & 2nd reads (only for paired ddRAD data)\n\nHISTOGRAMS\n
    """

    for i in histo:
        print >>outstats, "sample: "+i[0]+"\n"+i[1]
    
    outstats.close()
    for handle in dfiles:
        nothere = 0
        try: 
            submitted[handle]
        except KeyError:
            nothere = 1
        if not nothere:
            if submitted[handle]:
                if os.path.exists(outfolder+"/"+handle.split("/")[-1].\
                                  replace(".edit", ".clust.gz")):
                    os.remove(outfolder+"/"+handle.split("/")[-1].\
                              replace(".edit", ".clust.gz"))

if __name__ == "__main__":
    main()

