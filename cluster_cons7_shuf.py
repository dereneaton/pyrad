#!/usr/bin/env python2
import os
import sys
import itertools
import numpy
import random
import glob
import subprocess
import gzip
import copy
from consensdp import unhetero, uplow, breakalleles



def comp(seq):
    """ returns complement of sequence including ambiguity characters,
    and saves lower case info for multiple hetero sequences"""
    seq = seq.replace("A",'u')\
             .replace('T','v')\
             .replace('C','p')\
             .replace('G','z')\
             .replace('u','T')\
             .replace('v','A')\
             .replace('p','G')\
             .replace('z','C')
    seq = seq.replace('R','u')\
             .replace('Y','v')\
             .replace('K','p')\
             .replace('M','z')\
             .replace('u','Y')\
             .replace('v','R')\
             .replace('p','M')\
             .replace('z','K')
    seq = seq.replace('r','u')\
             .replace('y','v')\
             .replace('k','p')\
             .replace('m','z')\
             .replace('u','y')\
             .replace('v','r')\
             .replace('p','m')\
             .replace('z','k')
    return seq


def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True, 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0


def cluster(vsearch, handle, ID, datatype,
            quiet, WORK, gid, MASK):

    if datatype == 'pairddrad':
        " use first files for split clustering "
        if gid:
            "hierarchical clustering save temps "
            N = " -notmatched "+WORK+"prefix/"+handle.split("/")[-1].replace(".firsts_","._temp_")
            U = " -userout "+WORK+"prefix/"+handle.split("/")[-1].replace(".firsts_",".u_")
        else:
            N = ""
            U = " -userout "+handle.replace(".firsts_",".u")
    else:
        " use haplos files "
        if gid:
            "hierarchical clustering save temps "
            N = " -notmatched "+WORK+"prefix/"+handle.split("/")[-1].replace(".haplos_","._temp_")
            U = " -userout "+WORK+"prefix/"+handle.split("/")[-1].replace(".haplos_",".u_")
        else:
            N = ""
            U = " -userout "+handle.replace(".haplos_",".u")

    C = " -cluster_smallmem "+handle
    if datatype in ['gbs','pairgbs','merge']:
        P = " -strand both"
        COV = " -query_cov .90 "  ## this can vary 
    else:
        P = " -leftjust "
        COV = " -query_cov .90 "
    if 'vsearch' not in vsearch:
        Q = ""
        T = " -threads 1"
    else:
        Q = " -qmask "+MASK
        T = " -threads 6"
    cmd = vsearch+\
        C+\
        P+\
        Q+\
        T+\
        " -id "+ID+\
        U+\
        " -userfields query+target+id+gaps+qstrand+qcov"+\
        " -maxaccepts 1"+\
        " -maxrejects 0"+\
        " -fulldp"+\
        " -usersort"+\
        COV+\
        N
    os.system(cmd)


def makeclust(handle,datatype,gid,
              minmatch,WORK):

    " read in cluster hits and seeds files "
    if not gid:
        Userout = open(handle.replace(".haplos_",".u"),'r')
        outfile = gzip.open(handle.replace(".haplos_"+gid,".clust_"+gid+".gz"),'w')
    else:
        Userout = open(WORK+'prefix/'+handle.split("/")[-1].replace(".haplos_",".u_") ,'r')
        nomatch = open(WORK+'prefix/'+handle.split("/")[-1].replace(".haplos_","._temp_"),'r')
        outfile = open(WORK+'prefix/'+handle.split("/")[-1].replace(".haplos_",".seed_"),'w')
        outfilename = WORK+'prefix/'+handle.split("/")[-1].replace(".haplos_",".seed_")

    " load full fasta file into a Dic "
    D = {}
    if datatype == 'pairddrad':
        if gid:
            f = open(handle.replace(".haplos_"+gid,".firsts_"+gid))
        else:
            f = gzip.open(handle.replace(".haplos_"+gid,".consens_"+gid+".gz"))
    else:
        f = gzip.open(handle.replace(".haplos_"+gid,".consens_"+gid+".gz"))

    L = itertools.izip(*[iter(f)]*2)
    while 1:
        try: a,b = L.next()
        except StopIteration: break
        D[a.strip()] = b.strip()
    f.close()

    " load .u info into a Dic "
    U = {}
    for line in [line.split("\t") for line in Userout.readlines()]:
        if ">"+line[1] in U:
            U[">"+line[1]].append([">"+line[0],line[4]])
        else:
            U[">"+line[1]] = [[">"+line[0],line[4]]]

    " if tier 1 of hierarchical clustering "
    if gid:
        if int(minmatch) == 1:
            " no reduction, write seeds only "
            # if datatype == 'pairddrad':
            #     singles = itertools.izip(*[iter(open(handle.replace(".haplos_",".firsts_")))]*2)
            # else:
            #     singles = itertools.izip(*[iter(open(handle))]*2)
            singles = nomatch.read().split(">")[1:]
            for i in singles:
                i = i.split("\n")[0]+"\n"+"".join(i.split("\n")[1:]).upper()
                #print ">"+i+"\n//"
                print >>outfile, ">"+i+"\n//"
                #print "//\n".join(i)
                #outfile.write("//\n".join(i))
            #    i,j = i.split('\n')[0], "\n".join(i.split('\n')[1:])
            #    outfile.write("//\n".join(i+j))
            del singles
            #outfile.write("//\n".join(LLL))
            # LLL = []
            # while 1:
            #     try: a,b = singles.next()
            #     except StopIteration: break
            #     LLL.append(a+b)
            #outfile.write("//\n".join(LLL))
            #del LLL
        else:       
            for key,values in U.items():
                ## reduction, only write seed if minimum hits reached
                if (len(values)+1) >= int(minmatch):
                    ## fix for if short seqs are excluded during clustering
                    if D.get(key):
                        seq = key+"\n"+D[key]+"\n"
                        seq += "//\n"
                        outfile.write(seq)

    else:
        " map sequences to clust file in order "
        seq = ""
        for key,values in U.items():
            if D.get(key):   ## fix for if short seqs are excluded during clustering
                seq = key+"\n"+D[key]+'\n'
                S = [i[0] for i in values]
                R = [i[1] for i in values]
                for i in range(len(S)):
                    if R[i] == "+":
                        seq += S[i] + '\n' + D[S[i]] + "\n"
                    else:
                        seq += S[i] + '\n' + comp(D[S[i]][::-1]) + "\n"
                seq += "//\n"
                outfile.write(seq)
    outfile.close()
    Userout.close()
    if gid: nomatch.close()



def splitter(handle):
    infile = open(handle)
    if os.path.exists(handle.replace(".haplos",".firsts")):
        os.remove(handle.replace(".haplos",".firsts"))
        
    orderfirsts = open(handle.replace(".haplos",".firsts"),'w')
    dp = itertools.izip(*[iter(infile)]*2)
    ff = []
    cnts = 0
    for d in dp:
        n,s = d
        ## checking fix to pairddrad splitting problem...
        ## backwards compatible with pyrad v2
        s1 = s.replace("X","x").replace("x","n").split("nn")[0]
        ff.append(n+s1+"\n")
        cnts += 1
    orderfirsts.write("".join(ff))
    orderfirsts.close()
    return handle.replace(".haplos",".firsts")



def makecons(vsearch, ID, datatype, 
             outg, seed, gid, minmatch, inlist,
             WORK, quiet, outhandle):

    " find usearch"
    if not cmd_exists(vsearch):
        print "\tcannot find usearch (or vsearch), edit path in param file"
        sys.exit()

    " make list of consens files "
    FS = [i for i in inlist if "/cat.cons" not in i]
    FS = [i for i in FS if "/cat.group" not in i]
    if not FS:
        print "no consens files found"
        sys.exit()

    " and a list including outgroups "
    fs = copy.copy(inlist)
    
    " are files gzipped ? "
    if any(['.gz' in i[-4:] for i in FS]):
        gz = ".gz"
    else:
        gz = ""

    " remove previous files if present "
    if os.path.exists(WORK+'clust'+ID+'/cat.consens_'+gid+gz):
        os.remove(WORK+'clust'+ID+'/cat.consens_'+gid+gz)
    if os.path.exists(WORK+'clust'+ID+'/cat.group_'+gid+gz):
        os.remove(WORK+'clust'+ID+'/cat.group_'+gid+gz)


    " remove outgroup sequences, add back in later to bottom after shuffling "
    if outg:
        outgroup = outg.strip().split(",")
        if len(outgroup) > 1:
            for s in outgroup:
                if WORK+"clust"+ID+"/"+s+".consens"+gz in FS:
                    FS.remove(WORK+"clust"+ID+"/"+s+".consens"+gz)
        else:
            outgroup = WORK+"clust"+ID+"/"+outg+".consens"+gz
            if outgroup in FS:
                FS.remove(outgroup)
                
    " create file with consens seqs from all taxa in list "
    out = gzip.open(WORK+'clust'+ID+'/cat.group_'+gid+gz,'w')

    for qhandle in FS:
        if gz:
            f = gzip.open(qhandle)
        else:
            f = open(qhandle)
        k = itertools.izip(*[iter(f)]*2)
        while 1:
            try: a = k.next()
            except StopIteration: break
            print >>out, a[0].strip()+"    "+a[1].strip()
        f.close()
    out.close()

    " message to shell "
    if gid:
        sys.stderr.write('\n\tstep 6: clustering across '+str(len(FS))+' samples at '+`ID`+\
                         ' similarity \n\tfor group ('+str(gid)+') retaining seeds w/ minimum of '+str(minmatch)+' hits\n\n')
    else:
        sys.stderr.write('\n\tstep 6: clustering across '+str(len(FS))+' samples at '+`ID`+' similarity \n\n')

    " make list of random number and data "
    if seed:
        random.seed(seed)
    source = gzip.open(WORK+'clust'+ID+'/cat.group_'+gid+".gz",'r')
    data = [ (random.random(), line) for line in source ]
    source.close()
    " sort by random number "
    data.sort()

    " order by size while retaining randomization within size classes "
    D = [line.split('    ') for _, line in data]
    DD = ["".join([i[0]+" "*(50-len(i[0])),i[1]]) for i in D]
    DD.sort(key=len, reverse=True)
    k = iter(["**".join([i.split(" ")[0],i.split(" ")[-1]]) for i in DD])

    " write output to .consens_.gz file "
    out = gzip.open(WORK+'clust'+ID+'/cat.consens_'+gid+".gz",'w')
    while 1:
        try: a,b = k.next().split("**")
        except StopIteration: break
        print >>out, a+'\n'+b.strip()

    
    """ add outgroup taxa back onto end of file."""
    if outg:
        " append to existing consens_ file "
        outgroup = outg.strip().split(',')
        if len(outgroup) > 1:
            for s in outgroup:
                xoutg = WORK+"clust"+ID+"/"+s+".consens.gz"
                if xoutg in fs:
                    f = gzip.open(xoutg)
                    k = itertools.izip(*[iter(f)]*2)
                    while 1:
                        try: a = k.next()
                        except StopIteration: break
                        print >>out, a[0].strip()+"\n"+a[1].strip()
                    f.close()
        elif len(outgroup) == 1:
            xoutg = WORK+"clust"+ID+"/"+outgroup[0]+".consens.gz"
            if xoutg in fs:
                f = gzip.open(xoutg)
                k = itertools.izip(*[iter(f)]*2)
                while 1:
                    try: a = k.next()
                    except StopIteration: break
                    print >>out, a[0].strip()+"\n"+a[1].strip()
                f.close()
        else:
            None
    out.close()        


    """ convert ambiguity codes into a sampled haplotype for any sample
    to use for clustering, but save ambiguities for later """

    " output file"
    outhaplos = open(outhandle,'w')

    " input file "
    infile = gzip.open(WORK+"clust"+ID+"/cat.consens_"+gid+".gz")
    lines = iter(infile.readlines())
    infile.close()
    
    " write to haplo files in fasta format "
    writinghaplos = []

    for line in lines:
        if ">" in line:
            writinghaplos.append(line.strip())
        else:
            allele = breakalleles(line)[0]
            writinghaplos.append(allele.strip())
    outhaplos.write("\n".join(writinghaplos))
    outhaplos.close()


def main(vsearch, ID, datatype, 
         outg, seed, gid, minmatch, inlist,
         WORK, MASK, quiet):

    outhandle = WORK+"clust"+ID+"/cat.haplos_"+gid

    makecons(vsearch,ID,datatype,
             outg,seed,gid,minmatch,
             inlist,WORK,quiet,outhandle)

    if datatype == 'pairddrad':
        splithandle = splitter(outhandle)
        cluster(vsearch,splithandle,ID,datatype,quiet,WORK, gid, MASK)
    else:
        cluster(vsearch,outhandle,ID,datatype,quiet,WORK, gid, MASK)

    " remake clusters with .haplos, .u, and .temp files"
    makeclust(outhandle,datatype,gid,minmatch,WORK)





