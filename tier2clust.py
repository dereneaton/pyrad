#!/usr/bin/env python2

import os
import sys
import itertools
import numpy
import random
import glob
import subprocess
import pickle
import gzip
from cluster_cons7_shuf import comp


def cluster(UCLUST, ID, datatype, WORK, MASK):
    C = " -cluster_smallmem "+WORK+"prefix/cat.consens_"

    if datatype in ['gbs','pairgbs','mergegbs']:
        P = " -strand both"
        COV = ".90"
    else:
        P = " -leftjust "
        COV = ".90"
    if 'vsearch' not in UCLUST:
        Q = ""
        T = " -threads 1"
    else:
        Q = " -qmask "+MASK
        ## TODO: figure out optimized threads setting...
        T = " -threads 6"
    U = " -userout "+WORK+"prefix/cat.u"
    cmd = UCLUST+\
        C+\
        P+\
        " -id "+ID+\
        Q+\
        T+\
        U+\
        " -userfields query+target+id+gaps+qstrand+qcov"+\
        " -maxaccepts 1"+\
        " -maxrejects 0"+\
        " -fulldp"+\
        " -query_cov "+str(COV)+\
        " -notmatched "+WORK+"prefix/cat._tempU"
    os.system(cmd)
    #subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)


def flip(a):
    if a == "+":
        return "-"
    elif a == "-":
        return "+"


def makeclust(ID, datatype, WORK):

    " load tier 2 hits (names,direction) into a Dic with seeds as keys"
    Uin = open(WORK+"prefix/cat.u")
    Fseeds = {}    
    for line in [line.split("\t") for line in Uin.readlines()]:
        if line[1] not in Fseeds:
            Fseeds[line[1]] = [(line[0],line[4])]
        else:
            Fseeds[line[1]].append((line[0],line[4]))
    Uin.close()


    " load tier 1 hits (names,direction) into a Dictionary with seeds as keys"
    FS = glob.glob(WORK+"prefix/cat.u_*")
    Useeds = {}
    for f in FS:
        infile = open(f)
        for line in [line.split("\t") for line in infile.readlines()]:
            if line[1] not in Useeds:
                Useeds[line[1]] = [(line[0],line[4])]
            else:
                Useeds[line[1]].append((line[0],line[4]))
        infile.close()


    " Make one dictionary with combining Fseeds and Useeds matching to Fseeds"
    D = {}
    for seed in Fseeds:
        # add matches to seed to D[seed]
        Fhits = Useeds.get(seed)
        # add matches to hits to seed to D[seed]
        Mhits = []
        for hit in Fseeds[seed]:
            Mhits.append(hit)
            ugh = Useeds.get(hit[0])
            if ugh:
                if hit[1] == "-":
                    if len(ugh) == 1:
                        Mhits += [(ugh[0][0],flip(ugh[0][1]))]
                    elif len(ugh) > 1:
                        for child in ugh:
                            Mhits += [(child[0], flip(child[1]))]
                else:
                    Mhits += ugh
        if Fhits:
            D[(seed,'s')] = Fhits+Mhits
        else:
            D[(seed,'s')] = Mhits
    

    " load seeds of tier 2 into D and set its Useed hits"
    f = open(WORK+"prefix/cat._tempU")
    lines = f.readlines()
    for line in lines:
        if ">" in line:
            if (line.strip()[1:],'s') not in D:
                if Useeds.get(line.strip()[1:]):
                    D[(line.strip()[1:],'s')] = Useeds.get(line.strip()[1:])
    f.close()

    " load .consens files into Dics "
    FS = glob.glob(WORK+"clust"+ID+"/cat.consens_*.gz")
    Seqs = {}
    for f in FS:
        with gzip.open(f) as ff:
            k = itertools.izip(*[iter(ff)]*2)
            while 1:
                try: a = k.next()
                except StopIteration: break
                Seqs[a[0].strip()] = a[1].strip()
    

    " write clust file "
    outfile = gzip.open(WORK+"prefix/cat.clust_.gz", 'w')
    for i in D:
        thisclust = []
        outfile.write(">"+i[0]+'\n'+Seqs[">"+i[0]].upper()+'\n')
        thisclust.append(">"+i[0]+'\n'+Seqs[">"+i[0]].upper())
        for m in D[i]:
            if ">"+m[0]+'\n'+Seqs[">"+m[0]].upper() not in thisclust:
                if m[1] == "-":
                    outfile.write(">"+m[0]+'\n'+comp(Seqs[">"+m[0]].upper())[::-1]+'\n')
                    thisclust.append(">"+m[0]+'\n'+comp(Seqs[">"+m[0]].upper())[::-1])
                else:
                    outfile.write(">"+m[0]+'\n'+Seqs[">"+m[0]].upper()+'\n')
                    thisclust.append(">"+m[0]+'\n'+Seqs[">"+m[0]].upper())
        outfile.write("//\n")
    outfile.close()
    
    

def main(UCLUST, ID, datatype,
         gids, seed, WORK, MASK):
    
    sys.stderr.write('\n\tstep 6: clustering across cons-samples at '+`ID`+' similarity \n')

    " read in all seeds and hits "
    seeds = [WORK+"prefix/cat.seed_"+gid for gid in gids]
    temps = [WORK+"prefix/cat._temp_"+gid for gid in gids]

    #print seeds
    #print temps

    " read in all seeds and make same length for randomizing "
    out = gzip.open(WORK+'prefix/cat.group_.gz','wb')
    for handle in seeds:
        f = open(handle,'r')
        k = itertools.izip(*[iter(f)]*3)
        while 1:
            try: a = k.next()
            except StopIteration: break
            if len(a[0].strip()) < 100:
                " seriously, don't have names longer than 100 chars "
                out.write(a[0].strip()+" "*(100-len(a[0].strip()))+a[1])
            else:
                out.write(a[0].strip()+" "*((len(a[0].strip())+3)-len(a[0].strip()))+a[1])
                print "long name lengths may cause errors"
        f.close()
    out.close()

    """ randomize input order """
    if seed:
        random.seed(seed)
    with gzip.open(WORK+'prefix/cat.group_.gz','rb') as source:
        data = [ (random.random(), line) for line in source ]
    data.sort()

    """ sort by length while preserving randomization within size classes """
    D = [line for _,line in data]
    D.sort(key=len, reverse=True)
    k = iter(D)
    out = open(WORK+'prefix/cat.consens_','w')
    while 1:
        try: a = k.next().split(" ")
        except StopIteration: break
        ss = a[-1].replace("a","A").replace("g","G").replace("c","C").replace("t","T").strip()
        print >>out, a[0]+'\n'+ss
    out.close()

    cluster(UCLUST, ID, datatype, WORK, MASK)
    makeclust(ID, datatype, WORK)


    #if glob.glob(WORK+"prefix/*.seed_*") or glob.glob(WORK+"prefix/*._temp_*"):
    #    os.system("rm "+WORK+"prefix/*.seed_*")
    #    os.system("rm "+WORK+"prefix/*._temp_*")
    #    os.system("rm "+WORK+"prefix/*.u_*")

