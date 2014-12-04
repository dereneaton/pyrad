#! /usr/bin/env python2
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
from potpour import Worker


## reset multithread affinity if numpy reset it
## works only for linux machines
#try: os.system("taskset -p 0xff %d" % os.getpid())
#except ValueError: None


def makederepclust(outfolder,handle,w1,datatype):
    if os.path.exists(outfolder+"/"+handle.split("/")[-1].replace(".edit",".u")):
        Userout = open(outfolder+"/"+handle.split("/")[-1].replace(".edit",".u"), 'r').readlines()
    else:
        print "\n\tSkipping: no '.u' file available for sample",handle.split("/")[-1]
        sys.exit()
    outfile = gzip.open(outfolder+"/"+handle.split("/")[-1].replace(".edit",".clust.gz"),'w')

    " load reads into a Dictionary"
    D = {}
    f = open(handle.replace(".edit",".derep")).read()
    for line in f.split(">")[1:]:
        line = line.replace(";\n","\n",1) # Workaround to comply with both vsearch and usearch
        a,b,c = line.replace("\n",";",1).replace("\n","").split(";")
        D[">"+a+";"+b+";"] = [int(b.replace("size=","")),c.strip()]

    " create dictionary of .u file cluster hits info "
    U = {}
    for line in [line.split("\t") for line in Userout]:
        if line[1].endswith(";") == False: # Workaround to comply with both vsearch and usearch
                line[1] += ";"
                line[0] += ";"
                line[5] = re.sub("\..*\n","\n", line[5])
        if ">"+line[1] in U:
            U[">"+line[1]].append([">"+line[0],line[4],line[5].strip(),line[3]])
        else:
            U[">"+line[1]] = [[">"+line[0],line[4],line[5].strip(),line[3]]]

    " map sequences to clust file in order"
    seq = ""
    SEQS = []
    for key,values in U.items():
        seq = key+"\n"+D[key][1]+'\n'
        S    = [i[0] for i in values]       ## names of matches
        R    = [i[1] for i in values]       ## + or - for strands
        Cov  = [int(float(i[2])) for i in values]  ## query coverage (overlap)
        ins = [int(i[3]) for i in values]
        " allow only 'w1' indels in hits to seed"
        if not any([int(i) > int(w1) for i in ins]):
            for i in range(len(S)):
                if R[i] == "+":
                    " only match forward reads if high Cov"
                    if Cov[i] >= 90:
                        seq += S[i]+'+\n' + D[S[i]][1] + "\n"
                else:
                    " name change for reverse hits"
                    " allow low Cov for reverse hits"
                    #seq += S[i].replace("_r1;","_c1;")+'-\n' + comp(D[S[i]][1][::-1]) + "\n"
                    seq += S[i].replace("_r1;","_c1;")+'-\n' + comp(D[S[i]][1][::-1]) + "\n"
        SEQS.append(seq)
    outfile.write("//\n//\n".join(SEQS))

    " make Dict. from seeds (_temp files) "
    I = {}
    invar = open(outfolder+"/"+handle.split("/")[-1].replace(".edit","._temp"),'r')
    A = [(">"+i.replace("\n","")).split(';') for i in invar.read().split(">")[1:]]
    invar.close()
    for i in A:
        try: I[i[0]+';'+i[1]+';'] = i[2]
        except IndexError:
            None  ## skip binary errors 
    del A

    " create a set for keys in I not in U"
    set1 = set(I.keys())       ## temp file (no hits) seeds
    set2 = set(U.keys())       ## u file (with hits) seeds
    diff = set1.difference(set2)  ## seeds in 'temp not matched to in 'u
    if len(diff) > 1:
        for i in diff:
            outfile.write("//\n//\n"+i+"\n"+D[i][1].upper()+'\n')
    else:
        if diff:
            pp = diff.pop()
            outfile.write("\n//\n//"+pp+"\n"+D[i][1].upper()+'\n')
    outfile.write("//\n//\n")
    outfile.close()
    del f
    del Userout



def derep(UCLUST, handle, datatype, minuniq):
    """ dereplicates reads and write to .step file """
    if 'vsearch' in UCLUST:
        T = ""
    else:
        T = "-threads 1"
    C = " -derep_fulllength "+handle
    if datatype in ['pairgbs','gbs']:
        P = " -strand both "
    else:
        P = " "
    if minuniq:
        M = " -minuniquesize "+str(minuniq)
    else:
        M = " "
    cmd = UCLUST+\
        C+\
        P+\
        M+\
        " -output "+handle.replace(".edit",".step")+\
        " -sizeout "+\
        T
    subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)



def sortbysize(UCLUST, handle):
    """ sorts dereplicated file (.step) so reads that were highly
    replicated are at the top, and singletons at bottom, writes
    output to .derep file """
    cmd = UCLUST+\
          " -sortbysize "+handle.replace(".edit",".step")+\
          " -output "+handle.replace(".edit",".derep")
    subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
    cmd2 = "/bin/rm "+handle.replace(".edit",".step")
    subprocess.call(cmd2, shell=True)



## SKIP FOR VSEARCH
def splitbigfilesforderep(FS,UCLUST,datatype,minuniq):
    """ work around 4GB limit of 32-bit usearch
    by splitting files for derep then rejoining """
    for handle in FS:
        if not os.path.exists(handle.replace(".edit",".derep")):
            "break every 7.5M reads"
            bsize = 15000000
            statinfo = os.stat(handle)
            size = statinfo.st_size
            " if file size is over 3G"
            if size > 300000000:
                infile = open(handle,'r')
                L = infile.readlines()
                infile.close()
                breaks = len(L)/bsize
                l = 0
                if breaks > 1:
                    for b in range(breaks+1):
                        out = open(handle+"_piece_"+str(b),'w')
                        out.write("".join(L[l:l+bsize]))
                        out.close()
                        derep(UCLUST,handle+"_piece_"+str(b),datatype,minuniq)
                        l += bsize
                else:
                    b=0
                    out = open(handle+"_piece_"+str(b),'w')
                    out.write("".join(L[l:l+bsize]))
                    out.close()
                    derep(UCLUST,handle+"_piece_"+str(b),datatype,minuniq)
                    l += bsize


                out = open(handle+"_piece_"+str(b+1),'w')
                out.write("".join(L[l:]))
                out.close()
                del L
                derep(UCLUST, handle+"_piece_"+str(b+1), datatype, minuniq)
                cmd = "/bin/cat "+handle.replace(".edit",".step")+"_piece* > "+handle.replace(".edit",".derep")
                os.system(cmd)
                cmd = "/bin/rm "+handle.replace(".edit",".step")+"_piece*"
                os.system(cmd)
                if os.path.exists(handle+"_piece_0"):
                    cmd = "/bin/rm "+handle+"_piece*"
                    os.system(cmd)
        else:
            print 'skipping derep of ', handle.replace(".edit",".derep"), ', aleady exists'



def fullcluster(UCLUST, outfolder, handle, wclust, Parallel, datatype, fileno, MASK, threads):
    if datatype == 'pairddrad':
        C = " -cluster_smallmem "+handle.replace(".edit",".firsts")
    else:
        C = " -cluster_smallmem "+handle.replace(".edit",".derep")
    if datatype in ['gbs','mergegbs']:
        P = " -strand both "
        COV = " -query_cov .35 " 
    elif datatype == 'pairgbs':
        P = " -strand both "
        COV = " -query_cov .90 " 
    else:     ## rad, ddrad, ddradmerge
        P = " -leftjust "
        COV = " -query_cov .90"
    if 'vsearch' not in UCLUST:
        Q = ""
    else:
        Q = " -qmask "+MASK
    cmd = UCLUST+\
        C+\
        P+\
        COV+\
        Q+\
        " -id "+wclust+\
        " -userout "+outfolder+"/"+handle.split("/")[-1].replace(".edit",".u")+\
        " -userfields query+target+id+gaps+qstrand+qcov"+\
        " -maxaccepts 1"+\
        " -maxrejects 0"+\
        " -minsl 0.5"+\
        " -fulldp"+\
        " -threads "+`threads`+\
        " -usersort "+\
        " -notmatched "+outfolder+"/"+handle.split("/")[-1].replace(".edit","._temp")
    subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)



def stats(outfolder, handle, mindepth, multihits):
    temphandle = outfolder+"/"+handle.split("/")[-1].replace(".edit",".clustS.gz")
    infile = gzip.open(temphandle)
    L = itertools.izip(*[iter(infile)]*2)
    try: a = L.next()[0]
    except StopIteration: print "no clusters found in ",temphandle+"\n\n"; sys.exit()
    depth = []
    d = int(a.split(";")[1].replace("size=",""))
    while 1:
        try: a = L.next()[0]
        except StopIteration: break
        if a != "//\n":
            d += int(a.split(";")[1].replace("size=",""))
        else:
            depth.append(d)
            #keep = [i for i in depth if i>=mindepth]
            #print d, numpy.mean(depth), numpy.mean(keep), len(depth), len(keep)
            d = 0
    infile.close()
    keep = [i for i in depth if i>=mindepth]
    namecheck = temphandle.split("/")[-1].replace(".clustS.gz","")
    if depth:
        me = round(numpy.mean(depth),3)
        std = round(numpy.std(depth),3)
    else:
        me = std = 0.0
    if keep:
        mek = round(numpy.mean(keep),3)
        stdk = round(numpy.std(keep),3)
    else:
        mek = stdk = 0.0
    out = [namecheck, len(depth),
           me, std, len(keep), mek, stdk, multihits]

    bins = [0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 100, 250, 500, 99999]
    ohist, edges = numpy.histogram(depth,bins)
    hist = [float(i)/sum(ohist) for i in ohist]
    hist = [int(round(i*30)) for i in hist]

    sys.stderr.write("\tsample "+handle.split("/")[-1].split(".")[0]+" finished, "+str(len(depth))+"loci\n")
    del depth,keep
    return out,edges,hist,ohist


def comp(seq):
    return seq.replace("A",'t')\
           .replace('T','a')\
           .replace('C','g')\
           .replace('G','c')\
           .upper()


def sortalign(stringnames):
    """ parses muscle output from a string to two list """
    G = stringnames.split("\n>")
    GG = [i.split("\n")[0].replace(">","")+"\n"+"".join(i.split('\n')[1:]) for i in G]
    aligned = [i.split("\n") for i in GG]
    nn = [">"+i[0] for i in aligned]    #+["//"] #[">"+aligned[0][0]]
    seqs = [i[1] for i in aligned]      #+["//"]   #[aligned[0][1]]
    return nn,seqs



def alignfast(names,seqs,muscle):
    """ performs muscle alignments on cluster and returns output as string"""
    ST = "\n".join('>'+i+'\n'+j for i,j in zip(names,seqs))
    cmd = "/bin/echo '"+ST+"' | "+muscle+" -quiet -in -"
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
    (fin, fout) = (p.stdin, p.stdout)
    return fout.read()



def alignwrapPAIR(infile,mindepth,muscle,w2):
    """ same as alignwrap but for pairddrads,
    feeds in first and second reads separately """
    f = gzip.open(infile)
    k = itertools.izip(*[iter(f)]*2)
    OUT = []
    OUTP = []
    cnts = 0
    multihits = 0
    while 1:
        try: first = k.next()
        except StopIteration: break
        cnts += 1
        itera = [first[0],first[1]]
        names = []   
        seqs = []
        STACK = []
        BADPAIR = []
        nameiter = 0
        " read in all data for this stack "
        while itera[0] != "//\n":
            names.append(itera[0].strip()+"_"+str(nameiter))
            seqs.append(itera[1].strip().upper())
            itera = k.next()
            nameiter += 1

        " if longer than 1 it needs aligning "
        if len(names) > 1:
            firsts  = [i.split("X")[0] for i in seqs]
            seconds = [i.split("X")[-1] for i in seqs]
            
            " align first reads "
            stringnames = alignfast(names[0:200],firsts[0:200],muscle)
            nn, ss = sortalign(stringnames)
            D1 = {}
            for i in range(len(nn)):
                D1[nn[i]] = ss[i]
            " reorder keys by nameiter order "
            keys = D1.keys()
            #keys.sort(key=lambda x:int(x.split("_")[-1]),reverse=True)
            keys.sort(key=lambda x:int(x.split(";")[1].replace("size=","")), reverse=True)

            " align second reads "
            stringnames = alignfast(names[0:200],seconds[0:200],muscle)
            nn, ss = sortalign(stringnames)
            D2 = {}
            for i in range(len(nn)):
                D2[nn[i]] = ss[i]
            
            " check that second reads do not align poorly "
            badpair = any([i.count("-")>int(w2) for i in D2.values()])

            if not badpair:
                for key in keys:
                    STACK.append("_".join(key.split("_")[:-1])+'\n'+D1[key]+"XXXX"+D2[key])
            else:
                for key in keys:
                    BADPAIR.append("_".join(key.split("_")[:-1])+'\n'+D1[key]+"XXXX"+D2[key])
                multihits += 1
                
        else:
            STACK.append("_".join(names[0].split("_")[:-1])+'\n'+seqs[0])

        if STACK:
            OUT.append("\n".join(STACK))
        if BADPAIR:
            OUTP.append("\n".join(BADPAIR))

        if not cnts % 1000:
            if OUT:
                outfile = gzip.open(infile.replace(".clust",".clustS"),'a')
                outfile.write("\n//\n//\n".join(OUT)+"\n//\n//\n")
                outfile.close()
            OUT = []

    outfile = gzip.open(infile.replace(".clust",".clustS"),'a')
    if OUT:
        outfile.write("\n//\n//\n".join(OUT)+"\n//\n//\n")
    outfile.close()
    if OUTP:
        outbads = gzip.open(infile.replace(".clust",".badpairs"),'a')
        outbads.write("\n//\n//\n".join(OUTP)+"\n//\n//\n")
        outbads.close()

    return multihits


    
def alignwrap(infile,mindepth,muscle,w1):
    """ splits clusters and feeds them into alignfast function """
    f = gzip.open(infile)
    k = itertools.izip(*[iter(f)]*2)
    OUT = []
    cnts = 0
    while 1:
        try: first = k.next()
        except StopIteration: break
        cnts += 1
        itera = [first[0],first[1]]
        STACK = []
        names = []   
        seqs = []
        while itera[0] != "//\n":
            names.append(itera[0].strip())
            seqs.append(itera[1].strip().upper())
            itera = k.next()
        if len(names) > 1:
            " keep only the 200 most common dereps, aligning more is surely junk "
            stringnames = alignfast(names[0:200],seqs[0:200],muscle)
            nn, ss = sortalign(stringnames)
            D1 = {}

            leftlimit = 0
            for i in range(len(nn)):
                """ apply filter for number of indels again, post-alignment,
                this affects indels of sequences relative to each other, not
                just relative to the seed sequence """
                if ss[i].rstrip("-").lstrip("-").count("-") <= w1:
                    D1[nn[i]] = ss[i]
                #if 'gbs' in datatype: 
                " do not allow seqeuence to the left of the seed (may include adapter/barcodes)"
                if not nn[i].split(";")[-1]:
                    leftlimit = min([ss[i].index(j) for j in ss[i] if j!="-"])
                    
            " reorder keys by derep number "
            keys = D1.keys()
            keys.sort(key=lambda x:int(x.split(";")[1].replace("size=","")), reverse=True)
            for key in keys:
                STACK.append(key+'\n'+D1[key][leftlimit:])
        else:
            if names:
                STACK = [names[0]+"\n"+seqs[0]]

        if STACK:
            OUT.append("\n".join(STACK))

        if not cnts % 1000:
            if OUT:
                outfile = gzip.open(infile.replace(".clust",".clustS"),'a')
                outfile.write("\n//\n//\n".join(OUT)+"\n//\n//\n")
                outfile.close()
            OUT = []

    outfile = gzip.open(infile.replace(".clust",".clustS"),'a')
    if OUT:
        outfile.write("\n//\n//\n".join(OUT)+"\n//\n//\n")
    outfile.close()


def orderseqs(nn,seqs):
    """ reorders cluster by derep number because muscle output does
    not retain the order """
    try: dereps = [int(i.split(";")[1].replace("size=","")) for i in nn]
    except IndexError:
        print nn
    ordered = sorted(range(len(dereps)), key=lambda a: dereps[a], reverse=True)
    nnames = [nn[i] for i in ordered]
    sseqs = [seqs[i] for i in ordered]
    return nnames, sseqs



def splitter(handle):
    "splits paired reads and writes firsts to a file"
    dpairs = iter(open(handle.replace(".edit",".derep")).read().strip().split(">")[1:])
    ff = []
    cnts = 0
    for d in dpairs:
        dd = d.split('\n')
        n = dd[0]
        s = "".join(dd[1:])
        s1,s2 = s.split("xxxx")
        ff.append(">"+n+'\n'+s1+"\n")
        cnts += 1
        if not cnts % 10000:
            orderfirsts = open(handle.replace(".edit",".firsts"),'a')
            orderfirsts.write("\n".join(ff))
            orderfirsts.close()
            ff = []
    orderfirsts = open(handle.replace(".edit",".firsts"),'a')
    orderfirsts.write("\n".join(ff))
    orderfirsts.close()
    

def final(UCLUST, outfolder, handle, wclust, mindepth,
          Parallel, muscle, datatype, fileno,
          w1, w2, WORK, minuniq, MASK, threads, remake):

    multihits = 0
    if not remake:
        " de-replicate the reads if not done by big file method"
        if handle.replace(".edit",".derep") not in glob.glob(WORK+"edits/*"):
            derep(UCLUST, handle, datatype, minuniq)
            sortbysize(UCLUST,handle)

        if datatype == 'pairddrad':
            if handle.replace(".edit",".firsts") not in glob.glob(WORK+"edits/*"):
                splitter(handle)

        " cluster the reads "
        fullcluster(UCLUST, outfolder, handle, wclust, Parallel, datatype, fileno, MASK, threads)

    " build cluster files from .u & .temp files "
    makederepclust(outfolder, handle, w1, datatype)

    ## if using usearch align data right after clustering using the same processor
    ## but if using vsearch then align data later ...
    " align clusters w/ muscle "
    if 'pair' in datatype:
        multihits = alignwrapPAIR(outfolder+"/"+handle.split("/")[-1].replace(".edit",".clust.gz"),
                                  mindepth, muscle, w2)
    else:
        alignwrap(outfolder+"/"+handle.split("/")[-1].replace(".edit",".clust.gz"), mindepth, muscle, w1)

    out,edges,hist,ohist = stats(outfolder, handle, mindepth, multihits)
    end = handle.split("/")[-1].replace(".edit","")
    outwrite = open(outfolder+"/.temp."+end,'w')
    outwrite.write("\t".join([str(i) for i in out]))

    print >>outwrite, "\nbins\tdepth_histogram\tcnts"
    print >>outwrite, "   :\t0------------50-------------100%"

    for i,j,k in zip(edges,hist,ohist):
        firststar = " "
        if k>0:
            firststar = "*"
        print >>outwrite, i,'\t', firststar+"*"*j + " "*(34-j),k   ## HERE
    outwrite.close()
                   

def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True, 
        stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0


def main(WORK, Parallel, wclust, mindepth, UCLUST,
         subset, datatype, muscle, w1, w2, 
         minuniq, MASK, remake):


    " find usearch"
    if not cmd_exists(UCLUST):
        print "\tcannot find usearch, edit path in input file"
        sys.exit()

    " usearch or vsearch threading method"
    if 'vsearch' in UCLUST:
        ## trying current default of 6 threads per processor
        ## could do testing to find optimal partitioning...
        ## memory limitations become a concern eventually
        threads = 6
    else:
        threads = 1

    " find muscle"
    if not cmd_exists(muscle):
        print "\tcannot find muscle, edit path in input file"
        sys.exit()

    " find .edit files in edits/ directory "
    if not os.path.exists(WORK+'edits/'):
        print "\terror: could not find edits/ folder in working directory"
        sys.exit()

    " make output folder for clusters" 
    if not os.path.exists(WORK+'clust'+wclust):
        os.makedirs(WORK+'clust'+wclust)
    outfolder = WORK+'clust'+str(wclust)
    if not os.path.exists(WORK+'stats'):
        os.makedirs(WORK+'stats')

    " remake option... in development"
    if remake:
        for ufile in glob.glob(outfolder+"/*.u"):
            infile = open(ufile).readlines()
            cmd = "/bin/sed '$d' < " + ufile + " > tempfile"
            os.system(cmd)
            cmd = "/bin/mv "+ufile+" "+ufile+".backup"
            os.system(cmd)
            cmd = "/bin/mv tempfile "+ufile
            os.system(cmd)

    FS = []

    " if not only 1 sample "
    if len(glob.glob(WORK+"edits/"+subset+"*.edit*")) > 1:  
        for f in glob.glob(WORK+"edits/"+subset+"*.edit*"):
            " append files to list if not already clustered or empty"
            if not os.path.exists(outfolder+"/"+f.replace(".edit",".clustS.gz")):
                size = os.stat(f)
                if size.st_size > 0:
                    FS.append(f)
                else:
                    print "excluding "+str(f)+" file is empty"
            else:
                print f.replace(".edit",".clustS")+" already exists"
        " arranges files by decreasing size for fast clustering order"
        for i in range(len(FS)):
            statinfo = os.stat(FS[i])
            FS[i] = FS[i],statinfo.st_size
        FS.sort(key=operator.itemgetter(1), reverse = True)
        FS = [i[0] for i in FS]

    else:
        f = glob.glob(WORK+"edits/"+subset+"*.edit*")
        size = os.stat(f[0])
        if size.st_size > 0:
            FS = f
        else:
            print "excluding "+f[0]+" file is empty"

    sys.stderr.write("\n\tde-replicating files for clustering...\n")


    """ do not split big files if using 64-bit Usearch,
    or if using Vsearch, else do it to avoid 4GB limit of 32-bit usearch"""
    if not remake:
        if "64" not in UCLUST:
            if "vsearch" not in UCLUST:
                splitbigfilesforderep(FS, UCLUST, datatype, minuniq)

    " load work queue"
    work_queue = multiprocessing.Queue()
    result_queue = multiprocessing.Queue()

    " perform function 'final' on files in FS list "
    submitted = {}
    fileno = 1

    if not remake:
        sys.stderr.write("\n\tstep 3: within-sample clustering of "+\
                         `len(FS)`+" samples at \n\t        "+`wclust`+\
                         " similarity using up to "+`Parallel`+" processors\n")
    else:
        sys.stderr.write("\n\tstep 3: rebuilding clusters from unfinished step 3 files\n")

    for handle in FS:
        if outfolder+"/"+handle.split("/")[-1].replace(".edit",".clustS.gz") not in glob.glob(outfolder+"/*"):
            work_queue.put([UCLUST,outfolder,handle,wclust,mindepth,
                            Parallel,muscle,datatype,fileno, w1, w2, 
                            WORK, minuniq, MASK, threads, remake])
            submitted[handle] = 1
            fileno += 1
        else:
            print "\tskipping "+handle.split("/")[-1].replace(".edit",".clustS.gz")+\
                  ' already exists in '+WORK+outfolder.split("/")[-1]

    " create a queue to pass to workers to store the results"
    jobs = []
    for i in range( min(submitted,Parallel) ):
        worker = Worker(work_queue, result_queue, final)
        jobs.append(worker)
        worker.start()
    for j in jobs:
        j.join()

    " output statistics on depth of coverage"
    outstats = open(WORK+"stats/s3.clusters.txt",'a')
    print >>outstats, '\n'+'\t'.join(['taxa','total','dpt.me',
                                      'dpt.sd','d>'+`mindepth-1`+'.tot',
                                      'd>'+`mindepth-1`+'.me',
                                      'd>'+`mindepth-1`+'.sd',
                                      'badpairs'])

    RES = []
    HISTO = []
    #for ff in glob.glob(outfolder+"/.temp.*"):
    for ff in FS:
        end = ff.split("/")[-1].replace(".edit","")
        ff = outfolder+"/.temp."+end
        if os.path.exists(ff):
            line = open(ff).readlines()
            RES.append(line[0].strip().split("\t"))
            HISTO.append([line[0].split("\t")[0],"".join(line[1:])])
            os.remove(ff)
    RES.sort(key=lambda x:x[0])
    HISTO.sort(key=lambda x:x[0])
    
    for i in RES:
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

    for i in HISTO:
        print >>outstats, "sample: "+i[0]+"\n"+i[1]
    
    
    outstats.close()
    for handle in FS:
        nothere = 0
        try: submitted[handle]
        except KeyError:
            nothere = 1
        if not nothere:
            if submitted[handle]:
                if os.path.exists(outfolder+"/"+handle.split("/")[-1].replace(".edit",".clust.gz")):
                    cmd = "/bin/rm "+outfolder+"/"+handle.split("/")[-1].replace(".edit",".clust.gz")
                    subprocess.call(cmd, shell=True)


if __name__ == "__main__":
    main()

