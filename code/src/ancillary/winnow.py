#! /usr/bin/python

import re
import sys
import os
import shutil
from optparse import OptionParser
from datetime import date


fix_file = "/tmp/wcd.fix"


def InitStuff():
    usage = "usage: winnow.py [option] SEQUENCEFILE "
    parser = OptionParser(usage=usage)

    parser.add_option("-s","--start",dest="start",
                  default = "400",
                  metavar = "INT",
                  help = " start filter point [default: %default]")

    parser.add_option("-d","--delta",dest="delta",
                  metavar = "INT",
                  default = "50",
                  help = "step to reduce filter point")

    parser.add_option("-f","--final",dest="final",default=90,
                  metavar = "INT",
                  help = " final filter point  [default %default]")

    parser.add_option("-o","--output",dest="filename",default="",
                  metavar = "FILE",
                  help = " name of output file")

    parser.add_option("-w","--wcd",dest="wcd",default="/usr/local/bin/wcd",
                  metavar = "WCD",
                  help = " wcd binary ap")

    parser.add_option("-c","--cluster_threshold",dest="cthresh",default=1000,
                  metavar = "INT",
                  help = "threshold for cluster size to be CAP-ed")

    parser.add_option("-S","--SCHEDULE",dest="schedule",default="",
                  metavar = "INT,INT,..",
                  help = "which steps to take")


    parser.add_option("-m","--maxclustersize",dest="maxcluster",
                  default = "10000",
                  metavar = "INT",
                  help = " maximim cluster size [default: %default]")


    (options,args) = parser.parse_args()
    print options
    print args
    open(fix_file,"w").close()
    return (options,args)

(options,args)=InitStuff()


delta    = int(options.delta)
seqfname = args[0]
start = int(options.start)
stop  = int(options.final)
maxcluster = int(options.maxcluster)
cthresh    = int(options.cthresh)

if options.filename:
    outf = open(options.filename,"w")
else:
    outf = sys.stdout


if not os.path.exists(seqfname):
    sys.exit("File <%s> does not exist!"%seqfname);

if not (os.path.exists(seqfname+".nlc") and os.path.exists(seqfname+".nlc.ary")):
    sys.exit("File <%s> nlc or ary does not exist!"%seqfname);


if options.schedule:
    schedule = map(int,options.schedule.split(","))
else:
    schedule = range(start,stop-1,-delta)

schedule=schedule+[-1]

if stop > start:
    sys.exit("Stop is bigger than start")




suffix_cmd = options.wcd + (" -j %s "%fix_file)+" -c -F suffix -w %d %s -o %s"


def execute(cmd):
    print "  "+cmd[0:40]
    os.system(cmd)


def cluster_table_process(curr_clt_file,big_clusters):
    f = open(curr_clt_file)
    big=0
    for line in f:
        line = line.rstrip(".\n")
        nums = map(int, line.split(" "))
        size = len(nums)
        if size>big: big=size
        if size >= cthresh:
            big_clusters.append(nums)
    f.close()
    return big


def dump_clusters(cap_cluster):
    for cluster in cap_clusters:
       outf.write("%d"%cluster[0])
       for num in cluster[1:]:
           outf.write(" %d"%num)
       outf.write(".\n")



remove = []


cap_clusters = []

i=0


print "Schedule: ",schedule
curr = schedule.pop(0)
prev = curr*2

while curr != -1:
   print "Schedule: ",schedule
   maxlen = maxcluster+1
   while maxlen > maxcluster:
       print "Doing round %d. Current filter is %d. Delta is %d"%(i,curr,delta)
       curr_clt_file = "clt%04d.clt"%(curr)
       cmd = suffix_cmd%(curr,seqfname,curr_clt_file) 
       execute(cmd)
       big_clusters = []
       maxlen = cluster_table_process(curr_clt_file,big_clusters)
       print "  Biggest cluster is %d.  There are %d cap-clusters."%\
                 (maxlen,len(big_clusters))
       if maxlen >= maxcluster and curr<prev-1:
            print "  Backing up....."
	    schedule.insert(0,curr)
            curr = (curr+prev)/2
            os.remove(curr_clt_file)
       elif len(big_clusters)>0:
           print "  Updating fix file. "
           cap_clusters = cap_clusters+big_clusters
           fixf = open(fix_file,"a")
           fixf.write("fix")
           for cluster in big_clusters:
               for num in cluster:
                   fixf.write(" %s"%num)
           fixf.write(".\n")
           fixf.close()
           prev = curr
           curr = schedule.pop(0)
       else:
           print "  No big clusters this round!"
           prev=curr
           curr = schedule.pop(0)
           delta = int(options.delta)
       print
       i=i+1


dump_clusters(cap_clusters)

     
outf.close()           

