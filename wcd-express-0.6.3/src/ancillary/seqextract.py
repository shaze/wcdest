#!/sw/bin/python

import re
import sys
from optparse import OptionParser
import os


unit = { 'Mega': 1048576, 'Kilo':1024, 'bytes':1, 'Giga':1073741824}

usage = "usage: %prog [option] filename"
parser = OptionParser(usage=usage)



parser.add_option("-b","--breakup",dest="breakup",
	          action  = "store",
                  default = "",
                  metavar = "FILE",
                  help = " break up the file")


parser.add_option("-u","--unit",dest="unit",
                  default = "Mega",
                  metavar = "UNIT",
                  help = " unit (Mega, Kilo, bytes, Giga [bytes] [default: %default]")


parser.add_option("-t","--target",dest="target",
                  default = "",
                  metavar = "INT",
                  help = " target (seq num) [default: %default]")

parser.add_option("-f","--first",dest="first",
                  default = "",
                  metavar = "INT",
                  help = " first (seq num) [default: %default]")


parser.add_option("-l","--last",dest="last",
                  default = "",
                  metavar = "INT",
                  help = " last (number blocks) [default: %default]")



parser.add_option("-o","--outfile",dest="out",
                  default = "",
                  metavar = "FILENAME",
                  help = " baseout put [default: %default]")


parser.add_option("-s", "--subset", dest="subsetf",
                  default = "",
                  metavar = "FILE",
                  help = "subsetfilename [default: %default") 


parser.add_option("-S", "--split", dest="split", default="",
                  metavar = "FILENAME",
                  help  = "split file prefix")


(options,args) = parser.parse_args()


if options.breakup and not options.subsetf:
   parser.error("To use breakup you must use --subsetf")

if options.unit not in unit.keys():
    sys.exit("Unknown unit <%s>: Choose bytes, Kilo, Mega, Giga"%(options.unit))
size = unit[options.unit]

curr_hdr = ""
curr_seq = ""
bytes_used=0
first_seq = 0
big_num   = 1000*sys.maxint # some arb very big number
last_seq  = big_num 


if options.out:
    if options.split: sys.exit("Can't have out and split selected")
    outf = open(options.out,"w")
elif options.split:
   outf = open(options.split+"inc.fa","w")
   outf2= open(options.split+"exc.fa","w")
else:
    outf = sys.stdout

if options.first:    first_seq = int(options.first)
if options.last:    last_seq = int(options.last)


def output_cluster(curr_cluster, x):
    base = os.path.join(options.breakup, os.path.basename(options.breakup))
    sub = base + "%s%02d.fa"%(options.unit[0],x)
    curr_cluster.sort()
    clf = open(sub,"w")
    for n in curr_cluster:
       clf.write(seqcont[n])
    clf.close()


def breakUp():
    f = open(options.subsetf)
    x = curr_size = 0
    curr_cluster = []
    for line in f:
        line = line.rstrip(".\n ")
        line = line.rstrip(" ")
        nums = map(int,re.split("[, ;]",line))
	for n in nums:
            curr_size = curr_size + seqlen[n]
        curr_cluster = curr_cluster+nums
        if curr_size > target:
	   output_cluster(curr_cluster, x)
           curr_cluster = []
           x = x+1
           curr_size = 0
    if curr_size > target:
       output_cluster(curr_cluster, x)



def readFile(seqcont, seqlen):
    curr_hdr = curr_seq = ""
    print args[0]
    f = open(args[0])
    curr_hdr = f.readline()
    curr_seq = ""
    for line in f:
    	if line[0] == ">":
           seqcont.append(curr_hdr+curr_seq)
	   seqlen.append(len(curr_seq))
           curr_hdr = line
           curr_seq = ""
        else:
           curr_seq = curr_seq+line
    curr_seq=curr_seq+line
    seqcont.append(curr_hdr+curr_seq)
    seqlen.append(len(curr_seq))
    f.close()


if options.target:
    target = int(options.target)*size
else:
    target = big_num

    

if options.breakup: 
   seqcont = []
   seqlen = []
   readFile(seqcont, seqlen)
   breakUp()
   sys.exit(0)

if options.subsetf:
    subset = []
    first_seq = big_num
    last_seq = 0
    if "asis:" in options.subsetf:
       source = [options.subsetf[5:]]
    else:
       source = open(options.subsetf)
    for line in source:
        line = line.rstrip("\n.")
        nums = map(int,re.split("[, ;]",line))
        subset = subset + nums
        for n in nums:
            if n>last_seq: last_seq=n
            if n<first_seq: first_seq=n

oldlen=0



num=0
for line in open(args[0]):
    if line[0] == ">":
        if num > first_seq:
            if (not options.subsetf or num-1 in subset) :
               outf.write(curr_hdr+curr_seq)
            if (options.split and num-1 not in subset):
               outf2.write(curr_hdr+curr_seq)
        curr_hdr=line
        curr_seq=""
        if bytes_used >= target: break
        if num > last_seq and not options.split: break
        num = num+1
    else:
        curr_seq=curr_seq+line
        if num >= first_seq: bytes_used=bytes_used+len(line)-1
else:
    if last_seq+1 >= num >= first_seq:
        outf.write(curr_hdr+curr_seq)

outf.close()
if options.split: outf2.close()
