#!/usr/bin/python

import re
import sys
from optparse import OptionParser

usage = "usage: %prog [option] filename"
parser = OptionParser(usage=usage)



parser.add_option("-o","--outfile",dest="out",
                  default = "",
                  metavar = "FILENAME",
                  help = " baseout put [default: %default]")


parser.add_option("-s", "--subset", dest="subsetf",
                  default = "",
                  metavar = "FILE",
                  help = "subsetfilename [default: %default") 

parser.add_option("-S", "--skipsafe", dest="skipsafe",
                  action  = "store_true",
                  default = False,
                  help = "skip sequences that don't appear in file") 


parser.add_option("-i", "--index", dest="index",
                  action = "store_true",
                  default = False,
                  help  = "print sequence index rather than accid")

parser.add_option("-R", "--report", dest="report",
                  action  = "store_true",
                  default = False,
                  help = "report sequences that don't appear inf ile") 


(options,args) = parser.parse_args()


if options.index and not options.subsetf:
   parser.error("Must select --subsetf option if using --index")

if not options.subsetf and (options.report or options.skipsafe):
   parser.error("You must use the --subset option to use skip or report")

if options.report and options.skipsafe:
    parser.error("You can't use both report and skipsafe")


def getid(line):
    mm = re.search(">gi\|[^|]+\|gb\|([^|]+).*",line)
    if mm: return m.group(1)
    mm = re.search(">gi\|([^|]+)",line)
    if mm: return m.group(1)    
    mm = re.search(">([^| \t]+)",line)
    return mm.group(1).rstrip("\n")

if options.out:
    outf = open(options.out,"w")
else:
    outf = sys.stdout


id2num = {}
num2id = {}
lineseq= {}

num=0
linenum=1


def printlist(outf, convert, accs):
   if len(accs) == 0: return
   acc =convert[accs[0]]
   outf.write("%s"%acc)
   for acc in accs[1:]:
      acc = convert[acc]
      outf.write(" %s"%acc)
   outf.write(".\n")

for line in open(args[0]):
    if line[0] == ">":
        id = getid(line)
        if id2num.has_key(id):
            sys.exit("<%s> is a duplicate ID in the file"%id)
        id2num[id]=num
        num2id[str(num)]=id
        lineseq[num]=linenum
        num = num+1


all_seen = {}


if options.subsetf:
    tot=0
    if options.index:
       convert = id2num
    else:
       convert = num2id
    if "asis:" in options.subsetf:
       source = [options.subsetf[5:]]
    else:
       source = open(options.subsetf)
    for line in source:
        line = line.rstrip("\n.")
        accids = re.split("[, ;]",line)
        tot += len(accids)
        valid   = filter(convert.has_key, accids)
        for acc in accids: all_seen[acc] = True
        missing = filter(lambda x: not (convert.has_key(x)), accids)
        if options.report:
           printlist(outf,convertmissing)
        elif len(missing) != 0 and not options.skipsafe:
           sys.exit("From line %s: There are no IDs <%s>\n"%(line,str(missing)))           
        else:
           printlist(outf,convert,valid)
    if options.report: sys.exit(0)
    sys.exit(0)
        
    if "asis:" not in options.subsetf:
       for k in id2num.keys():
          if not all_seen.has_key(k): 
             if options.index:
                outf.write("%s.\n"%id2num[k])
             else:
                outf.write("%s.\n"%k)
       outf.close()
    sys.exit(0)

for i in range(num):
   outf.write("%d %10s %-5d\n"%(i,num2id[i],lineseq[i]))

outf.close()

    



