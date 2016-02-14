#! /usr/bin/python

from string import rstrip
import sys
import re
import random
from optparse import OptionParser
import string


compl = {"A":"T", "C":"G", "G":"C", "T":"A","N":"N","V":"X","W":"X","K":"X","S":"X","R":"X","Y":"X","B":"X","M":"X","D":"X","H":"X","Q":"X","-":"",
         "a":"T", "c":"G", "g":"C", "t":"a","n":"N", "x":"X", "X":"X"}

parser = OptionParser()
parser.add_option("-p","--nopos",action="store_true",dest = "nopos")
parser.add_option("-f","--fasta",action="store_true",dest = "fasta")
parser.add_option("-r","--rand", action="store",dest="random",default="0")
parser.add_option("-n","--norc",action="store_true",dest = "norc")
parser.add_option("-x","--Xaway",action="store_true",dest = "xaway")
parser.add_option("-d","--dust",action="store",dest = "dust")
parser.add_option("-o","--outfile",action="store",dest = "fname")


alt = {}
alt['A']="CGT"
alt['C']="AGT"
alt['T']="ACG"
alt['G']="ACT"
alt['a']="cgt"
alt['c']="agt"
alt['t']="acg"
alt['g']="act"

def randomify(rawseq):
    if options.random==0: return rawseq
    rawseq = list(rawseq)
    errors = float(len(rawseq))*options.random/100
    low = int(errors)
    hi  = low+1
    rnd = random.randrange(low,hi)
    if rnd < errors: 
      trg = low
    else:
      trg= hi
    sample = random.sample(range(len(rawseq)),trg)
    for pos in sample:
      rawseq[pos] = random.choice(alt[rawseq[pos]])
    return string.join(rawseq,"")
      
      



def pseq(id,seq):
  if options.fasta: 
     print id
     if len(seq)<1: seq="NNNN"
     prefix=""
  else:
     prefix=str(n)+" "
  if options.dust:
     seq = robj.sub("",seq)
  rseq = ""
  tmp=seq
  if options.fasta:
      tmp = ""
      tseq = randomify(seq)
      for i in range(len(seq)/80+1):
        tmp = tmp+tseq[i*80:min((i+1)*80,len(seq))]
        if (i+1)*80<len(seq): tmp=tmp+"\n"
  if options.nopos:
      pseq = ""
  else:
      pseq = prefix+tmp+"\n"
  if options.norc:
      rseq=""
  else:
    for ch in seq:
      rseq = compl[ch]+rseq
    rseq = prefix+rseq+"\n"
  wfile.write(pseq+rseq)




(options,args) = parser.parse_args()
options.random=float(options.random)

finp = file(args[0])

n = 0

if options.dust:
    repeat = r"(.)\1{"+options.dust+",}"
    robj = re.compile(repeat)
    options.dust = int(options.dust)


clean = re.compile(r"[^>aACcGgTt]")


if options.fname:
    wfile = open(options.fname,"w")
else:
    wfile = sys.stdout

line = finp.readline()
id = line.rstrip("\n")
line = finp.readline()
seq = ""


cseq =0
while line:
    while 0 < len(line) <5 and (len(line)==0 or line[0] != ">") : line=finp.readline()
    line=rstrip(line)
    if line:
       if line[0]==">":
           pseq(id, seq)
	   id = line
           seq = ""
           n=n+1
       else:
         if options.xaway:
           line = clean.sub("",line)
         seq = seq+line
    line = finp.readline()




pseq(id, seq)
wfile.close()

