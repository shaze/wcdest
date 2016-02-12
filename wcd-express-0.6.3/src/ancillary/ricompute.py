#!/usr/bin/python
# Scott Hazelhurst (C) 2011

import sys
from optparse import OptionParser
import math


from   string import split
import re



usage = "usage: %prog [option] fname1 fname2"
parser = OptionParser(usage=usage)

parser.add_option("-i","--index",dest="index",
                  default = "all",
                  metavar = "index-name",
                  help = " index type [default: %default]")

parser.add_option("-l","--label",dest="labels",
                  default = "",
                  metavar = "string,*",
                  help = " labels for table [default: %default]")


parser.add_option("-d","--diff",action="store_true",dest="diff",
                  default = False,
                  help = " show differences ")

parser.add_option("-a","--c1notc2only",action="store_true",
                  dest = "c1notc2only",
                  default = False,
                  help = " only diffs of pairs clustered in c1 but not c2")

parser.add_option("-b","--c2notc1only",action="store_true",
                  dest = "c2notc1only",
                  default = False,
                  help = " only diffs of pairs clustered in c2 but not c1")

parser.add_option("-e","--c1max",
                  dest = "c1maxsize",
                  default = "",
                  help = " only show pairs of sequences")

parser.add_option("-f","--c2max",
                  dest = "c2maxsize",
                  default = "",
                  help = " only diffs of pairs clustered in c2 but not c1")


parser.add_option("-t","--tabular",
                  dest = "tabular",
                  default = "",
                  help = "tabular")


parser.add_option("-x","--missing-as-singletons",
                  dest = "singletons",
                  action = "store_true",
                  default = "",
                  help = "missing seq is considered as a singleton")


parser.add_option("-o","--out",
                  dest = "outf",
                  default = "",
                  help = " OUTPUT")

(options,args) = parser.parse_args()

if options.outf:
   outf = open(options.outf,"w")
else:
   outf = sys.stdout


show_c1notc2 = show_c2notc1 = True
if options.c1notc2only and options.c2notc1only:
   parser.error("Doesn't make sense to use c1notc2only and c2notc1only")
if options.c1notc2only: show_c2notc1 = False
if options.c2notc1only: show_c1notc2 = False



if not options.diff and (options.c1maxsize or options.c2maxsize):
   parser.error("You can only set the maxsize if you are using the diff option")



c1maxsize = c2maxsize = sys.maxint
if options.c1maxsize:
   c1maxsize = int(options.c1maxsize)
if options.c2maxsize:
   c2maxsize = int(options.c2maxsize)



ind = {'ji':0, 'ri':1, 'se':2, 'sp':3, 'ppv':4, 'npv':5, 'cc':6 }

def do_nothing(a,b,c,d,e): return (a,b,c,d,e)


def get_indexes(tp,tn,fp,fn,n):
   ji = float(tp)/(tp+fn+fp)
   ri = float(tp+tn)/n
   se = float(tp)/(tp+fn)
   sp = float(tn)/(tn+fp)
   ppv= float(tp)/(tp+fp)
   npv= float(tn)/(fn+tn)
   cc = (tp*tn-fp*fn)/math.sqrt((tp+fp)*(tn+fn)*(tp+fn)*(tn+fp))
   return (ji,ri,se,sp,ppv,npv,cc)

def jaccard(tp,tn,fp,fn,n):
   (ji,ri,se,sp,ppv,npv,cc)=get_indexes(tp,tn,fp,fn,n)
   out = "JI=%7.5f\n"%(ji)
   return out


def se(tp,tn,fp,fn,n):
   (ji,ri,se,sp,ppv,npv,cc)=get_indexes(tp,tn,fp,fn,n)
   out = "SE=%7.5f\n"%(se)
   return out



def all(tp,tn,fp,fn,n):
   (ji,ri,se,sp,ppv,npv,cc)=get_indexes(tp,tn,fp,fn,n)
   out = "JI=%7.5f\nRI=%7.5f\nSE=%7.5f\nSP=%7.5f\nCC=%7.5f\nPPV=%7.5f\nNPV=%7.5f\n"%(ji,ri,se,sp,cc,ppv,npv)
   return out

def compReader(inp,clustering,lencluster):
    """ reads in cluster table from inp and produces a
        dictionary in clustering """
    cnum=0
    max = 0
    data = inp.readline()
    data = data.strip("\n.")
    while len(data) != 0:
        nums = split(data)
        rep  = nums[0]
        for n in nums:
            clustering[n]= rep
	    lencluster[n]= len(nums)
        data = inp.readline()
        data = data.strip("\n.")
	



# Following method works but is very inefficient -- it's n^2.
# indexCompute is the same in the worst case but is O(p^2) where
# p is the size of the largest cluster -- even in bad clusterings
# this will be at least a constant times faster
def cIndex(clustering1, clustering2,index):
    # computes the rand index between clustering1 and clustering2
    # these are 
    n=a=d=0
    for i in clustering1.keys():
        for j in clustering2.keys():
            if i != j:
                n=n+1
                if clustering1[i] == clustering1[j] and clustering2[i]==clustering2[j]: a=a+1
                if clustering1[i] != clustering1[j] and clustering2[i]!=clustering2[j]: d=d+1
    return index(a,d,n)

def get_cluster(data):
    data = data.strip("\n.")
    nums = re.split("[ ,;]",data)
    return nums

def indexCompute(c1,c2,f1,f2,index):
   n=tp = tn = fp = fn = 0
   for line in f1:
      cluster = get_cluster(line)
      n = n+len(cluster)
      for i in range(len(cluster)):
         for j in range(i+1,len(cluster)):
	     a = cluster[i]
	     b = cluster[j]
             if options.singletons and \
               not(c2.has_key(a) and c2.has_key(b)):
                  fp = fp+1
             elif c2[a]==c2[b]: tp = tp+1
             else: fp = fp + 1
   for line in f2:
      cluster = get_cluster(line)
      for i in range(len(cluster)):
         for j in range(i+1,len(cluster)):
	     a = cluster[i]
	     b = cluster[j]
             if options.singletons and not(c1.has_key(a) and c1.has_key(b)):
                fn = fn+1
             elif c1[a]!=c1[b]: 
                fn = fn+1
   for k in c2.keys():
      if not c1.has_key(k): n=n+1
   n=n*(n-1)/2
   tn = n-tp-fn-fp
   outf.write("TP=%d;  TN=%d;  FP=%d;  FN=%d\n"%(tp,tn,fp,fn))
   return index(tp,tn,fp,fn,n)


def diffCompute(f1,size1,c2,size2,c2len):
   tick = {}
   for line in f1:
      cluster = get_cluster(line)
      for i in range(len(cluster)):
	 a = cluster[i]
         for j in range(i+1,len(cluster)):
            b = cluster[j]
            if c2[a]!=c2[b]: 
	       if len(cluster)<=size1 and c2len[a]<=size2 and c2len[b]<=size2:
                     if not tick.has_key(a) and not tick.has_key(b):
		     	outf.write(a+" "+b+"\n")
		     	tick[a]=1
			tick[b]=1



def compareCompute(c1,c2,f1,f2,lenc1,lenc2):
   if show_c1notc2: 
      outf.write("Differences of pairs in c1 but NOT c2\n")
      diffCompute(f1,c1maxsize,c2,c2maxsize,lenc2)
   if show_c2notc1: 
      outf.write("Differences of pairs in c2 but NOT c1\n")
      diffCompute(f2,c2maxsize,c1,c1maxsize,lenc1)


if options.index == "all":
   index = all
elif options.index == "se":
   index = se
elif options.index == "jaccard":
   index = jaccard
else:
   sys.exit("No such index: %s"%options.index)


def InitClusterTables(f, c, lencl):
    for i in range(len(args)):
       f[i] = open(args[i])
       compReader(f[i],c[i],lencl[i])
       f[i].close()


def ProduceComparisonTables():
   table = [ [0]*len(args) for i in range (len(args)) ]
   chosen = options.tabular.split(',')
   lower = ind[chosen[1]]
   upper = ind[chosen[0]]
   symmetric = (lower==upper) and chosen[0] in ["ri","ji","cc"]
   if options.labels:
      labels = options.labels.split(',')
      if len(labels) != len(args):
         sys.exit("Number of labels and input files differs")
   else:
      labels=["      "]*len(args)
   for i in range(len(args)):
      if symmetric:
         jstart = i+1
      else:
         jstart = 0
      for j in range(jstart,len(args)):
         if i==j: continue
         f[i] = open(args[i])
         f[j] = open(args[j])
         all_indices = indexCompute(c[i],c[j],f[i],f[j],get_indexes)
         table[i][j]= all_indices[upper]
         if not symmetric:
            f[i] = open(args[i])
            f[j] = open(args[j])
            all_indices = indexCompute(c[j],c[i],f[j],f[i],get_indexes)
         table[j][i]= all_indices[lower]
   return(table,labels,symmetric)



def ProduceOutput(outf,table,labels,symmetric):
   outf.write("       ")
   for j in range(len(args)):
      outf.write("%8s"%labels[j])
   outf.write("\n")
   for i in range(len(args)):
      outf.write("%8s"%labels[i])
      for j in range(len(args)):
         if i==j or \
            i>j and symmetric: outf.write("        ")
         else:
            outf.write("%8.4f"%table[i][j])
      outf.write("\n")




f     = [0]*len(args)
c     = [ {} for i in range(len(args)) ]
lencl = [ {} for i in range(len(args)) ]


InitClusterTables(f, c, lencl)

if not options.tabular:
   f[0]=open(args[0])
   f[1]=open(args[1])
   outf.write(indexCompute(c[0],c[1],f[0],f[1],index))



if options.diff:
   f[0]=open(args[0])
   f[1]=open(args[1])
   compareCompute(c[0],c[1],f[0],f[1],lencl[0],lencl[1])
   sys.exit(0)


if options.tabular:
   (table,labels,symmetric)=ProduceComparisonTables()
   ProduceOutput(outf,table,labels,symmetric)


outf.close()
