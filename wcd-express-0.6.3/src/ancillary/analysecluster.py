#!/usr/bin/python


from optparse import OptionParser

import sys

parser = OptionParser()

parser.add_option("-x","--extract",action="store", dest="theta")
(options,args)= parser.parse_args();

def show_cluster(inds):
    inds = map(int,inds)
    inds.sort()
    for i in inds: print i,
    print "."


maxclustersize = 20000000

theta=0
if options.theta: theta = int(options.theta)




cluster = [0]*maxclustersize

f = open(args[0])

for line in f:
    line=line.rstrip()
    line=line.rstrip(".")
    inds = line.split()
    size = len(inds)
    if size==0:
       print line,inds
    cluster[size]=cluster[size]+1
    if theta !=0 and size > theta: show_cluster(inds)

if theta == 0:
    for size in range(maxclustersize):
        if cluster[size]:
            print ("%d\t%d")%(size,cluster[size])

        
