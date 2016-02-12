

'''
Takes two arguments --
1. A list of sequence indices
2. A cluster file

Produces the corresponding cluster file with all the indices in (1) removed.
The file 1 is a list of integers, separated by spaces and new lines (no semantic difference). Full stops (.) at the end of lines are ignored
'''

import sys

thinf = sys.argv[1]
clustf= sys.argv[2]

remove = []

for line in open(thinf):
    nums=line.rstrip("\n.").split()
    remove = remove+nums



for line in open(clustf):
    nums=line.rstrip("\n.").split()
    curr = []
    for n in nums:
        if n not in remove:
            curr.append(n)
    if len(curr)==0: continue
    sys.stdout.write("%s"%curr[0])
    for n in curr[1:]:
        sys.stdout.write(" %s"%n)
    sys.stdout.write("\n")
