import sys
import os


finp = os.popen(r"fgrep  '>' " + sys.argv[1] + r" | nl -n rn -s : | cut -f 1 -d_ | sed -e 's/ *\(.*\):>\(.*\)/\2 \1/' | sort ")
last = ""
c_id = -1
for line in finp:
    (cluster,index) = line.rstrip("\n").split()
    index = int(index)-1
    if cluster != c_id:
        print last
        c_id=cluster
        last="%d."%index
    else:
        last="%d "%index+last
print last
