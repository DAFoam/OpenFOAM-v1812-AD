import numpy as np
import sys

file1=sys.argv[1]
file2=sys.argv[2]

f1=open(file1,"r")
f2=open(file2,"r")
lines1=f1.readlines()
lines2=f2.readlines()

val1=[]
val2=[]

for line in lines1:
    cols=line.split()
    val1.append(float(cols[3]))

for line in lines2:
    cols=line.split()
    val2.append(float(cols[3]))

val1=np.asarray(val1)
val2=np.asarray(val2)

diff = abs(val1-val2)
maxIdx = np.argmax(diff)
print ("max:", max(diff), " maxIdx: ", maxIdx, " val1: ", val1[maxIdx], " val2: ", val2[maxIdx])

print ("mean:", sum(diff)/len(diff))
