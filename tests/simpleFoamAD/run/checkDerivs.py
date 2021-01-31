"""
Dot product tests for forward and reverse mode AD
We also compare matrix-vector products between forward AD and FD
"""

import numpy as np
from numpy import linalg as LA
import sys

nProcs = int(sys.argv[1])
mode = sys.argv[2]

if mode == "state":
    partName = "dRdW"
elif mode == "point":
    partName = "dRdXv"

dot1 = 0.0
dot2 = 0.0

for n in range(nProcs):

    print("Processing %d" % n)

    stateFADVal = open(partName + "Psi_%d_AD_Values.txt" % n, "r")
    stateFADSeed = open(partName + "Psi_%d_AD_Seeds.txt" % n, "r")
    stateFDVal = open(partName + "Psi_%d_FD_Values.txt" % n, "r")
    stateRADVal = open(partName + "TPsi_%d_AD_Values.txt" % n, "r")
    stateRADSeed = open(partName + "TPsi_%d_AD_Seeds.txt" % n, "r")

    stateFADValLines = stateFADVal.readlines()
    stateFADSeedLines = stateFADSeed.readlines()
    stateFDValLines = stateFDVal.readlines()
    stateRADValLines = stateRADVal.readlines()
    stateRADSeedLines = stateRADSeed.readlines()

    stateFADVal.close()
    stateFADSeed.close()
    stateFDVal.close()
    stateRADVal.close()
    stateRADSeed.close()

    stateFADValList = []
    stateFADSeedList = []
    stateFDValList = []
    stateRADValList = []
    stateRADSeedList = []
    for idx in range(len(stateFADValLines)):
        stateFADValList.append(float(stateFADValLines[idx]))
        stateFDValList.append(float(stateFDValLines[idx]))
        stateRADSeedList.append(float(stateRADSeedLines[idx]))

    for idx in range(len(stateFADSeedLines)):
        stateFADSeedList.append(float(stateFADSeedLines[idx]))
        stateRADValList.append(float(stateRADValLines[idx]))

    stateFADValList = np.asarray(stateFADValList)
    stateFADSeedList = np.asarray(stateFADSeedList)
    stateFDValList = np.asarray(stateFDValList)
    stateRADValList = np.asarray(stateRADValList)
    stateRADSeedList = np.asarray(stateRADSeedList)

    diff = stateFADValList - stateFDValList

    dot1 += np.dot(stateFDValList, stateRADSeedList)
    dot2 += np.dot(stateRADValList, stateFADSeedList)

    print("Max diff in dRdW*Psi: ", diff.max())
    print("Diff norm in dRdW*Psi: ", LA.norm(diff))

print("Dot-product 1: ", dot1)
print("Dot-product 2: ", dot2)

if (dot1-dot2)/dot1 > 0.0001:
    print("\n**********************************************")
    print("Test Failed!!!!!! Relative error > 0.0001")
    print("**********************************************")
    exit(1)
else:
    print("\n**********************************************")
    print("Test Passed!")
    print("**********************************************")

