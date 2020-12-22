"""
Dot product tests for forward and reverse mode AD
We also compare matrix-vector products between forward AD and FD
"""

import numpy as np
from numpy import linalg as LA
import sys

nProcs = sys.argv[1]

dot1 = 0.0
dot2 = 0.0

for n in range(nProcs):

    print("Processing %d" % n)

    stateFADVal = open("dRdWPsi_%d_AD_Values.txt" % n, "r")
    stateFADSeed = open("dRdWPsi_%d_AD_Seeds.txt" % n, "r")
    stateFDVal = open("dRdWPsi_%d_FD_Values.txt" % n, "r")
    stateRADVal = open("dRdWTPsi_%d_AD_Values.txt" % n, "r")
    stateRADSeed = open("dRdWTPsi_%d_AD_Seeds.txt" % n, "r")

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
        stateFADSeedList.append(float(stateFADSeedLines[idx]))
        stateFDValList.append(float(stateFDValLines[idx]))
        stateRADValList.append(float(stateRADValLines[idx]))
        stateRADSeedList.append(float(stateRADSeedLines[idx]))

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
