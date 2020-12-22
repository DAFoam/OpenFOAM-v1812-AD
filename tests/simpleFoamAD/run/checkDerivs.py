import numpy as np
from numpy import linalg as LA

stateFADVal = open("dRdWPsi_0_AD_Values.txt", "r")
stateFADSeed = open("dRdWPsi_0_AD_Seeds.txt", "r")
stateFDVal = open("dRdWPsi_0_FD_Values.txt", "r")
stateRADVal = open("dRdWTPsi_0_AD_Values.txt", "r")
stateRADSeed = open("dRdWTPsi_0_AD_Seeds.txt", "r")

stateFADValLines = stateFADVal.readlines()
stateFADSeedLines = stateFADSeed.readlines()
stateFDValLines = stateFDVal.readlines()
stateRADValLines = stateRADVal.readlines()
stateRADSeedLines = stateRADSeed.readlines()

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

diffADFD = stateFADValList - stateFDValList

dot1 = np.dot(stateFDValList, stateRADSeedList)
dot2 = np.dot(stateRADValList, stateFADSeedList)

print("maxDiff: ", diffADFD.max())
print("normDiff: ", LA.norm(diffADFD))
print("MVF", dot1)
print("MVR", dot2)
