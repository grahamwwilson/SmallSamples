import numpy as np
import math
import random
from scipy.special import comb
import statistics as st
from itertools import combinations
from ROOT import TFile,TH1D
from histUtils import histInfo
import analyze
import myPythonCheck

myPythonCheck.Check()                         # Enforce use of python3

# Check numpy version
print("numpy version",np.__version__)

f   = TFile("Permutations-New.root", "recreate")
hT1 = TH1D("hT1","T1; T1",100,0.0,2.0)
hT2 = TH1D("hT2","T2; T2",100,0.0,2.0)
hT3 = TH1D("hT3","T3; T3",150,0.0,7.5)

x, y, z = analyze.GetData()                   # Read data from file or generate our own
T1Obs, T2Obs, T3Obs = analyze.Analyze(x,y,z)  # Statistical analysis with scipy 
                                              # returning observed statistics for further 
                                              # stand-alone permutation tests

# Conduct our own permutation tests

nx = len(x)
ny = len(y)
N  = len(x) + len(y)

T1List = []
T2List = []
T3List = []
n1 = 0
n2 = 0
n3 = 0
ntot = comb(N,nx,exact=True)
EPS = 1.0e-12                  # Allow for machine precision issues 

debug = True

icombo=0
for i in combinations(z,nx):   # Enumerate all combinations of partitions into two groups of size nx and ny
    icombo +=1
    ilist = list(i)            # make the permuted version for x
    jlist = z.copy()           # make the complementary permuted list for y
    for element in ilist:
        if element in jlist:
            jlist.remove(element)
    T1 = st.mean(ilist)   - st.mean(jlist)
    T1List.append(T1)
# note we use abs below to do 2-sided tests    
    hT1.Fill(abs(T1))
    T2 = st.median(ilist) - st.median(jlist)
    T2List.append(T2)
    hT2.Fill(abs(T2))
    T3 = T1/math.sqrt( (st.variance(ilist)/nx) + (st.variance(jlist)/ny) )
    T3List.append(T3)
    hT3.Fill(abs(T3))
    if abs(T1) >= abs(T1Obs)-EPS:
       n1 += 1
       if debug:
          print('Significant T1 ',T1,'for icombo ',icombo,ilist,jlist)             
    if abs(T2) >= abs(T2Obs)-EPS: n2 += 1
    if abs(T3) >= abs(T3Obs)-EPS: 
       n3 += 1
       if debug:
          print('Significant T3 ',T3,'for icombo ',icombo,ilist,jlist)       

# Range checks
print('min values ',min(T1List),min(T2List),min(T3List))
print('max values ',max(T1List),max(T2List),max(T3List)) 
print('length     ',len(T1List),len(T2List),len(T3List)) 
print('p-value for T1 ',n1/ntot,n1,ntot)
print('p-value for T2 ',n2/ntot,n2,ntot)
print('p-value for T3 ',n3/ntot,n3,ntot)

histList = [ hT1, hT2, hT3 ]
for h in histList:
    histInfo(h)

f.Write()
