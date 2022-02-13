import numpy as np
import math
import random
from pylab import genfromtxt
from scipy.special import comb
from scipy import stats
import statistics as st
from itertools import combinations
from ROOT import TFile,TH1D
from histUtils import histInfo
import myPythonCheck

myPythonCheck.Check()                         # Enforce use of python3

# Check numpy version
print("numpy version",np.__version__)

f   = TFile("Permutations-New.root", "recreate")
hT1 = TH1D("hT1","T1; T1",100,0.0,2.0)
hT2 = TH1D("hT2","T2; T2",100,0.0,2.0)
hT3 = TH1D("hT3","T3; T3",150,0.0,7.5)

# Read data from two separate files.
xnp = genfromtxt('Cvalues.dat',usecols=0)
ynp = genfromtxt('Mvalues.dat',usecols=0)
znp = np.append(xnp,ynp)

fromfile = False

x = []; y = []; z = []
if fromfile==True:
    x = xnp.tolist()
    y = ynp.tolist()
    z = znp.tolist()
else:
    random.seed(400)
    for i in range(9):
        xi = round(random.gauss(-1.5,1.0),3)
        x.append(xi)
        z.append(xi)
    random.seed(401)
    for i in range(11):
        yi = round(random.gauss( 1.5,1.0),3)
        y.append(yi)
        z.append(yi)

print('x:',x)
print('y:',y)
print('z:',z)
print('x, N:',len(x),'mean:',st.mean(x),'median:',st.median(x),'var:',st.variance(x))
print('y, N:',len(y),'mean:',st.mean(y),'median:',st.median(y),'var:',st.variance(y))
print('z, N:',len(z),'mean:',st.mean(z),'median:',st.median(z),'var:',st.variance(z))

nx = len(x); ny = len(y); N = nx+ny
print('N,nx,ny:',N,nx,ny)
print('Combinations (N_choose_nx):',comb(N,nx,exact=True))
print('Permutations (N!):',math.factorial(N))
print('Permutations (nx!):',math.factorial(nx))
print('Permutations (ny!):',math.factorial(ny))
print('Note combinations = N!/(nx! ny!) ')
Nperms = comb(N,nx,exact=True)

# First let's apply classic tests. Both implemented as 2-sided.
print('Students t-test: ',stats.ttest_ind(x,y))                    # Student's t-test
print('Welchs t-test:   ',stats.ttest_ind(x,y,equal_var=False))    # Welch's t-test
print('Wilcoxon rank-sum test:  ',stats.ranksums(x,y))             # Wilcoxon rank-sum test
print('Mann-Whitney U-test:  ',stats.mannwhitneyu(x,y,method='asymptotic'))             # Mann-Whitney U-test

dp=1/Nperms           # Should really be 2/Nperms for double sided though
print('Single occurrence p-value quantization = ',dp)
# Redo these as exact tests using all combinations of potential partitions
print('Exact version of Students t-test: ',stats.ttest_ind(x,y,permutations=Nperms))        # Exact Student's t-test
print('Exact version of Welchs t-test  : ',stats.ttest_ind(x,y,equal_var=False,permutations=Nperms))  # Exact Welch's t-test
print('Exact version of Mann-Whitney U-test:  ',stats.mannwhitneyu(x,y,method='exact'))             # Exact Mann-Whitney U-test

# Calculate the 3 observed statistics that we will use in our own exact permutation test
T1Obs = st.mean(x)   - st.mean(y)
T2Obs = st.median(x) - st.median(y)
T3Obs = T1Obs/math.sqrt( (st.variance(x)/len(x)) + (st.variance(y)/len(y)) )

print('mean(x) - mean(y)  T1Obs:',T1Obs)
print(' med(x) -  med(y)  T2Obs:',T2Obs)
print(' Welch t           T3Obs:',T3Obs)

# Now let's look more carefully at these definitions.
# https://en.wikipedia.org/wiki/Student%27s_t-test
# They differ for sure when nx != ny.

spB = math.sqrt( ((nx-1)*st.variance(x) + (ny-1)*st.variance(y))/(nx+ny-2) )
sDelta = math.sqrt( ( st.variance(x)/nx ) + ( st.variance(y)/ny ) )

# Equal variance version (Student's t-test)
t1 = (st.mean(x) - st.mean(y)) / (spB*math.sqrt((1/nx) + (1/ny)))
print('nx,ny,spB,Student t1:',nx,ny,spB,t1)

# Unequal variance version (Welch's t-test)
t2 = (st.mean(x) - st.mean(y)) / sDelta
print('nx,ny,sDelta,Welch t2:',nx,ny,sDelta,t2)

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
