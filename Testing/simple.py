import numpy as np
import math
import random
from pylab import genfromtxt
from scipy.special import comb
from scipy import stats
import statistics as st
import myPythonCheck

myPythonCheck.Check()                         # Enforce use of python3

# Check numpy version
print("numpy version",np.__version__)

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
        xi = round(random.gauss(-1.5, 1.0),3)
        x.append(xi)
        z.append(xi)
    random.seed(401)
    for i in range(11):
        yi = round(random.gauss( 1.5, 1.0),3)
        y.append(yi)
        z.append(yi)

print('x:',x)
print('y:',y)
print('z:',z)
print('x, N:',len(x),'mean:',st.mean(x),'median:',st.median(x),'var:',st.variance(x))
print('y, N:',len(y),'mean:',st.mean(y),'median:',st.median(y),'var:',st.variance(y))
print('z, N:',len(z),'mean:',st.mean(z),'median:',st.median(z),'var:',st.variance(z))
zcopy = z.copy()
zcopy.sort()
print('z sorted ',zcopy)
nx = len(x); ny = len(y); N = nx+ny
print('N,nx,ny:',N,nx,ny)
print('Combinations (N_choose_nx):',comb(N,nx,exact=True))
Nperms = comb(N,nx,exact=True)

# First let's apply classic tests. Both implemented as 2-sided.
print('Students t-test: ',stats.ttest_ind(x,y))                    # Student's t-test
print('Welchs t-test:   ',stats.ttest_ind(x,y,equal_var=False))    # Welch's t-test
print('Wilcoxon rank-sum test:  ',stats.ranksums(x,y))             # Wilcoxon rank-sum test
print('Mann-Whitney U-test:  ',stats.mannwhitneyu(x,y,method='asymptotic'))             # Mann-Whitney U-test

# Redo these as exact tests using all combinations of potential partitions
dp=1/Nperms
print('Single occurrence p-value quantization = ',dp)
print('Exact version of Students t-test: ',stats.ttest_ind(x,y,permutations=Nperms))        # Exact Student's t-test
print('Exact version of Welchs t-test  : ',stats.ttest_ind(x,y,equal_var=False,permutations=Nperms))  # Exact Welch's t-test
print('Exact version of Mann-Whitney U-test:  ',stats.mannwhitneyu(x,y,method='exact'))             # Exact Mann-Whitney U-test

# Calculate the 3 observed statistics
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
