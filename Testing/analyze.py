# analyze.py
#
# Put the common code here for easier maintenance

import math
from scipy.special import comb
from scipy import stats
import statistics as st
import numpy as np
import random
from pylab import genfromtxt

def GetData(fromFile = False):

    x = []; y = []; z = []

    if fromFile==True:
# Read data from two separate files.
        xnp = genfromtxt('Cvalues.dat',usecols=0)
        ynp = genfromtxt('Mvalues.dat',usecols=0)
        znp = np.append(xnp,ynp)    
        x = xnp.tolist()
        y = ynp.tolist()
        z = znp.tolist()
    else:
# Generate appropriate random samples    
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
    
    return x,y,z

def Analyze(x,y,z):

    print('x:',x)
    print('y:',y)
    print('z:',z)
# Round these a little when printing   
    print('x, nx:',len(x),'mean:',round(st.mean(x),6),' +- ',round(math.sqrt(st.variance(x)/len(x)),6),
          'median:',round(st.median(x),6),'sx:',round(math.sqrt(st.variance(x)),6))
    print('y, ny:',len(y),'mean:',round(st.mean(y),6),' +- ',round(math.sqrt(st.variance(y)/len(y)),6),
          'median:',round(st.median(y),6),'sy:',round(math.sqrt(st.variance(y)),6))          
    print('z,  N:',len(z),'mean:',round(st.mean(z),6),' +- ',round(math.sqrt(st.variance(z)/len(z)),6),
          'median:',round(st.median(z),6),'sz:',round(math.sqrt(st.variance(z)),6))
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
    print('med(x)  -  med(y)  T2Obs:',T2Obs)
    print('Welch t            T3Obs:',T3Obs)

# Now let's look more carefully at these definitions.
# https://en.wikipedia.org/wiki/Student%27s_t-test 
# They differ for sure when nx != ny. 

    spB = math.sqrt( ((nx-1)*st.variance(x) + (ny-1)*st.variance(y))/(nx+ny-2) )
    pooledsE = spB*math.sqrt((1/nx) + (1/ny))   # on the difference
    sDelta = math.sqrt( ( st.variance(x)/nx ) + ( st.variance(y)/ny ) )

# Equal variance version (Student's t-test)
    t1 = (st.mean(x) - st.mean(y)) / pooledsE
    print(' ')
    print('nx,ny,spB,pooledsE,Student t1:',nx,ny,spB,pooledsE,t1)
    
# Maybe use code like
# https://stats.stackexchange.com/questions/475289/confidence-interval-for-2-sample-t-test-with-scipy 
# for confidence intervals

    pt = 2.0*stats.t.cdf(-abs(t1), nx+ny-2)
    print('2-sided equal variance t-test p-value ',pt,' (dof = ',nx+ny-2,')')
    
    lbound1 = (st.mean(x) - st.mean(y)) - stats.t.ppf(0.975, nx+ny-2)*pooledsE
    ubound1 = (st.mean(x) - st.mean(y)) + stats.t.ppf(0.975, nx+ny-2)*pooledsE
    print('95% CL interval on x - y : [',lbound1,',',ubound1,']')
    
# Unequal variance version (Welch's t-test)
    t2 = (st.mean(x) - st.mean(y)) / sDelta
    print(' ')
    print('nx,ny,sDelta,Welch t2:',nx,ny,sDelta,t2)
    vx = st.variance(x)
    vy = st.variance(y)
    df = ( (vx / nx) + (vy / ny))**2 / ( (vx**2 / (nx**2 * (nx - 1)) + vy**2 / (ny**2 * (ny - 1)))) 
    ptw = 2.0*stats.t.cdf(-abs(t2), df)
    print('2-sided non-equal variance t-test p-value ',ptw,' (dof =',df,')')
    lbound2 = (st.mean(x) - st.mean(y)) - stats.t.ppf(0.975, df)*sDelta
    ubound2 = (st.mean(x) - st.mean(y)) + stats.t.ppf(0.975, df)*sDelta
    print('95% CL interval on x - y : [',lbound2,',',ubound2,']')    
    
# Note I have checked these against an R implementation, they seem to be sound.          
    
    return T1Obs, T2Obs, T3Obs
