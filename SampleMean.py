# SampleMean.py
#
# Check small sample results and correct for the bias 
# in the estimated uncertainty on the sample mean and 
# the estimate of the standard deviation when the true 
# rms is unknown, but the distribution is assumed to be normal.
#
# Samples of 10 events are drawn from the normal distribution x ~ N(0,1)
#
# When trying to estimate the population parameters from each sample we'd prefer 
# an unbiased estimate.
# The bias correction for the variance is well known.
# Here we also apply the exact correction for the bias in the 
# standard deviation estimate following Gurland-Tripathi eqn 2, and 
# consequently for the estimate of the standard error on the mean (sd/sqrt(N)).
# Reference is: John Gurland and Ram Tripathi, 
# The American Statistician, October 1971, Vol 25. No. 4 pp 30--32.
#
# While this is good for quoting unbiased uncertainties on the mean from small 
# samples, it seems that Student's t test (2-sample one) is the more 
# appropriate statistic for comparing the means of two different 
# samples, although this also assumes the distributions are normal and 
# also that they have the same underlying variance.
#
# Code uses python3 assumptions on integers and will only work correctly 
# if invoked with python3. It throws an exception if called with python 2.
#
#                       Graham W. Wilson
#                            05-FEB-2022
#
import random
import math
from ROOT import TFile, TH1D, TMath
from histUtils import histInfo
import sys

if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required for proper execution")

# Histogram and histogram file initialization
hmean = TH1D("hmean","hmean; Sample mean",250,-2.5,2.5)
hsmean = TH1D("hsmean","hsmean; #sqrt{N} sample mean",300,-7.5,7.5)
hvar1 = TH1D("hvar1","hvar1; Population variance estimate I",300,0.0,6.0)
hvar2 = TH1D("hvar2","hvar2; Population variance estimate II",300,0.0,6.0)
hsd1 = TH1D("hsd1","hsd1; Population standard deviation estimate I",300,0.0,3.0)
hsd2 = TH1D("hsd2","hsd2; Population standard deviation estimate II",300,0.0,3.0)
hsem1 = TH1D("hsem1","hsem1; Sample mean error estimate I",200,0.0,1.0)
hsem2 = TH1D("hsem2","hsem2; Sample mean error estimate II",200,0.0,1.0)

f  = TFile("SmallSamples.root", "recreate")

Nexp = 1000000
N = 10

# Implement corrective constant a la Gurland and Tripathi equation 2
# to make standard deviation estimate unbiased (assumes normal).
cN = math.sqrt((N-1)/2) * TMath.Gamma((N-1)/2)/ TMath.Gamma(N/2)
print("Gurland and Tripathi exact correction for N =",N,":",cN)

# Exact method
random.seed(200)
for i in range(Nexp):                     # Nexp repetitions
    xsum = 0.0
    xxsum = 0.0
    for j in range(N):                    # Small sample of N from standardized normal
        x = random.gauss(0.0,1.0)
        xsum += x
        xxsum += x**2
    xave  = xsum/N
    xxave = xxsum/N
    var1 = xxave - xave**2                # Naive biased variance estimate.
    var2 = (N/(N-1))*var1                 # Variance estimate with the N/(N-1) bias correction.

    sdnaive = math.sqrt(var2)
    hmean.Fill(xave)                      # Sample mean
    hsmean.Fill(math.sqrt(N)*xave)        # Sample mean scaled by sqrt(N)
    
    hvar1.Fill(var1)                      # Variance estimate (biased)
    hsd1.Fill(sdnaive)                    # Standard deviation estimate (biased)
    hsem1.Fill(sdnaive/math.sqrt(N))      # Standard error on the mean estimate (biased)

    hvar2.Fill(var2)                      # Variance estimate (unbiased)
    sd = cN*sdnaive 
    hsd2.Fill(sd)                         # With G-T correction (unbiased)
    hsem2.Fill(sd/math.sqrt(N))           # With G-T correction (unbiased)
    
histList = [ hmean, hsmean, hvar1, hsd1, hsem1, hvar2, hsd2, hsem2 ]
for h in histList:
    histInfo(h)

f.Write()
