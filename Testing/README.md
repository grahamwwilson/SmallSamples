# SmallSamples/Testing
2-sample test statistics. Either read data from files or generate our own.

# simple.py
Apply various 2-sample tests using sample sizes in each sample of about 10.
Interestingly, scipy, has direct implementations of exact permutation tests.
This runs in about 3 seconds for me

# permute.py
Do our own exact permutation test with 3 statistics, also serves as 
a cross-check of the scipy implementation and provides histograms, and 
more potential insight. Runs in about a minute for problems of size 
around n1+n2=20.

# analyze.py
This contains the GetData and Analyze functions now used by 
simple.py and permute.py. (It is no longer a main program).
