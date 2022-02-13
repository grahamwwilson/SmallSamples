# histUtils.py
from ROOT import TH1D

def histInfo(h):
    h.SetLineWidth(2)
    h.SetLineColor(4)
    print(" ")
    h.Print()
    print("Mean: ",h.GetMean()," +- ",h.GetMeanError())
    print("RMS : ",h.GetRMS() ," +- ",h.GetRMSError())
#    h.Write()      # Make sure it gets written to output file too
