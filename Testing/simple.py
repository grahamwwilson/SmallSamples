import numpy as np
import myPythonCheck
import analyze

myPythonCheck.Check()                         # Enforce use of python3

# Check numpy version
print("numpy version",np.__version__)

x, y, z = analyze.GetData()
T1Obs, T2Obs, T3Obs = analyze.Analyze(x,y,z)
