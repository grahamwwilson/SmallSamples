# myPythonCheck.py
import sys

def Check():
    if sys.version_info[0] < 3:
        raise Exception("Python 3 or a more recent version is required for proper execution")
