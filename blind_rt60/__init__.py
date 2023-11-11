"""
blind_rt60 is a Python-based implementation of blind estimation of reverberation time (rt60)
================================================
Documentation is available in the docstrings and
online at https://github.com/nuniz/pyBlindRT/blob/main/README.md.

Contents
--------
blind_rt60 imports all the functions from numpy, and in addition provides:
 BlindRT60       --- A python module that estimates rt60 based on an input signal

Public API in the main TorchGating namespace
--------------------------------------
::
 __version__       --- blind_rt60 version string

References
--------------------------------------
The algorithm was originally proposed by Ratnam et al. [1]
[1] Ratnam, Rama & Jones, Douglas & Wheeler, Bruce & O'Brien, William & Lansing, Charissa & Feng, Albert. (2003).
Blind estimation of reverberation time. The Journal of the Acoustical Society of America. 114. 2877-92.
10.1121/1.1616578.

"""

from .estimation import BlindRT60
from .version import __version__
