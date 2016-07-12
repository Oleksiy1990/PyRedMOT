import numpy as np
from pyConstants import * 

mSr84 = 84*amu
mSr86 = 86*amu
mSr87 = 87*amu
mSr88 = 88*amu

blueGamma = 2*np.pi*30.5*10e6
redGamma = 2*np.pi*7.4*10e3

blueFreqInvCm = 21698.452
blueFreqInvM = 100*blueFreqInvCm
redFreqInvCm = 14504.334
redFreqInvM = 100*redFreqInvCm

blueIsat = 407 # [W/m^2]
redIsat = 3/100 # [W/m^2 = 3 microW/cm^2]

kVecBlue = 2*np.pi*blueFreqInvM
kVecRed = 2*np.pi*redFreqInvM
