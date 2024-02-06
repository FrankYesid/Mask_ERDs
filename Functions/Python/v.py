import numpy as np
from matplotlib import mlab

data = range(1,1000) #Dummy data. Just for testing

Fs = 8000
tWindow = 64e-3
NWindow = Fs*tWindow
window = np.hamming(NWindow)

NFFT = 512
NOverlap = NWindow/2

[s, f, t] = mlab.specgram(data, NFFT = NFFT, Fs = Fs, window = window, noverlap = NOverlap)

print s.shape
print s
print f.shape
print f
print t.shape
print t