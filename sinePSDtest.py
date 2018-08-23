### BEWARE: specifying the file doesn't seem to always work... unclear why it does in some cases but not others --> specify the directory only!
sys.path.insert(0, '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Analysis/') 
import OpenEphys as OE

import numpy as np
import scipy.io as sio
import pandas as pd
import sys
import os.path
import glob
import pdb
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import seaborn as sns
from matplotlib import style
from matplotlib import pyplot as plt 
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path #, PureWindowsPath # only works on 3.4+ (not 2.7)
from scipy.signal import welch, hanning, windows

### use this to get a traceback for warnings!
#import warnings
#warnings.simplefilter("error")
#warnings.simplefilter("ignore", DeprecationWarning)

###### user parameters
SAMPLE_RATE = 30000 # in hz
RECORDING_DURATION_MINS = 10
NOISE_FACTOR = 10
MAX_FREQ_PLOTTED = 20 # hz

plt.close()
recording_duration = RECORDING_DURATION_MINS * 60 # in secs
numDataPts = recording_duration * SAMPLE_RATE
timePts = np.linspace(0,recording_duration, numDataPts)
freq1 = 2 # hz
sine1 = np.sin(2 * np.pi * freq1 * timePts)
freq2 = 0.1
sine2 = np.sin(2 * np.pi * freq2 * timePts)
freq3 = 8
sine3 = np.sin(2 * np.pi * freq3 * timePts)
freq4 = 200
sine4 = np.sin(2 * np.pi * freq4 * timePts)
freq5 = 15
sine5 = np.sin(2 * np.pi * freq5 * timePts)
freq6 = 4
sine6 = np.sin(2 * np.pi * freq6 * timePts)
noise = NOISE_FACTOR * np.random.normal(0,1,timePts.shape)

#signal = sine1 + sine2 + sine3 + sine4 + sine5 + sine6 + noise
#signal = sine1 + noise

        
### sine stim's bin'd sample rate
sineStimSampleRate = SAMPLE_RATE
print('sineStimSampleRate = {0}'.format(sineStimSampleRate))
### neuron's bin'd sample rate
#numBins, remainderSecs = np.divmod(totalStimDurationSecs, BIN_SIZE_SEC)
#bins = np.linspace(FIRST_SAMPLE, FIRST_SAMPLE + round(float(numBins) * float(BIN_SIZE_SEC)), int(numBins), endpoint=True) 
#binSampleRate = 1./np.mean(np.diff(bins))

### Welch PSD
### sine stimulus
welchFreqs_signal, Pxx_den = welch(signal, nperseg = SAMPLE_RATE * 10, fs = sineStimSampleRate) ######## WARNING: Are the timePts }{ sample rate reasonable

NUM_ROWS = 2
NUM_COLUMNS = 1 
### gridspec.GridSpec(NUM_ROWS,NUM_COLUMNS) # state-based versions for subplot
### plt.subplot(611) ### doesn't require gridspec but is kinda clunky for dynamic changes        
FIG_SIZE = (NUM_ROWS, NUM_COLUMNS)
fig = plt.figure(figsize = FIG_SIZE)
rowPosition = 0
columnPosition = 0

### raw sine stimulus
ax2 = plt.subplot2grid((NUM_ROWS, NUM_COLUMNS), (rowPosition,columnPosition))
rowPosition += 1
ax2.set_title('raw signal')
ax2.plot(timePts, signal)

### sine stimulus psd
ax2 = plt.subplot2grid((NUM_ROWS, NUM_COLUMNS), (rowPosition,columnPosition))
rowPosition += 1
ax2.set_title('sine stimulus Welch PSD')
ax2.set_xlabel('frequency [Hz]')
ax2.set_ylabel('PSD [V**2/Hz?]')
ax2.semilogy(welchFreqs_signal,Pxx_den, '-o')
ax2.set_xlim(0, MAX_FREQ_PLOTTED)

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.show()
#ax1.set_xlim(0, MAX_FREQ_PLOTTED)