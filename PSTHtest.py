import numpy as np
import scipy.io as sio
import pandas as pd
import sys
import os.path
import glob
import pdb

### use this to get a traceback for warnings!
#import warnings
#warnings.simplefilter("error")
#warnings.simplefilter("ignore", DeprecationWarning)

### BEWARE: specifying the file doesn't seem to always work... unclear why it does in some cases but not others --> specify the directory only!
sys.path.insert(0, '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Analysis/') 
import OpenEphys as OE

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import seaborn as sns
from matplotlib import style
from matplotlib import pyplot as plt 
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path #, PureWindowsPath # only works on 3.4+ (not 2.7)
from scipy.signal import welch, hanning, windows

mat = np.full((10,100), np.nan)
mat = np.random.random((10,100)) *10

plt.figure()

for unitInd, unitSpikeTime in enumerate(mat[:,0]):
    print("{0} {1}".format(unitInd,unitSpikeTime))
    plt.plot(mat[unitInd],label=str(unitInd))

plt.legend()