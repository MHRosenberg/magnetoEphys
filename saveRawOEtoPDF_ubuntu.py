import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import style
from matplotlib.backends.backend_pdf import PdfPages
import math
import sys
import os
import glob
import time
import datetime
import webbrowser

### BEWARE: specifying the file doesn't seem to always work... unclear why it does in some cases but not others --> specify the directory only!
sys.path.insert(0, '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Analysis/') 
import OpenEphys as OE

#WIP: some day make a slider but it sounds like a pain from matplotlib.widgets import Slider

CHROME_PATH = 'C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s'

DIR_OF_RECORDINGS = 'C:\\Users\\yatang\\Documents\\Ephys\\Open Ephys'
LOAD_MOST_RECENT_REC = 0 # use 1 to load the most recent and 0 if you wish to specify the path
LOAD_HERE = 'yes' # use yes to load data in cwd or no to not. NOTE: NO EFFECT IF LOAD_MOST_RECENT_REC = 1!
LOAD_THIS = '' #specific the raw data path here. NOTE: NO EFFECT IF EITHER LOAD_MOST_RECENT OR LOAD_HERE IS SELECTED!

NUM_CH = 64
SAMPLE_RATE_HZ = 30000 # TO DO: load this from the events file
FIRST_CH = 1 # 1 indexed
LAST_CH = 3

###### plotting parameters
X_TICK_STEP = 0.050 # in same units as time
Y_TICK_STEP = 100 # in microvolts (Mu agrees that it should be in microV)
### currently only works for a second at a time
FIRST_SEC_PLOTTED = 1.5 # set to 0 to begin at the first sample
LAST_SEC_PLOTTED = 2 # set to -1 if you want to plot to the last sample 

if LOAD_MOST_RECENT_REC ==1:
	# all_subdirs = [d for d in os.listdir(DIR_OF_RECORDINGS) if os.path.isdir(d)]
	# print(all_subdirs)
	# latest_subdir = max(all_subdirs, key=os.path.getmtime)
	latest_subdir = max(glob.glob(os.path.join(DIR_OF_RECORDINGS, '*\\')), key=os.path.getmtime) ## TO DO: make this platform independent
	pathToFile = latest_subdir + '\\100_CH1.continuous'
	RAW_DATA_PATH = latest_subdir
	print('loading most recent dir: {0}\nChange parameters at the top of the script to change this behavior.'.format(RAW_DATA_PATH))
else:
	if LOAD_HERE.lower() == 'yes':
		RAW_DATA_PATH = os.getcwd()
		pathToFile = RAW_DATA_PATH + '/100_CH1.continuous'
		print('loading recordings from cwd')
		
	else:
		if LOAD_THIS != '':
			### laptop paths
			# RAW_DATA_PATH = '/home/orthogonull/a_MHR/a_Research/a_WIP/2018-04-29_23-20-48'
			# pathToFile = '/home/orthogonull/a_MHR/a_Research/a_WIP/2018-04-29_23-20-48/100_CH1.continuous'
			# figSaveDir = '/home/orthogonull/a_MHR/a_Research/a_WIP'

			### ephys computer paths
			pathToFile = 'C:\\Users\\yatang\\Documents\\Ephys\\Open Ephys\\2018-05-04_13-55-49_MR1_r2_pre\\100_CH1.continuous'
			RAW_DATA_PATH = 'C:\\Users\\yatang\\Documents\\Ephys\\Open Ephys\\2018-05-04_13-55-49_MR1_r2_pre\\'
			figSaveDir = 'C:\\Users\\yatang\\Documents\\Ephys\\Open Ephys\\'
		else:
			sys.exit('ERROR: user parameters in this script failed to find a valid directory of recordings!')

print('loading: ' + RAW_DATA_PATH)

data_dict = OE.load(pathToFile) # returns a dict with data, timestamps, etc.
print(data_dict)

# open ephys's BROKEN (float given but expects int error) downsample code
# def downsample(trace,down):
#     downsampled = scipy.signal.resample(trace,np.shape(trace)[0]/down)
#     return downsampled
data = OE.loadFolderToArray(RAW_DATA_PATH, channels = range(FIRST_CH, LAST_CH+1), chprefix = 'CH', dtype = float, session = '0', source = '100')
print(np.shape(data))
totalNumSamples = np.shape(data)[0] #no longer used but available for later
print(data)

### FIGURE STYLING:
titleFontSize = 120
LINE_WIDTH = 0.5
style.use('fivethirtyeight')
sns.set_style('white')
axes_font_size = 100
axes_font_weight = 'demi'
x_axes_font = {'fontsize': axes_font_size ,
               'fontweight': axes_font_weight,
               'verticalalignment': 'top',
               'horizontalalignment': 'center'}
y_axes_font = {'fontsize': axes_font_size ,
               'fontweight': axes_font_weight,
               'verticalalignment': 'bottom',
               'horizontalalignment': 'center'}
figTitle = 'num samples, num channels = ' + str(np.shape(data))
plt.rc('xtick', labelsize=32)
plt.rc('ytick',labelsize=32)
###

if LAST_SEC_PLOTTED == -1:
	sample_range_to_plot = [FIRST_SEC_PLOTTED*SAMPLE_RATE_HZ, totalNumSamples] # specify as a 2 int list [first sample included, last sample included]
elif LAST_SEC_PLOTTED > 0:
	sample_range_to_plot = [round(FIRST_SEC_PLOTTED*SAMPLE_RATE_HZ), round(SAMPLE_RATE_HZ*(LAST_SEC_PLOTTED))] # specify as a 2 int list [first sample included, last sample included]
else:
	print('ERROR: invalid LAST_SEC_PLOTTED parameter (must be -1 to run until the end or greater than 0')
	sys.exit() 

numSamples = sample_range_to_plot[-1] - sample_range_to_plot[0]

# SAMPLES_PER_PAGE = 30000 #this gives proper xticks and shows 1 sec per page
SAMPLES_PER_PAGE = 5000 

X_LABEL = 'seconds'
Y_LABEL = 'microvolts'
# xs = xs = np.linspace(0,1,np.shape(data)[0]) 

############## WIP don't forget to rescale to get correct units

print('sample range to plot: ' + str(sample_range_to_plot))

sampleInds = np.asarray(range(sample_range_to_plot[0],sample_range_to_plot[1]))
# sampleTimePts = sampleInds/SAMPLE_RATE_HZ
sampleTimePts = sampleInds/SAMPLE_RATE_HZ
print('plotted time points: ' + str(sampleTimePts[0]) + ' to ' + str(sampleTimePts[-1]))

numPages = math.ceil(numSamples/SAMPLES_PER_PAGE)
subplot_dim = (LAST_CH-FIRST_CH+1, 1)
# numXticks = -999 # temporary value designed to create an error if not handled properly within the loop 

ts = time.time()
timeStr = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d__%H_%M_%S')
pdfName = 'rawTraces' + timeStr + '.pdf'

with PdfPages(pdfName) as pdf_out:

	for pageInd in range(1,numPages+1):
	# for pageInd in range(1,numPages+1):
	# for pageInd in range(1,numPages-1):
		f = plt.figure(figsize=(200,200)) # REDUCE THESE NUMBERS IF YOU GET SEGMENTATION/MEMORY ERRORS!
		
		firstInd = (pageInd-1)*SAMPLES_PER_PAGE
		lastInd = pageInd*SAMPLES_PER_PAGE

		# print('firstInd: ' + str(firstInd))
		# print('lastInd: ' + str(lastInd))

		xs = sampleTimePts[firstInd:lastInd]

		min_xs = xs[0]
		max_xs = xs[-1]

		if pageInd == 1:
			numXticks = int(round((max_xs - min_xs)/X_TICK_STEP))

		for chInd in range(FIRST_CH-1,LAST_CH):
			print('plotting channel:' + str(chInd+1) + '\npage num: ' + str(pageInd) + ' of ' + str(numPages) + ' total')

			ys = data[firstInd:lastInd,chInd]

			sp = plt.subplot2grid(subplot_dim, (chInd, 0), aspect="auto",adjustable='box-forced',  colspan=1, rowspan=1)
			sp.autoscale(True)	
			
			# ax = plt.gca() # to get current axis

			try:
				sp.plot(xs,ys,linewidth=LINE_WIDTH)
			except:
				print('ERROR: the number of samples can not be less than the value specified in SAMPLES_PER_PAGE')
				print('WIP: this error can probably be handled by explicitly processing the remainder of samples')
				### 
				# sys.exit() # to halt execution upon error
				continue # to disregard and continue to plot
				###

			plt.xticks(np.linspace(min_xs,max_xs,num=numXticks)) 

		f.tight_layout() # makes the figure fill the whole window
		f.subplots_adjust(hspace=0.2, left=0.1, right = 0.9, bottom = 0.05, top = 0.95)
		f.suptitle(RAW_DATA_PATH, fontsize= titleFontSize)

		plt.xlabel(X_LABEL, x_axes_font)
		plt.ylabel(Y_LABEL, y_axes_font)

		print('saving page: ' + str(pageInd))
		pdf_out.savefig(f)
		
		f.clf() # might help with out of memory error

print('finished: rawTraces.pdf saved to cwd')

webbrowser.get(CHROME_PATH).open(pdfName)
# plt.subplot_tool() # creates a gui for window border spacing etc

# create subplots
# sp1 = plt.subplot2grid(subplot_dim, (0, 0), aspect="auto",adjustable='box-forced', colspan=3, rowspan=2)
# sp1 = plt.subplot2grid(subplot_dim, (0, 0), aspect="auto",adjustable='box-forced', colspan=1, rowspan=1)
# sp = plt.subplot2grid(subplot_dim, (0, chInd), aspect="auto",adjustable='box-forced', sharex='col', sharey='row')
# spn = plt.subplot2grid(subplot_dim, (chInd, 0), aspect="auto",adjustable='box-forced', sharex=sp1, sharey=sp1)

# ###works ish
# fig, ax = plt.subplots(64,1, figsize=(10,100),sharex='col', sharey='row', squeeze=True)

# fig.tight_layout()

# ax[chInd].plot(xs, ys, linewidth = 0.7, linestyle = '--')

# ax = fig.add_subplot(65,1,chInd+1)	

# ax = fig.add_subplot(111)
# plt.show()
# # plt.savefig(figSaveDir + figTitle, bbox_inches='tight', pad_inches=0.003)
# plt.close(fig)
