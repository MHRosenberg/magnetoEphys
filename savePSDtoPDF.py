
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import style
from matplotlib import pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
import math
import sys
sys.path.insert(0, '/media/matthew/Data/gitResearch/OpenEphysAnalysisTools/Python3/OpenEphys.py') # TO DO: get the import folder thing to work with a path to local open ephys github folder
import OpenEphys as OE
import os
import glob
import time
import datetime
import webbrowser
from pathlib import Path, PureWindowsPath # only works on 3.4+ (not 2.7)
import platform

### CHECK USER PARAMETERS THAT FOLLLOW THESE FUNCTIONS!

def genericInvalidParamError(userParamStr):
	return print('ERROR: You provided an invalid user parameter! Correct this param: ' + userParamStr)

def loadRawOEdir(rawDataDir):
	# print('loading: ' + pathToFile)
	# data_dict = OE.load(pathToFile) # returns a dict with data, timestamps, etc.
	# print(data_dict)

	# # open ephys's BROKEN (float given but expects int error) downsample code
	# # def downsample(trace,down):
	# #     downsampled = scipy.signal.resample(trace,np.shape(trace)[0]/down)
	# #     return downsampled
	data = OE.loadFolderToArray(rawDataDir, channels = range(FIRST_CH, LAST_CH+1), chprefix = 'CH', dtype = float, session = '0', source = '100')
	print(np.shape(data))
	totalNumSamples = np.shape(data)[0] #no longer used but available for later
	print(data)
	return data

def saveRawDataPathsToTxt(pathList):
	ts = time.time()
	timeStr = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d__%H_%M_%S')
	inputDataPathsTextFileName = 'inputDataPaths' + timeStr + '.txt'
	with open(inputDataPathsTextFileName, 'w') as pathFile:
		for path in pathList:
			# cwd = os.getcwd()
			# path = path.replace(path[:2],cwd + '/')
			path = path.replace('/Continuous_Data.openephys','')
			pathFile.write("%s\n" % path)
	print('these paths are saved rawDataPaths.txt in the cwd')

def generateRecordingIDsToRun():
	recPathlist = []
	with open('rawDataPaths.txt') as paths:
		for recording in paths:
			recording = recording.strip() # remove new line
			slashInd = recording.rfind('/')
			recID = recording[slashInd+1:]
			recPathlist.append(recID)

	recIDlistFile = open('recordingIDsToRun.txt', 'w')
	for recID in recPathlist:
		recIDlistFile.write("%s\n" % recID)
	recIDlistFile.close()
	return recPathlist

### INPUT DATA PARAMETERS:
### choose among: 'mostRecent' # 'singleDir' 'singleFile' 'allInChildDirs' 
DATA_SELECTION_SCHEME = 'allInChildDirs' # 'singleDir' 'singleFile' 'allInChildDirs' 'fromTextFileInCwd' 
EPHYS_PARENT_PATH = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Data/raw/noiseTests/' #use '' for current work directory
SUMMARY_OUTPUT_PATH = './NoiseTests/PSDs' # or modify to custom dir 
PATH_a_Analysis = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Analysis'
PATH_rawParent = os.getcwd()
PATH_a_Analysis = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Analysis'
PATH_rawParent = os.getcwd()
BASHRC_DIR = '/home/matthew'

### ANALYSIS PARAMETERS:
NUM_CH = 64
SAMPLE_RATE_HZ = 30000
FIRST_CH = 1 # 1 indexed
LAST_CH = 3
INSPECT_INDIVIDUAL_CHAN = 0 # 0 for no; 1 for yes
OPEN_EACH_PDF = 1 # 1 to open each psd plot pdf before proceeding to the next plot
SAVE_ALL_TO_PDF = 1
## for the zoomed in plot (there's a plot of all of the frequencies as well)
HIGH_FREQ_CUTOFF = 20 # potentially NOT in Hz!

NFFT = 1024

###data selection 
print('\nDATA_SELECTION_SCHEME: ' + DATA_SELECTION_SCHEME)
if DATA_SELECTION_SCHEME == 'mostRecent':
	latest_subdir = max(glob.glob(Path(EPHYS_PARENT_PATH).joinpath('*')), key=os.path.getmtime)
	pathToFile = Path(latest_subdir).joinpath('/100_CH1.continuous')
	RAW_DATA_PATH = latest_subdir
elif DATA_SELECTION_SCHEME == 'fromTextFileInCwd':
	CHOSEN_REC_DIR = 'populateThisAppropriately' # 
elif DATA_SELECTION_SCHEME == 'singleFile':
	print()
elif DATA_SELECTION_SCHEME == 'allInChildDirs':
	print('processing all raw Open Ephys data in child directories of EPHYS_PARENT_PATH (use \'\' for cwd)\n')
	pathList = []
	if EPHYS_PARENT_PATH == '':
		EPHYS_PARENT_PATH = os.getcwd()
	print('processing data in: ' + EPHYS_PARENT_PATH)    
	recPathlist = list(glob.iglob(str(Path(EPHYS_PARENT_PATH).joinpath('**/*.openephys')), recursive = True)) ### recursive file selection
	print('\n'.join(recPathlist))
	saveRawDataPathsToTxt(recPathlist)
	print('\nremove some of the above from recordingIDsToRun.txt and change DATA_SELECTION_SCHEME to fromTextFileInCwd to run on a subset of the child dir data\n')
else:
	genericInvalidParamError()

### directiory handling
if platform.system() == 'Linux':
	CHROME_PATH = '/usr/bin/google-chrome' # ubuntu tower path
	print('Linux system detected\nChrome path set to: ' + CHROME_PATH)
	if not os.path.exists(SUMMARY_OUTPUT_PATH):
		os.makedirs(SUMMARY_OUTPUT_PATH)
elif platform.system() == 'Windows':
	CHROME_PATH = 'C:/Program Files (x86)/Google/Chrome/Application/chrome.exe %s' # windows stim/ephys computer path	
	DIR_OF_RECORDINGS = 'C:\\Users\\yatang\\Documents\\Ephys\\Open Ephys' # windows ephys/stim path
	### ephys computer paths
	pathToRawDataFile  = Path(EPHYS_PARENT_PATH).joinpath(CHOSEN_REC_DIR).with_suffix('100_CH1.continuous')
	RAW_DATA_PATH = 'C:\\Users\\yatang\\Documents\\Ephys\\Open Ephys\\2018-05-04_13-55-49_MR1_r2_pre\\'
	figSaveDir = 'C:\\Users\\yatang\\Documents\\Ephys\\Open Ephys\\'
else:
	print('WARNING: code has been optimized for Windows and Mac operating systems only.\nSet your own system parameters at the top of .py file yourself >:p')


print('skipping individual channel inspection\n set INSPECT_INDIVIDUAL_CHAN to 1 for inspection and choose CHAN_NUM')
for rawDataPath in recPathlist:
	rawDataPath = rawDataPath.replace('/Continuous_Data.openephys','')
	recordingName = str(os.path.basename(os.path.normpath(rawDataPath)))
	recordingName = recordingName.replace('/','')
	data = loadRawOEdir(rawDataPath)
	ts = time.time()
	timeStr = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d__%H_%M_%S')
	pdfName = Path(SUMMARY_OUTPUT_PATH).joinpath(recordingName + timeStr + '.pdf')
	numChannelsPlotted = LAST_CH - FIRST_CH + 1
	with PdfPages(pdfName) as pdf_out:
		for chInd in range(FIRST_CH-1, LAST_CH):
			# f = plt.figure(figsize=(200,200)) # REDUCE THESE NUMBERS IF YOU GET SEGMENTATION/MEMORY ERRORS!
			f = plt.figure() # first page 
			ax1 = plt.subplot(311)
			plt.plot(data[:,chInd])
			plt.subplot(312)
			Pxx, freqs = plt.psd(data[:,chInd], 2048, SAMPLE_RATE_HZ)
			
			# freqs, Pxx  = plt.psd(data[:,chInd], 2048, SAMPLE_RATE_HZ)
			plt.subplot(313)
			# extraTicks = [-50, -25, -10, -5, -2, -1, 0,1]
			# plt.semilogy(freqs[0:HIGH_FREQ_CUTOFF],Pxx[0:HIGH_FREQ_CUTOFF], yticks = list(plt.yticks()[0]) + extraTicks)
			plt.semilogy(freqs[0:HIGH_FREQ_CUTOFF],Pxx[0:HIGH_FREQ_CUTOFF])

			if SAVE_ALL_TO_PDF == 1:
				pdf_out.savefig(f)

			f = plt.figure() # second page
			
			### testing
			plt.subplot(311)
			plt.magnitude_spectrum(data[:,chInd], Fs = SAMPLE_RATE_HZ, scale = 'dB')

			plt.subplot(312)
			
			Pxx, freqs, bins, im = plt.specgram(data[:,chInd], NFFT=NFFT, Fs=SAMPLE_RATE_HZ, noverlap=900)


			if INSPECT_INDIVIDUAL_CHAN ==1:
				print('close the plot window to advance to the next channel\nOR set INSPECT_ONE_BY_ONE to 0')
				plt.show()
			
			if chInd == LAST_CH:
				f.suptitle(RAW_DATA_PATH, fontsize= titleFontSize)
			print('saving PSD for channel: ' + str(chInd+1) + ' of ' + str(numChannelsPlotted))
			if SAVE_ALL_TO_PDF == 1:
				pdf_out.savefig(f)
			plt.close(f)
			f.clf() # might help with out of memory error
			
	print('finished: ' + str(pdfName) + ' saved to cwd')
	if OPEN_EACH_PDF == 1:
		webbrowser.get(str(CHROME_PATH)).open(str(pdfName))
	else:
		print('PSD results saved to cwd')	




# elif SAVE_ALL_TO_PDF == 0:
# 	print('skipping pdf output saving step \nset SAVE_ALL_TO_PDF to 1 or output pdf')
# else:
# 	print('ERROR: SAVE_ALL_TO_PDF must be 0 or 1')



# ### FIGURE STYLING:
# titleFontSize = 120
# LINE_WIDTH = 0.5
# style.use('fivethirtyeight')
# sns.set_style('white')
# axes_font_size = 100
# axes_font_weight = 'demi'
# x_axes_font = {'fontsize': axes_font_size ,
#                'fontweight': axes_font_weight,
#                'verticalalignment': 'top',
#                'horizontalalignment': 'center'}
# y_axes_font = {'fontsize': axes_font_size ,
#                'fontweight': axes_font_weight,
#                'verticalalignment': 'bottom',
#                'horizontalalignment': 'center'}
# figTitle = 'num samples, num channels = ' + str(np.shape(data))
# plt.rc('xtick', labelsize=32)
# plt.rc('ytick',labelsize=32)
# ###

# X_TICK_STEP = 0.050 # in same units as time
# Y_TICK_STEP = 100 # in microvolts (Mu agrees that it should be in microV)

# ### currently only works for a second at a time
# FIRST_SEC_PLOTTED = 1 # set to 0 to begin at the first sample
# LAST_SEC_PLOTTED = 1.25 # set to -1 if you want to plot to the last sample 

# if LAST_SEC_PLOTTED == -1:
# 	sample_range_to_plot = [FIRST_SEC_PLOTTED*SAMPLE_RATE_HZ, totalNumSamples] # specify as a 2 int list [first sample included, last sample included]
# elif LAST_SEC_PLOTTED > 0:
# 	sample_range_to_plot = [round(FIRST_SEC_PLOTTED*SAMPLE_RATE_HZ), round(SAMPLE_RATE_HZ*(LAST_SEC_PLOTTED))] # specify as a 2 int list [first sample included, last sample included]
# else:
# 	print('ERROR: invalid LAST_SEC_PLOTTED parameter (must be -1 to run until the end or greater than 0')
# 	sys.exit() 

# numSamples = sample_range_to_plot[-1] - sample_range_to_plot[0]

# # SAMPLES_PER_PAGE = 30000 #this gives proper xticks and shows 1 sec per page
# SAMPLES_PER_PAGE = 5000 

# X_LABEL = 'seconds'
# Y_LABEL = 'microvolts'
# # xs = xs = np.linspace(0,1,np.shape(data)[0]) 

# ############## WIP don't forget to rescale to get correct units

# print('sample range to plot: ' + str(sample_range_to_plot))

# sampleInds = np.asarray(range(sample_range_to_plot[0],sample_range_to_plot[1]))
# # sampleTimePts = sampleInds/SAMPLE_RATE_HZ
# sampleTimePts = sampleInds/SAMPLE_RATE_HZ
# print('plotted time points: ' + str(sampleTimePts[0]) + ' to ' + str(sampleTimePts[-1]))

# numPages = math.ceil(numSamples/SAMPLES_PER_PAGE)
# subplot_dim = (LAST_CH-FIRST_CH+1, 1)
# # numXticks = -999 # temporary value designed to create an error if not handled properly within the loop 

# ts = time.time()
# timeStr = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d__%H_%M_%S')
# pdfName = 'rawTraces' + timeStr + '.pdf'

# with PdfPages(pdfName) as pdf_out:

# 	for pageInd in range(1,numPages+1):
# 	# for pageInd in range(1,numPages+1):
# 	# for pageInd in range(1,numPages-1):
# 		f = plt.figure(figsize=(200,200)) # REDUCE THESE NUMBERS IF YOU GET SEGMENTATION/MEMORY ERRORS!
		
# 		firstInd = (pageInd-1)*SAMPLES_PER_PAGE
# 		lastInd = pageInd*SAMPLES_PER_PAGE

# 		# print('firstInd: ' + str(firstInd))
# 		# print('lastInd: ' + str(lastInd))

# 		xs = sampleTimePts[firstInd:lastInd]

# 		min_xs = xs[0]
# 		max_xs = xs[-1]

# 		if pageInd == 1:
# 			numXticks = int(round((max_xs - min_xs)/X_TICK_STEP))

# 		for chInd in range(FIRST_CH-1,LAST_CH):
# 			print('plotting channel:' + str(chInd+1) + '\npage num: ' + str(pageInd) + ' of ' + str(numPages) + ' total')

# 			ys = data[firstInd:lastInd,chInd]

# 			sp = plt.subplot2grid(subplot_dim, (chInd, 0), aspect="auto",adjustable='box-forced',  colspan=1, rowspan=1)
# 			sp.autoscale(True)	
			
# 			# ax = plt.gca() # to get current axis

# 			try:
# 				sp.plot(xs,ys,linewidth=LINE_WIDTH)
# 			except:
# 				print('ERROR: the number of samples can not be less than the value specified in SAMPLES_PER_PAGE')
# 				print('WIP: this error can probably be handled by explicitly processing the remainder of samples')
# 				### 
# 				# sys.exit() # to halt execution upon error
# 				continue # to disregard and continue to plot
# 				###

# 			plt.xticks(np.linspace(min_xs,max_xs,num=numXticks)) 

# 		f.tight_layout() # makes the figure fill the whole window
# 		f.subplots_adjust(hspace=0.2, left=0.1, right = 0.9, bottom = 0.05, top = 0.95)
# 		f.suptitle(RAW_DATA_PATH, fontsize= titleFontSize)

# 		plt.xlabel(X_LABEL, x_axes_font)
# 		plt.ylabel(Y_LABEL, y_axes_font)

# 		print('saving page: ' + str(pageInd))
# 		pdf_out.savefig(f)
		
# 		f.clf() # might help with out of memory error

# print('finished: rawTraces.pdf saved to cwd')

# webbrowser.get(CHROME_PATH).open(pdfName)

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