
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
import pdb


### CHECK USER PARAMETERS THAT FOLLLOW THESE FUNCTIONS!

def genericInvalidParamError(userParamStr):
	return print('ERROR: You provided an invalid user parameter! Correct this param: ' + userParamStr)

def loadRawOEdir(rawDataDir):
	data = OE.loadFolderToArray(rawDataDir, channels = range(FIRST_CH, LAST_CH+1), chprefix = 'CH', dtype = float, session = '0', source = '100')
	print(np.shape(data))
	totalNumSamples = np.shape(data)[0] #no longer used but available for later
	print(data)
	return data

def loadTTLdata(rawDataDir): # from the .events file
	event_path = str(Path(rawDataDir).joinpath('all_channels.events'))
	print('loading: ' + event_path)
	eventsDict = OE.load(event_path) ### load TTL dict file
	print(eventsDict)

	print(eventsDict['eventType'])
	print(eventsDict['eventId'])


	print('\nall the events that are saved have type TTL = 3 ; Network Event = 5')
	eventIDs = eventsDict['eventType']	 
	print('event IDs: ' + str(eventIDs))
	TTL_values = eventsDict['eventId']
	print('shape of TTL values' + str(TTL_values.shape))
	print('TTL values: ' + str(TTL_values))
	TTL_timeStamps = eventsDict['timestamps']
	print('shape of timestamps: ' + str(TTL_timeStamps.shape))
	print('time stamps: ' + str(TTL_timeStamps))
	return TTL_timeStamps, TTL_values

def saveRawDataPathsToTxt(pathList):
	ts = time.time()
	timeStr = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d__%H_%M')
	inputDataPathsTextFileName = 'inputDataPaths' + timeStr + '.txt'
	with open(inputDataPathsTextFileName, 'w') as pathFile:
		for path in pathList:
			# cwd = os.getcwd()
			# path = path.replace(path[:2],cwd + '/')
			path = path.replace('/Continuous_Data.openephys','')
			pathFile.write("%s\n" % path)
	print('these paths are saved inputDataPaths.txt in the cwd')
	return inputDataPathsTextFileName

def getRecordingNameFromPath(recordingPath):
	recordingPath = recordingPath.strip() # remove new line
	slashInd = recordingPath.rfind('/')
	recordingName = recordingPath[slashInd+1:]
	return recordingName

def generateRecordingIDsToRun(inputDataPathsTextFileName):
	recPathlist = []
	with open(inputDataPathsTextFileName) as paths:
		for recording in paths:
			recording = recording.strip() # remove new line
			# slashInd = recording.rfind('/')
			# recID = recording[slashInd+1:]
			recording = recording.replace('/Continuous_Data.openephys','')
			recPathlist.append(recording)
	recIDlistFile = open('recordingIDsToRun.txt', 'w')
	for recID in recPathlist:
		recIDlistFile.write("%s\n" % recID)
	recIDlistFile.close()
	return recPathlist

### INPUT DATA PARAMETERS:
### choose among: 'mostRecent' # 'singleDir' 'singleFile' 'allInChildDirs' 
DATA_SELECTION_SCHEME = 'allInChildDirs' # 'singleDir' 'singleFile' 'allInChildDirs' 'fromTextFileInCwd' 
EPHYS_PARENT_PATH = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Data/raw/MR1_r1_5_2_18/' #use '' for current work directory
SUMMARY_OUTPUT_PATH = './stimReconstruction/TTLtoSine' # or modify to custom dir 
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
INSPECT_EACH_PLOT = 1 # 0 for no; 1 for yes
OPEN_EACH_PDF = 1 # 1 to open each psd plot pdf before proceeding to the next plot
SAVE_ALL_TO_PDF = 1


## TTL plotting paramete
FIRST_TTL_IND = 0 # skip the first one because this is presumably truncated wrt to the normal period
LAST_TTL_IND = -1 # use -1 to go until the end
NUM_PLATEU_PTS = 4 # number of points to add per high/low TTL pulse to recreate the square wave from transistion points
PHASE_OFFSET = 2 # phase offset (in radians?) to align the sine wave to the squarewave from TTLs

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

	inputDataPathsTextFileName = saveRawDataPathsToTxt(recPathlist)
	recPathlist = generateRecordingIDsToRun(inputDataPathsTextFileName) # pulls from the inputDataPaths
	print('\n'.join(recPathlist))
	
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
	print('WARNING: code has been optimized for Windows and Linux operating systems only.\nYou may need to set your system parameters at the top of .py file yourself.')


for rawDataPath in recPathlist:
	# data = loadRawOEdir(rawDataPath)
	print(rawDataPath)
	TTL_timeStamps, TTL_values = loadTTLdata(str(rawDataPath))

	if not np.argwhere(TTL_values ==1).shape == (0,1):
		firstHighTTLind = np.argwhere(TTL_values ==1)[0]
		lastHighTTLind = np.argwhere(TTL_values ==1)[-1] 
		TTLsFound = True
	else:
		print('\nWARNING: skipping this events file due to no TTL values found: ' + str(rawDataPath))
		TTLsFound = False

	if TTLsFound == True:

		if TTL_values[FIRST_TTL_IND] == 1:
			currentValue = 1
			FIRST_TTL_IND = FIRST_TTL_IND + 2 # skip two cuz the next one is the low value and we want the next high value
			print('skipping the first TTL high pulse cuz it is probably truncated')
		elif TTL_values[FIRST_TTL_IND] == 0:
			currentValue = 0

		TTL_values = TTL_values[FIRST_TTL_IND:LAST_TTL_IND] # truncate data
		TTL_timeStamps = np.divide(TTL_timeStamps[FIRST_TTL_IND:LAST_TTL_IND],SAMPLE_RATE_HZ)  # yeilds all time stamps in secs
		numTimeStamps_TTL = TTL_values.shape[0]
		

		

		print('num time stamps: ' + str(numTimeStamps_TTL))
		print('TTL_values shape' + str(TTL_values.shape))
		print('current value: ' +  str(currentValue))

		TTL_squ_ts = np.full((numTimeStamps_TTL * NUM_PLATEU_PTS), np.nan) 	
		TTL_squ_vals = np.full((numTimeStamps_TTL * NUM_PLATEU_PTS), np.nan)
		sine_ts = np.full(TTL_squ_ts.shape, np.nan)
		sine_vals = np.full(TTL_squ_ts.shape, np.nan) 

		
		numTimeStamps_squ = sine_ts.shape[0] # TEST FIX

		# print(np.isnan(sine_ts).any())
		# print(np.isnan(sine_vals).any())

		for ttlInd in range(0,numTimeStamps_TTL): # TEST FIX
			currentTimeStamp = TTL_timeStamps[ttlInd]
			
			if ttlInd < numTimeStamps_TTL-1:
			# if ttlInd < numTimeStamps: # TEST FIX
				nextTimeStamp = TTL_timeStamps[ttlInd+1]
			newTimeStamps = np.linspace(currentTimeStamp,nextTimeStamp,NUM_PLATEU_PTS, endpoint=False)
			firstNewInd = ttlInd*NUM_PLATEU_PTS
			lastNewInd = (ttlInd+1)*(NUM_PLATEU_PTS)


			TTL_squ_ts[firstNewInd:lastNewInd] = newTimeStamps 
			TTL_squ_vals[firstNewInd:lastNewInd] = currentValue
			if currentValue == 0:
				currentValue = 1
			elif currentValue ==1:
				currentValue = 0

			### calculate sine wave from ttl
			# TTL_values_slice = TTL_values[int(firstHighTTLind):int(lastHighTTLind)]
			# TTL_values_slice = TTL_values[ttlInd:ttlInd + NUM_PLATEU_PTS]
			TTL_values_slice = TTL_squ_vals[firstNewInd:lastNewInd]

			# if ttlInd + NUM_PLATEU_PTS < numTimeStamps:
			if lastNewInd <= numTimeStamps_squ:
				# freq = (np.argwhere(TTL_values_slice.shape == 1).shape[0] -1) / (TTL_timeStamps[lastHighTTLind-2] - TTL_timeStamps[firstHighTTLind]) ### TEST FIX
				# freq = (np.argwhere(TTL_values_slice.shape == 1).shape[0] -1) / (TTL_timeStamps[ttlInd] - TTL_timeStamps[ttlInd+NUM_PLATEU_PTS])
				freq = (np.argwhere(TTL_values_slice.shape == 1).shape[0] -1) / (TTL_squ_ts[lastNewInd-1] - TTL_squ_ts[firstNewInd])
				# sine_vals = 0.5 + 0.5 * np.sin(2 * np.pi * freq * TTL_squ_ts +  PHASE_OFFSET / np.pi)
				# sine_ts = np.linspace(np.min(TTL_timeStamps), np.max(TTL_timeStamps), TTL_timeStamps.size * 10) # from just TTLs 

				sine_ts[firstNewInd:lastNewInd] = np.linspace(TTL_squ_ts[firstNewInd],TTL_squ_ts[lastNewInd-1], NUM_PLATEU_PTS) # from reconstructed square wave
				sine_vals[firstNewInd:lastNewInd] = 0.5 + 0.5 * np.sin((2 * np.pi * freq * sine_ts[firstNewInd:lastNewInd]) +  PHASE_OFFSET / np.pi)

			

		pdb.set_trace()
			# print('calculated freq from TTL peaks to be: ' + str(freq)) # WIP replace with mean and variance measures

		plt.figure()
		plt.plot(TTL_timeStamps, TTL_values, 'o')
		plt.plot(TTL_squ_ts, TTL_squ_vals, 'o')
		plt.plot(sine_ts, sine_vals)
		plt.legend(['TTL transitions', 'TTL square', 'reconstructed sine stim'], loc='best')
		plt.show()

		ts = time.time()
		timeStr = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d__%H_%M')
		recordingName = getRecordingNameFromPath(rawDataPath)
		print('naming the file:' + recordingName + timeStr + '.pdf')

		pdfName = Path(SUMMARY_OUTPUT_PATH).joinpath(recordingName + timeStr + '.pdf')
		# with PdfPages(pdfName) as pdf_out:
		# 	for chInd in range(FIRST_CH-1, LAST_CH):
		# 		# f = plt.figure(figsize=(200,200)) # REDUCE THESE NUMBERS IF YOU GET SEGMENTATION/MEMORY ERRORS!
		# 		f = plt.figure() 
		# 		plt.subplot(411)
		# 		plt.plot(data[:,chInd])
		# 		plt.subplot(412)
		# 		Pxx, freqs = plt.psd(data[:,chInd], 2048, SAMPLE_RATE_HZ)
				
		# 		# freqs, Pxx  = plt.psd(data[:,chInd], 2048, SAMPLE_RATE_HZ)
		# 		plt.subplot(413)
		# 		# extraTicks = [-50, -25, -10, -5, -2, -1, 0,1]
		# 		# plt.semilogy(freqs[0:HIGH_FREQ_CUTOFF],Pxx[0:HIGH_FREQ_CUTOFF], yticks = list(plt.yticks()[0]) + extraTicks)
		# 		plt.semilogy(freqs[0:HIGH_FREQ_CUTOFF],Pxx[0:HIGH_FREQ_CUTOFF])

		# 		### testing
		# 		plt.subplot(414)
		# 		plt.magnitude_spectrum(data[:,chInd], Fs = SAMPLE_RATE_HZ, scale = 'dB')

		# 		if INSPECT_EACH_PLOT ==1:
		# 			print('close the plot window to advance to the next channel\nOR set INSPECT_EACH_PLOT to 0')
		# 			plt.show()
				
		# 		if chInd == LAST_CH:
		# 			f.suptitle(RAW_DATA_PATH, fontsize= titleFontSize)
		# 		print('saving PSD for channel: ' + str(chInd+1) + ' of ' + str(numChannelsPlotted))
		# 		if SAVE_ALL_TO_PDF == 1:
		# 			pdf_out.savefig(f)
		# 		plt.close(f)
		# 		f.clf() # might help with out of memory error
				
		# print('finished: ' + str(pdfName) + ' saved to cwd')
		# if OPEN_EACH_PDF == 1:
		# 	webbrowser.get(CHROME_PATH).open(pdfName)
		# else:
		# 	print('PSD results saved to cwd')	
	



	

	




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