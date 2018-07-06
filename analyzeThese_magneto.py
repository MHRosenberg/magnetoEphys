import matlab.engine # conflicts with conda therefore imported below
import glob

import os.path
import os
import bashrcMod # don't forget that importing a file executes it!
# import reinstateCondaPathToBashrc
# import removeCondaPathFromBashrc

### temp parameters DO NOT POPULATE!
# rawMDApath = ''
# prePath = ''
# PATH_a_Analysis = ''
# rawMDApaths = ''
###

def removeCondaPathFromBashrc(BASHRC_DIR):
	commentOutCondaPathInBashrc = "sed -i 's/export PATH=\"\/home\/matthew\/miniconda3\/bin:\$PATH\"/\#export PATH=\"\/home\/matthew\/miniconda3\/bin:\$PATH\"/g' ./.bashrc"
	os.chdir(BASHRC_DIR)
	os.system(commentOutCondaPathInBashrc)
	os.system('. ~/.bashrc')
	print(BASHRC_DIR + '.bashrc modified to remove conda path as follows: ')
	print(commentOutCondaPathInBashrc)

def reinstateCondaPathInBashrc(BASHRC_DIR):
	reinstateCondaPathInBashrc = "sed -i 's/\#export PATH=\"\/home\/matthew\/miniconda3\/bin:\$PATH\"/export PATH=\"\/home\/matthew\/miniconda3\/bin:\$PATH\"/g' ./.bashrc"
	os.system(reinstateCondaPathInBashrc)
	os.system('. ~/.bashrc')
	print(BASHRC_DIR + '.bashrc modified to reinstate conda path as follows: ')
	print(reinstateCondaPathInBashrc)

### recursive search 
### stores results in pathList (don't forget to clear results in pathList before reusing)
def recursiveSearch(exten): 
	# The top argument for walk. Th
	topdir = '.'
	# Start the walk
	os.path.walk(topdir, step, exten)
	return pathList

def step(ext, dirname, names):
	ext = ext.lower()
	for name in names:
		if name.lower().endswith(ext):
			tmpPath = os.path.join(dirname, name) 
			pathList.append(tmpPath)
### end method group

def saveRawDataPathsToTxt(pathList):
	pathFile = open('rawDataPaths.txt', 'w')
	for path in pathList:
		cwd = os.getcwd()
		# path = path[]
		path = path.replace(path[:2],cwd + '/')
		path = path.replace('/Continuous_Data.openephys','')
		pathFile.write("%s\n" % path)
	pathFile.close()     

def generateRecordingIDsToRun():
	recIDlist = []
	with open('rawDataPaths.txt') as paths:
		for recording in paths:
			recording = recording.strip() # remove new line
			slashInd = recording.rfind('/')
			recID = recording[slashInd+1:]
			recIDlist.append(recID)

	recIDlistFile = open('recordingIDsToRun.txt', 'w')
	for recID in recIDlist:
		recIDlistFile.write("%s\n" % recID)
	recIDlistFile.close()
	return recIDlist

def convertOEtoMDA():
	### mountain sort loop
	eng = matlab.engine.start_matlab()
	eng.addpath(PATH_a_Analysis)
	os.chdir(PATH_rawParent)
	with open('recordingIDsToRun.txt') as recIDs, open('rawDataPaths.txt') as rawPaths:
		recIDs = recIDs.read().splitlines()
		rawPaths = rawPaths.read().splitlines()
		for ind in range(0,len(recIDs)):
			saveLocation = '../../spikesorting/mountainSorted/' + recIDs[ind].strip() #+
			if not os.path.exists(saveLocation):
				os.makedirs(saveLocation)
				
				print('created mountainSort dir: ' + saveLocation)
				print('converting open ephys to MDA format')

				eng.workspace['BASE_FILE_PATH'] = rawPaths[ind].strip() 
				eng.workspace['RAW_MDA_OUTPUT_NAME'] = recIDs[ind].strip()
				eng.workspace['RAW_MDA_SAVE_PATH'] = saveLocation
				eng.workspace['CHAN_NUM'] = CHAN_NUMS
				eng.convertToMDA_1_30_18(eng.workspace['BASE_FILE_PATH'], 
					eng.workspace['RAW_MDA_OUTPUT_NAME'], eng.workspace['RAW_MDA_SAVE_PATH'], 
					eng.workspace['CHAN_NUM'])
			else:
				print(saveLocation + ' already exists --> skipping this conversion')

def getRawMDApaths():
	os.chdir(PATH_rawParent)
	os.chdir('../../spikesorting/mountainSorted/')
	# print('cd' + os.getcwd())
	MDA_parentPath = os.getcwd() + '/'

	pathList = recursiveSearch('raw.mda')
	rawMDApaths = []
	for MDApath in pathList:
		MDApath = MDApath.replace('./', MDA_parentPath)
		MDApath = MDApath.strip()
		rawMDApaths.append(MDApath)
	return rawMDApaths

def runMntSrt(rawMDApaths):
	# removeCondaPathFromBashrc() # anaconda path in bashrc breaks mountainSort but is necessary for phy (kilosort visualization)
	os.chdir(PATH_rawParent)
	print('cd into ' + os.getcwd())
	with open('recordingIDsToRun.txt') as recIDs:
		recIDs = recIDs.read().splitlines()
		os.chdir('../../spikesorting/a_paramsNconfigs/mountainSorting')
		print('running mountainSort commands from' + os.getcwd() + '\n' + 'running mountain sort now')
		print('reminder: are you passing the correct probe geom.csv file based on the probe type?')
		for recID in recIDs:
			for rawMDApath in rawMDApaths:
				if recID in rawMDApath: 
					print('\nspike sorting: \n' + rawMDApath) 
					rawMDAdir = os.path.dirname(rawMDApath)

					mntSrt_run_options = 'mlp-run mountainsort3.mlp sort --raw=' + rawMDApath + ' --geom=' + PROBE_GEOM_CSV_NAME + '.csv --firings_out=firings_curated.mda --firings_original_out=firings_original.mda --filt_out=filt.mda --pre_out=pre.mda --cluster_metrics_out=metrics.json --_params=params.json --curate=true' 
					os.system(mntSrt_run_options)
					mntSort_saveDir = os.path.dirname(rawMDApath) + '/'
					os.system('mv *.mda metrics.json ' + mntSort_saveDir)
					os.system('cp *.csv *.mlp *params.json ' + mntSort_saveDir)

	# reinstateCondaPathInBashrc()    

## WIP: 
# def viewMntSrt(recIDs, mntSrt_run_options):
#     os.chdir(PATH_rawParent)
#     with open('recordingIDsToRun.txt') as recIDs: # WIP: change to recordingIDsToView.txt
#         recIDs = recIDs.read().splitlines()
#         for recID in recIDs:
#             for rawMDApath in rawMDApaths:
#                 if recID in rawMDApath: 
#                     print('\nspike sorting: \n' + rawMDApath) 
#                     os.chdir(PATH_a_Analysis)
#                     print('cd into' + os.getcwd() + '\n' + 'running mountain sort now') 
#                     rawMDAdir = os.path.dirname(rawMDApath)
#                     prePath = rawMDAdir + '/' + recID
#                     print(prePath)
#                     os.chdir(prePath)
#                     print('cd into ' + os.getcwd())
#                     print(rawMDApath)
#                     # '.mlp sort --raw=' + rawMDApath + ' --geom=geom.csv --firings_out=' + prePath + ' _firings_curated.mda --firings_original_out=' + prePath + '_firings_original.mda --filt_out=' + prePath + '_filt.mda --pre_out=' + prePath + '_pre.mda --cluster_metrics_out=' + prePath + '_metrics.json --_params=' + PATH_a_Analysis + '/params.json --curate=true' 
#                     mntSrt_view_options = 'mountainview --raw=' + rawMDApath + ' --filt=filt.mda --pre=pre.mda --geom=geom.csv --firings=firings_curated.mda --cluster_metrics=metrics.json'                    
#                     os.system(mntSrt_view_options)

def runKiloSort():
	os.chdir(PATH_rawParent)
	eng = matlab.engine.start_matlab()
	with open('recordingIDsToRun.txt') as recIDs, open('rawDataPaths.txt') as rawPaths:
		recIDs = recIDs.read().splitlines()
		rawPaths = rawPaths.read().splitlines()
		print('\nPreparing to sort: ')
		print('\n'.join(recIDs))
		print('\namong the following (possible) raw data in subdirs of the cwd: ')
		print('\n'.join(rawPaths))
		
		for recID in recIDs:
			for rawPath in rawPaths:
				if recID in rawPath:
					os.chdir(PATH_rawParent)
					try:
						os.chdir('../../spikesorting/kiloSorted')
					except:
						print('\n dir not found, therefore creating dir: ../../spikesorting/kiloSorted')
						os.chdir('../../')
						os.system('mkdir spikesorting/kiloSorted') # UNTESTED: can this make a nested directory?
					kiloDir = os.getcwd()
					kiloSaveDir = kiloDir + '/' + recID
					if not os.path.exists(kiloSaveDir):
						os.chdir(rawPath)

						if not os.path.exists('./100_CH1.continuous'):
							os.system('cp 100_ch1.continuous ./100_CH1.continuous') # alternatively: os.system("rename 's/ch/CH/' *.continuous")
							print('copied continuous OE file from ch1 to CH1')
						else:
							print('data already named correctly; no need to rename OE file from ch1 to CH1')
						os.chdir(PATH_rawParent)
						os.chdir('../../spikesorting/a_paramsNconfigs/kiloSorting/')
						eng.addpath(os.getcwd())

						if os.path.exists(PROBE_CHAN_MAP_PATH):
							eng.workspace['CHAN_MAP_MAT_PATH'] = PROBE_CHAN_MAP_PATH
						else:
							print('ERROR: ' + PROBE_CHAN_MAP_PATH + 'does not exist')
						print('reminder: have you checked that the parameters in the separate matlab files are correct? \nLook in: ' + PROBE_CHAN_MAP_PATH)

						eng.workspace['RAW_DATA_ROOT_PATH'] = rawPath.strip()

						if POSTHOC_MERGE == 'yes':
							eng.master_file_2_9_18posthocMerged(nargout = 0)
							kiloSaveDir = kiloSaveDir + 'posthocMerged'
						elif POSTHOC_MERGE == 'no':
							eng.master_file_2_5_18(nargout = 0)
						
						os.makedirs(kiloSaveDir)
						os.system('mv ' + rawPath + '/kilosort.dat ' + kiloSaveDir)
						os.system('mv ' + rawPath + '/*.npy ' + kiloSaveDir)
						os.system('mv ' + rawPath + '/params.py ' + kiloSaveDir)
						os.system('mv ' + rawPath + '/rez.mat ' + kiloSaveDir)
					else:
						print('\nsorting results already exist for: ' + recID)
						print('rename, move, or delete this dir to rerun the analysis of this dataset')
	print('\nfinished with kiloSorting')

### MAIN CODE: 

## common parameters:
CHAN_NUMS = 64
SPIKESORTER = 'mountainSort' # use 'kiloSort', mountainSort', 'none', or 'both'
SORT_ALL_IN_SUBDIRS = 'yes' # 'yes' or 'no' # no requires a pre-populated recordingIDsToRun.txt file in the cwd of this script 
PATH_a_Analysis = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Analysis'
PATH_rawParent = os.getcwd()
BASHRC_DIR = '/home/matthew'

## kiloSort specific parameters:
PROBE_CHAN_MAP_PATH = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Analysis/spikeSorting/kiloSorting/probe64F_chanMap.mat' #make this file using createChannelMapFile.m       
POSTHOC_MERGE = 'no' # 'yes' or 'no' yes merges together clusters kilosort suspects as being the same etc

## mountainSort specific parameters
PROBE_GEOM_CSV_NAME = 'geom64F' # DO NOT ADD THE .csv file ext

### CHOOSING DATA: output = recordingIDsToRun.txt
## get all data paths for subdirs in ./
pathList = []    
if SORT_ALL_IN_SUBDIRS == 'yes':
	print('sorting all in subdirs from the location of this script')
	exten = '.openephys' # make this a UNIQUE FILE within the raw open ephys dir

	recursiveSearch(exten)
	saveRawDataPathsToTxt(pathList)
	### comment out and modify recordingIDsToRun.txt manually to run on a subset of the raw data in the child dirs of this directory
	recIDlist = generateRecordingIDsToRun() 
	print('\n'.join(recIDlist))
	print('\nremove some of the above from recordingIDsToRun.txt to run on a subset of the child dir data\n')
elif SORT_ALL_IN_SUBDIRS == 'no':
  print('sorting based on contents of recordingIDsToRun.txt')
#   ## WIP check for file existence
else:
	print('error: invalid selection for SORT_ALL_IN_SUBDIRS param')


### SPIKE SORTING
print('using spikesorter: ' + SPIKESORTER) ########## WIP: programmatically solve conda bashrc conflict b/w kilosort and mntSort
if SPIKESORTER == 'kiloSort':
	print('spikesorting with kilosort')
	runKiloSort()
elif SPIKESORTER == 'mountainSort':
	### FILE TYPE CONVERSION: only needs to be run the first time
	convertOEtoMDA()
	## searches through all MDA files in this project dir tree for matches with recordingIDsToRun.txt
	pathList = [] # REQUIRED FOR getRawMDApaths() !!!
	rawMDApaths = getRawMDApaths()
	print('\nraw MDA paths: \n')
	print('\n'.join(rawMDApaths))

	runMntSrt(rawMDApaths) # WIP: make param for mntSrt_run_options 
elif SPIKESORTER == 'both':
	#### mountainSort
	### FILE TYPE CONVERSION: only needs to be run the first time
	convertOEtoMDA()
	## searches through all MDA files in this project dir tree for matches with recordingIDsToRun.txt
	pathList = []
	rawMDApaths = getRawMDApaths()
	print(rawMDApaths)
	runMntSrt(rawMDApaths) # WIP: make param for mntSrt_run_options 
	#### kiloSort
	print('spikesorting with kilosort')
	runKiloSort()
elif SPIKESORTER == 'none':
	print('you chose not to spikesort')
else:
	print('invalid sorter specified as parameter. use kiloSort or mountainSort')

### VIEW RESULTS
# reinstateCondaPathInBashrc(BASHRC_DIR)
