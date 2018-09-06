import matlab.engine 
import glob
import os.path
import os
# import bashrcMod # don't forget that importing a file executes it!
import pdb
import time
import datetime
import sys
# import reinstateCondaPathToBashrc
# import removeCondaPathFromBashrc

### CHECK THESE PARAMETERS!!!!!!
SPIKESORTER = 'mountainSort' # use 'kiloSort', mountainSort', 'none', or 'both'
SORT_ALL_IN_SUBDIRS = 'no' # 'yes' or 'no' # yes generates rawDataPaths.txt and recordingIDsToRun.txt; no requires a pre-populated recordingIDsToRun.txt file and a rawDataPaths to match it in the cwd of this script 
PATH_a_Analysis = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Analysis/'
PATH_rawParent = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Data/raw' #use os.getcwd() to process data in child directories of the cwd
PATH_spikeSortedParent = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Data/spikesorted'

### kiloSort specific parameters:
KILOSORT_PARAMS_DIR = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Analysis/spikeSorting/kiloSorting/' 
KILOSORT_RESULTS_DIR = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Data/spikesorting/kiloSorted'
POSTHOC_MERGE = 'no' # 'yes' or 'no' yes merges together clusters kilosort suspects as being the same etc
RUN_PHY = 'no' # 'yes' or 'no' for opening the phy results in a loop

### mountainSort specific parameters
MNTSRT_PARAMS_DIR = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Analysis/spikeSorting/mountainSorting/'
VIEW_MTNSRT_RESULTS = 'no' # 'yes' or 'no'

BASHRC_DIR = '/home/matthew'
callingDir = os.getcwd()
### TO DO: optional integration with Spikalyze to plot summary data for validation
#### visualize results
#PLOT_RASTERS = 'no' # 'yes' or 'no'
#
#### save results
SAVED_FIG_TYPE = 'pdf' #TO DO: jpg

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

def moveKiloResultsOutOfRawDir():
    FILES_TO_MOVE_PATH = os.getcwd() # WARNING: probably doesn't support changes to this at present!
    ts = time.time()
    timeStr = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d__%H_%M')
    newSortingDirName = os.path.basename(os.getcwd()) + timeStr
    fullDestinationPath = KILOSORT_RESULTS_DIR + '/' + newSortingDirName
    if not os.path.exists(fullDestinationPath):
        os.makedirs(fullDestinationPath)
    os.system('mv ' + FILES_TO_MOVE_PATH + '/kilosort.dat ' + fullDestinationPath)
    os.system('mv ' + FILES_TO_MOVE_PATH + '/*.npy ' + fullDestinationPath)
    os.system('mv ' + FILES_TO_MOVE_PATH + '/params.py ' + fullDestinationPath)
    os.system('mv ' + FILES_TO_MOVE_PATH + '/rez.mat ' + fullDestinationPath)
    
def runKiloSort():
    os.chdir(PATH_rawParent)
    print('starting matlab engine for kilosorting')
    eng = matlab.engine.start_matlab()
    print('engine running')
    with open('recordingIDsToRun.txt') as recIDs, open('rawDataPaths.txt') as rawPaths:
        recIDs = recIDs.read().splitlines()
        rawPaths = rawPaths.read().splitlines()
        print('\nPreparing to sort: ')
        print('\n'.join(recIDs))
        print('\namong the following (possible) raw data in subdirs of the cwd: ')
        print('\n'.join(rawPaths))
        
        for recID in recIDs:
            for rawPath in rawPaths:
                if recID in rawPath: # nice implicit string search here

                    ### format output results directory
                    try:
                        os.chdir(KILOSORT_RESULTS_DIR)
                    except:
                        print('\n dir not found, therefore creating dir: ' + KILOSORT_RESULTS_DIR)
                        os.system('mkdir ' + KILOSORT_RESULTS_DIR) # UNTESTED: can this make a nested directory?                    
                    newKiloSaveDir = KILOSORT_RESULTS_DIR + '/' + recID
                    
                    
                    if POSTHOC_MERGE.lower() == 'yes':
                        newKiloSaveDir = newKiloSaveDir + 'posthocMerged'
                    elif POSTHOC_MERGE.lower() != ('yes' and 'no'):
                        sys.exit("\n\nPOSTHOC_MERGE param invalid! must be yes or no\n\n")

                    if not os.path.exists(newKiloSaveDir):
                        os.chdir(rawPath)

                        ### rename first OE file if it has the weird capitalization
                        if not os.path.exists('./100_CH1.continuous'):
                            os.system('cp 100_ch1.continuous ./100_CH1.continuous') # alternatively: os.system("rename 's/ch/CH/' *.continuous")
                            print('copied continuous OE file from ch1 to CH1')
                        else:
                            print('data already named correctly; no need to rename OE file from ch1 to CH1')
                        

                        ########## TO DO: replace this with a functional call after it's tested
                        # read probe type from probeType.txt
                        with open('probeType.txt') as f:
                            probeType = f.readline()
                            probeType = probeType.strip()
                        print('using chanmap for probe: ' + probeType)

                        ### run the sorting
                        print('now kilosorting: ' + rawPath + '\n\nJust be patient now... :)')
                        os.chdir(KILOSORT_PARAMS_DIR)
                        eng.addpath(os.getcwd())
                        eng.workspace['RAW_DATA_ROOT_PATH'] = rawPath.strip()

                        if probeType.lower() == '64f':
                            if POSTHOC_MERGE.lower() == 'yes':
                                eng.master_file_2_9_18posthocMerged_64F(nargout = 0)
                            elif POSTHOC_MERGE.lower() == 'no':
                                eng.master_file_6_12_18_64F(nargout = 0)
                        elif probeType.lower() == '64g':
                            if POSTHOC_MERGE.lower() == 'yes':
                                eng.master_file_2_9_18posthocMerged_64G(nargout = 0)
                            elif POSTHOC_MERGE.lower() == 'no':
                                eng.master_file_6_12_18_64G(nargout = 0)                                
                        else:
                            sys.exit('ERROR: chanmap unavailable for that probe type or probeType.txt missing from raw data dir')

                        os.chdir(rawPath)
                        moveKiloResultsOutOfRawDir()

                    else:
                        print('\nsorting results already exist for: ' + newKiloSaveDir)
                        print('rename, move, or delete this dir to rerun the analysis of this dataset')
    print('\nfinished with kiloSorting: ' + recID)

def usePhyForKilosortManualValidation():
    os.chdir(KILOSORT_RESULTS_DIR)
    listOfResultDirs = glob.glob('./*')
    listOfResultDirs.sort()
    listOfResultDirs = [s.replace('./', '') for s in listOfResultDirs] 
    os.chdir(PATH_a_Analysis)

    with open('recordingIDsToRun.txt') as recIDs:
        recIDs = recIDs.read().splitlines()
        recIDs.sort()
        print('\nPreparing to open Phy for manual stage of sorting for the following: ')
        print('\n'.join(recIDs))
        print('\namong the following raw paths in rawDataPaths.txt: ')
        print('\n'.join(listOfResultDirs))
        print('\nuse SORT_ALL_IN_SUBDIRS = yes to populate rawDataPaths.txt with all raw data in subdirs')
        
        for recID in recIDs:
            for kiloOutputDir in listOfResultDirs:
                
                if recID in kiloOutputDir: # nice implicit (sub)string search here
                    print('recID:' + recID)
                    print('kilosort output dir: ' + kiloOutputDir)
                    outputDirFullPath = KILOSORT_RESULTS_DIR + '/' + kiloOutputDir
                    os.chdir(outputDirFullPath)

                    ### run phy command
                    try:
                        os.system("bash -c 'source activate phy; phy template-gui ./params.py'")
                    except:
                        sys.exit('ERROR: conda phy environment activation failed... (omit "source" in "source activate phy" for windows)')

def saveRawDataPathsToTxt(pathList):
    with open('rawDataPaths.txt', 'w') as pathFile:
        for path in pathList:
            # cwd = os.getcwd()
            # path = path.replace(path[:2],cwd + '/')
            path = path.replace('/Continuous_Data.openephys','')
            pathFile.write("%s\n" % path)

def generateRecordingIDsToRun():
    recIDlist = []
    with open('rawDataPaths.txt') as paths:
        for recording in paths:
            recording = recording.strip() # remove new line
            slashInd = recording.rfind('/')
            recID = recording[slashInd+1:]
            recIDlist.append(recID)

    with open('recordingIDsToRun.txt', 'w') as recIDlistFile:
        for recID in recIDlist:
            recIDlistFile.write("%s\n" % recID)
    return recIDlist

### recursive search 
def recursiveSearch(exten): 
    pathList = []
#    for path in glob.iglob(PATH_rawParent + '/**/*' + exten, recursive=True):
    for path in glob.iglob(os.getcwd() + '/**/*' + exten, recursive=True):
        # print(path)
        pathList.append(path)
    return pathList

def convertOEtoMDA():
    os.chdir(callingDir)
    with open('recordingIDsToRun.txt') as recIDs, open('rawDataPaths.txt') as rawPaths:
        recIDs = recIDs.read().splitlines()
        rawPaths = rawPaths.read().splitlines()
        engineRunning = False 
        for recID in recIDs:
            for rawPath in rawPaths:
                if recID in rawPath:
#        for ind in range(0,len(recIDs)):
                    os.chdir(PATH_spikeSortedParent)
                    saveLocation = PATH_spikeSortedParent +'/mountainSorted/' + recID.strip() 
                    if not os.path.exists(saveLocation):
                        if engineRunning == False:
                            eng = matlab.engine.start_matlab() ### only start the matlab engine if it's not already running
                            engineRunning = True
                        eng.addpath(MNTSRT_PARAMS_DIR)
                        
                        os.makedirs(saveLocation)
                        print('created mountainSort dir: ' + saveLocation)
                        print('converting open ephys to MDA format')
                        os.chdir(MNTSRT_PARAMS_DIR) ### not sure if this is necessary ##########################
                        eng.workspace['BASE_FILE_PATH'] = rawPath.strip()
                        eng.workspace['RAW_MDA_OUTPUT_NAME'] = recID.strip()
                        eng.workspace['RAW_MDA_SAVE_PATH'] = saveLocation
                        proteType, numChannels = getProbeType(rawPath)
                        eng.workspace['CHAN_NUM'] = numChannels
                        eng.convertToMDA_9_4_18(eng.workspace['BASE_FILE_PATH'], 
                            eng.workspace['RAW_MDA_OUTPUT_NAME'], eng.workspace['RAW_MDA_SAVE_PATH'], 
                            eng.workspace['CHAN_NUM'])
                    else:
                        print(saveLocation + ' already exists --> skipping this conversion')
                    print('to reconvert raw data to MDA format: delete/move/rename the dirs in [project]a_Data/spikesorted/mountainSorted')
    os.chdir(callingDir) #### may not be necessary ###############################

def getRawMDApaths():
    os.chdir(PATH_spikeSortedParent)
    os.chdir('mountainSorted')
    # print('cd' + os.getcwd())
    MDA_parentPath = os.getcwd() + '/'

    pathList = recursiveSearch('raw.mda')
    rawMDApaths = []
    for MDApath in pathList:
        MDApath = MDApath.replace('./', MDA_parentPath)
        MDApath = MDApath.strip()
        rawMDApaths.append(MDApath)
    return rawMDApaths

def getProbeGeometryCsv(probeType):
    if probeType.lower() == '64f':
        probeGeometryCsvFullpath = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Analysis/spikeSorting/mountainSorting/geom64F.csv'
    elif probeType.lower() == '64g':
        probeGeometryCsvFullpath = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Analysis/spikeSorting/mountainSorting/geom64G.csv'
    elif probeType.lower() == 'ADD NEW PROBE TYPES HERE!': ### <-- add new ones here
        probeGeometryCsvFullpath = 'ADD NEW PROBE FULL PATH (INCLUDING geomX.csv) HERE'
    else:
        sys.exit('ERROR: probeType.txt does not match a geometry file in spikesorting/mountainSorting')
    print('using probe geometry path:\n{0}'.format(probeGeometryCsvFullpath))
    return probeGeometryCsvFullpath

### read probe type from user populated text file in raw data dir
def getProbeType(path):
    print('path: {0}'.format(path))
    cwd = os.getcwd()      
    os.chdir(path)
    if not os.path.exists('probeType.txt'):
        os.chdir(callingDir)
        with open('rawDataPaths.txt') as rawPathsFile:
            rawPaths = rawPathsFile.read().splitlines()
            correctedProbeTypePath = [rawPath for rawPath in rawPaths if os.path.basename(path) in rawPath]
            if correctedProbeTypePath:
                print('WARNING: probeType.txt not found in cwd\nChecking the raw data dir with the same recording ID corrected: {0}\n'.format(correctedProbeTypePath))
                os.chdir(correctedProbeTypePath[0])
            elif not correctedProbeTypePath:
                sys.exit('ERROR: unable to find a probeType.txt giving the type of the probe used for this recording')
    with open('probeType.txt') as f:
        probeType = f.readline()
        probeType = probeType.strip()
    print('using chanmap for probe: ' + probeType)

    if '64' in probeType: ## WIP TEST THIS:
        numChannels = 64
    elif '128' in probeType:
        numChannels = 128
    elif '256' in probeType:
        numChannels = 256
    else:
        sys.exit('probeType does not contain an expected number of channels\ncheck the probeType.txt in the raw data dir')
    print('number of channels: {0}'.format(numChannels))    
    os.chdir(cwd)
    return probeType, numChannels

def runMntSrt(rawMDApaths):

    removeCondaPathFromBashrc(BASHRC_DIR) # anaconda path in bashrc breaks mountainSort but is necessary for phy (kilosort visualization)
    os.chdir(callingDir)
    with open('recordingIDsToRun.txt') as recIDs, open('rawDataPaths.txt') as rawPaths:
        recIDs = recIDs.read().splitlines()
        os.chdir(MNTSRT_PARAMS_DIR)
        print('running mountainSort commands from' + os.getcwd() + '\n' + 'running mountain sort now')
        print('reminder: are you passing the correct probe geom.csv file based on the probe type in probeType.txt in the raw data dir?')
        for recID in recIDs:
            for rawMDApath in rawMDApaths:
                if recID in rawMDApath: 
                    print('\nspike sorting: \n' + rawMDApath)
                    rawMDAdir = os.path.dirname(rawMDApath)
                    
                    # PROBE_GEOM_CSV_NAME = 'geom64F' # DO NOT ADD THE .csv file ext  # <- old version
                    probeType, numChannels = getProbeType(rawMDAdir) ########################################## WRONG DIRECTORY!!!!!!!!!!!!!1
                    probeGeometryCsvFullpath = getProbeGeometryCsv(probeType)
                    mntSrt_run_options = 'mlp-run mountainsort3.mlp sort --raw=' + rawMDApath + ' --geom=' + probeGeometryCsvFullpath + ' --firings_out=firings_curated.mda --firings_original_out=firings_original.mda --filt_out=filt.mda --pre_out=pre.mda --cluster_metrics_out=metrics.json --_params=params.json --curate=true' 
                    os.system(mntSrt_run_options)
                    mntSort_saveDir = os.path.dirname(rawMDApath) + '/'
                    os.system('mv *.mda metrics.json ' + mntSort_saveDir)
                    os.system('cp *.csv *.mlp *params.json ' + mntSort_saveDir)
    reinstateCondaPathInBashrc(BASHRC_DIR)    

# WIP: 
def viewMntSrt():
    os.chdir(callingDir)
    rawMDApaths = getRawMDApaths()
    with open('recordingIDsToRun.txt') as recIDs: # maybe not? --> WIP: change to recordingIDsToView.txt
        recIDs = recIDs.read().splitlines()
        for recID in recIDs:
            for rawMDApath in rawMDApaths:
                if recID in rawMDApath: 
                    print('\nspike sorting: \n' + rawMDApath) 
                    os.chdir(PATH_a_Analysis)
                    print('cd into' + os.getcwd() + '\n' + 'running mountain sort now') 
                    rawMDAdir = os.path.dirname(rawMDApath)
                    prePath = rawMDAdir + '/' + recID
                    print(prePath)
                    os.chdir(prePath)
                    print('cd into ' + os.getcwd())
                    print(rawMDApath)
                    # '.mlp sort --raw=' + rawMDApath + ' --geom=geom.csv --firings_out=' + prePath + ' _firings_curated.mda --firings_original_out=' + prePath + '_firings_original.mda --filt_out=' + prePath + '_filt.mda --pre_out=' + prePath + '_pre.mda --cluster_metrics_out=' + prePath + '_metrics.json --_params=' + PATH_a_Analysis + '/params.json --curate=true' 
                    mntSrt_view_options = 'mountainview --raw=' + rawMDApath + ' --filt=filt.mda --pre=pre.mda --geom=geom.csv --firings=firings_curated.mda --cluster_metrics=metrics.json'                    
                    os.system(mntSrt_view_options)

### MAIN CODE: 

### CHOOSING DATA: output = recordingIDsToRun.txt
## get all data paths for subdirs in ./
# pathList = []    
if SORT_ALL_IN_SUBDIRS == 'yes':
    print('sorting all in subdirs from the location of this script')
    exten = '.openephys' # make this a UNIQUE FILE within the raw open ephys dir
    pathList = recursiveSearch(exten)
    print('\nraw paths found:')
    print('\n'.join(pathList))
    saveRawDataPathsToTxt(pathList)
    ### comment out and modify recordingIDsToRun.txt manually to run on a subset of the raw data in the child dirs of this directory
    recIDlist = generateRecordingIDsToRun()
    print('\nraw data directory names:')
    print('\n'.join(recIDlist))
    print('\nremove some of the above from recordingIDsToRun.txt to run on a subset of the child dir data\n')
elif SORT_ALL_IN_SUBDIRS == 'no':
    if os.path.exists('recordingIDsToRun.txt'):
        print('sorting based on contents of recordingIDsToRun.txt')
    else:
        print('error: invalid selection for SORT_ALL_IN_SUBDIRS param')


### SPIKE SORTING
print('using spikesorter: ' + SPIKESORTER) ########## WIP: programmatically solve conda bashrc conflict b/w kilosort and mntSort
if SPIKESORTER == 'kiloSort':
    print('spikesorting with kilosort')
    runKiloSort()
elif SPIKESORTER == 'mountainSort':
    convertOEtoMDA() # note: will only run if the mountainsorted result dir does not exist
    ## searches through all MDA files in this project dir tree for matches with recordingIDsToRun.txt
    pathList = [] # REQUIRED FOR getRawMDApaths() !!!
    rawMDApaths = getRawMDApaths()
    print('\nraw MDA paths: \n')
    print('\n'.join(rawMDApaths))
    runMntSrt(rawMDApaths) # WIP: make param for mntSrt_run_options 
    print('finished spikesorting')
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

### only applies for kilosort
if RUN_PHY.lower() == 'yes':
    print('looping over all kilosort outputs in phy for manual stage:')
    usePhyForKilosortManualValidation()

if VIEW_MTNSRT_RESULTS.lower() == 'yes':
    viewMntSrt()


### To Do: spikalyze integration to give summary plots