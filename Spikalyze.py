import numpy as np
import scipy.io as sio
import pandas as pd
import sys
import os.path
import glob

### BEWARE: specifying the file doesn't seem to always work... unclear why it does in some cases but not others --> specify the directory only!
sys.path.insert(0, '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Analysis/') 
import OpenEphys as OE

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import style
from matplotlib import pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path #, PureWindowsPath # only works on 3.4+ (not 2.7)

### the majority of the __init__ code is rom Kyu Lee's analyzeSortedSpikes.py
# path to sorting folder
# path = 'J:/ksort/20180618'

class Spikalyze:
    def __init__(self, spikeDir):
        # load stuff
        spike_times = np.load(os.path.join(spikeDir,'spike_times.npy'))
        spike_clusters = np.load(os.path.join(spikeDir,'spike_clusters.npy'))
        spike_templates = np.load(os.path.join(spikeDir,'spike_templates.npy'))
        templates = np.load(os.path.join(spikeDir,'templates.npy'))
        cluster_groups = pd.read_csv(os.path.join(spikeDir,'cluster_groups.csv'),sep='\t')
        clusters = np.array(cluster_groups['cluster_id'])
        group = np.array(cluster_groups['group'])
        channel_map = np.load(os.path.join(spikeDir,'channel_map.npy'))

        # find good clusters
        good_clusters = clusters[group=='good']

        every_channel = {}
        max_channel = {} # will hold channel w/ max amplitude
        self.output = {} # will hold spikes
        for cluster in good_clusters:
            self.output[cluster] = spike_times[spike_clusters==cluster]
            # for clusters that were merged, compute the average waveform for every channel
            # by taking a weighted average across the clusters that were merged;    
            template_list = spike_templates[spike_clusters==cluster]
            unique_template_list = np.unique(template_list)
            preaverage = np.empty([len(unique_template_list),templates.shape[1],templates.shape[2]])
            for ind in range(len(unique_template_list)):
                preaverage[ind] = templates[unique_template_list[ind]]*len(template_list==unique_template_list[ind])/len(template_list)
            every_channel[cluster] = np.average(preaverage, axis=0)
            
            # then identify the channel with the largest peak to peak amplitude (max - min)
            max_p2p_every_channel=np.max(every_channel[cluster],0)-np.min(every_channel[cluster],0)
            max_channel[cluster] = np.where(max_p2p_every_channel==max(max_p2p_every_channel))[0][0]
            
        ### output
        cellid = np.array(list(max_channel.keys()))
        masks=channel_map[list(max_channel.values())].flatten()

        sp = np.array(list(self.output.values()))
        
        self.units_by_spikesLst = []
        self.units_by_spikesLst = self.spikeDictToList()
        
        ### declare stimulus numpy arrays
        self.numTimeStamps_TTL = []
        self.TTL_squ_timesInSamples = []
        self.TTL_squ_vals = []
        self.sine_timesInSamples = []
        self.sine_vals = []
        self.frequencies = []
        self.TTL_timesInSecs = []
        self.sine_timesInSecs = []

        ### save as mat
        sio.savemat(os.path.join(spikeDir,'sp.mat'),{'cellid':cellid,'masks':masks,'sp':sp})
        ### end of Kyu's code

        ### save as numpy
        np.save('spikeData', {'cellid':cellid,'masks':masks,'sp':sp})
        
    def recursiveSearch(self, exten): 
            pathList = []
            cwd = os.getcwd()
            for path in glob.iglob(cwd + '/**/*' + exten, recursive=True):
                # print(path)
                pathList.append(path)
            return pathList
        
    
    def spikeDictToList(self):
        print('loading the following units classified as good by the spikesorting user:')
        for unit in self.output.keys():    
            print(unit)
            unitSpikeTimeLst = []
            for spikeTime in self.output[unit][:][:]:
                unitSpikeTimeLst.append(int(spikeTime))
            self.units_by_spikesLst.append(unitSpikeTimeLst)
        return self.units_by_spikesLst
    
    def getRawDataDir(self):
        RAW_DATA_SUBDIR = 'a_Data/raw'
        print('\nsearching **/a_Data/raw/ for raw data dir with a similar name to the spikesorting dir')
        print('\ncwd of script: ' +  os.getcwd())
        try:
            if not RAW_DATA_SUBDIR in os.getcwd():
                try:
                    print('searching for raw data in: ../../../raw/')
                    os.chdir('../../../raw/')
                except:
                    try:
                        print('searching for raw data in: ../../../../raw/')
                        os.chdir('../../../../raw/')
                    except:
                        try:
                            print('searching for raw data in: ../../raw/')
                            os.chdir('../../raw/')
                        except:
                            sys.exit('ERROR: unable to locate raw data dir in reasonable number of parent directories')
            SEARCH_FOR_THIS = '.events' 
            pathLst = self.recursiveSearch(SEARCH_FOR_THIS)
            for pathInd, path in enumerate(pathLst): 
                pathLst[pathInd] = os.path.dirname(path) ### strip away the file name
            pathLst = list(set(pathLst)) ### remove duplicates
            pathLst.sort() ### alphabetize
#                [print(path) for path in pathLst] ### display found raw data dirs
            matchFound = False
            for path in pathLst:
                if os.path.basename(path) in SPIKE_DIR:
                    if matchFound == True:
                        sys.exit('ERROR: more than one .events file found for this data set... there should only be one\nduplicate: ' + path)
                    matchFound = True
                    rawDataDir = path
                    print('\nUsing TTLs from .events file in: ' + rawDataDir)
        except:
            sys.exit('ERROR: unable to load ttl data\nCheck that the data exists in **/a_data/raw/ or specific the directory as a parameter to loadTTLdata')
        return str(rawDataDir)

    def loadTTLdata(self, *RAW_DATA_DIR): # from the .events file
        if RAW_DATA_DIR:
            print('raw data dir specified by user: ' + RAW_DATA_DIR)
            rawDataDir = RAW_DATA_DIR
        else:
            rawDataDir = self.getRawDataDir()
        event_path = str(Path(rawDataDir).joinpath('all_channels.events'))
        print('\nloading: ' + event_path)
        eventsDict = OE.load(event_path) ### load TTL dict file
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
    
    def reconstructSineFromTTLs(self, TTL_values, TTL_times):
        ### I: remove incomplete periods from the front and back of the recording
        ###        i. truncate start
        ### burn extra values cuz there's still something weird...
        TTL_values = TTL_values[2::] # skip two cuz the next one is the low value and we want the next high value
        TTL_times = TTL_times[2::]

        ###            a. TTL values begin high
        if TTL_values[0] == 1:
            TTL_values = TTL_values[2::] # skip two cuz the next one is the low value and we want the next high value
            TTL_times = TTL_times[2::]
        ###            b.TTL values begin low
        elif TTL_values[FIRST_TTL_IND] == 0:
            TTL_values = TTL_values[1::] # skip one to get next high TTL
            TTL_times = TTL_times[1::]
        ###        ii. truncate end
        ###            a. TTL values end high --> no action necessary
        ###            b. TTL values end low
        if TTL_values[-1] == 0:
            TTL_values = TTL_values[:-1]
            TTL_times = TTL_times[:-1]
        
        ### preallocate numpy arrays
        self.numTimeStamps_TTL = TTL_values.shape[0]
        self.TTL_squ_timesInSamples = np.full((self.numTimeStamps_TTL * NUM_PLATEU_PTS), np.nan)     
        self.TTL_squ_vals = np.full((self.numTimeStamps_TTL * NUM_PLATEU_PTS), np.nan)
        self.sine_timesInSamples = np.full((self.numTimeStamps_TTL-2)*NUM_SINE_VALS_PER_TTL//2, np.nan)
        self.sine_vals = np.full(self.sine_timesInSamples.shape, np.nan)
        self.frequencies = np.full(self.numTimeStamps_TTL//2, np.nan) 
        print('num time stamps: ' + str(self.numTimeStamps_TTL))
        print('TTL_values shape' + str(TTL_values.shape))
        self.TTL_timesInsecs = np.divide(TTL_times,SAMPLE_RATE_HZ)  # convert times to units of seconds

        ### II: make the sine wave
        NUM_TTLs_IN_WINDOW = 3 #number of points to average over when calculating the freq of the sine wave from TTL inputs (WIP  CODE DOES NOT SUPPORT OTHER VALUES AT PRESENT! )
        freqInd = 0
        sineInd = 0
        for ttlInd in range(0, self.numTimeStamps_TTL - NUM_TTLs_IN_WINDOW-1, NUM_TTLs_IN_WINDOW-1): # (need 3 TTLs to define 2 intervals; 2 intervals = 1 full sine period)
            ### i: select window
            TTL_times_in_window = TTL_times[ttlInd:ttlInd+NUM_TTLs_IN_WINDOW] # 
            ### ii: calculate the freq/period
            self.frequencies[freqInd] = 1 / (TTL_times_in_window[-1] - TTL_times_in_window[0])
            ### iii: calculate times for sine
            times = np.linspace(TTL_times_in_window[0], TTL_times_in_window[-1], NUM_SINE_VALS_PER_TTL, endpoint = False)
            self.sine_timesInSamples[sineInd:sineInd+NUM_SINE_VALS_PER_TTL] = times
            ### iv: calculate vals for sine
            try:
                self.sine_vals[sineInd:sineInd+NUM_SINE_VALS_PER_TTL] = MEAN_SINE_VAL + SINE_AMPLITUDE * np.sin(2 * np.pi * self.frequencies[freqInd] * (self.sine_timesInSamples[sineInd:sineInd+NUM_SINE_VALS_PER_TTL] - self.sine_timesInSamples[sineInd]))
            except:
                print('ttlInd:' + str(ttlInd))
                print('final incomplete sine period discarded --> TO DO: add code to handle this')            
            ### advance indices
            freqInd += 1
            sineInd += NUM_SINE_VALS_PER_TTL

        ### III: generate square wave from TTL highs and lows
        currentValue = TTL_values[0]
        for ttlInd in range(0,self.numTimeStamps_TTL): # TEST FIX
            currentTimeStamp = TTL_times[ttlInd]            
            if ttlInd < self.numTimeStamps_TTL-1: 
                nextTimeStamp = TTL_times[ttlInd+1]
            newTimeStamps = np.linspace(currentTimeStamp,nextTimeStamp,NUM_PLATEU_PTS, endpoint=False)
            firstNewInd = ttlInd*NUM_PLATEU_PTS
            lastNewInd = (ttlInd+1)*(NUM_PLATEU_PTS)
            self.TTL_squ_timesInSamples[firstNewInd:lastNewInd] = newTimeStamps 
            self.TTL_squ_vals[firstNewInd:lastNewInd] = currentValue
            if currentValue == 0:
                currentValue = 1
            elif currentValue ==1:
                currentValue = 0    
            
        ### remove trailing nans 
        self.sine_timesInSamples = self.sine_timesInSamples[~np.isnan(self.sine_timesInSamples)] 
        self.sine_vals = self.sine_vals[~np.isnan(self.sine_vals)]
        
        ### uncomment this to plot only the sine wave for debugging purposes
#        f = plt.figure()
#        plt.plot(self.TTL_squ_timesInSamples, self.TTL_squ_vals, 'o', color = 'b')
##            plt.plot(TTL_times, TTL_values, 'o', color = 'r')
#        plt.plot(self.sine_timesInSamples, self.sine_vals, 'o', color = 'g')
#        f.show()

    def plotRasters(self, *args): 
        print('plotting spike rasters:')
        if 'sine' in args: 
            
            for ttlTime in self.TTL_timesInsecs:
                
                ### get all spikes that are nearest to this ttl
                
                ## subtract ttl high time from spike times
                print('to do: plot by ttl phase')

            ### plot values
            f = plt.figure()
            plt.plot(self.TTL_squ_times, self.TTL_squ_vals, 'o', color = 'b')
#            plt.plot(TTL_times, TTL_values, 'o', color = 'r')
            plt.plot(self.sine_times, self.sine_vals, 'o', color = 'g')
    
            ### plt.plot( sineAlignedSpikes)
        
            ### vertical line for debuggin
            # FIRST_VLINE_IND = 243650
            # LAST_VLINE_IND = 243700
            # for i in range(FIRST_VLINE_IND,LAST_VLINE_IND):
            #     ind = int(i)
            #     print(ind)
                
            #     val = sine_times[ind]
            #     print(val)
            #     plt.axvline(x = val)

        plt.legend(['TTL square', 'TTL transitions', 'reconstructed sine stim'], loc='best')
        if INSPECT_EACH_PLOT ==1:
            print('close the plot window to advance to the next channel\nOR set INSPECT_EACH_PLOT to 0')
            plt.show()
            
            print('ploting rasters aligned to the sine period.')
            

            ### load TTL pulses

            ### recreate square wave

            ### recreate sine wave

        elif not args:
            print('no argument provided --> plotting generic rasters')
        

        ### plot spikes over recording length with TTL sine wave
        f = plt.figure()
        numUnits = len(self.units_by_spikesLst)
        ###### TO DO: auto generate appropriate line lengths and colorCodes for the given num of neurons
        lineLengths = np.ones(numUnits)
        colorCodes = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 0, 1]])
        plt.eventplot(self.units_by_spikesLst[0:8][:], color=colorCodes) #, linelengths = lineLengths)
        
        Y_OFFSET_FOR_TTLs = 2
        sine_vals_offset = self.sine_vals - Y_OFFSET_FOR_TTLs
        TTL_squ_vals_offset = self.TTL_squ_vals - Y_OFFSET_FOR_TTLs
        try:        
            plt.plot(self.TTL_squ_timesInSamples, TTL_squ_vals_offset, 'o', color = 'b')
            plt.plot(self.sine_timesInSamples, sine_vals_offset, 'o', color = 'g')
    ##            plt.plot(TTL_times, TTL_values, 'o', color = 'r')        
        except:
            print('WARNING: plotting TTL stimuli failed either because TTL were not present or due to an error')        
        
        plt.show()
        
        ### plot spikes color coded by unit ID and depth on a single sine period

        ### save plots to output dir

    def plotSpectra(self):
        print('TO DO: plotting spectra for units and stimulus')

    ### Welch PSD

    ### periodogram

    ### magnitude spectrum

    ### save plots to output dir
    
### ANALYSIS PARAMETERS:
NUM_CH = 64
SAMPLE_RATE_HZ = 30000
FIRST_CH = 1 # 1 indexed
LAST_CH = 3
INSPECT_EACH_PLOT = 0 # 0 for no; 1 for yes
OPEN_EACH_PDF = 0 # 1 to open each psd plot pdf before proceeding to the next plot
SAVE_TO_PDF = 0


## TTL plotting paramete
FIRST_TTL_IND = 2 # skip the first one because this is presumably truncated wrt to the normal period
LAST_TTL_IND = -1 # use -1 to go until the end
NUM_PLATEU_PTS = 20 # number of points to add per high/low TTL pulse to recreate the square wave from transistion points
MEAN_SINE_VAL = 0.5
SINE_AMPLITUDE = 0.6
NUM_SINE_VALS_PER_TTL = 512

### main code execution (call desired )

SPIKE_DIR = os.getcwd() ### directory of spike sorted data
data = Spikalyze(SPIKE_DIR) ### load sorted data


#rawDataDir = data.getRawDataDir() ### find the raw data dir from which the sorted data originated
TTL_times, TTL_values = data.loadTTLdata()
data.reconstructSineFromTTLs(TTL_values, TTL_times)


data.plotRasters()
# data.plotRasters('sine') ### TO DO

data.plotSpectra()

print('finished!')

