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

### the majority of the __init__ code is rom Kyu Lee's analyzeSortedSpikes.py
# path to sorting folder
# path = 'J:/ksort/20180618'

### NOTE: assume all times to be in secs/Hz unless specified otherwise!

class Spikalyze:
    def __init__(self, spikeDir):
        # load stuff
        spike_times = np.load(os.path.join(spikeDir,'spike_times.npy'))  ### in samples
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
        
        ### save as mat
        sio.savemat(os.path.join(spikeDir,'sp.mat'),{'cellid':cellid,'masks':masks,'sp':sp})
        ### end of Kyu's code
        
        self.unitsXspikesLst = []
        self.unitsXspikesLst = self.spikeDictToList() ### spike times converted from sample numbers to secs here!
        
        ### declare stimulus numpy arrays
        self.numTimeStamps_TTL = []
        self.TTL_squ_timesInSecs = [] 
        self.TTL_squ_vals = []
        self.sine_timesInSecs = [] 
        self.sine_vals = []
        self.frequencies = []
        self.sampleRate = []

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
        print('loading the following units classified as good by the spikesorting user from Kyu\'s dict to my list of lists:')
        SAMPLE_RATE = 30000
        for unit in self.output.keys():    
            print(unit)
            unitSpikeTimeLst = []
            for spikeTime in self.output[unit][:][:]:
                spikeTimeInSecs = float(spikeTime / SAMPLE_RATE) ### TO DO: make this flexible based on the found sample rate
#                unitSpikeTimeLst.append(int(spikeTime)) # int req'd to strip the number out of the 1 element array
                unitSpikeTimeLst.append(spikeTimeInSecs) # int req'd to strip the number out of the 1 element array
            self.unitsXspikesLst.append(unitSpikeTimeLst)
        return self.unitsXspikesLst
    
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
        eventsDict = OE.load(event_path) ### load TTL dict file
        eventIDs = eventsDict['eventType']   

        print('\nloading: ' + event_path)
        print(eventsDict['eventType'])
        print(eventsDict['eventId'])
        print('\nall the events that are saved have type TTL = 3 ; Network Event = 5')
        print('event IDs: ' + str(eventIDs))
        
        self.sampleRate = float(eventsDict['header']['sampleRate'])
        TTL_timesInSamples = eventsDict['timestamps']
        self.TTL_values = eventsDict['eventId']
        self.TTL_timesInSecs = TTL_timesInSamples / self.sampleRate
        
        print('shape of TTL values' + str(self.TTL_values.shape))
        print('TTL values: ' + str(self.TTL_values))
        print('shape of timestamps: ' + str(TTL_timesInSamples.shape))
        print('time stamps: ' + str(TTL_timesInSamples))
        print('sample rate (from all_channels.events): ' + str(self.sampleRate))
        print('converting times from sample number to secs')
        
        ### skip weird bad values at the beginning... TO DO: find out why these are so rapid yet still have event id 3. check the other data set
#        TTL_timeStamps = TTL_times[4::]
#        self.TTL_values = self.TTL_values[4::]
        
        return self.TTL_timesInSecs, self.TTL_values
    
    def reconstructSineFromTTLs(self, TTL_values, TTL_times, *args):
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
        elif TTL_values[0] == 0:
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
        self.TTL_squ_timesInSecs = np.full((self.numTimeStamps_TTL * NUM_PLATEU_PTS), np.nan)     
        self.TTL_squ_vals = np.full((self.numTimeStamps_TTL * NUM_PLATEU_PTS), np.nan)
        
#        self.sine_timesInSecs = np.full((self.numTimeStamps_TTL-2)*NUM_SINE_VALS_PER_TTL//2, np.nan) # original 
        self.sine_timesInSecs = np.full((self.numTimeStamps_TTL-2)*NUM_SINE_VALS_PER_TTL, np.nan) # DEBUGGING 8/8/18
         
        self.sine_vals = np.full(self.sine_timesInSecs.shape, np.nan)
        self.frequencies = np.full(self.numTimeStamps_TTL//2, np.nan) 
        print('num time stamps: ' + str(self.numTimeStamps_TTL))
        print('TTL_values shape' + str(TTL_values.shape))
#        self.TTL_timesInsecs = np.divide(TTL_times,SAMPLE_RATE_HZ)  # convert times to units of seconds

        ### II: make the sine wave
        NUM_TTLs_IN_WINDOW = 3 #number of points to average over when calculating the freq of the sine wave from TTL inputs (WIP  CODE DOES NOT SUPPORT OTHER VALUES AT PRESENT! )
        freqInd = 0
        sineInd = 0
        for ttlInd in range(0, self.numTimeStamps_TTL - NUM_TTLs_IN_WINDOW-1, NUM_TTLs_IN_WINDOW-1): # (need 3 TTLs to define 2 intervals; 2 intervals = 1 full sine period)
            ### i: select window
            TTL_times_in_window = TTL_times[ttlInd:ttlInd+NUM_TTLs_IN_WINDOW]  
            ### ii: calculate the freq/period
            self.frequencies[freqInd] = 1 / (TTL_times_in_window[-1] - TTL_times_in_window[0])
            ### iii: calculate times for sine
            
            times = np.linspace(TTL_times_in_window[0], TTL_times_in_window[-1], NUM_SINE_VALS_PER_TTL, endpoint = False) 
            
            self.sine_timesInSecs[sineInd:sineInd+NUM_SINE_VALS_PER_TTL] = times
            self.sine_timesInSecs[sineInd:sineInd+NUM_SINE_VALS_PER_TTL] = times
            ### iv: calculate vals for sine
            try:
                self.sine_vals[sineInd:sineInd+NUM_SINE_VALS_PER_TTL] = MEAN_SINE_VAL + SINE_AMPLITUDE * np.sin(2 * np.pi * self.frequencies[freqInd] * (self.sine_timesInSecs[sineInd:sineInd+NUM_SINE_VALS_PER_TTL] - self.sine_timesInSecs[sineInd]))
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
            self.TTL_squ_timesInSecs[firstNewInd:lastNewInd] = newTimeStamps 
            self.TTL_squ_vals[firstNewInd:lastNewInd] = currentValue
            if currentValue == 0:
                currentValue = 1
            elif currentValue ==1:
                currentValue = 0    
            
        ### remove trailing nans 
        self.sine_timesInSecs = self.sine_timesInSecs[~np.isnan(self.sine_timesInSecs)] ### TO DO: explicitly remove }{ calculate appropriately from the beginning
        self.sine_vals = self.sine_vals[~np.isnan(self.sine_vals)]
        
        if 'debug' in args:
            print('plotting stimulus alone for debugging purposes (remove debug from function call to skip this)') 
            f = plt.figure()
            plt.ylabel('TTL (eg sine) freq in Hz:')
            plt.xlabel('TTL num')
#            freqsInSecs = np.multiply(self.frequencies, self.sampleRate)
#            plt.plot(freqsInSecs)
            plt.plot(self.frequencies)
            
            g = plt.figure()
            plt.ylabel('AU')
            plt.xlabel('seconds')
            plt.plot(self.TTL_squ_timesInSecs, self.TTL_squ_vals, 'o', color = 'b')
            plt.plot(TTL_times, TTL_values, 'o', color = 'r')
#            sine_timesSecs = self.sine_timesInSecs * self.sampleRate ### outdated cuz conversion to secs happens above now
            plt.plot(self.sine_timesInSecs, self.sine_vals, 'o', color = 'g')
            
            h = plt.figure()
            plt.ylabel('interval between TTL pulses in secs')
            plt.xlabel('TTL num')
            plt.plot(np.diff(TTL_times))
            
            plt.show()
            ### these create a weird terminating warning in IPython interactive 
#            f.show()
#            g.show()
#            h.show()

    def plotRasters(self, *args): 
        print('plotting spike rasters:')
        if 'sine' in args: 
            print('plotting ttl aligned spikes')
            
            ### TO DO: make this robust to low TTL first without dropping data
             
            self.units_x_SpikeSineTimes = []
            totalNumMatches = 0
            unitNumMatches = 0            
            ### IIa: for loop over all the units
            for unit in range(0,len(self.unitsXspikesLst)):
                print('aligning spikes to high TTLs for unit: ' + str(unit))
                ttlSubtractedUnitSpikeTimes = []
                ### IIb: for loop over high TTLs Oo. skip every other starting with high TTL
                for highTTL in range(0, (len(TTL_times)//2)-2, 2):
                    ### III: get earliest and latest time for this sine/TTL period
                    earliestTime = TTL_times[highTTL]
                    latestTime = TTL_times[highTTL+2] ## skip the low TTL; ################# WARNING: does this need to have something added to avoid double counting the boundaries?
                    ### IV: find all values in spikes in III's time window and append unit's sine aligned spike list
                    for spikeTime in self.unitsXspikesLst[unit]:
                        if spikeTime > earliestTime and spikeTime < latestTime:
                            ### V: subtract high TTL times from IV
                            ttlSubtractedUnitSpikeTimes.append(spikeTime - earliestTime) 
                            print('match added to sine aligned 2d list')
                            print('unit: ' + str(unit) + '; time win: ' + str(earliestTime) + ' - ' + str(latestTime) + '; spikeTime: ' + str(spikeTime) + '; total matches: ' + str(totalNumMatches) + '; unit matches: ' + str(unitNumMatches))
                            totalNumMatches += 1
                            unitNumMatches += 1
                    ### VII: append unit spike list to unit list
                self.units_x_SpikeSineTimes.append(ttlSubtractedUnitSpikeTimes)
                ttlSubtractedUnitSpikeTimes = []
                unitNumMatches = 0
            
            print('plotting spikes aligned to TTL (sine) phase:')
            plt.figure()
            LINE_LENGTHS = 0.8
            plt.eventplot(self.units_x_SpikeSineTimes, linelengths = LINE_LENGTHS) #, linelengths = lineLengths)
            plt.show()
            print('done')

        elif not args:
            print('no argument provided --> plotting generic rasters')

        ### plot spikes over recording length (with TTL stimuli if present)
        plt.figure()
        numUnits = len(self.unitsXspikesLst)
        
        ####################### TO DO: auto generate appropriate line lengths and colorCodes for the given num of neurons
        lineLengths = np.ones(numUnits)
        colorCodes = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 0, 1]])
#        plt.eventplot(self.unitsXspikesLst[0:8][:], color=colorCodes) #, linelengths = lineLengths)
        plt.eventplot(self.unitsXspikesLst[:][:]) #, linelengths = lineLengths)
        
        
        try:       
            Y_OFFSET_FOR_TTLs = 2
            sine_vals_offset = self.sine_vals - Y_OFFSET_FOR_TTLs
            TTL_squ_vals_offset = self.TTL_squ_vals - Y_OFFSET_FOR_TTLs
            TTL_vals_offset = TTL_values - Y_OFFSET_FOR_TTLs
            plt.plot(self.TTL_squ_timesInSecs, TTL_squ_vals_offset, 'o', color = 'b')
            plt.plot(self.sine_timesInSecs, sine_vals_offset, 'o', color = 'g')
            plt.plot(TTL_times, TTL_vals_offset, 'o', color = 'r')
            print('stimuli plotted as negative y values')        
        except:
            print('WARNING: plotting TTL stimuli failed either because TTL were not present or due to an error')        
        plt.show()
        

        ### TO DO: save plots to output dir

    def plotSpectra(self, FIRST_SAMPLE, LAST_SAMPLE):
        print('\nplotting spectra for units and stimulus')
        
        totalStimDurationSecs = LAST_SAMPLE - FIRST_SAMPLE ### in secs
        if totalStimDurationSecs > 10: ## stimulus duration must be longer than 10 secs to plot it
            plotSine = True
        else:
            plotSine = False
            FIRST_SAMPLE = 0
            LAST_SAMPLE = np.max(np.max(self.unitsXspikesLst))
            totalStimDurationSecs = LAST_SAMPLE - FIRST_SAMPLE ### in secs
        
        ### neuron's bin'd sample rate
        numBins, remainderSecs = np.divmod(totalStimDurationSecs, BIN_SIZE_SEC) ### totalStimDuration is from 0 to the last spike time if no TTLs are present
        bins = np.linspace(FIRST_SAMPLE, FIRST_SAMPLE + round(float(numBins) * float(BIN_SIZE_SEC)), int(numBins), endpoint=True) 
        binSampleRate = 1./np.mean(np.diff(bins))
        
        
        print('binSampleRate = {0}'.format(binSampleRate))
        
        ### Welch PSD
        ### sine stimulus
        if plotSine == True:
            ### sine stim's bin'd sample rate
            try:
                sineStimSampleRate = self.sine_vals.shape[0] / totalStimDurationSecs
            except:
                sys.exit('\n\nERROR: You probably need to uncomment data.reconstructSineFromTTLs(TTL_values, TTL_times) below to get the TTLs...\n\n')
            print('sineStimSampleRate = {0}'.format(sineStimSampleRate))
            welchFreqs_sine, Pxx_den = welch(self.sine_vals, nperseg = sineStimSampleRate * 100, fs = sineStimSampleRate) ######## TO DO: Fix the sample rate }{ time units
            
            ### set up figure assuming the sine is present
            NUM_ROWS = 3
            NUM_COLUMNS = 1 
            ### gridspec.GridSpec(NUM_ROWS,NUM_COLUMNS) # state-based versions for subplot
            ### plt.subplot(611) ### doesn't require gridspec but is kinda clunky for dynamic changes        
            FIG_SIZE = (NUM_ROWS, NUM_COLUMNS)
            fig = plt.figure(figsize = FIG_SIZE)
            rowPosition = 0
            columnPosition = 0
            
            ### sine stimulus
            ax1 = plt.subplot2grid((NUM_ROWS, NUM_COLUMNS), (rowPosition,columnPosition))
            rowPosition += 1
            ax1.set_title('sine stimulus Welch PSD')
            ax1.set_xlabel('frequency [Hz]')
            ax1.set_ylabel('PSD [V**2/Hz?]')
            ax1.semilogy(welchFreqs_sine,Pxx_den, '-o')
            ax1.set_xlim(0, MAX_FREQ_PLOTTED)
        
        elif plotSine == False:
            print('sine stimulus sample rate calculation failed (this is expected if there are no TTLs for this recording)')
            ### set up the figure if there's no sine stimulus
            NUM_ROWS = 2
            NUM_COLUMNS = 1 
            FIG_SIZE = (NUM_ROWS, NUM_COLUMNS)
            fig = plt.figure(figsize = FIG_SIZE)
            rowPosition = 0
            columnPosition = 0
        
        
        
        ### neurons
        ax2 = plt.subplot2grid((NUM_ROWS, NUM_COLUMNS), (rowPosition,columnPosition))
        rowPosition += 1
        ax2.set_title('neurons\' Welch PSD (sci py ver)')
        ax2.set_xlabel('frequency [Hz]')
        ax2.set_ylabel('PSD [V**2/Hz?]')
        ax2.set_xlim(0, MAX_FREQ_PLOTTED)
        
        ax3 = plt.subplot2grid((NUM_ROWS, NUM_COLUMNS), (rowPosition,columnPosition))
        rowPosition += 1
        ax3.set_title('neurons\' PSD (matplotlib ver)')
        ax3.set_xlim(0, MAX_FREQ_PLOTTED)
        
#        ax4 = plt.subplot2grid((NUM_ROWS, NUM_COLUMNS), (rowPosition,columnPosition))
#        rowPosition += 1
#        ax4.set_title('neurons\' magnitude spectrum')
#        ax4.set_xlim(0, MAX_FREQ_PLOTTED)
        
#        ax5 = plt.subplot2grid((NUM_ROWS, NUM_COLUMNS), (rowPosition,columnPosition))
#        rowPosition += 1
#        ax5.set_title('neurons\' spectrogram')
#        ax5.set_ylabel('Freq in Hz')
#        ax5.set_xlabel('time in Secs')
#        ax5.set_ylim(0,100)
        
        if UNITS_PLOTTED == [0, 0]:
            firstUnit = 0
            lastUnit = len(self.unitsXspikesLst)-1 ### return to original after debugging
        else:
            firstUnit = UNITS_PLOTTED[0]
            lastUnit = UNITS_PLOTTED[-1]
        for unit in range(firstUnit,lastUnit+1):
            t = "unit: {0}".format(unit)
            print(t)
            
#            unit_hist, binEdges = np.histogram(self.unitsXspikesLst[unit], bins, (FIRST_SAMPLE, LAST_SAMPLE))
            unit_hist, binEdges = np.histogram(self.unitsXspikesLst[unit], bins)
            
            ### Welch PSD            
            welchFreqs_neuron, Pxx_den = welch(unit_hist,
               fs= binSampleRate, # sample rate
               window='hanning',   # apply a Hanning window before taking the DFT
               nperseg= self.sampleRate,        # compute periodograms of 1024-long segments of x
               detrend='constant')
            ax2.semilogy(welchFreqs_neuron,Pxx_den, '-o', label = str(unit))
#            ax2.plot(welchFreqs_neuron,Pxx_den, '-o')
#            ax2.legend(['unit: {0}'.format(unit)])
            
            ### Welch PSD
#            ax3.psd(unit_hist, binSampleRate, pad_to=1024, detrend = 'mean') # try pad_to to increase the number of points... 
            ax3.psd(unit_hist, Fs= binSampleRate, NFFT= self.sampleRate, marker='o') 
            
            ### magnitude spectrum
#            ax4.magnitude_spectrum(unit_hist, Fs = binSampleRate, scale = 'dB', marker = 'o')
#            try:
#                ax4.magnitude_spectrum(unit_hist, Fs = binSampleRate, marker = 'o')
#            except:
#                print('magnitude spectrum threw an error for unit {0};\n error: noverlap must be less than n'.format(unit))
            
            ### Spectrogram
#                Pxx, freqs, bins, im = ax5.specgram(self.unitsXspikesLst[unit], Fs=self.sampleRate, noverlap=100)                
#                Pxx, freqs, bins, im = ax5.specgram(unit_hist, Fs=self.sampleRate, noverlap=100)
#            Pxx, freqs, bins, im = ax5.specgram(unit_hist, Fs= binSampleRate) ### TO DO: FIX RUNTIME WARNING: DIVIDE BY ZERO

#        fig.colorbar(im).set_label('Intensity [dB]')
#        fig.tight_layout()
        fig.subplots_adjust(hspace = 1.0)
        fig.set_size_inches(w=15,h=8)
        fig_name = 'plot.png'
        fig.savefig(fig_name)
        
        plt.show()
               
    
        ### save plots to output dir

### TTL plotting parameters
NUM_PLATEU_PTS = 20 # number of points to add per high/low TTL pulse to recreate the square wave from transistion points
MEAN_SINE_VAL = 0.5
SINE_AMPLITUDE = 0.6
NUM_SINE_VALS_PER_TTL = 1000

####### main code execution (call desired )
SPIKE_DIR = os.getcwd() ### directory of spike sorted data
data = Spikalyze(SPIKE_DIR) ### load sorted data


#rawDataDir = data.getRawDataDir() ### find the raw data dir from which the sorted data originated
TTL_times, TTL_values = data.loadTTLdata()


#data.reconstructSineFromTTLs(TTL_values, TTL_times, 'debug') # optinal debug flab at end plots the stimulus alone
data.reconstructSineFromTTLs(TTL_values, TTL_times) # req'd for plot spectra TO DO: remove this dependence

### plot rasters
#data.plotRasters()
data.plotRasters('sine') # flag aligns all spikes to the TTL pulses

### plot PSDs
### neural PSD parameters
FIRST_SAMPLE = TTL_times[0] 
LAST_SAMPLE = TTL_times[-2] 
BIN_SIZE_SEC = 5/1000 ### ms bin size
MAX_FREQ_PLOTTED = 15
UNITS_PLOTTED = [0, 0] # use [0, 0] to plot all of them
data.plotSpectra(FIRST_SAMPLE, LAST_SAMPLE)



print('finished!')

