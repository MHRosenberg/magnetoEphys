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


#        pdb.set_trace()

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
        for unit in self.output.keys():
            print(unit)
            unitSpikeTimeLst = []
            for spikeTime in self.output[unit][:][:]:
                unitSpikeTimeLst.append(int(spikeTime))
            self.units_by_spikesLst.append(unitSpikeTimeLst)
        return self.units_by_spikesLst

    def loadTTLdata(self, *RAW_DATA_DIR): # from the .events file
        RAW_DATA_SUBDIR = 'a_Data/raw'
        if not RAW_DATA_DIR:
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
                                sys.exit('ERROR: unable to locate raw data dir in reasonable parent directories')
            
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
        
        else:
            
            rawDataDir = RAW_DATA_DIR
        
        print()
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

    def plotRasters(self, *args): 
        print('plotting spike rasters:')
        if 'sine' in args: 
            print('TO DO: plot rasters aligned to the sine period.')

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
        plt.show()
        
        ### plot spikes color coded by unit ID and depth on a single sine period

        ### save plots to output dir

    def plotSpectra():
        print('TO DO: plotting spectra for units and stimulus')

    ### Welch PSD

    ### periodogram

    ### magnitude spectrum

    ### save plots to output dir
SPIKE_DIR = os.getcwd()
data = Spikalyze(SPIKE_DIR)


TTL_timeStamps, TTL_values = data.loadTTLdata()


data.plotRasters()
# data.plotRasters('sine') ### TO DO

print('finished!')

