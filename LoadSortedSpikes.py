import numpy as np
import scipy.io as sio
import pandas as pd
import os.path
import pdb
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import style
from matplotlib import pyplot as plt 
from matplotlib.backends.backend_pdf import PdfPages

### code from Kyu Lee
# path to sorting folder
# path = 'J:/ksort/20180618'


# path = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Data/spikesorting/kiloSorted/2018-05-04_13-55-49_MR1_r2_pre2018_07_10__12_55'

class LoadSortedSpikes:

    def __init__(self, path):

        # load stuff
        spike_times = np.load(os.path.join(path,'spike_times.npy'))
        spike_clusters = np.load(os.path.join(path,'spike_clusters.npy'))
        spike_templates = np.load(os.path.join(path,'spike_templates.npy'))
        templates = np.load(os.path.join(path,'templates.npy'))
        cluster_groups = pd.read_csv(os.path.join(path,'cluster_groups.csv'),sep='\t')
        clusters = np.array(cluster_groups['cluster_id'])
        group = np.array(cluster_groups['group'])
        channel_map = np.load(os.path.join(path,'channel_map.npy'))


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

        ### save as mat
        sio.savemat(os.path.join(path,'sp.mat'),{'cellid':cellid,'masks':masks,'sp':sp})
        ### end of Kyu's code

        ### save as numpy
        np.save('spikeData', {'cellid':cellid,'masks':masks,'sp':sp})
        
    
    def spikeDictToList(self):
        for unit in self.output.keys():
            print(unit)
            unitSpikeTimeLst = []
            for spikeTime in self.output[unit][:][:]:
                unitSpikeTimeLst.append(int(spikeTime))
            self.units_by_spikesLst.append(unitSpikeTimeLst)
        return self.units_by_spikesLst

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
path = os.getcwd()
data = LoadSortedSpikes(path)
spikes = data.spikeDictToList()
data.plotRasters()
# data.plotRasters('sine') ### TO DO

