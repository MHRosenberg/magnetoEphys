import numpy as np
import scipy.io as sio
import pandas as pd
import os.path
import matplotlib.pyplot as plt
from scipy.signal import welch, hanning

### adapted from Kyu Lee's code
KILO_RESULTS_DIR = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Data/spikesorting/kiloSorted/MR1_2018-05-02_14-56-40_baseline2018_07_03__01_24'

# load stuff
spike_times = np.load(os.path.join(KILO_RESULTS_DIR,'spike_times.npy'))
spike_clusters = np.load(os.path.join(KILO_RESULTS_DIR,'spike_clusters.npy'))
spike_templates = np.load(os.path.join(KILO_RESULTS_DIR,'spike_templates.npy'))
templates = np.load(os.path.join(KILO_RESULTS_DIR,'templates.npy'))
cluster_groups = pd.read_csv(os.path.join(KILO_RESULTS_DIR,'cluster_groups.csv'),sep='\t')
clusters = np.array(cluster_groups['cluster_id'])
group = np.array(cluster_groups['group'])
channel_map = np.load(os.path.join(KILO_RESULTS_DIR,'channel_map.npy'))

# find good clusters
good_clusters = clusters[group=='good']

every_channel = {}
max_channel = {} # will hold channel w/ max amplitude
output = {} # will hold spikes
for cluster in good_clusters:
    output[cluster] = spike_times[spike_clusters==cluster]
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
    
# output
cellid = np.array(list(max_channel.keys()))
masks=channel_map[list(max_channel.values())].flatten()
spikes = np.array(list(output.values())) ### index with: [unitNum][spikeInd]

# save as mat
if not os.path.exists(KILO_RESULTS_DIR+'/spikes.mat'):
    sio.savemat(os.path.join(KILO_RESULTS_DIR,'spikes.mat'),{'cellid':cellid,'masks':masks,'sp':spikes})
###


### user parameters 
SAMPLE_RATE = 30000 # in Hz
FIRST_SAMPLE = 0
LAST_SAMPLE = 8000000007 # TO DO: calculate a reasonable value here from the last spike, last TTL pulse, etc
BIN_SIZE_IN_MS = 10


binSizeInSamples = BIN_SIZE_IN_MS * SAMPLE_RATE / 1000   
numBins = LAST_SAMPLE // binSizeInSamples # note: this discards the remainder
lastSampleIncluded = numBins * binSizeInSamples 


bins = np.linspace(FIRST_SAMPLE,lastSampleIncluded, numBins)
hist, binEdges = np.histogram(spikes[0],np.linspace(0,lastSampleIncluded,numBins,endpoint=False)) ################ TO DO: replace w/ unit ind


nblock = 1024
overlap = 128
win = hanning(nblock, True)
f, Pxx_den = welch(hist, SAMPLE_RATE, window=win, noverlap=overlap, nfft=nblock, return_onesided=True)

print(hist)
print(binEdges)
print('plotting (kinda slow at present)')


plt.xlabel('frequency [Hz]')
plt.ylabel('PSD [V**2/Hz?]')

# plt.plot(bins[:-1], hist)#, width = binSizeInSamples)
# plt.xlim(min(binEdges), max(binEdges))

plt.semilogy(f,Pxx_den, '-o')
plt.show()


for unit in spikes:
    print(unit.shape)
