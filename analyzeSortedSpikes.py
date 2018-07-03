import numpy as np
import scipy.io as sio
import pandas as pd
import os.path

### code from Kyu Lee
# path to sorting folder
# path = 'J:/ksort/20180618'
path = os.getcwd()

# load stuff
spike_times = np.load(os.path.join(path,'spike_times.npy'))
spike_clusters = np.load(os.path.join(path,'spike_clusters.npy'))
spike_templates = np.load(os.path.join(path,'spike_templates.npy'))
templates = np.load(os.path.join(path,'templates.npy'))
cluster_groups = pd.read_csv(os.path.join(path,'cluster_groups.csv'),sep='\t')
clusters = np.array(cluster_groups['cluster_id'])
group = np.array(cluster_groups['group'])
channel_map = np.load(os.path.join(path,'channel_map.npy'))

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
sp = np.array(list(output.values()))

# save as mat
sio.savemat(os.path.join(path,'sp.mat'),{'cellid':cellid,'masks':masks,'sp':sp})
### end of Kyu's code