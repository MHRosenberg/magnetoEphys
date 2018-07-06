%  create a channel map file

Nchannels = 64;
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;

%%%% MHR: probe sensor geometry goes here! ONLY RUN THIS FIRST BLOCK! 5_7_18
xcoords   = [0
0
0
0
0
0
0
0
0
4
20
20
20
20
20
20
20
20
20
20
20
20
36
40
40
40
40
40
40
40
40
40
300.1
300.1
300.1
300.1
300.1
300.1
300.1
300.1
300.1
304.1
320.1
320.1
320.1
320.1
320.1
320.1
320.1
320.1
320.1
320.1
320.1
320.1
336.1
340.1
340.1
340.1
340.1
340.1
340.1
340.1
340.1
340.1
]';
ycoords   = [50
0
99.9
150
200
250
300
350
400
450
475
375
275
175
75
-25
25
125
225
325
425
500
450
400
350
300
250
200
150
0
50
100
50
0
99.9
150
200
250
300
350
400
450
475
375
275
175
75
-25
25
125
225
325
425
500
450
400
350
300
250
200
150
0
50
100
]';
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

fs = 30000; % sampling frequency
save('/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Analysis/spikeSorting/kiloSorting/probe64G_chanMap.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')

%%
%%%%% MHR 5_7_18: THIS DISTORTS THE X AND Y COORDINATES; I DON'T KNOW WHAT THE PURPOSE OF THIS BLOCK IS... 


Nchannels = 64; %%% MHR changed from 32 2/4/18
connected = true(Nchannels, 1);
chanMap   = 1:Nchannels;
chanMap0ind = chanMap - 1;

xcoords   = repmat([1 2 3 4]', 1, Nchannels/4);
xcoords   = xcoords(:);
ycoords   = repmat(1:Nchannels/4, 4, 1);
ycoords   = ycoords(:);
kcoords   = ones(Nchannels,1); % grouping of channels (i.e. tetrode groups)

fs = 30000; % sampling frequency %%% MHR changed 25000 from 32 2/4/18

save('/media/matthew/Data/a_Ephys/a_Magnetoreception/a_Data/spikesorting/a_paramsNconfigs/kiloSorting/probe64F_chanMap_unclearPurpose.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
%%

% kcoords is used to forcefully restrict templates to channels in the same
% channel group. An option can be set in the master_file to allow a fraction 
% of all templates to span more channel groups, so that they can capture shared 
% noise across all channels. This option is

% ops.criterionNoiseChannels = 0.2; 

% if this number is less than 1, it will be treated as a fraction of the total number of clusters

% if this number is larger than 1, it will be treated as the "effective
% number" of channel groups at which to set the threshold. So if a template
% occupies more than this many channel groups, it will not be restricted to
% a single channel group. 