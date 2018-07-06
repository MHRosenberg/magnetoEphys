import os.path
import os
import time
import datetime

FILES_TO_MOVE_PATH = os.getcwd() # WARNING: probably doesn't support changes to this at present!
KILO_SAVE_PARENT_DIR = '/media/matthew/Data/a_Ephys/a_Projects/a_Magnetoreception/a_Data/spikesorting/kiloSorted'

ts = time.time()
timeStr = datetime.datetime.fromtimestamp(ts).strftime('%Y_%m_%d__%H_%M')
newSortingDirName = os.path.basename(os.getcwd()) + timeStr

fullDestinationPath = KILO_SAVE_PARENT_DIR + '/' + newSortingDirName

if not os.path.exists(fullDestinationPath):
	os.makedirs(fullDestinationPath)
os.system('mv ' + FILES_TO_MOVE_PATH + '/kilosort.dat ' + fullDestinationPath)
os.system('mv ' + FILES_TO_MOVE_PATH + '/*.npy ' + fullDestinationPath)
os.system('mv ' + FILES_TO_MOVE_PATH + '/params.py ' + fullDestinationPath)
os.system('mv ' + FILES_TO_MOVE_PATH + '/rez.mat ' + fullDestinationPath)