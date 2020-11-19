#!/usr/bin/env python3

#  dirs and files
BaseDir = '/Users/wbhawley/Research/Seismology/Blanco-SW/'
ConfigDir = BaseDir+'config/'
DataDir = BaseDir+'data/'
EventsDataDir = DataDir+'SAC_Events/'
NoiseDataDir = DataDir+'SAC_Noise/'
#EventsFileName = ConfigDir+'BlancoEventTest_M6.5.txt'
EventsFileName = ConfigDir+'BlancoProblematicEvt.txt'
DayFileName = ConfigDir+'BlancoProblematicDay.txt'
ANTDayFileName = ConfigDir+'BlancoANTDaysTest_M6.5.txt'
# two station lists, one for BXH, one for all stations
XStafn = ConfigDir+'X9_stations_X.txt'
Stafn = ConfigDir+'X9_stations.txt'
#Stafn = ConfigDir+'7D_stations.txt'
#Stafn = ConfigDir+'SC_stations.txt'

#  Event download
minMag = 6.5
webservice = "IRIS"
#webservice = "SCEDC"
network = "X9"  # X9 = Blanco
#network = "7D"  # 7D = Cascadia Initiative
#network = "CI"
isCMT_params = 1  # use GCMT parameters for SAC header; 0 = use IRIS
isCentroid = 1  # if isCMT_params = 1, use centroid; 0 = epicentral

#  Noise Download
trLen = 60 * 60 * 24  # seconds 
noDays = 4  # number of days prior to event to use
isCalDay = 1  # 0 to start each day at 00:00; 0 to use 24h segments prior to eq

tstart = '2012-10-09T00:00:00'
tend = '2012-10-12T23:00:00'

tstart = '2013-07-10T00:00:00'
tend = '2013-08-17T23:59:59'

tstart = '2013-02-01T00:00:00'
tend = '2013-02-28T23:59:59'

tstart = '2013-07-21T00:00:00'
tend = '2013-07-21T10:00:00'

# Stations to download
# need one for low-pass filtered stations
inFile = open(XStafn)
XStaList = []
for line in inFile:
    XStaList.append(line.rstrip('\n'))
XChanList = ['LHZ', 'LH1', 'LH2', 'BXH']

# and for all stations...
inFile = open(Stafn)
StaList = []
for line in inFile:
    StaList.append(line.rstrip('\n'))
ChanList = ['LHZ', 'LH1', 'LH2', 'BDH']
#ANTChanList = ['LHZ']

# Event trace info
isDownsamp = 1  # downsample option
isRemoveResp = 1  # remove response option
srNew = 1  # new sample rate, samples/sec
trLen = 6000  # len of traces (sec)

# Frequency Limits for Response Removal
# Hi corner deternimed by Nyquist
# Lo corner defined by
LoFreq1 = 0.001
LoFreq2 = 0.005