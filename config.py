#!/usr/bin/env python3

## dirs and files
BaseDir = '/Users/wbhawley/Research/Seismology/Blanco_SurfaceWave/'
ConfigDir = BaseDir+'config/'
DataDir = BaseDir+'data/'
EventsFileName = ConfigDir+'EventsTest.txt'
# two station lists, one for BXH, one for BDH
XStafn = ConfigDir+'X9_stations_X.txt'
DStafn = ConfigDir+'X9_stations_D.txt'
Stafn = ConfigDir+'X9_stations.txt'

## Event download
minMag = 6.5
webservice = "IRIS"
network = "X9" # X9 = Blanco
isCMT_params = 1 # use GCMT parameters for SAC header; 0 = use IRIS
isCentroid = 1 # if isCMT_params = 1, use centroid; 0 = epicentral

tstart = '2012-11-16T00:00:00'
tend = '2012-11-16T23:00:00'

## Stations to download
# need one for low-pass filtered stations
inFile = open(XStafn)
XStaList = []
for line in inFile:
    XStaList.append(line.rstrip('\n'))
XChanList = ['LHZ','LH1','LH2','BXH']

# and for regular stations...
inFile = open(DStafn)
DStaList = []
for line in inFile:
    DStaList.append(line.rstrip('\n'))
DChanList = ['LHZ','LH1','LH2','BDH']

# and for all stations...
inFile = open(Stafn)
StaList = []
for line in inFile:
    StaList.append(line.rstrip('\n'))
ChanList = ['LHZ','LH1','LH2','BDH']
    
# trace info
isDownsamp = 1 # downsample option
srNew = 1 # new sample rate, samples/sec
trLen = 6000 # len of traces (sec)
