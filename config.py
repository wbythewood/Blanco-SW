#!/usr/bin/env python3

#  dirs and files
#BaseDir = '/Users/whawley/Research/github/Blanco-SW/'

#minMag = 6.5
#Label = 'GordaAddition'
minMag = 5.0
Label = 'CoverAzimuth'

# for the extra events to download based on azimuth
NAziBin = 8 # number of bins in azimuth we want to make sure there are enough events from
LargeEqCutoff = 6.3 # count no. of earthquakes in each bin down to this mag; use that no. to find smaller eqs in other bins

# old mac
#BaseDir = '/Users/wbhawley/Research/Seismology/Blanco-SW/'
# gaherty mac
BaseDir = '/Users/whawley/Research/github/Blanco-SW/'
# new mac
#BaseDir = '/Users/whawley/Research/Blanco-SW/'

ConfigDir = BaseDir+'config/'
DataDir = BaseDir+'data/'
EventsDataDir = DataDir+'SAC_Events/'
EventsDataDir = DataDir+'SAC_Events_'+Label+'/'
NoiseDataDir = DataDir+'SAC_Noise_'+Label+'/'
#EventsFileName = ConfigDir+'BlancoEventTest_M6.5.txt'
#EventsFileName = ConfigDir+'BlancoEvents_'+Label+'_M'+str(minMag)+'.txt'
#EventsFileName = ConfigDir+'BlancoANT_Test_12h.txt'
#for azi part:
EventsFileName = ConfigDir+'BlancoEvents_'+Label+'_M'+str(LargeEqCutoff)+'-'+str(minMag)+'.txt'

#DayFileName = ConfigDir+'BlancoANT_Test_12h.txt'
#DayFileName = ConfigDir+'BlancoNoise_'+Label+'_M'+str(minMag)+'.txt'
#for azi part:
DayFileName = ConfigDir+'BlancoNoise_'+Label+'_M'+str(LargeEqCutoff)+'-'+str(minMag)+'.txt'

#ANTDayFileName = ConfigDir+'BlancoANT_Test_12h.txt'
# two station lists, one for BXH, one for all stations
XStafn = ConfigDir+'X9_stations_X.txt'
#Stafn = ConfigDir+'X9_stations.txt'# all blanco stations
#Stafn = ConfigDir+'7D_stations_Gorda.txt' # only new gorda stations
Stafn = ConfigDir+'7D_stations_withGorda.txt' # all ci stations, including above gorda
#Stafn = ConfigDir+'G02BTest.txt' # test this station data
#Stafn = ConfigDir+'SC_stations.txt' # stations in so cal for testing

#  Event download
webservice = "IRIS"
#webservice = "SCEDC"
#network = "X9"  # X9 = Blanco
network = "7D"  # 7D = Cascadia Initiative
#network = "CI" # socal station code - confusing, this is not cascadia initiative
isCMT_params = 1  # use GCMT parameters for SAC header; 0 = use IRIS
isCentroid = 1  # if isCMT_params = 1, use centroid; 0 = epicentral

ArrayLoc = [43.75,-128.5]

#  Noise Download
trLen = 60 * 60 * 24  # seconds
#trLen = 60 * 60 * 12  # seconds
noDays = 4  # number of days prior to event to use
isCalDay = 1  # 0 to start each day at 00:00; 0 to use 24h segments prior to eq

tstart = '2012-10-09T00:00:00'
tend = '2012-10-12T23:00:00'

tstart = '2013-07-10T00:00:00'
tend = '2013-08-17T23:59:59'

tstart = '2013-02-01T00:00:00'
tend = '2013-02-28T23:59:59'

#tstart = '2013-07-21T00:00:00'
#tend = '2013-07-21T10:00:00'

# for new mac test - only one m6.5+
tstart = '2013-02-13T00:00:00'
tend = '2013-02-14T23:59:59'

# G02B test... why no seismic data?
#tstart = '2012-09-26T00:00:00'
#tend = '2012-09-26T23:59:59'

# the entire x9 experiment time
tstart = '2012-09-18T00:00:00'
tend = '2013-10-05T23:59:59'


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
