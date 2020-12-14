#!/urs/bin/env python 3.7

from subprocess import call
import sys
import re
import os

import config
dataDir = 'data/CORRSEIS_SAC/'
eventDirs = os.listdir(dataDir) 
eventDirs = [x for x in eventDirs if not x.startswith('.')] # get rid of hidden dirs
MatlabConfigFile = '03_ASWMS/setup_parameters.m'

# first get the periods

LineName = 'parameters.periods'
periods = []

MCF = open(MatlabConfigFile,'r')

for line in MCF:
    if re.match(LineName,line):
        
        # get rid of everything after comment sign
        periodsStr = line.split(';')[0]
        # get rid of everything before equals sign
        periodsStr = periodsStr.split('=')[-1]
        # get rid of brackets
        periodsStr = periodsStr.strip("[ ]")
        for i in periodsStr.split(' '):
            periods.append(int(i))

        # let's just take the first one
        break

# Now we want to set up the filters, which have a width ~10% of center freq.

# higher frequency = shorter period
HiCorner = []
LoCorner = []
for i in periods:
    fac = 0.05 * i
    hi = i - fac
    lo = i + fac
    HiCorner.append(1/hi)
    LoCorner.append(1/lo)

# Now we want to run the sac macro for each of the event dirs
for event in eventDirs:
    os.chdir(dataDir+event+'/')
    print(periods)
    print(len(periods))
    for i in range(0,len(periods)):
        T = periods[i]
        print(T)
        H = HiCorner[i]
        L = LoCorner[i]
        ofn = 'rec_section_'+str(T)+'s'

        call(['../../../sac_rec_sections.csh',str(L),str(H),ofn])

    os.chdir('../../..')
