# Programmer: Zhehao Lu Tony
# 2AP2, International Curriculum Center,
# Shenzhen Foreign Languages School, China
# Date: June 28, 2016
# For making software injections on raw LIGO signal data.
# Recommend make a copy of the original file and then do the injection,
# because the injections cannot be deleted easily.
# One can make injections on as many files one time as he/she wants
# by altering the 'fileName[]' array.

import h5py
import numpy as np
import matplotlib.pyplot as plt
import readligo as rl
from random import randint

#-- Sampling rate equals to 4096 Hz
fs = 4096

#-- Read in the waveform
temp_time, temp = np.genfromtxt('GW150914_4_NR_waveform.txt')\
                  .transpose()

#-- Set the amplification
print 'Please input the amplification of the signal.'
amplification = input()
temp *= amplification

#-- Define the function to segment the segments into suitable segments
def slice_filter(dq, hw):
    #-- Cut the whole segment in to small ones of 60s, with spacing of 1s
    i = 0
    while i < len(dq):
        if dq[i].stop - dq[i].start < 60*fs:
            del dq[i]
            continue
        dq.insert(i, slice(dq[i].start, dq[i].start + 60*fs))
        dq[i+1] = slice(dq[i].stop + 1*fs, dq[i+1].stop)
        i += 1
    #-- Delete the segments with injection activated
    i = 0
    j = 0
    if len(hw) == 0:
        return dq
    while i < len(dq) and j < len(hw):
        #-- dq is 100% before hw
        if dq[i].stop < hw[j].start:
            i += 1
        #-- dq is 100% after hw
        elif dq[i].start > hw[j].stop:
            j += 1
        #-- dq and hw overlap
        else:
            del dq[i]

    return dq

#-- Prepare the files to load
fileName = []
fileName.append('H-H1_LOSC_4_V1-931127296-4096')
fileName.append('H-H1_LOSC_4_V1-934846464-4096')
fileName.append('H-H1_LOSC_4_V1-941707264-4096')
fileName.append('H-H1_LOSC_4_V1-941785088-4096')
fileName.append('H-H1_LOSC_4_V1-947154944-4096')
fileName.append('H-H1_LOSC_4_V1-952623104-4096')
fileName.append('H-H1_LOSC_4_V1-959344640-4096')
fileName.append('H-H1_LOSC_4_V1-963629056-4096')
fileName.append('H-H1_LOSC_4_V1-967442432-4096')
fileName.append('H-H1_LOSC_4_V1-971407360-4096')

#-- Prepare the document to record the injections
f = open('Injection_M_{0}.txt'.format(amplification), 'w')

#===========================#
#========Main Part==========#
#===========================#
for name in fileName:

    inj_samp = []
    
    #-- Load the original file
    rawFile = h5py.File(name + '.hdf5', 'r')
    strain_raw = rawFile['strain/Strain'].value
    dq_raw = rawFile['quality/injections/Injmask'].value
    CBCHIGH_CAT4 = (rawFile['quality/simple/DQmask'].value >> 4) & 1
    HW_CBC = (rawFile['quality/injections/Injmask'].value >> 0) & 1
    GPSstart = rawFile['meta/GPSstart'].value
    rawFile.close()

    #-- Load the file again which will be injected
    dataFile = h5py.File(name + '.hdf5', 'r+')
    strain = dataFile['strain/Strain']
    dqInj = dataFile['quality/injections/Injmask']

    #-- Getting the segement lists
    segList = rl.dq_channel_to_seglist(CBCHIGH_CAT4)
    segList_HW = rl.dq_channel_to_seglist(HW_CBC)
    segList = slice_filter(segList, segList_HW)

    #-- Pick 10 (or if not more than 10, all) segments randomly
    inj_num = []
    if len(segList) <= 10:
        inj_num = range(0, len(segList))
        #-- If less than 10, throw a caution
        print 'Caution: there are less than ten suitable spots in %s!'\
                      % name
        print '{0} only!'.format(len(segList))
    else:
        while len(inj_num) < 10:
            i = randint(0, len(segList)-1)
            if i in inj_num:
                continue
            inj_num.append(i)
        inj_num.sort()

    #-- Make the injection to every piece
    for i in inj_num:
        seg = segList[i]
        #-- Get a random starting point
        inj_sample = randint(seg.start, seg.stop - temp.size)
        inj_samp.append(inj_sample)
        #-- Superpose the waveform to the signal
        for j in range(0, temp.size):
            strain[inj_sample + j] += temp[j]
        #-- Write onto dq flag: HW_CBC
        Start = seg.start / fs
        End = seg.stop / fs
        for j in range(Start,End):
            dqInj[j] = (1 << 0) | dqInj[j]
            dqInj[j] = (1 << 1) | dqInj[j]

    #-- Record the injections
    for sample in inj_samp:
        time = sample * 1. / fs
        f.write('{0} {1} {2} {3}\n'.format(GPSstart, GPSstart+time,\
                                           time, sample))
    print 'File %s injection at M={0} succeed.'.format(amplification)\
                                          % name

    #-- Close the file
    dataFile.flush()
    dataFile.close()

f.close()

    
