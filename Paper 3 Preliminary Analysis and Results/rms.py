# Programmer: Zhehao Lu Tony
# 2AP2, International Curriculum Center,
# Shenzhen Foreign Languages School, China
# Date: July 6, 2016
# For calculating the RMS

from math import sqrt
import h5py
import readligo as rl

#-- Read in the names
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

#-- Define the function for RMS
def qmean(a):
    return sqrt(sum(n**2 for n in a) / len(a))

#-- Prepare the document to record the RMSs
f = open('rms.txt','w')

for Name in fileName:
    
    #-- Read in the data files
    dataFile = h5py.File(Name + '.hdf5', 'r')
    GPSstart = dataFile['meta/GPSstart'].value
    strain = dataFile['strain/Strain'].value

    RMS = qmean(strain)
    print '%s\'s RMS = {0}'.format(RMS) % Name
    f.write('{0} {1}\n'.format(GPSstart, RMS))

f.close()
