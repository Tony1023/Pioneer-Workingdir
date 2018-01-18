# Programmer: Tony Lu Zhehao
# From: 2AP2 International Cirriculum Center,
# Shenzhen Foreign Languages School, Shenzhen, China
# This script is to make and recover injections on LIGO's raw data files.
# It generally randomly picks 5 points in each data file to work on.
# The recovery algorithm is adapted from the tutorial of GW150914 on LOSC.
# It has been very slightly altered, but 99% of it is the same.
# URL: https://losc.ligo.org/s/events/GW150914/LOSC_Event_tutorial_GW150914.html
# NOTICE:
# Highly recommend backing up the files, because the injections will not be
# easily removed.
# The output is the average SNR of each injection.

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.signal as sig
import readligo as rl
from random import randint


#-- Sampling rate equals to 4096 Hz
fs = 4096

#-- Read in the waveform
temp_time, temp_strain = np.genfromtxt('GW150914_4_NR_waveform.txt').transpose()

#-- Read in the template
f_template = h5py.File('GW150914_4_template.hdf5', "r")
template_p, template_c = f_template['template'].value
template = template_p + template_c * 1.j

# to remove effects at the beginning and end of the data stretch, window the data
# https://en.wikipedia.org/wiki/Window_function#Tukey_window
try:   dwindow2 = sig.tukey(template.size, alpha=1./8)  # Tukey window preferred, but requires recent scipy version 
except: dwindow2 = sig.blackman(template.size)          # Blackman window OK if Tukey is not available

# the length and sampling rate of the template MUST match that of the data.
datafreq = np.fft.fftfreq(template.size)*fs
df = np.abs(datafreq[1] - datafreq[0])
# prepare the template fft.
template_fft = np.fft.fft(template*dwindow2) / fs

# -- To calculate the PSD of the data, choose an overlap and a window (common to all detectors)
#   that minimizes "spectral leakage" https://en.wikipedia.org/wiki/Spectral_leakage
NFFT = 4*fs
psd_window = np.blackman(NFFT)
# and a 50% overlap:
NOVL = NFFT/2

#-- Set the amplifications of signal to be injected
A = [0,1,2] + range(4,31,2)

#-- Define the function to segment the segments into suitable segments
def slice_filter(dq, hw):
    #-- Cut the whole segment in to small ones of 60s, with spacing of 1s
    i = 0
    while i < len(dq):
        if dq[i].stop - dq[i].start < 32*fs:
            del dq[i]
            continue
        dq.insert(i, slice(dq[i].start, dq[i].start + 32*fs))
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

#-- Prepare the files to work on
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


#===========================#
#========Main Part==========#
#===========================#
for name in fileName:

    #-- Prepare the document to record the table (SNR vs Ampli.)
    f = open(name + '.txt', 'w')

    #-- Load data from the original file
    rawFile = h5py.File(name + '.hdf5', 'r')
    strain_raw = rawFile['strain/Strain'].value
    CBCHIGH_CAT4 = (rawFile['quality/simple/DQmask'].value >> 4) & 1
    HW_CBC = (rawFile['quality/injections/Injmask'].value >> 0) & 1
    GPSstart = rawFile['meta/GPSstart'].value
    rawFile.close()

    #-- Load the file again which will be injected
    dataFile = h5py.File(name + '.hdf5', 'r+')
    strain = dataFile['strain/Strain']
    dqInj = dataFile['quality/injections/Injmask']

    #-- Getting the suitable segement lists
    segList = rl.dq_channel_to_seglist(CBCHIGH_CAT4)
    segList_HW = rl.dq_channel_to_seglist(HW_CBC)
    segList = slice_filter(segList, segList_HW)

    #-- Pick 10 (or if not more than 10, all) segments randomly
    inj_num = []
    if len(segList) <= 5:
        inj_num = range(0, len(segList))
        #-- If less than 5, throw a caution
        print 'Caution: there are less than five suitable spots in %s!' % name
        print '{0} only!'.format(len(segList))
    else:
        while len(inj_num) < 5:
            i = randint(0, len(segList)-1)
            if i in inj_num:
                continue
            inj_num.append(i)
    inj_num.sort()

    List = []
    for i in inj_num:
        List.append(segList[i])

    #-- Make and recover injections at each amplitude
    for k in range(0, len(A)):
        f.write('{0} '.format(A[k]))
        
        #-- Amplify the template with the difference to next amplification
        if k == 0:
            temp = temp_strain * (A[k])
        else:
            temp = temp_strain * (A[k] - A[k-1])

        #-- Make the injection to every piece
        SNRsum = 0
        
        for seg in List:
            #-- Injection starts at the middle
            inj_sample = seg.start + 16*fs
            #-- Superpose the waveform to the signal
            for j in range(0, temp.size):
                strain[inj_sample + j] += temp[j]
                    
            #-- Injection Recovery --#
            #using the segments in the above section
            data = strain[seg]

            # to remove effects at the beginning and end of the data stretch, window the data
            # https://en.wikipedia.org/wiki/Window_function#Tukey_window
            try:   dwindow = sig.tukey(data.size, alpha=1./8)  # Tukey window preferred, but requires recent scipy version 
            except: dwindow = sig.blackman(data.size)          # Blackman window OK if Tukey is not available

            #-- Calculate the PSD of the data.  Also use an overlap, and window:
            # This is where the only change is made, I replaced the PSD of the 32-s segment with
            # that of the whole file, because when the file is small the PSD will be greatly affected by the
            # signal in the file, whereas we only want the PSD of the background noise. So taking the PSD
            # of a much larger interval can help eliminating the effect.
            data_psd, freqs = mlab.psd(strain_raw, Fs = fs, NFFT = NFFT, window=psd_window, noverlap=NOVL)

            # Take the Fourier Transform (FFT) of the data and the template (with dwindow)
            data_fft = np.fft.fft(data*dwindow) / fs

            #-- Interpolate to get the PSD values at the needed frequencies
            power_vec = np.interp(np.abs(datafreq), freqs, data_psd)

            #-- Calculate the matched filter output in the time domain:
            # Multiply the Fourier Space template and data, and divide by the noise power in each frequency bin.
            # Taking the Inverse Fourier Transform (IFFT) of the filter output puts it back in the time domain,
            # so the result will be plotted as a function of time off-set between the template and the data:
            optimal = data_fft * template_fft.conjugate() / power_vec
            optimal_time = 2*np.fft.ifft(optimal)*fs

            #-- Normalize the matched filter output:
            # Normalize the matched filter output so that we expect a value of 1 at times of just noise.
            # Then, the peak of the matched filter output will tell us the signal-to-noise ratio (SNR) of the signal.
            sigmasq = 1*(template_fft * template_fft.conjugate() / power_vec).sum() * df
            sigma = np.sqrt(np.abs(sigmasq))
            SNR_complex = optimal_time/sigma

            # shift the SNR vector by the template length so that the peak is at the END of the template
            peaksample = int(data.size / 2)  # location of peak in the template
            SNR_complex = np.roll(SNR_complex,peaksample)
            SNR = abs(SNR_complex)

            #-- Find the time and SNR value at maximum:
            indmax = np.argmax(SNR)
            SNRmax = SNR[indmax]
            SNRsum += SNRmax

        SNR_avr = SNRsum / len(inj_num)
        f.write('{0}\n'.format(SNR_avr))
        dataFile.flush()

        print '%s at amplification = {0} finished.'.format(A[k]) % name

    dataFile.close()
    f.close()



            
            
        

    


