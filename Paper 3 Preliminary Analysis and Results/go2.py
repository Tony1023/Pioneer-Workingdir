# Programmer: Zhehao Lu Tony
# 2AP2, International Curriculum Center,
# Shenzhen Foreign Languages School(Shenzhen, China)
# Date: June 30, 2016
# For recovering the self-made injections in raw LIGO data


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import scipy.signal as sig
import readligo as rl
import h5py 

#-- Read in the template
f_template = h5py.File('GW150914_4_template.hdf5', "r")
template_p, template_c = f_template['template'].value
temp_strain = template_p + template_c * 1.j

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

#-- Prepare the document to record the SNRs
print 'For what level of amplification are you recovering?'
amplification = input()
f = open('SNR_M={0}.txt'.format(amplification), 'w')

for Name in fileName:
    print Name + ':'
    #-- Load the data
    fileName = Name + '.hdf5'
    raw_strain, raw_time, dq = rl.loaddata(fileName)
    fs = int(1. / (raw_time[1] - raw_time[0]))

    #-- Prepare the segment list
    segList = rl.dq_channel_to_seglist(dq['HW_CBC'])

    #-- Prepare the txt file to record the result
    filename = Name + '.txt'

    for seg in segList:
        
        #-- If not self-injection, skip it
        if seg.stop - seg.start == 30:
            continue
        strain = raw_strain[seg]
        time = raw_time[seg]

        #-- Compute the PSD
        NFFT = fs * 4
        NOVL = NFFT / 2
        psd_window = np.blackman(NFFT)
        Pxx, fpsd = mlab.psd(raw_strain, Fs=fs, NFFT=NFFT, \
                             window=psd_window, noverlap = NOVL)
        Pxx_temp, fpsd_temp = mlab.psd(temp_strain, Fs=fs, NFFT=NFFT)

        #-- Optimal matched filter
        dwindow = sig.tukey(strain.size, alpha = 1./8)
        strain_fft = np.fft.fft(strain * dwindow) / fs
        temp_padded = np.append(temp_strain,\
            np.zeros(strain.size - temp_strain.size))
        dwindow2 = sig.tukey(temp_padded.size, alpha = 1./8)
        temp_fft = np.fft.fft(temp_padded * dwindow2)
        seg_freq = np.fft.fftfreq(strain.size) * fs
        #Interpolate to get the PSD values at the needed frequencies
        power_vec = np.interp(np.abs(seg_freq), fpsd, Pxx)
        #Op.Ma.Fi.:
        optimal = strain_fft * temp_fft.conjugate() / power_vec
        optimal_time = 2 * np.fft.ifft(optimal) * fs
        #Normalize:
        df = np.abs(seg_freq[1] - seg_freq[0])
        sigmasq = 1 * (temp_fft * temp_fft.conjugate() / \
                       power_vec).sum() * df
        sigma = np.sqrt(np.abs(sigmasq))
        SNR = abs(optimal_time / sigma)

        #-- Print the result
        Max = SNR.max()
        found_time = time[np.where(Max == SNR)] + 15
        print 'Recovered SNR: {0} at {1}'.format(Max,found_time[0])

        #-- Record the result
        f.write('{0} {1} {2}\n'.format(int(raw_time[0]),\
                                       found_time[0], Max))
    
f.close()
