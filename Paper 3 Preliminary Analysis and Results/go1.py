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

#-- Read in the data
Name = 'H-H1_LOSC_4_V1-959344640-4096'
fileName = Name + '.hdf5'
raw_strain, raw_time, dq = rl.loaddata(fileName)
fs = int(1. / (raw_time[1] - raw_time[0]))

#-- Read in the template
f_template = h5py.File('GW150914_4_template.hdf5', "r")
template_p, template_c = f_template['template'].value
temp_strain = template_p + template_c * 1.j


#-- Find the injection segment
segList = rl.dq_channel_to_seglist(dq['HW_CBC'])
print 'Please choose the segment you are goint to investigate.'
print 'Please input an integer.\nRefer to: "injection info.txt"'

#Plotting the dq flag
plt.figure()
plt.plot(np.arange(0, raw_time.size / fs), dq['HW_CBC'],\
                 label = 'HW_CBC')
plt.ylim(-1,2)
plt.legend(loc = 1)
plt.show()
#Making the choice
num = input()
try:
    seg = segList[num - 1]
except:
    raise ValueError('Not within one of the segments')
strain = raw_strain[seg]
time = raw_time[seg]

#-- Plot the ASD
    #PSD first:
NFFT = fs * 4
NOVL = NFFT / 2
psd_window = np.blackman(NFFT)
Pxx, fpsd = mlab.psd(raw_strain, Fs=fs, NFFT=NFFT, window=psd_window\
                     ,noverlap = NOVL)
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
optimal = strain_fft * temp_fft.conjugate() / power_vec #vs. frequency
optimal_time = 2 * np.fft.ifft(optimal) * fs #vs. time
#Normalize:
df = np.abs(seg_freq[1] - seg_freq[0])
sigmasq = 1 * (temp_fft * temp_fft.conjugate() / power_vec).sum() * df
sigma = np.sqrt(np.abs(sigmasq))
SNR = abs(optimal_time / sigma)

#-- Plot it:
plt.figure()
time_shifted = time - time[0]
plt.plot(time_shifted, SNR)
plt.xlabel('Seconds after {0}'.format(time[0]))
plt.title('Optimal Matched Filter')
plt.show()

#-- Print the result
Max = SNR.max()
found_time = time[np.where(Max == SNR)] + 15 #The peak in template
print 'Recovered SNR: {0}.'.format(Max)
print 'Recovered at GPS time {0}.'.format(found_time[0])

#-- Record the result
print 'Do you want to record the result? Yes[y] No, anykey'
decision = raw_input()
if decision == 'y':   
    filename = Name + ' {0}.txt'.format(num + 1)
    f = open(filename, 'w')
    f.write('Starting GPS time    SNR\n')
    f.write('{0}    {1}\n'.format(found_time[0], Max))
    f.close()

