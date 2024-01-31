import sys

from baseband import vdif
import numpy as np
from matplotlib import pyplot as plt

from tqdm import tqdm

import multiprocessing as mp

channel_order = {"8": "CH01", "12": "CH02", "0": "CH03", "4": "CH04", "9": "CH05", "13": "CH06", "1": "CH07",
                 "5": "CH08", "10": "CH09", "14": "CH010", "2": "CH11", "6": "CH12", "11": "CH13", "15": "CH14",
                 "3": "CH15", "7": "CH16"}


channel_polarization = {"CH01": "RCP", "CH02": "LCP", "CH03": "RCP", "CH04": "LCP", "CH05": "RCP", "CH06": "LCP",
                        "CH07": "RCP", "CH08": "LCP", "CH09": "RCP", "CH010": "LCP", "CH11": "RCP", "CH12": "LCP",
                        "CH13": "RCP", "CH14": "LCP", "CH15": "RCP", "CH16": "LCP"}

def order_chan(psd_in):
    psd_out = np.zeros((Nf+1, 16))
    for chan in range(0, 16):
        channale_output_index = int(channel_order[str(chan)].replace('CH', '')) - 1
        psd_out[:, channale_output_index] = psd_in[:, chan]
    return psd_out

    
def combine_polarazations(psd_in):
    chan_indices = np.arange(0, 8)
    a_indices = 2 * chan_indices
    b_indices = a_indices + 1

    return (psd_in[:, a_indices] + psd_in[:, b_indices]) / 2

def run_fft(index):
    Psd = np.zeros((16, Nf+1))
    Sxx = np.zeros((16, 129))
                      
    fh = vdif.open(file_name, 'rs')
    fh.seek(index * sample_in_second, 1)
    one_second_data = fh.read(sample_in_second)
    
    a = 0
    b = fs
    for n in range(0, ns):
        one_second_data_n = one_second_data[a:b]
    
        for chan in range(0, 16):
             
             Psd[chan, :] += plt.psd(one_second_data_n[:, chan], NFFT=Nfft, Fs=Fs, Fc=fc-bw/2, detrend='mean', window=np.hanning(Nfft), noverlap=Nf)[0]/ns
             
        a += fs
        b += fs
    
    Psd = Psd.T                    # shape(Nf+1,Nrt)
    Psd = np.flipud(Psd)           # Psd(Nf+1,Nrt), band is mirrorred due to mixing scheme
    plt.close()        
    fh.close()

    return Psd



file_name = "/mnt/VLBI/data/stefk1/stefk1_ir_no0001.m5a"
fh = vdif.open(file_name, 'rs')

print(fh.info)
print("\n\n")

nr_samples = fh.info.shape[0]
start_time = fh.start_time
stop_time = fh.stop_time
duration = (stop_time.mjd - start_time.mjd) * 24 * 60 * 60
sample_in_second = int(nr_samples / duration)

fs = 1024  * 128 
ns = sample_in_second // fs    


print("duration [s]:", duration)
print("start time:", start_time)
print("stop_time:", stop_time)

Nfft = 2**8                # number of FFT bins    
Nf = int(Nfft/2)           # number of frequency bins
Fs = sample_in_second      # sample rate, samples per second (Nyquist, band is flipped)
bw = Fs/2                  # bandwidth in Hz
df = bw/Nf                 # channel width in Hz
fc = 6644000000

fh.close()


seconds = duration
inc = 20

min_time = 0
max_time = inc

results = []
indexs = np.arange(min_time, max_time, 1, dtype=int)

for i in tqdm(range(0, int(seconds / inc))):
    with mp.Pool(46) as pool:
        results.extend(pool.map(run_fft, indexs))
        min_time += inc
        max_time += inc
    
pool.close()
pool.join()


results_array = np.array(results)

results_array = np.array([order_chan(result) for result in results_array])
results_array = np.array([combine_polarazations(result) for result in results_array])

fig, ax = plt.subplots(nrows=2, ncols=4, figsize=(16, 16), dpi=100)
fig2, ax2 = plt.subplots(nrows=2, ncols=4, figsize=(16, 16), dpi=100)
plt_index = np.arange(0, 8)

freq = [np.linspace(6636.49, 6644.49, Nf+1), np.linspace(6644.49, 6652.49, Nf+1), np.linspace(6652.49, 6660.49, Nf+1), np.linspace(6660.49, 6668.49, Nf+1), np.linspace(6668.49, 6676.49, Nf+1), np.linspace(6676.49, 6684.49, Nf+1), np.linspace(6684.49, 6692.49, Nf+1), np.linspace(6692.49, 6700.49, Nf+1)]


for i, (ax1, ax2) in enumerate(zip(ax.flat, ax2.flat)):
    data_tmp = results_array[:, :, plt_index[i]].T
    
    ax1.plot(freq[plt_index[i]], np.median(data_tmp, axis=1))
    ax2.imshow(data_tmp, aspect="auto", extent=[0, seconds, freq[plt_index[i]][0], freq[plt_index[i]][1]], vmin=np.percentile(data_tmp, 1), vmax=np.percentile(data_tmp, 95))


plt.show()

sys.exit(0)

