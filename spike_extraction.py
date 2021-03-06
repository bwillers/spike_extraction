#!/usr/bin/python
"""Simple spike extraction. Assumes input .mat files contain variable LFPVoltage,
   with each channel a separate file. Does simple threshold crossing detection,
   followed by alignment on centre of mass.

Usage:
    spike_extraction.py [options] <inputfiles>...
    spike_extraction.py [options] --inputdir=<inputdir>

Options:
    --output=<outfilename>   Output filename [default: output.spike]
    --samplingfreq=<fs>      Sampling frequency in Hz [default: 25000]
    --segmentlength=<mins>   Size of data chunks (in min) to process at a time [default: 10]
    --waveformlength=<ms>    Length of waveforms (in ms) to extract [default: 1]
    --threshold=<sigmamult>  Multiple of std to use as threshold [default: 5]
    --filter_low=<lowcut>    Filter range low cut (Hz) [default: 600]
    --filter_high=<highcut>  Filter range high cut (Hz) [default: 6000]
    --invert                 Spikes are downward deflections (negate for peak).

"""

from __future__ import print_function
from __future__ import division

import os
import struct
import sys

from docopt import docopt

import numpy as np
import scipy.io
import scipy.signal
import datetime

def amplitudeThresholdSpikeExtraction(filtered_data, sigma_multiplier, width, padding, fs):
    signal = filtered_data.T
    signal = signal - np.mean(signal, axis=0)
    sigma = 1.4826 * np.median(np.abs(signal), axis=0)
    threshold = sigma_multiplier * sigma
    crossing, crossing_chan = np.nonzero((signal[1:,:] > threshold) & (signal[:-1,:] <= threshold))
    valid = np.ones(crossing.shape, dtype=np.bool)
    # crudely enforce the censor period
    temp = [crossing[0]]
    for c in crossing[1:]:
        if (c - temp[-1] > width) :
            temp.append(c)
    crossing = np.array(temp)
    # check that we dont over or underflow the sample length
    crossing = crossing[(crossing >= 2 * padding) & ( crossing <  signal.shape[0] - width)]
    # we'll extract and extra 5 samples on each side to do the alignment
    wv = np.zeros((crossing.size, filtered_data.shape[0], width + 2 * padding))
    for i, ind in enumerate(crossing):
        wv[i,:,:] = filtered_data[:,(ind-2*padding):(ind+width)]
    return wv, crossing

def alignWaveformsCoM(wv, pad):
    poss_shifts = np.arange(-pad, pad+1)
    p = np.sum(wv * wv, axis=1)
    p = (p.T / np.sum(p, axis=1)).T
    shifts = np.round(np.sum(p * np.arange(wv.shape[2]), axis=1)) - 3 * pad
    shifts[shifts < -pad] = -pad
    shifts[shifts > pad] = pad
    sp = np.zeros((wv.shape[0], wv.shape[1], wv.shape[2] - 2 * pad))
    for i in np.unique(shifts):
        w = shifts == i
        sp[w,:,:] = wv[w,:,(i + pad):(i + wv.shape[2] - pad)]
    return sp

def alignWaveformsPeak(wv, pad):
    poss_shifts = np.arange(-pad, pad+1)
    mu = np.mean(wv, axis=1)  # take mean across channels, left with NxL
    shifts = np.argmax(mu[:, pad:-pad], axis=1) -pad
    shifts[shifts < -pad] = -pad
    shifts[shifts > pad] = pad
    sp = np.zeros((wv.shape[0], wv.shape[1], wv.shape[2] - 2 * pad))
    for i in np.unique(shifts):
        w = shifts == i
        sp[w,:,:] = wv[w,:,(i + pad):(i + wv.shape[2] - pad)]
    return sp

if __name__ == '__main__':
    args = docopt(__doc__, version='Simple spike extractor 0.0.1')
    fs = int(args['--samplingfreq'])
    outf = args['--output']
    segment_length = float(args['--segmentlength'])
    wv_length = float(args['--waveformlength'])
    sigma_mult = float(args['--threshold'])
    filter_low_cut = float(args['--filter_low'])
    filter_high_cut = float(args['--filter_high'])
    wv_buffer_frac = 0.2
    inverted = args['--invert']

    width = int(np.ceil(wv_length * fs / 1.0e3))
    pad = int(np.ceil(wv_buffer_frac * width))

    print("Parameters:")
    print("\tSampling freq:", fs, "Hz")
    print("\tSegment length:", segment_length, "min")
    print("\tWaveform length:", wv_length, "ms (", width, " samples )")
    print("\tThreshold:", sigma_mult, "sigma")
    print("\tFilter range:", filter_low_cut, "-", filter_high_cut, "Hz")
    print("\tInvert:", inverted)
    print("\tInput files:", " ".join(args['<inputfiles>']))
    print("\tOutput file:", outf)

    print("")

    if args['<inputfiles>']:
        input_files = args['<inputfiles>']
    elif args['--inputdir']:
        print("Scanning {} for .mat input files".format(args['--inputdir']))
        input_files = sorted([os.path.join(os.path.expanduser(args['--inputdir']), f) for
            f in os.listdir(os.path.expanduser(args['--inputdir'])) if f.endswith('.mat')])
    else:
        raise ValueError("No input files specified.")

    print("Loading data")
    data = []
    for f in input_files:
        print("\t",f)
        data.append(scipy.io.loadmat(f)['LFPvoltage'][0,:])
    data = np.array(data)
    if inverted:
        data = -data
    data_split = np.array_split(data, np.ceil(data.shape[1] / (60 * fs * segment_length)), axis=1)

    fn = fs / 2
    b,a = scipy.signal.butter(4, [filter_low_cut / fn, filter_high_cut / fn], 'bandpass')

    print("")
    wv = []
    ts = []
    print("Processing...")
    for count, segment in enumerate(data_split):
        print("\tsegment ", count+1, "of", len(data_split), end=" ")
        print("...filtering", end=" ")
        segment_ft = scipy.signal.filtfilt(b,a, segment, axis=1)
        print("...extracting spike waveforms", end=" ")
        _w, _t = amplitudeThresholdSpikeExtraction(segment_ft, sigma_mult, width, pad, fs)
        print("...extracted", _w.shape[0], "waveforms")
        wv.append(_w)
        offset = 0
        for pre in data_split[:count]: offset = offset + pre.shape[1]
        ts.append(_t + offset)
    wv = np.concatenate(wv, axis=0)
    ts = ((np.concatenate(ts, axis=0) / fs) * 1e6).astype(np.uint64)

    print("")
    print("Aligning on peak")
    #sp = alignWaveforms(wv, pad)
    sp = alignWaveformsPeak(wv, pad)

    int_conversion = np.max(np.abs(sp)) / (2**15 - 2)
    sp_conv = (sp / int_conversion).astype(np.int16)

    def packStringBinary(s):
        rv = struct.pack('<I', len(s))
        if sys.version_info[0] == 3:
            return rv + s.encode('ascii')
        else:
            return rv + s

    print("")
    print("Writing output")
    f = open(outf, 'wb')
    # write version no
    f.write(struct.pack('<H', 1))
    # write num spikes
    f.write(struct.pack('<Q', sp.shape[0]))
    # write num channels
    f.write(struct.pack('<H', sp.shape[1]))
    # write num samples per waveform
    f.write(struct.pack('<H', sp.shape[2]))
    # write sampling frequency
    f.write(struct.pack('<I', int(fs)))
    # write align point
    f.write(struct.pack('<H', 6)) # hax
    # write a2d conversion factor (spikes are stored as int16 for space reasons)
    f.write(struct.pack('<' + ('d' * sp.shape[1]), *([int_conversion] * sp.shape[1])))
    # write datetime string - well jsut write todays date since we dont know when data was recorded
    f.write(packStringBinary(str(datetime.datetime.today())))
    # write subject string
    f.write(packStringBinary("unknown"))
    # write filter descrtipn
    f.write(packStringBinary("Butterworth4, %3.0f-%3.0f" % (filter_low_cut, filter_high_cut)))
    # write the spikes
    dt = np.dtype([('time', '<Q'), ('spikes', np.dtype('<h'), (sp.shape[2], sp.shape[1]))])
    temp = np.zeros((sp.shape[0],), dtype=dt)
    temp['time'] = ts
    temp['spikes'] = np.transpose(sp_conv, (0,2,1))
    temp.tofile(f)
    f.close()

    print("Extracted a total of ", ts.size, "spikes, stored in ", outf)
