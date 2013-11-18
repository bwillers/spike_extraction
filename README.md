spike_extraction
================

Simple spike extraction script in python.

Takes as input a series of matlab files, each of which contains a variable 'LFPVoltage' representing
the voltage on a single recording channel. Saves files to .spike format to be read in pyclust.

Usage:                                                                                        
    spike_extraction.py [options] \<inputfiles\>...
                                                                                              
Options:                                                                                      
    --output=\<outfilename\>   Output filename [default: output.spike]                          
    --samplingfreq=\<fs\>      Sampling frequency in Hz [default: 25000]                        
    --segmentlength=\<mins\>   Size of data chunks (in min) to process at a time [default: 10]  
    --waveformlength=\<ms\>    Length of waveforms (in ms) to extract [default: 1]              
    --threshold=\<sigmamult\>  Multiple of std to use as threshold [default: 5]                 
    --filter_low=\<lowcut\>    Filter range low cut (Hz) [default: 600]                         
    --filter_high=\<highcut\>  Filter range high cut (Hz) [default: 6000]                       

Sample usage:

   ./spike_extraction.py --output=foo.spike LFPvoltage_ch1.mat LFPvoltage_ch2.mat LFPvoltage_ch3.mat LFPvoltage_ch4.mat
