#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 26 12:53:10 2022

@author: ayush
"""

import numpy as np
# import matplotlib.pyplot as plt
from scipy.io import wavfile
from am_analysis import am_analysis as ama

    

def modSpec(x,fs):
    win_size_sec = 0.04  # window length for the STFFT (seconds)
    win_shft_sec = 0.01  # shift between consecutive windows (seconds)


    stft_modulation_spectrogram = ama.strfft_modulation_spectrogram(x, fs, win_size = round(win_size_sec*fs), win_shift = round(win_shft_sec*fs), channel_names = x_name)
    X_plot=ama.plot_modulation_spectrogram_data(stft_modulation_spectrogram, 0 , modf_range = np.array([0,20]), c_range =  np.array([-90, -50]))


    return X_plot


if __name__ == "__main__":
    
    # speech signal
    # The speech signal p234_004.wav is one sample from the: 
    # CSTR VCTK Corpus: English Multi-speaker Corpus for CSTR Voice Cloning Toolkit
    # avialable in: https://datashare.is.ed.ac.uk/handle/10283/2651
    fs, x = wavfile.read('/home/ayush/00001.wav')
    y=x
    x_name = ['speech']
    x = x / np.max(x)
    # 1s segment to analyze
    x = x[int(fs*1.6) : int(fs*3.6)]
    
    X_plot=modSpec(x,fs)