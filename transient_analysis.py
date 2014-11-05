import numpy as np

def get_transient_amplitude_simple(signal, trans_window=(20,80)):
    '''Given a baseline-subtraced, calibrated transient signal,
       return the amplitude of the maximum deviation from the baseline. The 
       search for signal deviations is limited to the signal region given 
       by trans_window'''
    sig_max = signal[trans_window[0]:trans_window[1]].max()
    sig_min = np.abs(signal[trans_window[0]:trans_window[1]].min())
    if sig_max > sig_min: trans_amp = sig_max
    else: trans_amp = sig_min
    return trans_amp

def get_transient_noise_simple(signal, trans_window=(20,80)):
    '''Return the sample standard deviation of the signal outside of the 
       transient window.'''
    non_trans_signal = np.concatenate((signal[0:trans_window[0]], signal[trans_window[1]:]))
    return non_trans_signal.std()
