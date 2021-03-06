from __future__ import division
import numpy as np

# dtypes
interactionType = np.dtype({'names':['energy', 'x', 'y', 'z', 'dT', 'dE', 'det'],\
                            'formats':[np.float32, np.float32, np.float32, np.float32,\
                                       np.float32, np.float32, np.uint8]}, align=True )
readoutType = np.dtype( [('timestamp', np.uint64),('ADC_value', '<f4'), \
              ('detector', np.float), ('trigger', '<u2'), ('pileup','<u2'), \
              ('retrigger', '<u2'), ('rid', '<u4')] )
# Detector strip mappings
PITCH = 2 # mm
DETECTOR_SEPARATION = 10 # mm
NUMCH = 38

# Functions for splitting event data by detector and side
def inge1( ev ):
    dets = ev['detector']
    return ev[(dets < NUMCH) | ((dets >= NUMCH*2)&(dets < NUMCH*3)) ]

def inge2( ev ):
    dets = ev['detector']
    return ev[((dets >= 1*NUMCH)&(dets < 2*NUMCH)) | (dets >= NUMCH*3) ]

def onAC( ev ):
    return ev[ev['detector'] >= 2*NUMCH]

def onDC( ev ):
    return ev[ev['detector'] < 2*NUMCH]

# Energy matching
def determineEnergyDifference(en1, en2):
    return np.abs(en1 - en2)

def checkForEnergyMatch(en1, en2, sigma=2):
    '''Check to see if the two energies agree to within sigma*sqrt of the 
       maximum input energy.'''
    maxEnergy = max( (en1, en2) )
    other = min( (en1, en2) )
    if other >= maxEnergy - sigma*np.sqrt( maxEnergy ): retVal = True
    else: retVal = False
    return retVal

# Trigger correction
def correct_trigger_flags(ev, transient_en_cut=10.0):
    '''Use results of slow filter (ADC_value) in conjunction with the existing
       trigger value to better determine whether charge was collected or not 
       on each strip in the ev. 
       
       transient_en_cut determines the energy limit at which a signal is 
       considered a transient or not'''
    trig_mask = ev['ADC_value'] >= transient_en_cut
    # Values below the cut are likely transients
    ev['trigger'][np.invert(trig_mask)] = 0
    ev['trigger'][trig_mask] = 1
    return ev
