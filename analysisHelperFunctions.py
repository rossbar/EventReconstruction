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
  dets = ev.detector
  return ev[(dets < NUMCH) | ((dets >= NUMCH*2)&(dets < NUMCH*3)) ]

def inge2( ev ):
  dets = ev.detector
  return ev[((dets >= 1*NUMCH)&(dets < 2*NUMCH)) | (dets >= NUMCH*3) ]

def onAC( ev ):
  return ev[ev.detector >= 2*NUMCH]

def onDC( ev ):
  return ev[ev.detector < 2*NUMCH]

# Check for energy matching
def checkForEnergyMatch( en1, en2, sigma=2 ):
  '''Check to see if the two energies agree to within sigma*sqrt of the 
     maximum input energy.'''
  maxEnergy = max( (en1, en2) )
  other = min( (en1, en2) )
  if other >= maxEnergy - sigma*np.sqrt( maxEnergy ): retVal = True
  else: retVal = False
  return retVal

def convertToInteraction( pod, interpFns ):
  '''Converts a 1-1 readout pod to an interaction.'''
  if type(pod) == np.ndarray: pod = pod.view(np.recarray)
  ac, dc = onAC( pod ), onDC( pod )
  # Energy first
  energy = pod.ADC_value.max()
  dE = abs( ac.ADC_value - dc.ADC_value )
  # Figure out which detector the interaction was in
  isge1 = inge1( pod )
  if len(isge1) == 0:
    det = 2
    fn = interpFns[-1]
    detSep = DETECTOR_SEPARATION + 15
  else: 
    det = 1
    fn = interpFns[0]
    detSep = 0
  # Convert strips to x-y coords
  if det == 1:
    x = (ac.detector - 2*NUMCH) * PITCH - 1
    y = (dc.detector - 0*NUMCH) * PITCH - 1
  else:
    x = (ac.detector - 3*NUMCH) * PITCH - 1
    y = (dc.detector - 1*NUMCH) * PITCH - 1
  # Time differencing
  dT = float(ac.timestamp) - float(dc.timestamp)
  # Z interpolation 
  z = fn( dT ) + detSep
  return np.array([(energy, x, y, z, dT, dE, det)], dtype=interactionType )

def filterBadStrips( edata, stripFile ):
  '''
    Usage: filterBadStrips( edata, stripFile )
           edata = recarray of event data
           stripFile = string: path to file containing the numbers of the bad strips
    Given the name of a file that contains a list of the known bad strips, load the information and trim the given
    edata. Return the edata minus readouts from the strips specified in the stripFile.
  '''
  # Load the strips to filter from the specified file
  badStrips = np.loadtxt(stripFile)
  # Construct the string to be executed from the badStrips data
  cmdStr = 'edata = edata['
  for s in badStrips: cmdStr += '(edata.detector != %i)&' %(s)
  cmdStr = cmdStr[:-1]+']'
  exec(cmdStr)
  return edata
