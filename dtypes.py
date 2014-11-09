import numpy as np

interactionType = np.dtype({'names':['energy', 'x', 'y', 'z', 'dT', 'dE', 'det'],\
                            'formats':[np.float32, np.float32, np.float32, np.float32,\
                                       np.float32, np.float32, np.uint8]}, align=True )
readoutType = np.dtype( [('timestamp', np.uint64),('ADC_value', '<f4'),\
              ('detector', np.float), ('trigger', '<u2'), ('pileup','<u2'),\
              ('retrigger', '<u2'), ('rid', '<u4')] )
full_readout_type = np.dtype( [('timestamp', '<u8'), ('ADC_value', '<f4'), ('detector', '<u2'), ('trigger', '<u2'), ('signal', (np.float, 256))], align=True )

refined_edata_type = np.dtype({'names':['timestamp', 'ADC_value', 'detector',\
                                          'trigger', 't50', 'R', 'sigma_R'],\
                                 'formats':[np.uint64, np.float32, np.uint16,\
                                            np.uint16, np.float32, np.float32,\
                                            np.float32]}, align=True )

def merge_edata_with_signals(edata, rdata):                                     
    out = np.zeros(edata.shape[0], dtype=full_readout_type)                     
    out['timestamp'] = edata['timestamp']                                       
    out['ADC_value'] = edata['ADC_value']                                       
    out['detector'] = edata['detector']                                         
    out['trigger'] = edata['trigger']                                           
    out['signal'] = rdata                                                       
    return out

def create_refined_edata(edata):
    '''Convert an array with edata dtype to a refined edata by copying all
       shared fields.'''
    out = np.zeros(edata.shape[0], dtype=refined_edata_type)
    out['timestamp'] = edata['timestamp']                                       
    out['ADC_value'] = edata['ADC_value']                                       
    out['detector'] = edata['detector']                                         
    out['trigger'] = edata['trigger']                                           
    return out
