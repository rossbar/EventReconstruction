import numpy as np
from analysisHelperFunctions import onAC, onDC, inge1, inge2,\
                                    checkForEnergyMatch

#class PotentialEvent(object):
#    '''Class for encapsulating a potential gamma ray event. Consists only of
#       time-correlated readouts before any reconstruction takes place.'''
#
#    def __init__(self, ev):
#        '''Take time-correlated readouts and break them up by detector and
#           side.'''
#        ge1, ge2 = inge1(ev), inge2(ev)

class PotentialSingleDetectorEvent(object):
    '''Class encapsulating all of the time-correlated strip fires from in a
       single detector'''

    def __init__(self, ev):
        '''Usage:
             potential_event = PotentialSingleDetectorEvent(ev)
           input:
             ev is a time-correlated group of readout strips.
             ev = ndarray with edata dtype; has at least the following columns:
                 timestamp
                 ADC_value
                 detector
                 trigger
           Description:
             Creates an object that splits the time-correlated input by detector
             side (AC/DC) and determines whether or not the potential event is
             reconstructible by checking that a) There is at least one readout
             per side and b) That the total energy collected by each side is
             commensurate.'''
        # Initial values
        self.ac = None
        self.dc = None
        self.can_reconstruct = None
        self.type = None

        # Split up readouts based on detector side. If there are no readouts on
        # one of the sides, set the error case appropriately
        ac, dc = onAC(ev), onDC(ev)
        if len(ac) > 0:
            self.ac = ac
        else:
            self.can_reconstruct = False
            self.type = 'DC_ONLY'
        if len(dc) > 0:
            self.dc = dc
        else:
            self.can_reconstruct = False
            self.type = 'AC_ONLY'

        # Check the total readout energies to see if the event can be
        # reconstructed
        if self.can_reconstruct is True or self.can_reconstruct is None:
            tot_en_AC, tot_en_DC = self.ac['ADC_value'].sum(),\
                                   self.dc['ADC_value'].sum()
            if checkForEnergyMatch(tot_en_AC, tot_en_DC):
                self.can_reconstruct = True
            else:
                self.can_reconstruct = False
                self.type = 'ENERGY_MISMATCH'
            # Total energy difference between the two electrods
            self.dE = np.abs(tot_en_AC - tot_en_DC)

    def determine_event_type(self):
        '''Simple pattern matching to determine the type of the event. Can 
           implement more complicated matching patterns as needed.'''
        # If there is only one readout per side, it is a 1-1 interaction
        if len(self.ac) == 1 and len(self.dc) == 1:
            self.type = '1-1'
            self.can_reconstruct = True
        # If there is one readout on one side and two on the other, it is a 1-2
        # interaction. These can arise from several scenarios
        elif len(self.ac) == 1 and len(self.dc) == 2:
            self.type = '1-2'
            self.double_side = self.dc
            self.single_side = self.ac
        elif len(self.ac) == 2 and len(self.dc) == 1:
            self.type = '1-2'
            self.double_side = self.ac
            self.single_side = self.dc
        # If there are two readouts on each side, 2-2 interaction
        elif len(self.ac) == 2 and len(self.dc) == 2:
            self.type = '2-2'
        # If there are more than two readouts per side, it is considered a high
        # order event, for now
        else:
            self.type = 'Event multiplicity > 2'
