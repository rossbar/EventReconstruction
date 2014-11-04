import numpy as np
from analysisHelperFunctions import *
from error_codes import *

class PotentialEvent(object):
    def __init__(self, ev, refine_triggers=True):
        # Correct the trigger flags using info from the slow filter
        if refine_triggers: ev = correct_trigger_flags(ev)
        # Break data up by detector
        ge1, ge2 = inge1(ev), inge2(ev)
        self.ge1 = Detector(ge1)
        self.ge2 = Detector(ge2)

    def __str__(self):
        outstr = 'GeI\n'
        outstr += self.ge1.__str__()
        outstr += 'GeII\n'
        outstr += self.ge2.__str__()
        return outstr

    def has_orphans(self):
        return self.ge1.has_orphans() or self.ge2.has_orphans()

    def passes_energy_match(self):
        return self.ge1.passes_energy_match() and self.ge2.passes_energy_match()

    def clusterize(self):
        '''Call side clusterization.'''
        self.ge1.ac.clusterize_data()
        self.ge1.dc.clusterize_data()
        self.ge2.ac.clusterize_data()
        self.ge2.dc.clusterize_data()

class Detector(object):
    errors = []
    def __init__(self, ev):
        ac, dc = onAC(ev), onDC(ev)
        self.ac = Side(ac)
        self.dc = Side(dc)

    def __str__(self):
        outstr = '  AC\n'
        outstr += self.ac.__str__()
        outstr += '  DC\n'
        outstr += self.dc.__str__()
        return outstr

    def has_orphans(self):
        if (len(self.ac.data) == 0) and (len(self.dc.data) > 0):
            self.errors.append(DCONLY)
            return True
        elif (len(self.dc.data) == 0) and (len(self.ac.data) > 0):
            self.errors.append(ACONLY)
            return True
        return False

    def passes_energy_match(self):
        if checkForEnergyMatch(self.ac.total_energy, self.dc.total_energy):
            return True
        else: 
            self.errors.append(ENERGY_MATCH_FAILURE)
            return False


class Side(object):
    def __init__(self, ev):
        # Sort data by detector upon receiving
        ev.sort(order='detector')
        # Set recarray to data
        self.data = ev
        # Compute total energy collected on electrode
        self.total_energy = ev[ev['trigger'] == 1]['ADC_value'].sum()
        # No clusterization yet
        self.cluster_inds = None

    def __str__(self):
        outstr = ''
        for i, row in enumerate(self.data):
            if self.cluster_inds is not None:
                if i in self.cluster_inds: outstr += '\n'
            outstr += '    ' + str(row) + '\n'
        return outstr

    def clusterize_data(self):
        '''Given all the readouts on the side, try to break them down into
           strips that correlate with each other. Use strip adjacency as a 
           starting point.'''
        det_diffs = self.data['detector'][1:] - self.data['detector'][:-1]
        ind_mask = det_diffs != 1
        self.cluster_inds = (np.arange(len(self.data))+1)[ind_mask]
