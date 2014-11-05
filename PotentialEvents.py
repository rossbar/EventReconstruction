import numpy as np
from analysisHelperFunctions import *
from error_codes import *
from dtypes import create_refined_edata

# Framewor imports
import sys
import os
sys.path.append(os.path.expanduser('~')+r'/PythonFramework')
from src.CAnalysis.canalysis import search_t50

class PotentialEvent(object):
    total_energy = 0
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
        ge1match = self.ge1.passes_energy_match()
        ge2match = self.ge2.passes_energy_match()
        if ge1match and ge2match:
            self.total_energy = self.ge1.total_energy + self.ge2.total_energy
            return True
        else:
            return False

    def clusterize(self):
        '''Call side clusterization.'''
        self.ge1.clusterize()
        self.ge2.clusterize()

    def trim_unused_transients(self):
        '''Call side transient removal'''
        self.ge1.trim_unused_transients()
        self.ge2.trim_unused_transients()

class Detector(object):
    errors = []
    total_energy = 0
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
            self.total_energy = np.max((self.ac.total_energy, self.dc.total_energy))
            return True
        else: 
            self.errors.append(ENERGY_MATCH_FAILURE)
            return False

    def clusterize(self):
        '''Call clusterization for each side'''
        self.ac.clusterize_data()
        self.dc.clusterize_data()

    def trim_unused_transients(self):
        '''Call transient trimming for each side'''
        self.ac.remove_redundant_transients()
        self.dc.remove_redundant_transients()

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
        # Reset the clustering
        if self.cluster_inds is not None: self.cluster_inds = None
        # Clustering by trigger
        for i in range(len(self.data)-1):
            next_trig = self.data['trigger'][i+1]
            cur_trig = self.data['trigger'][i]
            if next_trig == 0 and cur_trig == 0 and (1 in self.data['trigger'][i+1:]):
                new_ind = np.array([i+1], dtype=int)
                if self.cluster_inds is None:
                    self.cluster_inds = new_ind
                else:
                    self.cluster_inds = np.concatenate((self.cluster_inds, np.array([i+1], dtype=int)))

    def remove_redundant_transients(self):
        '''Large charge deposits may cause non-adjacent neighbors to fire
           (n+2, etc.). These signals are not currently used for anything.
           This function removes them from the data.'''
        # Only attempt if there is data to operate on
        if len(self.data) > 1:
            self.keep_mask = np.ones(len(self.data), dtype=bool)
            # If there are three adjacent transient signals, remove the center one
            for i in np.arange(1, len(self.data)-1):
                prev_trig = self.data['trigger'][i-1]
                cur_trig = self.data['trigger'][i]
                next_trig = self.data['trigger'][i+1]
                # trig 0 readout between two other trig 0 readouts is redundant
                if prev_trig == 0 and cur_trig == 0 and next_trig == 0:
                    self.keep_mask[i] = False
                # Another neighbor trig but subsequent readout isn't adjacent
                elif prev_trig == 0 and cur_trig == 0 and next_trig == 1:
                    diff = self.data['detector'][i+1] - self.data['detector'][i]
                    if diff > 1:
                        self.keep_mask[i] = False
            # Handle edge cases
            if self.data['trigger'][0] == 0 and self.data['trigger'][1] == 0:
                self.keep_mask[0] = False
            if self.data['trigger'][-2] == 0 and self.data['trigger'][-1] == 0:
                self.keep_mask[-1] = False
            # Remove the unused transients
            self.data = self.data[self.keep_mask]

class ReadoutCluster(object):
    def __init__(self, ev):
        # Sort by detector
        ev.sort(order='detector')
        # Store data
        self.data = ev
        self.num_trigs = (ev['trigger'] == 1).sum()
        self.num_strips = len(ev)
        self.trigger_pattern = ev['trigger']

    def condense_to_edata(self, rdata):
        '''Take the data associated with the readout cluster, and condense it
           into a single edata-type output which includes stuff from the
           raw signals, like transient ratios, t50, etc.'''
        # Charge collection on only one strip
        if self.num_trigs == 1:
            # Setup output
            cc_strip = self.data[self.data['trigger'] == 1]
            edata_out = create_refined_edata(cc_strip)
            # Fill in t50
            signal = rdata[cc_strip['rid']]
            t50 = search_t50(signal)
            # If it has a transient on each side, do transient analysis
            if np.array_equal(self.trigger_pattern, np.array([0,1,0])):
                # Get signals
                lsig = rdata[self.data[0]['rid'],:]
                rsig = rdata[self.data[2]['rid'],:]
                # Get transient amplitudes
                lta = get_transient_amplitude_simple(lsig)
                rta = get_transient_amplitude_simple(rsig)
                # Get signal noise
                lsig_noise = get_transient_noise_simple(lsig)
                rsig_noise = get_transient_noise_simple(rsig)
