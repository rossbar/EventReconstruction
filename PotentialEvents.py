import numpy as np
from matplotlib.pyplot import *
from analysisHelperFunctions import *
from error_codes import *
from dtypes import create_refined_edata, refined_edata_type
from transient_analysis import *

# Framework imports
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

    def visualize_signals(self, rdata):
        '''If there is only one cluster per side for any detector containing
           clusters, plot the signals associated with that cluster.'''
        if (self.ge1.ac.cluster_inds is None) or (self.ge2.ac.cluster_inds is None):
            print "Data has not been clustered!"
            return
        cluster_lens = np.array([len(self.ge1.ac.clusters), len(self.ge1.dc.clusters), len(self.ge2.ac.clusters), len(self.ge2.dc.clusters)], dtype=int)
        for cl in cluster_lens:
            if cl > 1:
                print "Multiple clusters detected, no simple visualization"
                return
        fig, ax = subplots(2, 2)
        for a, c, t in zip(ax.ravel(), (self.ge1.ac.clusters, self.ge1.dc.clusters, self.ge2.ac.clusters, self.ge2.dc.clusters), ('GeI AC', 'GeI DC', 'GeII AC', 'GeII DC')):
            c = c[0]
            c.visualize_signals(rdata, ax=a, title=t)
        fig.canvas.draw()

    def condense(self, rdata):
        '''Call self.condense on all segmented readout clusters'''
        out = np.zeros(0, dtype=refined_edata_type)
        for i,cluster_list in enumerate((self.ge1.ac.clusters, self.ge1.dc.clusters, self.ge2.ac.clusters, self.ge2.dc.clusters)):
            if len(cluster_list) > 0:
                for cluster in cluster_list:
                    ret = cluster.condense_to_edata(rdata)
                    if ret is not None:
                        out = np.concatenate((out, cluster.edata))
                    else:
                        print 'Encountered cluster for which analysis not implemented'
                        return None
        return out

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
        self.clusters = []

    def __str__(self):
        outstr = ''
        if len(self.clusters) == 0:
            for i, row in enumerate(self.data):
                if self.cluster_inds is not None:
                    if i in self.cluster_inds: outstr += '\n'
                outstr += '    ' + str(row) + '\n'
        else:
            for i, cluster in enumerate(self.clusters):
                for row in cluster.data:
                    outstr += '    ' + str(row) + '\n'
                if i != len(self.clusters)-1: outstr += '\n'
        return outstr

    def clusterize_data(self):
        '''Given all the readouts on the side, try to break them down into
           strips that correlate with each other. Use strip adjacency as a 
           starting point.'''
        # If there is no data to cluster, ignore
        if len(self.data) == 0: return
        # Reset the clustering
        if self.cluster_inds is not None:
            self.cluster_inds = None
            self.clusters = []
        self.cluster_inds = np.array([0], dtype=int)
        # Clustering by trigger
        for i in range(len(self.data)-1):
            next_trig = self.data['trigger'][i+1]
            cur_trig = self.data['trigger'][i]
            if next_trig == 0 and cur_trig == 0 and (1 in self.data['trigger'][i+1:]):
                new_ind = np.array([i+1], dtype=int)
                self.cluster_inds = np.concatenate((self.cluster_inds, new_ind))
        if self.cluster_inds is None: self.cluster_inds = np.array([], dtype=int)
        # Apply cluster indices to data to make readout cluster objects
        if len(self.cluster_inds) == 1:
            self.clusters.append(ReadoutCluster(self.data))
        else:
            i = 0
            while i < len(self.cluster_inds)-1:
                cluster = ReadoutCluster(self.data[self.cluster_inds[i]:self.cluster_inds[i+1]])
                self.clusters.append(cluster)
                i += 1
            # Handle last cluster
            cluster = ReadoutCluster(self.data[self.cluster_inds[i]:])
            self.clusters.append(cluster)

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
        self.edata = None

    def __str__(self):
        return str(self.data)

    def visualize_signals(self, rdata, ax=None, title=None):
        '''Plot the signals corresponding to the readout cluster'''
        for rid in self.data['rid']: 
            if ax is None:
                plot(rdata[rid,:])
            else:
                ax.plot(rdata[rid,:])
                if title is not None:
                    ax.set_title(title)

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
            if np.isfinite(t50): edata_out['t50'] = t50
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
                # If signals are above noise, compute transient ratio
                if lta > lsig_noise and rta > rsig_noise:
                    R = (rta - lta)/(rta + lta)
                    # Compute uncertainty
                    SNR_l = lta/lsig_noise
                    SNR_r = rta/rsig_noise
                    R_sig = 2*((lta*rta)/(lta+rta)**2)*\
                            ((1/SNR_l)**2 + (1/SNR_r)**2)**(0.5)
                    edata_out['R'] = R
                    edata_out['sigma_R'] = R_sig
            self.edata = edata_out
            return 0
        else: return None
