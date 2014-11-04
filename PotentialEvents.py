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
        # Set recarray to data
        self.data = ev
        # Compute total energy collected on electrode
        self.total_energy = ev[ev['trigger'] == 1]['ADC_value'].sum()

    def __str__(self):
        self.data.sort(order='detector')
        outstr = ''
        for row in self.data:
            outstr += '    ' + str(row) + '\n'
        return outstr
