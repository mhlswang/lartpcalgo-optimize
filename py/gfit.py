#!/usr/bin/env python

plot = True 
report = True
guess = False

import mmap, contextlib, shelve, os
from io import BytesIO
import numpy as np
from lmfit import Minimizer, Parameter, Parameters, report_fit
from lmfit.models import *   # https://lmfit.github.io/lmfit-py/builtin_models.html
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


class Hit:
    def __init__(self, hitdata):
        self.info = {}
        self.model = None 
        nl = hitdata.find(b'\n')
        header = hitdata[:nl]
        l0 = header.split()
        for i in l0[1:]: 
            k,v= i.split(b'=')
            self.info[k.decode('utf-8')]=v.decode('utf-8')
        self.name = 'w' + self.info['wire'] + '_v' + self.info['view']
        self.x, self.y = np.genfromtxt(BytesIO(hitdata[nl+1:]), delimiter=',', unpack=True)
        self.peaks_in_interval = []
        self.center = []
        self.amplitude = []
        self.width = []
        self.num_peaks = 0
        self.pars = Parameters()
        
    def fit_model(self, method='leastsq'):
        self.__init_model()
        for i in range(self.num_peaks):
            this_model = self.__make_local_model(i) 
            if not self.model: self.model = this_model
            else: self.model = self.model + this_model
        return self.model.fit(self.y, x=self.x, method=method)
        
                
    def __init_model(self):
        self.peaks_in_interval = np.array(argrelextrema(self.y, np.greater)[0])
        self.num_peaks = len(self.peaks_in_interval)
        self.center = self.x[self.peaks_in_interval]
        self.amplitude = self.y[self.peaks_in_interval]
        self.width = np.zeros(self.num_peaks) + 2.0

    def __make_local_model(self, num):
        global guess
        pref = "n{0}_".format(num)
        #model = GaussianModel(prefix = pref)  # nonlinear model for peaks
        model = VoigtModel(prefix = pref)  # nonlinear model for peaks
        if guess:
            self.pars  = model.guess(self.y, x=self.x)
        else:
            model.set_param_hint(pref+'amplitude', value=self.amplitude[num], min=1, max=self.amplitude[num]*7)
            model.set_param_hint(pref+'center', value=self.center[num], min=self.center[num]-1., max=self.center[num]+1.)
            model.set_param_hint(pref+'sigma', value=self.width[num], min=0.5, max=3)
        # Add linear model for background
        model += LinearModel(prefix = "l{0}_".format(num)) 
        return model


def ReadHitData(path):
    global hits
    prev_loc=0
    
    prefix = os.path.basename(path).split('.')[0]
    if os.path.exists(prefix+'.shelve.dat'):
        print("Reading hit data from binary shelve file: %s" %prefix+'.shelve.dat')
        d = shelve.open(prefix+'.shelve',flag='r')
        hits = d['hits']
        d.close()
        return hits
    
    if not os.path.exists(path):
        sys.exit("Error: invalid data path: %s" % path)
    hits = []
    
    with open(path, 'r') as f:
        print("Reading hit data from ASCII file: %s" % path)
        with contextlib.closing(mmap.mmap(f.fileno(), 0, prot=mmap.PROT_READ)) as m:
            while prev_loc >= 0:
                m.seek(prev_loc+4)
                hit_loc = m.find(b'hit view')
                hits.append(Hit(m[prev_loc:hit_loc]))
                prev_loc = hit_loc
    try:
        d = shelve.open(prefix+'.shelve',flag='n')
        d['hits'] = hits
        d.close()
    except:
        print("Could not cache data with shelve")
    return hits


#hits = ReadHitData('all_hits.txt')
hits = ReadHitData('some_hits.txt')

if plot:
    pdf = PdfPages('results.pdf')

# Process hits
for hit in hits:
    print (hit.info)
    result = hit.fit_model(method='leastsq')
    if report: print(result.fit_report(min_correl=0.25))
    if not plot: continue
    fig = plt.figure(figsize=(4, 3))
    result.plot_fit()
    plt.title(hit.name)
    #plt.show()
    pdf.savefig(fig)

if plot: pdf.close()



