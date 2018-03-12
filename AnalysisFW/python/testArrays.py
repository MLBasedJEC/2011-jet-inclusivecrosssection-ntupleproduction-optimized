#!/bin/env python

import numpy as np
import pandas as pd
pd.options.display.max_rows = 60
for jc in ['ak5PFJets','ak7PFJets']:
    print "Loading the",jc,"array\n"
    arr = np.load("/uscms_data/d2/aperloff/YOURWORKINGAREA/MLJEC/CMSSW_5_3_32/src/2011-jet-inclusivecrosssection-ntupleproduction-optimized/AnalysisFW/python/%sParams0.npy" % jc)
    df = pd.DataFrame(arr,columns=['run', 'lumi', 'event', 'met', 'sumet', 'rho', 'pthat', 'mcweight', 'njet', 'jet_pt', 'jet_eta', 'jet_phi', 'jet_E', 'jet_area', 'jet_jes', 'chf', 'nhf', 'phf', 'elf', 'muf', 'hf_hf', 'hf_phf', 'hf_hm', 'hf_phm', 'chm', 'nhm', 'phm', 'elm', 'mum', 'beta', 'bstar','cand_ijet'])
    df = df.drop_duplicates()
    #df[['event','met','sumet']].head(60)
    #df[['event','njet','mcweight']].head(60)
    print df[['event','njet','mcweight','met','sumet']].head(60)
    print '\n\n'
