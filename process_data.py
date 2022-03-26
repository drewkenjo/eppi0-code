#!/usr/bin/env python

import ROOT, sys
ROOT.gROOT.SetBatch(True) 
ROOT.gErrorIgnoreLevel = ROOT.kWarning


from eppi0_columns import define_eppi0_columns_using_proton
from eppi0_binning_scheme_v3 import define_eppi0_bins
from eppi0_binning_scheme_v3 import save


for fname in ["data/lvl2_eppi0.inb.qa.ecorr.pcorr.root", "data/lvl2_eppi0.outb.qa.ecorr.pcorr.root"]:
    ist = 0
    df = ROOT.RDataFrame("h22", fname)
    df = define_eppi0_columns_using_proton(df)
    df = df.Filter("abs(dpx)<0.3 && abs(dpy)<0.3 && abs(dphi)<4 && dpz>-0.5 && dpz<0.9 && mm2>-0.3 && mm2<0.4")
    for status in ["true", " && ".join(f"(((int)status)&(1<<{i}))" for i in [6])]:
        rdf = df.Filter(status)
        rdf = define_eppi0_bins(rdf, fname)
        save(rdf, fname+f".{ist}.root")
        ist+=1


