#!/usr/bin/env python

import ROOT, sys
ROOT.gROOT.SetBatch(True) 
ROOT.gErrorIgnoreLevel = ROOT.kWarning


from eppi0_columns import define_eppi0_columns_using_proton
from eppi0_binning_scheme_v3 import define_eppi0_bins
from eppi0_binning_scheme_v3 import process_eppi0_bins
from eppi0_binning_scheme_v3 import save

ecloose, ectight, prochi2, prodc, gloosest, gloose, gbeta = range(7)

dpt,dphi = [float(vv) for vv in sys.argv[1:]]

for fname in ["data/lvl2_eppi0.inb.qa.ecorr.pcorr.root", "data/lvl2_eppi0.outb.qa.ecorr.pcorr.root"]:
    dfs = []
    df = ROOT.RDataFrame("h22", fname)
    df = define_eppi0_columns_using_proton(df)

    cut = "&&".join(f"(((int)status)&(1<<{i}))" for i in [prodc,gloosest])
    cut += f"&& abs(dpt)<{dpt} && abs(dphi)<{dphi}"
    cut += "&& thp<44.106+-6.625*pp+1.438*pp*pp"
    cut += "&& tt>tmin"

    df = df.Filter(cut)

    for dpz in [0.6, 0.7, 0.8]:
        for dmm2 in [0.3, 0.35, 0.4]:
            rdf = df.Filter(f"abs(dpz-0.2)<{dpz} && abs(mm2-0.05)<{dmm2}")
            rdf = define_eppi0_bins(rdf, fname)
            rdf = process_eppi0_bins(rdf, fname)
            dfs.append((str(dpz),str(dmm2),rdf))


    for dpz,dmm2,rdf in dfs:
        save(rdf, fname.replace(".root", f".{dpt}_{dphi}_{dpz}_{dmm2}.root").replace("data/","data/syst/"))


