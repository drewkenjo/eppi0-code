#!/usr/bin/env python

import ROOT, sys

fname = sys.argv[1]
rdf = ROOT.ROOT.RDataFrame("h22", fname)
names0 = list(rdf.GetColumnNames())

ecorrfile = '../code/utsav.inb.ele3phis.code' if 'inb' in fname else '../code/utsav.outb.ele3phis.code'
with open(ecorrfile) as ff:
    rdf = rdf.Define("fe",ff.read()+"return fe;")

rdf = rdf.Define("ex1","fe*ex").Define("ey1","fe*ey").Define("ez1","fe*ez")

names1 = list(rdf.GetColumnNames())
names1.remove("ex")
names1.remove("ey")
names1.remove("ez")

rdf.Snapshot("h22","/tmp/ecorr.tmp.root", names1)

rdf = ROOT.ROOT.RDataFrame("h22","/tmp/ecorr.tmp.root")
rdf = rdf.Define("ex","(float)ex1").Define("ey","(float)ey1").Define("ez","(float)ez1")
rdf.Snapshot("h22", fname.replace(".root",".ecorr.root"), names0)
