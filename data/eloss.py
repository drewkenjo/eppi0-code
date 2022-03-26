#!/usr/bin/env python

import ROOT, sys

fname = sys.argv[1]
rdf = ROOT.ROOT.RDataFrame("h22", fname)
names0 = list(rdf.GetColumnNames())


#_________________________________

lowbandcode = """
double dc1th = atan2(sqrt(dcx1*dcx1 + dcy1*dcy1), dcz1)*TMath::RadToDeg();
bool lowband = dc1th < (-53.14680163254601 + 79.61307254040804*pow(pp-0.3, 0.05739232362022314));
"""

inblosscode = lowbandcode + """
double dploss = exp(-2.739 - 3.932*pp) + 0.002907;
if(!lowband) dploss = exp(-1.2 - 4.228*pp) + 0.007502;
"""
outblosscode = lowbandcode + """
double dploss = exp(-2.739 - 3.932*pp) + 0.002907;
if(!lowband) dploss = exp(-1.871 - 3.063*pp) + 0.007517;
"""

losscode = inblosscode if 'inb' in fname else outblosscode

#_________________________________


rdf = rdf.Define("vals", """
double th = atan2(sqrt(px*px+py*py), pz)*TMath::RadToDeg();
double pp = sqrt(px*px + py*py + pz*pz);
""" + losscode + """
double feloss = (pp + dploss)/pp;
vector<double> vals = {feloss*px, feloss*py, feloss*pz};
return vals;
""")
rdf = rdf.Define("px1","vals[0]").Define("py1","vals[1]").Define("pz1","vals[2]")

names1 = list(rdf.GetColumnNames())
names1.remove("px")
names1.remove("py")
names1.remove("pz")

rdf.Snapshot("h22","/tmp/eloss.tmp.root", names1)

rdf = ROOT.ROOT.RDataFrame("h22","/tmp/eloss.tmp.root")
rdf = rdf.Define("px","(float)px1").Define("py","(float)py1").Define("pz","(float)pz1")
rdf.Snapshot("h22", fname.replace(".root",".eloss.root"), names0)

