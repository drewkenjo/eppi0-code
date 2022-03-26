#!/usr/bin/env python

import ROOT, sys
ROOT.gROOT.SetBatch(True) 

from eppi0_bsa import make_allbsas

outfile = ROOT.TFile('systematics.graphs.root', 'recreate')

for fname in sys.argv[1:]:
  print(fname)

  rdf = type('rdf',(),{})()
  rdf.shards = {}
  rdf.qdfs = [type('qdf',(),{})() for i in range(5)]
  for qdf in rdf.qdfs:
    qdf.tdfs = [type('tdf',(),{})() for i in range(3)]
    for tdf in qdf.tdfs:
      tdf.fdfs = [type('fdf',(),{})() for i in range(9)]

  ff = ROOT.TFile(fname)

  for kk in ff.GetListOfKeys():
    tpl = eval(kk.GetName())
    obj = kk.ReadObj()
    if len(tpl)==2:
      iqx,name = tpl
      df = rdf.qdfs[iqx]
    elif len(tpl)==3:
      iqx,itt,name = tpl
      df = rdf.qdfs[iqx].tdfs[itt]
    elif len(tpl)==4:
      iqx,itt,ifi,name = tpl
      df = rdf.qdfs[iqx].tdfs[itt].fdfs[ifi]

    rdf.shards[tpl[:-1]] = df
    setattr(df, 'fname', fname)
    setattr(df, name, obj)

  make_allbsas(rdf)

  outfile.cd()
  for kk in rdf.shards:
    df = rdf.shards[kk]
    for name in ['grbsa','grp1']:
      if hasattr(df,name):
        getattr(df,name).Write(str((fname,name)+kk))

outfile.Close()
