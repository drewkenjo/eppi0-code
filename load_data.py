#!/usr/bin/env python

import ROOT, os, sys
ROOT.gROOT.SetBatch(True) 

def load_rdf(fname):
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

    df.binid = tpl[:-1]
    rdf.shards[tpl[:-1]] = df
    setattr(df, 'fname', fname)
    setattr(df, name, obj)

  return rdf
