def define_eppi0_bins(rdf, fname):
	q2s = [
		   "xb<x1 && q2>2 && q2<(y0 + (xb-x0)/(x1-x0)*(y1-y0)) && xb<x01",
		   "xb<x1 && q2>2 && q2<(y0 + (xb-x0)/(x1-x0)*(y1-y0)) && xb>x01",
		   "xb<x1 && q2>2",
		   "xb>x1 && q2>2 && q2< (y4 + (xb-x4)/(x5-x4)*(y5-y4))",
		   "xb>x1 && q2>2",
		  ]

	rdf = rdf.Define("iqx","""
		double x0=0.15, y0=2.63;
		double x1=0.45, y1=4.5;
		double x01 = 0.34;

		double x4=0.45, y4=4.8;
		double x5=0.666, y5=6.2;

		std::vector<bool> q2s = {"""+",".join(q2s)+"""};
		for(int ii=0;ii<q2s.size();ii++)
			if(q2s[ii])
				return ii;
		return -1;
	""")

	if "inb" in fname:
		ttgrid = [
			[0.42,0.9,2],
			[0.42,0.9,2],
			[0.46,0.9,2],
			[0.62,0.98,2],
			[0.82,1.26,2]
		]
	else:
		ttgrid = [
			[0.62,1.09,2],
			[0.60,0.94,2],
			[0.62,1.02,2],
			[0.70,1.08,2],
			[0.94,1.36,2]
		]

	tarrstr = "},{".join(", ".join(map(str,tts)) for tts in ttgrid)

	rdf = rdf.Define("itt","""
		if(iqx<0) return -1;
		std::vector<std::vector<float>> tts = {{"""+tarrstr+"""}};
		for(int ii=0; ii<tts[iqx].size(); ii++)
			if( tt < tts[iqx][ii] )
				return ii;
		return -1;
	""")


	rdf = rdf.Define("ifi","(int) (phistar/40)")
	rdf.q2s = q2s
	rdf.ttgrid = ttgrid
	return rdf




def process_eppi0_bins(rdf, fname):
	rdf.qdfs = []
	rdf.shards = {}

	for iqx in range(5):
		qdf = rdf.Filter(f"iqx=={iqx}")
		rdf.qdfs.append(qdf)
		qdf.binid = (iqx,)
		qdf.tdfs = []
		rdf.shards[qdf.binid] = qdf
		
		for itt in range(3):
			tdf = qdf.Filter(f"itt=={itt}")
			qdf.tdfs.append(tdf)
			tdf.binid = (iqx,itt)
			tdf.fdfs = []
			rdf.shards[tdf.binid] = tdf
			
			for ifi in range(9):
				fdf = tdf.Filter(f"ifi=={ifi}")
				tdf.fdfs.append(fdf)
				fdf.binid = (iqx,itt,ifi)
				rdf.shards[fdf.binid] = fdf

				hneg = fdf.Filter("ihel<0").Histo1D(("hmgg","M_{#gamma#gamma};M_{#gamma#gamma} [GeV]", 195, 0.07, 0.2), "mgg")
				hpos = fdf.Filter("ihel>0").Histo1D(("hmgg","M_{#gamma#gamma};M_{#gamma#gamma} [GeV]", 195, 0.07, 0.2), "mgg")
				fdf.hneg = hneg
				fdf.hpos = hpos

	for binid in rdf.shards:
		shrd = rdf.shards[binid]
		shrd.fname = fname
		shrd.hqx = shrd.Histo2D(("hqx", "Q^{2} vs x_{B};x_{B};Q^{2} [GeV^{2}]",200,0.06,0.84,200,1.8,12), "xb", "q2")
		shrd.httfi = shrd.Histo2D(("httfi", "tt vs #phi;#phi [#circ];-tt [GeV^{2}]",180,0,360,200,0,2), "phistar", "tt")
		shrd.hmgg = shrd.Histo1D(("hmgg", "M_{#gamma#gamma};M_{#gamma#gamma} [GeV]",180,0.07,0.2), "mgg")

	return rdf





def save(rdf, outname):
	import ROOT
	fout = ROOT.TFile(outname, "recreate")
	for binid in rdf.shards:
		shrd = rdf.shards[binid]
		shrd.hqx.Write(str(binid+('hqx',)))
		shrd.httfi.Write(str(binid+('httfi',)))
		shrd.hmgg.Write(str(binid+('hmgg',)))
		if hasattr(shrd, 'hneg'):
			shrd.hneg.Write(str(binid+('hneg',)))
			shrd.hpos.Write(str(binid+('hpos',)))
	fout.Close()

