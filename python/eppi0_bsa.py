import math, ROOT

def make_fits(qdf, option):
	fpar1 = ROOT.TF1("p1","[0]*sin(x*TMath::DegToRad())",0,360)
	fpar3 = ROOT.TF1("p3","[0]*sin(x*TMath::DegToRad())/(1 + [1]*cos(x*TMath::DegToRad()) + [2]*cos(2*x*TMath::DegToRad()))",0,360)
	fpar2 = ROOT.TF1("p2","[0]*sin(x*TMath::DegToRad())/(1 + [1]*cos(2*x*TMath::DegToRad()))",0,360)
	fpar3f = ROOT.TF1("p3f","[0]*sin(x*TMath::DegToRad())/(1 + [1]*cos(x*TMath::DegToRad()) + [2]*cos(2*x*TMath::DegToRad()))",0,360)

	if 'fit' in option:
		qdf.fit = type("fit", (), {})()

	for fpar in [fpar1, fpar2, fpar3, fpar3f]:
		grt = ROOT.TGraphErrors()
		for tdf in qdf.tdfs:
			gr = tdf.grbsa.Clone()
			if 'fit' in option:
				gr = tdf.grf.Clone()

			f0 = fpar.Clone()

			if f0.GetName()=="p3":
				f0.SetParLimits(1,0,0.3)
				f0.SetParLimits(2,-1,0)
			elif f0.GetName()=="p3f":
				f0.FixParameter(1,0)
				f0.FixParameter(2,0)

			gr.Fit(f0, 'qr')

			fitbn = type(fpar.GetName(), (), {})()
			fitbn.grbsa = gr
			fitbn.fbsa = f0
			setattr(tdf, fpar.GetName(), fitbn)

			grt.SetPoint(grt.GetN(), tdf.tt, f0.GetParameter(0))
			grt.SetPointError(grt.GetN()-1, 0, f0.GetParError(0))

		if 'fit' in option:
			setattr(qdf.fit, 'gr'+fpar.GetName(), grt)
		else:
			setattr(qdf, 'gr'+fpar.GetName(), grt)




def make_allbsas(rdf, option='num'):
	for binid in rdf.shards:
		shrd = rdf.shards[binid]
		shrd.xb = shrd.hqx.GetMean(1)
		shrd.q2 = shrd.hqx.GetMean(2)
		shrd.fi = shrd.httfi.GetMean(1)
		shrd.tt = shrd.httfi.GetMean(2)
		shrd.pb = 0.86 if 'inb' in shrd.fname else 0.89

	for qdf in rdf.qdfs:
		for tdf in qdf.tdfs:
			make_bsas(tdf)
		make_fits(qdf, option)




def make_bsas(tdf):
	grn,grf = ROOT.TGraphErrors(),ROOT.TGraphErrors()
	tdf.totbg = []
	tdf.grn = grn
	tdf.grbsa = grn
	tdf.grf = grf

	for fdf in tdf.fdfs:
		nsigs, fsigs = [], []
		if fdf.hneg.GetEntries()>100 and fdf.hpos.GetEntries()>100:
			for h1 in [fdf.hneg, fdf.hpos]:
				nevs,binw = h1.GetEntries(),h1.GetBinWidth(1)
				mu,sig = 0.135, 0.012
				
				f1 = ROOT.TF1("f1","[0]*exp(-0.5*((x-[1])/[2])**2)/(sqrt(2*pi)*[2]) + (1.53846*[3] - 0.538462*[4]) + (7.69231*[4]-7.69231*[3])*x", 0.07, 0.2)
				#f1 = ROOT.TF1("f1","gausn(0) + pol0(3)", 0.07, 0.2)
				f1.SetLineWidth(1)
				f1.SetParameters(nevs*binw,mu,sig,1,1)

				f1.SetParLimits(0,0.3*nevs*binw,nevs*binw)
				f1.SetParLimits(1,0.1,0.17)
				f1.SetParLimits(2,0.005,0.02)
				f1.SetParLimits(3,0,h1.GetMaximum())
				f1.SetParLimits(4,0,h1.GetMaximum())
				h1.Fit(f1,"LRQ")
				mu,sig = f1.GetParameter(1), abs(f1.GetParameter(2))

				mggs = [mu-5*sig, mu-3*sig, mu+3*sig, mu+5*sig]
				il5sig,il3sig,ir3sig,ir5sig = [h1.FindBin(mgg) for mgg in mggs]
				tot,bg = h1.Integral(il3sig,ir3sig), h1.Integral(il5sig,il3sig-1) + h1.Integral(ir3sig+1,ir5sig)

				nsig = tot-bg
				dnsig = math.sqrt(math.sqrt(tot)**2 + 1.5*math.sqrt(bg)**2)

				fsig = f1.GetParameter(0)/binw
				dfsig = f1.GetParError(0)/binw
			
				nsigs.append((nsig, dnsig, tot, bg))
				fsigs.append((fsig, dfsig))

			a,da,tot1,bg1,b,db,tot2,bg2 = nsigs[0] + nsigs[1]
			nbsa = (a-b)/(a+b)/fdf.pb
			dnbsa = 2*math.sqrt(b*b*da*da + a*a*db*db)/fdf.pb/(a+b)/(a+b)
			
			a,da,b,db = fsigs[0] + fsigs[1]
			fbsa = (a-b)/(a+b)/fdf.pb
			dfbsa = 2*math.sqrt(b*b*da*da + a*a*db*db)/fdf.pb/(a+b)/(a+b)

			if dnbsa>0:
				grn.SetPoint(grn.GetN(), fdf.fi, nbsa)
				grn.SetPointError(grn.GetN()-1, 0, dnbsa)
				tdf.totbg.append((tot1, bg1))
				tdf.totbg.append((tot2, bg2))

				grf.SetPoint(grf.GetN(), fdf.fi, fbsa)
				grf.SetPointError(grf.GetN()-1, 0, dfbsa)
