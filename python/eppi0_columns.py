ecloose, ectight, prochi2, prodc, gloosest, gloose, gbeta = range(7)

def define_eppi0_columns_using_proton(df):
	vals = 'pe,the,fie,pp,thp,fip,pg1,thg1,fig1,pg2,thg2,fig2,mm2,mgg,yy,q2,ww,xb,tt,phistar,ftheta,misse,dpx,dpy,dpz,dphi,dtheta,mmp,dpt'

	rdf = df.Define("vals","""
	double E0 = 10.6041, Mpro = 0.938272;

	TLorentzVector beam(0,0,E0, E0), targ(0,0,0,Mpro);

	TLorentzVector ele, pro, g1, g2;
	ele.SetXYZM(ex,ey,ez,0);
	pro.SetXYZM(px,py,pz,Mpro);
	g1.SetXYZM(g1x,g1y,g1z,0);
	g2.SetXYZM(g2x,g2y,g2z,0);

	auto vqq = beam-ele;
	auto gg = g1+g2;
	auto epx = beam+targ-ele-pro;

	//auto nrm1 = vqq.Vect().Cross(pro.Vect());
	//auto nrm2 = vqq.Vect().Cross(gg.Vect());
	//double dfi = nrm2.Dot(pro.Vect()) > 0 ? 360 - nrm1.Angle(nrm2)*TMath::RadToDeg() : nrm1.Angle(nrm2)*TMath::RadToDeg();

	double dphi = (gg.Phi()-epx.Phi())*TMath::RadToDeg();
	double dtheta = (gg.Theta()-epx.Theta())*TMath::RadToDeg();

	double pe = ele.P();
	double the = ele.Theta()*TMath::RadToDeg();
	double fie = ele.Phi()*TMath::RadToDeg();
	fie += (fie<0) ? 360 : 0;

	double pp = pro.P();
	double thp = pro.Theta()*TMath::RadToDeg();
	double fip = pro.Phi()*TMath::RadToDeg();
	fip += (fip<0) ? 360 : 0;

	double pg1 = g1.P();
	double thg1 = g1.Theta()*TMath::RadToDeg();
	double fig1 = g1.Phi()*TMath::RadToDeg();
	fig1 += (fig1<0) ? 360 : 0;

	double pg2 = g2.P();
	double thg2 = g2.Theta()*TMath::RadToDeg();
	double fig2 = g2.Phi()*TMath::RadToDeg();
	fig2 += (fig2<0) ? 360 : 0;

	double ftheta = gg.Angle((beam+targ-ele-pro).Vect())*TMath::RadToDeg();

	auto epggx = beam+targ-ele-pro-gg;
	double misse = epggx.E();
	double dpx = epggx.Px();
	double dpy = epggx.Py();
	double dpz = epggx.Pz();
	double dpt = sqrt(dpx*dpx + dpy*dpy);

	double mgg = gg.M();
	double mm2 = (beam+targ-ele-pro).M2();
	double mmp = (beam+targ-ele-gg).M();
	double q2 = -vqq.M2();
	double ww = (beam+targ-ele).M();
	double xb = q2/(ww*ww - targ.M2() + q2);
	double tt = -(pro-targ).M2();

	double yy = (E0-pe)/E0;

	auto lnorm = vqq.Vect().Cross(beam.Vect());
	auto hnorm = pro.Vect().Cross(vqq.Vect());
	double phistar = lnorm.Dot(pro.Vect()) > 0 ? 360 - lnorm.Angle(hnorm)*TMath::RadToDeg() : lnorm.Angle(hnorm)*TMath::RadToDeg();

	return std::vector<double> {""" + vals+"""};
	""")

	vals = vals.split(',')
	for iv in range(len(vals)):
		rdf = rdf.Define(vals[iv], "vals[{}]".format(iv))

	if 'ex0' in list(rdf.GetColumnNames()):
		rdf = rdf.Define("phistar0", """
		double E0 = 10.6041, Mpro = 0.938272;
		TLorentzVector beam(0,0,E0, E0), targ(0,0,0,Mpro);
		TLorentzVector ele0, pro0;
		ele0.SetXYZM(ex0,ey0,ez0,0);
		pro0.SetXYZM(px0,py0,pz0,0);
		auto vqq = beam - ele0;
		auto lnorm = vqq.Vect().Cross(beam.Vect());
		auto hnorm = pro0.Vect().Cross(vqq.Vect());
		double phistar0 = lnorm.Dot(pro0.Vect()) > 0 ? 360 - lnorm.Angle(hnorm)*TMath::RadToDeg() : lnorm.Angle(hnorm)*TMath::RadToDeg();
		return phistar0;
		""")
	return rdf
