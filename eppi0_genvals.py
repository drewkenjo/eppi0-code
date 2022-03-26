def define_true_eppi0_columns(df):
	vals = 'q20,xb0,mt0'

	rdf = df.Define("genvals","""
	double E0 = 10.6041, Mpro = 0.938272;

	TLorentzVector beam(0,0,E0, E0), targ(0,0,0,Mpro);

	TLorentzVector ele0, pro0;
	ele0.SetXYZM(ex0,ey0,ez0,0);
	pro0.SetXYZM(px0,py0,pz0,Mpro);

	auto vqq0 = beam-ele0;

	double q20 = -vqq0.M2();
	double ww0 = (beam+targ-ele0).M();
	double xb0 = q20/(ww0*ww0 - targ.M2() + q20);
	double mt0 = -(pro0-targ).M2();

	return std::vector<double> {""" + vals+"""};
	""")

	vals = vals.split(',')
	for iv in range(len(vals)):
		rdf = rdf.Define(vals[iv], "genvals[{}]".format(iv))

	q2s = [
		"xb0<x1 && q20>2 && q20<(y0 + (xb0-x0)/(x1-x0)*(y1-y0)) && xb0<x01",
		"xb0<x1 && q20>2 && q20<(y0 + (xb0-x0)/(x1-x0)*(y1-y0)) && xb0>x01",
		"xb0<x1 && q20>2",
		"xb0>x1 && q20>2 && q20<(y4 + (xb0-x4)/(x5-x4)*(y5-y4))",
		"xb0>x1 && q20>2",
		]

	rdf = rdf.Define("iqx0","""
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
		mtgrid = [
			[0.42,0.9,2],
			[0.42,0.9,2],
			[0.46,0.9,2],
			[0.62,0.98,2],
			[0.82,1.26,2]
		]
	else:
		mtgrid = [
			[0.62,1.09,2],
			[0.60,0.94,2],
			[0.62,1.02,2],
			[0.70,1.08,2],
			[0.94,1.36,2]
		]

	tarrstr = "},{".join(", ".join(map(str,mts)) for mts in mtgrid)

	rdf = rdf.Define("imt0","""
		if(iqx0<0) return -1;
		std::vector<std::vector<float>> mts = {{"""+tarrstr+"""}};
		for(int ii=0; ii<mts[iqx0].size(); ii++)
			if( mt0 < mts[iqx0][ii] )
				return ii;
		return -1;
	""")

	return rdf


