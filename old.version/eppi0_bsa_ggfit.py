import math
from collections import defaultdict

ROOT.gStyle.SetOptFit(1)

def get_bsas(shards):
    hggms={}
    grbsa = defaultdict(lambda: ROOT.TGraphErrors())
    for kk in [kk for kk in shards if len(kk)==3]:
        sigs = []
        for ipb,pb in pbs+[(-ipb,-pb) for ipb,pb in pbs]:
            h2 = shards[kk]
            ipbbin = h2.GetYaxis().FindBin(ipb)
            h1 = h2.ProjectionX('hh'+str(ipb),ipbbin,ipbbin)
            f1 = ROOT.TF1("f1","gausn(0)+pol1(3)",mggmu-5*mggsig, mggmu+5*mggsig)
            f1.SetParameters(1,0.135,0.01,1,-1)
            h1.Fit(f1,"LRQ")
            sig = f1.GetParameter(0)/h1.GetBinWidth(1)
            dsig = f1.GetParError(0)/h1.GetBinWidth(1)
            if h1.GetEntries()>100:
                sigs.append((sig,dsig,pb))

        top = sum(sig/pb for sig,dsig,pb in sigs)
        btm = sum(sig for sig,dsig,pb in sigs)
            
        if btm>0:
            bsa = top/btm
            dbsa = math.sqrt(sum((dsig*(btm/pb-top)/btm/btm)**2 for sig,dsig,pb in sigs))

            iqx,itt,ifi = kk
            gr1 = grbsa[(iqx,itt)]
            setattr(gr1, 'kin', (shards[kk[:1]], shards[kk[:2]]))

            gr1.SetPoint(gr1.GetN(), ifi*40+20, bsa)
            gr1.SetPointError(gr1.GetN()-1, 0, dbsa)
    return grbsa
