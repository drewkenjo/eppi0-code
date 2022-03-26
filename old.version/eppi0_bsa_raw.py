import math
from collections import defaultdict

def get_bsas(shards):
    grbsa = defaultdict(lambda: ROOT.TGraphErrors())
    for kk in [kk for kk in shards if len(kk)==3]:
        sigs = []
        fsigs = []
        for ipb,pb in pbs+[(-ipb,-pb) for ipb,pb in pbs]:
            h2 = shards[kk]
            ipbbin = h2.GetYaxis().FindBin(ipb)
            h1 = h2.ProjectionX('hh'+str(ipb),ipbbin,ipbbin)
            iggs = [h1.FindBin(lm) for lm in gglims]

            f1 = ROOT.TF1("f1","gausn(0)+pol1(3)",mggmu-5*mggsig, mggmu+5*mggsig)
            f1.SetParameters(1,0.135,0.01,1,1)
            f1.SetLineWidth(3)
            h1.Fit(f1,"LRQ")

            mu,sig = f1.GetParameter(1), abs(f1.GetParameter(2))
            lms = [mu-5*sig, mu-3*sig, mu+3*sig, mu+5*sig]
            iggs = [h1.FindBin(lm) for lm in lms]
            
            tot,bg = h1.Integral(iggs[1]+1,iggs[2]-1), h1.Integral(iggs[0],iggs[1]) + h1.Integral(iggs[2],iggs[3])
            sig = tot
            dsig = math.sqrt(math.sqrt(tot)**2)

            if h1.GetEntries()>100:
                sigs.append((sig,dsig,pb))

        #if not any(sg[0]>0 and sg[0]<30 for sg in sigs):
        top = sum(sig/pb for sig,dsig,pb in sigs)
        btm = sum(sig for sig,dsig,pb in sigs)
        ftop = sum(sig/pb for sig,dsig,pb in fsigs)
        fbtm = sum(sig for sig,dsig,pb in fsigs)
            
        if btm>0:
            bsa = top/btm
            dbsa = math.sqrt(sum((dsig*(btm/pb-top)/btm/btm)**2 for sig,dsig,pb in sigs))

            if dbsa>0:
                iqx,itt,ifi = kk
                gr1 = grbsa[(iqx,itt)]
                setattr(gr1, 'kin', (shards[kk[:1]], shards[kk[:2]]))


                gr1.SetPoint(gr1.GetN(), ifi*40+20, bsa)
                gr1.SetPointError(gr1.GetN()-1, 0, dbsa)
    return grbsa
