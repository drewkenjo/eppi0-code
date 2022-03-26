def eps(q2,xb,E):
    y = q2/2/0.938/xb/E
    g2 = (2*xb*0.938)**2/q2
    eps = (1-y-0.25*y*y*g2)/(1-y+0.5*y*y+0.25*y*y*g2)
    return eps

def kinf(q2,xb,E):
    return math.sqrt(2*eps(q2,xb,E)*(1-eps(q2,xb,E)))


def get_alus(grbsa,shards):
    gralu = defaultdict(lambda: ROOT.TGraphErrors())
    for kk in grbsa:
        f1=ROOT.TF1("f1","[0]*sin(x*TMath::DegToRad())",0,360)
        f1.SetParameters(0.1,0)
        grbsa[kk].Fit(f1,"QR")
        iqx,itt = kk
        tmean = shards[(iqx,itt)][2].GetValue()
        gr1 = gralu[(iqx,)]
        gr1.SetPoint(gr1.GetN(), tmean, f1.GetParameter(0))
        gr1.SetPointError(gr1.GetN()-1, 0, f1.GetParError(0))
    return gralu


def get_sigs(grbsa, tshift=0):
    grsig = defaultdict(lambda: ROOT.TGraphErrors())
    for kk in grbsa:
        f1=ROOT.TF1("f1","[0]*sin(x*TMath::DegToRad())/(1+[1]*cos(2*x*TMath::DegToRad()))",0,360)
        f1.SetParLimits(1,-0.9,0)
        f1.SetParameters(0.1,-0.5)
        grbsa[kk].Fit(f1,"QR")
        iqx,itt = kk

        kin0,kin1 = grbsa[kk].kin
        qmean,xmean,tmean = kin1[0].GetValue(), kin1[1].GetValue(), kin1[2].GetValue()

        gr1 = grsig[(iqx,)]
        kf = kinf(qmean,xmean,10.604)
        gr1.SetPoint(gr1.GetN(), tmean+tshift, f1.GetParameter(0)/kf)
        gr1.SetPointError(gr1.GetN()-1, 0, f1.GetParError(0)/kf)

        qmean,xmean = kin0[0].GetValue(), kin0[1].GetValue()
        setattr(gr1, 'qmean', qmean)
        setattr(gr1, 'xmean', xmean)
    return grsig
