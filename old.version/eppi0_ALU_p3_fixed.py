def eps(q2,xb,E):
    y = q2/2/0.938/xb/E
    g2 = (2*xb*0.938)**2/q2
    eps = (1-y-0.25*y*y*g2)/(1-y+0.5*y*y+0.25*y*y*g2)
    return eps

def kinf(q2,xb,E):
    return math.sqrt(2*eps(q2,xb,E)*(1-eps(q2,xb,E)))



atts = [

(2.59, 0.29, 0.29, 0.08747182840987841, -0.2646635971329326),
(2.57, 0.29, 0.63, 0.05412491365991858, -0.27884956042056336),
(2.62, 0.28, 1.34, 0.013126068646483545, -0.09189728063288845),
(2.91, 0.39, 0.31, 0.08943073990845882, -0.2046850818931999),
(2.98, 0.4, 0.62, 0.05381806239684255, -0.348682768276539),
(3.18, 0.39, 1.28, 0.012260518885145961, -0.263346986336739),
(4.62, 0.36, 0.33, 0.04017884539501501, -0.08186086258265353),
(4.76, 0.38, 0.65, 0.025571255122763415, -0.12450887535213255),
(4.67, 0.37, 1.36, 0.005433899301106981, -0.08270846687572068),
(3.79, 0.49, 0.47, 0.05380278539788884, -0.15292811478736587),
(4.15, 0.52, 0.78, 0.028624326998027794, -0.2966218569907229),
(4.36, 0.53, 1.37, 0.005919843512053976, -0.4003962631127233),
(6.11, 0.51, 0.62, 0.021208075349101192, -0.08458313101545631),
(6.58, 0.56, 1.03, 0.008569837770477989, -0.16546293018702068),
(6.88, 0.59, 1.6, 0.0017521206122438396, -0.26271947636468534),
(2.55, 0.29, 0.46, 0.07117000133681459, -0.3009976233373564),
(2.54, 0.29, 0.83, 0.03905995772726458, -0.2398371118952545),
(2.56, 0.28, 1.48, 0.009921500475184489, -0.07326091267030761),
(2.9, 0.39, 0.45, 0.07580775061480176, -0.3087395982559496),
(2.94, 0.39, 0.76, 0.04062951659446603, -0.3543553325072214),
(3.1, 0.39, 1.33, 0.01129260180158205, -0.26675417392130074),
(4.54, 0.36, 0.46, 0.036362595721446976, -0.11352065900753575),
(4.64, 0.37, 0.8, 0.01926374142377073, -0.1235075722824796),
(4.57, 0.36, 1.44, 0.004556897069497473, -0.06964839201659874),
(3.82, 0.5, 0.54, 0.05012711794042936, -0.20242134818278218),
(4.17, 0.52, 0.88, 0.022371957571740573, -0.3289950153835684),
(4.29, 0.52, 1.46, 0.00482279254829135, -0.39516206253624436),
(6.03, 0.52, 0.72, 0.018899966325202668, -0.11657955108904135),
(6.46, 0.56, 1.14, 0.006733263519699129, -0.19460386378685166),
(6.63, 0.58, 1.64, 0.0016673399167083654, -0.2724815394730434),

]


def get_sigs(grbsa, tshift=0):
    grsig = defaultdict(lambda: ROOT.TGraphErrors())
    for kk in grbsa:
        iqx,itt = kk

        kin0,kin1 = grbsa[kk].kin
        qmean,xmean,tmean = kin1[0].GetValue(), kin1[1].GetValue(), kin1[2].GetValue()
        
        #s0,att,alt = [(st+eps*sl,eps*stt,math.sqrt(2*eps*(1+eps))*slt) for q2,xb,tt,eps,st,sl,stt,slt in atts if abs(q2-qmean)<0.05 and abs(xb-xmean)<0.05 and abs(tt-tmean)<0.05][0]
        alt,att = [(alt,att) for q2,xb,tt,alt,att in atts if abs(q2-qmean)<0.05 and abs(xb-xmean)<0.05 and abs(tt-tmean)<0.05][0]

        f1=ROOT.TF1("f1","[0]*sin(x*TMath::DegToRad())/(1+[1]*cos(x*TMath::DegToRad())+[2]*cos(2*x*TMath::DegToRad()))",0,360)
        f1.SetParameters(0.1,alt,att)
        f1.FixParameter(1,alt)
        f1.FixParameter(2,att)
        grbsa[kk].Fit(f1,"QR")

        gr1 = grsig[(iqx,)]
        kf = kinf(qmean,xmean,10.604)
        gr1.SetPoint(gr1.GetN(), tmean+tshift, f1.GetParameter(0)/kf)
        gr1.SetPointError(gr1.GetN()-1, 0, f1.GetParError(0)/kf)

        qmean,xmean = kin0[0].GetValue(), kin0[1].GetValue()
        setattr(gr1, 'qmean', qmean)
        setattr(gr1, 'xmean', xmean)
    return grsig
