#!/usr/bin/env python

import sys, ROOT, math
from collections import defaultdict


def eps(q2,xb,E):
  y = q2/2/0.938/xb/E
  g2 = (2*xb*0.938)**2/q2
  eps = (1-y-0.25*y*y*g2)/(1-y+0.5*y*y+0.25*y*y*g2)
  return eps

def kinf(q2,xb,E):                                                                                                         return math.sqrt(2*eps(q2,xb,E)*(1-eps(q2,xb,E)))


rootnames = [fn for fn in sys.argv if fn.endswith('.root')]

systmap = {}
for fname in [fn for fn in sys.argv if fn.endswith('.syst')]:
  with open(fname) as ff:
    for line in ff:
      vvs = line.split()
      kk = tuple(int(vv) for vv in vvs[1:-1])
      systmap[kk] = float(vvs[-1])



datasets = []

for fname in [fn for fn in sys.argv if fn.endswith('.dat')]:
  name,ext = fname.split("_")
  igr=len(datasets)
  with open(fname) as ff:
    datasets.append((int(igr),name,ff.read()))


grmap = defaultdict(lambda: ROOT.TGraphErrors())
for igr,name,dataset in datasets:
  for line in dataset.split('\n'):
    if line:
      xb,q2,tt,bsa,dbsa,kf = [float(v) for v in line.split()[:6]]
      gr = grmap[(name,igr,q2,xb)]
      if kf!=1:
        kf = kinf(q2,xb,5.776)
      gr.SetPoint(gr.GetN(), tt, bsa/kf)
      gr.SetPointError(gr.GetN()-1, 0, dbsa/kf)




phis = [0,40,80,120,160,200,240,280,320,360]


###########################################
###########################################
q2s = [2, 3, 4, 5, 6, 11]
tts = [
[0.0, 0.42, 1.0, 2.0],
[0.0, 0.5, 0.9, 2.0],
[0.0, 0.6, 1.02, 2.0],
[0.0, 1, 1.3, 2.0],
[0.0, 0.92, 1.38, 2.0],
]


###########################################
###########################################
if 'outb' in rootnames[0]:
  q2s = [2, 3, 4, 5, 6, 11]
  tts = [
[0.0, 0.62, 1.02, 2.0],
[0.0, 0.64, 1.08, 2.0],
[0.0, 0.72, 1.16, 2.0],
[0.0, 0.84, 1.28, 2.0],
[0.0, 1.04, 1.46, 2.0],
  ]



###########################################
###########################################
if '2019' in rootnames[0] or '2018' in rootnames[0]:
  q2s = [2,11]
  tts = [[0.0, 2.0]]


###########################################
###########################################
if len(rootnames)>1:
  q2s = [2, 3, 4, 5, 6, 11]
  tts = [
[0.0, 0.5, 0.88, 2.0],
[0.0, 0.56, 0.98, 2.0],
[0.0, 0.64, 1.08, 2.0],
[0.0, 0.74, 1.18, 2.0],
[0.0, 0.96, 1.4, 2.0],
  ]


statusarg=[arg for arg in sys.argv if arg.startswith('status=')]

ROOT.gStyle.SetLabelFont(42,'xy')
ROOT.gStyle.SetLabelSize(0.05,'xy')
ROOT.gStyle.SetTitleSize(0.07,'xy')
ROOT.gStyle.SetTitleOffset(0.75,'xy')



h22 = ROOT.TChain("h22")
for fname in rootnames:
  h22.Add(fname)
rdf = ROOT.RDataFrame(h22)
if h22.GetBranchStatus("status") and statusarg:
  stats = statusarg[0].split("=")[1]
  statline = " && ".join("(st&{:d})".format(2**int(stid)) for stid in stats)
  rdf = rdf.Filter("""
  int st = (int) status;
  return """+statline+";")

vals = ['mm2','mgg','q2','ww','xb','tt','phistar','iqx','itt','iphi']

rdf = rdf.Define("vals","""
TLorentzVector beam(0,0,10.6041,10.6041), targ(0,0,0,0.938);

TLorentzVector ele, pro, g1, g2;
ele.SetXYZM(ex,ey,ez,0);
pro.SetXYZM(px,py,pz,0.938);
g1.SetXYZM(g1x,g1y,g1z,0);
g2.SetXYZM(g2x,g2y,g2z,0);

auto vqq = beam-ele;
auto gg = g1+g2;

double mgg = gg.M();
double mm2 = (beam+targ-ele-pro).M2();
double q2 = -vqq.M2();
double ww = (beam+targ-ele).M();
double xb = q2/(ww*ww - targ.M2() + q2);
double tt = -(pro-targ).M2();

auto lnorm = vqq.Vect().Cross(beam.Vect());
auto hnorm = pro.Vect().Cross(vqq.Vect());
double phistar = lnorm.Dot(pro.Vect()) > 0 ? 360 - lnorm.Angle(hnorm)*TMath::RadToDeg() : lnorm.Angle(hnorm)*TMath::RadToDeg();

std::vector<double> q2s = {"""+",".join([str(q2) for q2 in q2s])+"""};
double iqx = -1;
for(int ii=0;ii<q2s.size();ii++)
  if(q2s[ii]<q2) iqx=ii;
  else break;

std::vector<std::vector<double>> tts = {"""+
",".join([("{"+",".join([str(tt) for tt in tlist])+"}") for tlist in tts])
+"""};

double itt = -1;
if(iqx>=0 && iqx<tts.size()){
  for(int ii=0;ii<tts[iqx].size();ii++)
    if(tts[iqx][ii]<tt) itt=ii;
    else break;
}

std::vector<double> phis = {"""+",".join(str(fi) for fi in phis)+"""};
double iphi = -1;
for(int ii=0;ii<phis.size();ii++)
  if(phis[ii]<phistar) iphi=ii;
  else break;


std::vector<double> vvs = {""" + ",".join(vals)+"""};
return vvs;
""")

for vv in vals:
  rdf = rdf.Define(vv, "vals[{}]".format(vals.index(vv)))

rdf = rdf.Define("ipb","-ihel*(run<5420 ? 1 : 2)")

hevs,htphis,hggs,hmggs,htts = [],[],[],[],[]
for iqx in range(len(q2s)-1):
  qdf = rdf.Filter('iqx=='+str(iqx))
  htt = qdf.Histo1D(('htt','',100,0,2), "tt")
  htts.append(htt)
  htphis.append(qdf.Histo2D(('htphi',';#phi;-t',100,0,360,100,0,2), "phistar", "tt"))
  hevs.append([])
  for itt in range(len(tts[iqx])-1):
    hevs[-1].append([])
    tdf = qdf.Filter("itt=={}".format(itt))
    hggs.append(tdf.Histo1D(("hmgg",";M_{#gamma#gamma} [GeV]",100,0.07,0.2), "mgg"))
    for iphi in range(len(phis)-1):
       bdf = tdf.Filter("iphi=={}".format(iphi))
       #hmggs.append(bdf.Histo1D(("hmgg","{} {} {}".format(iqx,itt,iphi),100,0.07,0.2), "mgg"))
       hmggs.append(bdf.Histo1D(("hmgg",";M_{#gamma#gamma} [GeV]",100,0.07,0.2), "mgg"))
       hevs[-1][-1].append((iqx,itt,iphi, qdf.Mean('q2'), qdf.Mean('xb'), tdf.Mean('tt'), tdf.Mean('xb'), bdf.Histo2D(("hgg","",100,0.07,0.2,5,-2,3),"mgg","ipb")))


hqx = rdf.Histo2D(('hqx',';x_{B};Q^{2} [GeV^{2}]',300,0,1,300,0,11), "xb", "q2")
#htx = rdf.Filter('iqx==4').Histo2D(('htx',';x_{B};-t [GeV^{2}]',300,0.3,0.8,300,0,2), "xb", "tt")


ROOT.gStyle.SetOptStat(0)
pdfname = rootnames[0]+'.bsa.pdf'
c1 = ROOT.TCanvas("c1","c1",1100,800)
c1.SetGrid(1)
c1.SetMargin(0.12,0.0,0.12,0.0)

c2 = ROOT.TCanvas("c2","c2",1100,470)
c2.SetGrid(1)
c2.SetMargin(0.12,0.0,0.12,0.0)

c1.Print(pdfname+'[')
c2.Print(pdfname.replace('.pdf','.phi.pdf['))


ll = ROOT.TLine()
lat = ROOT.TLatex();
lat.SetTextSize(0.06)

ROOT.gStyle.SetOptFit(1)
lms=[0.1325-5*0.0125, 0.1325-3*0.0125, 0.1325+3*0.0125, 0.1325+5*0.0125]
for hqele in hevs:
  grt = ROOT.TGraphErrors()
  grtsys = ROOT.TGraphErrors()
  xbis = []
  for htele in hqele:
    grbsa = ROOT.TGraphErrors()
    grsys = ROOT.TGraphErrors()
    for iqx,itt,iphi,q2m,xbm,ttm,xbi,h2 in htele:
      q2m = q2m.GetValue()
      xbm = xbm.GetValue()
      ttm = ttm.GetValue()
      xbis.append(xbi.GetValue())

      sigs = []
      goodbin = True
 
      yax = h2.GetYaxis()
      pbs = [(1,0.86), (2,0.89)]
      for ipb,pb in pbs+[(-ipb,-pb) for ipb,pb in pbs]:
        ipbbin = yax.FindBin(ipb)
        h1 = h2.ProjectionX('hh'+str(ipb),ipbbin,ipbbin)
        bns = [h1.FindBin(lm) for lm in lms]
        tot,bg = h1.Integral(bns[1],bns[2]), h1.Integral(bns[0],bns[1]) + h1.Integral(bns[2],bns[3])
        sig = tot-bg
        dsig = math.sqrt(math.sqrt(tot)**2 + 1.5*math.sqrt(bg)**2)
        if tot>30:
          sigs.append((sig,dsig,pb))
        elif tot>0:
          goodbin = False

      if goodbin:
        top = sum(sig/pb for sig,dsig,pb in sigs)
        btm = sum(sig for sig,dsig,pb in sigs)

        if btm>0:
          bsa = top/btm
          dbsa = math.sqrt(sum((dsig*(btm/pb-top)/btm/btm)**2 for sig,dsig,pb in sigs))

          grbsa.SetPoint(grbsa.GetN(), iphi*40+20, bsa)
          grbsa.SetPointError(grbsa.GetN()-1, 0, dbsa)
          kk = (iqx,itt,iphi*40+20)
          print('kim: '+str(kk)+str(sigs) + str([bsa,dbsa]))
          if kk in systmap:
            grsys.SetPoint(grsys.GetN(), iphi*40, -0.2)
            grsys.SetPointError(grsys.GetN()-1, 0, systmap[kk])
            grsys.SetPoint(grsys.GetN(), iphi*40+40, -0.2)
            grsys.SetPointError(grsys.GetN()-1, 0, systmap[kk])
          print('systbsa: {} {} {} {} {}'.format(iqx,itt,iphi*40+20,bsa,dbsa))

    grbsa.Draw("AP*")
    #fas = ROOT.TF1("fas","[0]*sin(x*TMath::DegToRad())",0,360)
    fas = ROOT.TF1("fas","[0]*sin(x*TMath::DegToRad())+[1]",0,360)
    fas.SetParameter(0,0.1)
    grbsa.Fit(fas,'Q')
    c2.cd()
    frame = c2.DrawFrame(0,-0.2,360,0.2)
    frame.SetTitle(";#phi [#circ];BSA")
    fraxis = frame.Clone()
    frame.SetAxisColor(16,'xy')
    fraxis.Draw("axissame")

    grsys.SetFillColor(17)
    grsys.Draw('3')

    grbsa.Draw("P0")
    grbsa.SetLineWidth(2)
    grbsa.SetMarkerSize(2)
    grbsa.SetMarkerStyle(20)

    c2.Print(pdfname.replace('.pdf','.phi.pdf'))

    kf = kinf(q2m,xbm,10.604)
    #kf=1
    grt.SetPoint(grt.GetN(), ttm, fas.GetParameter(0)/kf)
    grt.SetPointError(grt.GetN()-1, 0, fas.GetParError(0)/kf)


    kk = (iqx,itt)
    if kk in systmap:
      grtsys.SetPoint(grtsys.GetN(), tts[iqx][itt], -0.04)
      grtsys.SetPointError(grtsys.GetN()-1, 0, systmap[kk])
      grtsys.SetPoint(grtsys.GetN(), tts[iqx][itt+1], -0.04)
      grtsys.SetPointError(grtsys.GetN()-1, 0, systmap[kk])
 

    print('kim: {} {} {} {} {} 1'.format(xbm,q2m,ttm, fas.GetParameter(0)/kf, fas.GetParError(0)/kf))
    print('systalu: {} {} {} {} {} {} {} 1'.format(iqx,itt,xbm,q2m,ttm, fas.GetParameter(0)/kf, fas.GetParError(0)/kf))



  c1.cd()
  frame = c1.DrawFrame(0,-0.04,2,0.34)
  frame.SetTitle(";-t [GeV^{2}];#sigma_{LT'} / #sigma_{0}")
  #frame = c1.DrawFrame(0,-0.04,2,0.15)
  #frame.SetTitle(";-t [GeV^{2}];BSA")
  fraxis = frame.Clone()
  frame.SetAxisColor(16,'xy')
  fraxis.Draw("axissame")

  grtsys.SetFillColor(17)
  grtsys.Draw('3')

  grt.SetLineWidth(3)
  grt.SetMarkerSize(3)
  grt.SetMarkerStyle(20)
  grt.Draw("P")

  if len(rootnames)>1:
    datasetname = 'CLAS12:'
  elif 'inb.' in rootnames[0]:
    datasetname = 'INB:'
  elif '2019' in rootnames[0]:
    datasetname = '2019:'
  else:
    datasetname = 'OUTB:'
  lat.DrawLatexNDC(0.45,0.935,'{}: Q^{{2}}={:.3f}, x_{{B}}={:.3f}'.format(datasetname,q2m,xbm))


  cols = [1,2,4,6,3,7,8,9]
  ipos = 0
  for name,igr,q2,xb in sorted(grmap, key=lambda it:it[1]):
    gr6 = grmap[(name,igr,q2,xb)]
    if abs(xb-xbm)<0.1:
      gr6.SetLineWidth(3)
      gr6.SetMarkerSize(3)
      gr6.SetMarkerStyle(21)
      ipos += 1
      gr6.SetMarkerColor(cols[igr+1])
      gr6.SetLineColor(cols[igr+1])
      gr6.Draw("L" if 'gk' in name.lower() else 'P')
      lat.DrawLatexNDC(0.45,0.935-ipos*0.075,'{}: Q^{{2}}={:.3f}, x_{{B}}={:.3f}'.format(name,q2,xb)).SetTextColor(cols[igr+1])

  for xbi in set(xbis[:1]):
    xi = xbm / (2-xbm)
    tmin = -4*0.938**2 * xi**2/(1-xi**2)
    #ll.DrawLine(-tmin,-0.04,-tmin,0.34).SetLineColor(4)
  
  c1.Print(pdfname)




for iqx in range(len(q2s)-1):
  htt = htts[iqx]
  tot = htt.Integral()
  itts = [0]
  for itt in range(1,htt.GetNbinsX()+1):
    if htt.Integral(itts[-1]+1, itt) > tot/3.0 or itt==htt.GetNbinsX():
      itts.append(itt)
  print(str([round(htt.GetBinLowEdge(ii+1),2) for ii in itts])+",")



fl = ROOT.TF1("fl","-0.4836 + 21.3*x - 8.465*x*x",0,1)
fr = ROOT.TF1("fr","3.12*x/(1-x)",0,1)

#y = x(w2-M)/(1-x)


c0 = ROOT.TCanvas("c0","c0",800,1100)
c0.SetMargin(0.12,0.06,0.12,0.01)
c0.cd()
hqx.Draw("colz")

for obj in [ll,fl,fr]:
  obj.SetLineColor(2)
  obj.SetLineWidth(3)
lls = []
for y0,y1 in zip(q2s[:-1], q2s[1:]):
  xxs = []
  for f1 in [fl, fr]:
    x1 = f1.GetX(y0)
    x2 = f1.GetX(y1)
    xxs.append((x1,x2))
    fx = f1.Clone()
    fx.SetRange(x1,x2)
    fx.Draw("same")
    lls.append(fx)
  ll.DrawLine(xxs[0][0],y0, xxs[1][0],y0)
  ll.DrawLine(xxs[0][1],y1, xxs[1][1],y1)
c0.Print(rootnames[0]+'.bins.pdf[')
c0.Print(rootnames[0]+'.bins.pdf')
#htx.Draw("colz")
#c0.Print(rootnames[0]+'.bins.pdf')
c0.Print(rootnames[0]+'.bins.pdf]')




c1.cd()
for htphi in htphis:
  htphi.Draw("colz")
  c1.Print(pdfname)

ROOT.gStyle.SetOptFit(1)


hl = ROOT.TH1F("hl","left gg edge",100,0.07,0.2)
hm = ROOT.TH1F("hm","gg mean",100,0.07,0.2)
hr = ROOT.TH1F("hr","right gg edge",100,0.07,0.2)

for ipage in range(0,len(hmggs),9):
  c1.Clear()
  c1.Divide(3,3,0.001,0.001)
  for ipad in range(len(phis)-1):
    c1.cd(ipad+1).SetMargin(0.1,0.01,0.1,0.01)
    hmggs[ipage+ipad].Draw()
    f1 = ROOT.TF1("f"+str(ipad),"gaus(0)+pol1(3)",0.07,0.2)
    f1.SetParameters(hmggs[ipage+ipad].GetMaximum(),0.135,0.01,1)

    if hmggs[ipage+ipad].GetEntries()>100:
      hmggs[ipage+ipad].Fit(f1,'Q')
      hl.Fill(f1.GetParameter(1)-3*f1.GetParameter(2))
      hm.Fill(f1.GetParameter(1))
      hr.Fill(f1.GetParameter(1)+3*f1.GetParameter(2))

  c1.Print(pdfname)


c1.Clear()

for hgg in hggs:
  fgg = ROOT.TF1("fgg","gaus(0)+pol1(3)",0.07,0.2)
  fgg.SetParameters(hgg.GetMaximum(), 0.135, 0.013, 1, 1)
  hgg.Fit(fgg,'Q')
  c1.Print(pdfname)


hl.Draw()
c1.Print(pdfname)
hm.Draw()
c1.Print(pdfname)
hr.Draw()
c1.Print(pdfname)


#gr = ROOT.TGraph()
#xax = h2.GetXaxis()
#qax = h2.GetYaxis()
#for iby in range(qax.FindBin(2), qax.FindBin(10)):
#  ibx = h2.ProjectionX('hx',iby,iby).FindFirstBinAbove(0)
#  gr.SetPoint(gr.GetN(), xax.GetBinCenter(ibx), qax.GetBinCenter(iby))
#ROOT.gStyle.SetOptFit(1)
#fpol = ROOT.TF1("fpol", "pol2", 0.1,0.7)
#gr.Draw("A*")
#gr.Fit(fpol,'Q')
#c1.Print(pdfname)

c1.Print(pdfname+']')
c2.Print(pdfname.replace('.pdf', '.phi.pdf]'))

print('used:')
for tt in tts:
  print(tt)
