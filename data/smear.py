#!/usr/bin/env python

import ROOT, sys

fname = sys.argv[1]
rdf = ROOT.ROOT.RDataFrame("h22", fname)
names0 = list(rdf.GetColumnNames())


#_________________________________


rdf = rdf.Define("vals", """
auto smear = [&](float xx, float yy, float zz, double mass) {
  auto V4 = new TLorentzVector();
  V4->SetXYZM(xx, yy, zz, mass);

  double sP  = V4->P();
  double sTh = V4->Theta();
  double sPh = V4->Phi();

  //calculate resolutions
  double sThD = TMath::RadToDeg()*sTh;
  double momS1 = 0.0184291 -0.0110083*sThD + 0.00227667*sThD*sThD -0.000140152*sThD*sThD*sThD + 3.07424e-06*sThD*sThD*sThD*sThD;
  double momS2 = 0.02*sThD;
  double momR  = 0.01 * TMath::Sqrt( TMath::Power(momS1*sP,2) + TMath::Power(momS2,2) );
  momR *= 2.0;

  double theS1 = 0.004*sThD + 0.1;
  double theS2 = 0;
  double theR  = TMath::Sqrt(TMath::Power(theS1*TMath::Sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + TMath::Power(theS2,2) );
  theR *= 2.5;

  double phiS1 = 0.85-0.015*sThD;
  double phiS2 = 0.17-0.003*sThD;
  double phiR  = TMath::Sqrt(TMath::Power(phiS1*TMath::Sqrt(sP*sP+0.13957*0.13957)/(sP*sP),2) + TMath::Power(phiS2,2) );
  phiR *= 3.5;

  //overwrite EB
  sPh += TMath::DegToRad() * phiR * gRandom->Gaus();
  sTh += TMath::DegToRad() * theR * gRandom->Gaus();
  sP  += momR  * gRandom->Gaus() *  V4->P() ;
  // EB_rec_mom = GEN_mom + resolution_momentum x gaussian x GEN_mom
  // EB_rec_ang = GEN_ang + resolution_angle x gaussian 

  V4->SetE( TMath::Sqrt( sP*sP + mass*mass)  );
  V4->SetRho( sP );
  V4->SetTheta( sTh );
  V4->SetPhi( sPh );
  return V4;
};

auto newele = smear(ex,ey,ez,0.000511);
auto newpro = smear(px,py,pz,0.938);

return array<double, 6>{newele->X(), newele->Y(), newele->Z(), newpro->X(), newpro->Y(), newpro->Z()};
""")

names1 = list(rdf.GetColumnNames())
names1.remove("ex")
names1.remove("ey")
names1.remove("ez")
names1.remove("px")
names1.remove("py")
names1.remove("pz")

rdf.Snapshot("h22","/tmp/smeared.tmp.root", names1)

rdf = ROOT.ROOT.RDataFrame("h22","/tmp/smeared.tmp.root")
rdf = rdf.Define("ex","(float)vals[0]").Define("ey","(float)vals[1]").Define("ez","(float)vals[2]")
rdf = rdf.Define("px","(float)vals[3]").Define("py","(float)vals[4]").Define("pz","(float)vals[5]")
rdf.Snapshot("h22", fname.replace(".root",".smear.root"), names0)

