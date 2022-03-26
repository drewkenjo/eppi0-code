auto smear = [&](float xx, float yy, float zz, double mass) {
  TLorentzVector V4;
  V4.SetXYZM(xx, yy, zz, mass);

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
  sPh += TMath::DegToRad() * phiR * myMC->Gaus();
  sTh += TMath::DegToRad() * theR * myMC->Gaus();
  sP  += momR  * myMC->Gaus() *  V4->P() ;
  // EB_rec_mom = GEN_mom + resolution_momentum x gaussian x GEN_mom
  // EB_rec_ang = GEN_ang + resolution_angle x gaussian 

  V4->SetE( TMath::Sqrt( sP*sP + inM*inM )  );
  V4->SetRho( sP );
  V4->SetTheta( sTh );
  V4->SetPhi( sPh );
  return V4;
}

auto newele = smear(ex,ey,ez,0.000511);
auto newpro = smear(px,py,pz,0.938);

return array<double, 6>{newele.X(), newele.Y(), newele.Z(), newpro.X(), newpro.Y(), newpro.Z()};
