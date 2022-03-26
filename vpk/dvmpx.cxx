#include<iostream>
#include<TMath.h>

/**********************************************************************************/
/*                             pi0/eta 2019 Model begins here                              */
/**********************************************************************************/
/* https://github.com/vkubarovsky/Model-for-the-pi0-eta-exclusive-cross-section   */
/*

     REAL FUNCTION DVMPX(del2,xb,Q2, Phi_g,E,heli,MESONMASS)
C
c  dsigma/dQ2 dX dt dphi for ep-->ep pi0/eta
C
C exc pi0/eta x-section
c
c input:
c del2=t (negative GeV^2)           NEGATIVE !!!!
c xb,Q2 x and Q^2 (GeV^2)
c Phi_g angle in the photon frame (radians)
c E energy of the electron in GeV
c heli electron helicity -1.0 or +1.0
c MESONMASS is the mass of the pi0 or eta

PARAMETERS pi0/eta

       DATA P/
     *   6.2624,  0.8032, -0.1152,  1.7431,  0.0,
     *  92.9727,  3.4071, -1.9099, -0.1199,  0.0, 16.8603,  2.1722,
     *   6.9439,  1.7523, -1.2875,  0.6822,  0.0,
     *  17.0423,  1.1264,  0.0491,  1.6594,  0.0, 21.6379,  3.9873/
*/

Double_t EBeam = 7.535;
const Double_t ProtonMass = 0.93827;
const Double_t Pi0mass = 0.134976;
const Double_t Etamass = 0.5473;
const Double_t me2 = TMath::Power(0.000511,2.);
const Double_t alpha_em = 1./137.036;


double dvmpx(double t, double xb, double Q2, double PHI_G , double E, int heli=0) {
 //  dsigma/dQ2 dX dt dphi for ep-->ep pi0

  /* pi0 */  	
double p[12]={6.2624, 0.8032, -0.1152, 1.7431, 0.0, 92.9727, 3.4071, -1.9099, -0.1199,  0.0, 16.8603, 2.1722};

/* eta 
double p[12]={6.9439, 1.7523, -1.2875, 0.6822,  0.0, 17.0423, 1.1264,  0.0491, 1.6594,  0.0, 21.6379,  3.9873};
*/
 
 if(xb<Q2/(2.0*ProtonMass*E)||xb>1.0)return 0.;
 double nu  = Q2/(2.*ProtonMass*xb);
 double y = nu/E;
 double e1 = TMath::Power(y*xb*ProtonMass,2)/Q2;
 double EPS = (1.0-y-e1)/(1-y+y*y/2+e1);
 if(EPS<0.||EPS>1.)return 0.;
 double W2  = ProtonMass*ProtonMass + 2.0*ProtonMass*nu - Q2;
 if(W2<TMath::Power(ProtonMass+Pi0mass,2.))return 0.;
 double W = TMath::Sqrt(W2);
 double E1cm = ProtonMass*(ProtonMass + nu)/W;
 double P1cm = ProtonMass*TMath::Sqrt(nu*nu+Q2)/W;
 double E2cm = (W2 + ProtonMass*ProtonMass-Pi0mass*Pi0mass)/(2.*W);
 if(E2cm<ProtonMass)return 0.;
 double P2cm = TMath::Sqrt(E2cm*E2cm-ProtonMass*ProtonMass);

 double E3cm = (W2-Pi0mass*Pi0mass+ProtonMass*ProtonMass)/(2*W);
 double P3cm = TMath::Sqrt(E3cm*E3cm-ProtonMass*ProtonMass);

 double tmax = 2.0*(ProtonMass*ProtonMass - E1cm*E2cm + P1cm*P2cm);
 double tmin = 2.0*(ProtonMass*ProtonMass - E1cm*E2cm - P1cm*P2cm);
 if(t<tmin||t>tmax)return 0.;
 double FLUXW = alpha_em/(2*TMath::Pi()) * y*y/(1-EPS)*(1-xb)/xb/Q2;

 tmin = -(TMath::Power(Q2+Pi0mass*Pi0mass,2)/4./W2-TMath::Power(P1cm-P3cm,2));
 double SLOPE = 2.*1.1*TMath::Log(xb);

 double T  = TMath::Abs(t);
 double T0 = TMath::Abs(tmin);
// cout << "T-T0=" << T-T0 << " , T=" << T << " , T0=" << T0 << endl;

 double HT = p[0]*TMath::Exp(-(p[1]+p[2]*(TMath::Log(xb)-TMath::Log(0.15)))*T)*TMath::Power(Q2,p[3]/2.);
 double ET = p[5]*TMath::Exp(-(p[6]+p[7]*(TMath::Log(xb)-TMath::Log(0.15)))*T)*TMath::Power(Q2,p[8]/2.);
 double HTEBAR = p[10]*TMath::Exp(-p[11]*T);
 
 double pi = TMath::Pi();  
 double hc2= 389379.36;
 double ProtonMass2 = ProtonMass*ProtonMass;
 double ksi = xb/(2-xb)*(1.+ProtonMass2/Q2);
 double phase = 16.*pi*(W2-ProtonMass2)*TMath::Sqrt(W2*W2+Q2*Q2+ProtonMass2*ProtonMass2+2.*W2*Q2-2.*W2*ProtonMass2+2.*Q2*ProtonMass2);

 double S_T  = hc2*4.*pi*alpha_em/(2.*phase*Q2*Q2)*((1.-ksi*ksi)*HT*HT+(T-T0)/(8.*ProtonMass2) * ET*ET);
 double S_L  = 0.0;
 double S_LT = hc2*4.*pi*alpha_em/(TMath::Sqrt(2.)*phase*TMath::Power(Q2,1.5)) * ksi*TMath::Sqrt(1.-ksi*ksi)*TMath::Sqrt(T-T0)/(2.*ProtonMass)*HTEBAR*HTEBAR;
 double S_TT = -hc2*4.*pi*alpha_em/(2.*phase*Q2*Q2)*(T-T0)/(8.*ProtonMass2) * ET*ET;
 double S_LTP = 0.;

 /*
 double S_T =  (814.59900+ 0.0*TMath::Sqrt(T-T0))*TMath::Exp(SLOPE*T*0.44703799)*    1./TMath::Power(Q2+1.0876700,1.6934299);
 double S_L = Q2*(1404.0500 + 0.0*(T-T0))*TMath::Exp(SLOPE*T*0.69298601)*     1./TMath::Power(Q2+1.0876700,1.6934299);
 double S_LT =         608.48901*TMath::Sqrt(T-T0)*TMath::Exp(SLOPE*T*1.0290900)*    1./TMath::Power(Q2+1.0876700,1.6934299);
 double S_TT =              -5205.8999*(T-T0)*TMath::Exp(SLOPE*T*0.98690498)* 1./TMath::Power(Q2+1.0876700,1.6934299);
 double S_LTP = 0.;
 */
 //cout << FLUXW << " S_T=" << S_T << " " << EPS << " " << S_L << " " << S_TT << " S_LT=" << S_LT << endl;
 std::cout << EPS<<" "<< S_T << " " << S_L << " " << S_TT << " " << S_LT << std::endl;


 double DVMPX = FLUXW/(2.*TMath::Pi())*( S_T + EPS*S_L + EPS * S_TT  * TMath::Cos(2*PHI_G) + TMath::Sqrt(2.*EPS*(1.+EPS))*S_LT * TMath::Cos(PHI_G) + heli*TMath::Sqrt(2.*EPS*(1.-EPS))*S_LTP * TMath::Sin(2*PHI_G) ) ;
      if(DVMPX<0.) DVMPX=0.;
 return DVMPX;
}


/**********************************************************************************/

int main(int argc, char** argv) {
  double tt = atof(argv[1]);
  double xb = atof(argv[2]);
  double Q2 = atof(argv[3]);
  double phi = atof(argv[4]);
  double ee = EBeam;
  if(argc==6) ee = atof(argv[5]);

  dvmpx(tt, xb, Q2, TMath::DegToRad()*phi, ee);

  return 0;
}
