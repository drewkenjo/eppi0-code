#include "ROOT/RDataFrame.hxx"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/VectorUtil.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TDatabasePDG.h"
#include "TRandom2.h"
#include "TError.h"
#include "TChain.h"
#include <iostream>

double Chi2(const double *xx) {
  double pars[6][2];
  int ipar=0;
  for(int isec=0;isec<6;isec++) {
    double dp05=xx[ipar++], dp3=xx[ipar++];

    pars[isec][0] = (6*dp05 - dp3)/5.0;
    pars[isec][1] = (dp3 - dp05)/2.5;
  }

  double mele = 0.00051099891, mpip = 0.13957, mpim = 0.13957, mpro = 0.938272;

  auto beam = ROOT::Math::PxPyPzMVector(0,0,10.604,mele);
  auto targ = ROOT::Math::PxPyPzMVector(0,0,0,mpro);

  auto dpp = [&](float px, float py, float pz, int sec) {
    double pp = sqrt(px*px + py*py + pz*pz);

    double p0=pars[sec-1][0],
           p1=pars[sec-1][1];

    double dp = p0 + p1*pp
    return dp/pp;
  };

  ROOT::RDataFrame df("h22", "esec*.epippimp.root");

  auto rdf = df.Define("dmePipPimX",[&](float ex,float ey,float ez,float esec,float px,float py,float pz,float psec) {
      double fe = dpp(ex,ey,ez,esec,0) + 1;
      double fpip = dpp(pipx,pipy,pipz,pipsec,1) + 1;
      double fpim = dpp(pimx,pimy,pimz,pimsec,2) + 1;

      auto ele = ROOT::Math::PxPyPzMVector(ex*fe,ey*fe,ez*fe,mele);
      auto pip = ROOT::Math::PxPyPzMVector(pipx*fpip,pipy*fpip,pipz*fpip,mpip);
      auto pim = ROOT::Math::PxPyPzMVector(pimx*fpim,pimy*fpim,pimz*fpim,mpim);
 
      auto ePipPimX = beam+targ - ele-pip-pim;
      double mePipPimX = ePipPimX.mass();
      return pow(mePipPimX-0.02-0.015, 2.0);
    }, {"ex","ey","ez","esec","px","py","pz","psec"});

  return rdf.Sum("dmePipPimX").GetValue();
}



int main(int argc, char** argv)
{
  ROOT::EnableImplicitMT(18);

  ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  minimum->SetMaxFunctionCalls(100000);
  minimum->SetMaxIterations(10000);
  minimum->SetTolerance(0.01);
  minimum->SetPrintLevel(100);

  int npars = 12;
  ROOT::Math::Functor f(&Chi2, npars);
  double step[npars];
  for(int ipar=0;ipar<npars;ipar++)
    step[ipar]=0.005;

  minimum->SetFunction(f);
  for(int ipar=0;ipar<npars;ipar++) {
    minimum->SetVariable(ipar, Form("par%02d",ipar), 0.01, step[ipar]);
    minimum->SetVariableLimits(ipar, -0.1,0.3);
  }

  minimum->Minimize();

  const double *xs = minimum->X();

  std::cout << "Minimum: f() = " << minimum->MinValue()  << std::endl;
  for(int ipar=0;ipar<npars;ipar++)
    std::cout<<xs[ipar]<<std::endl;

  return 0;
}



