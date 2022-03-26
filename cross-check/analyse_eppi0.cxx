/// //////////////////////////////////////////////////////////////////////////////
///
///  ROOT Macro to analyse filtered data from CLAS12
///
///  To execute run ROOT and compile macro with .L analyse_eppi0.cxx++
///
///  Then call function:  analyse_eppi0("output_7/skim4_inb_pid_qa_andrey_selected", "output/analysed_eppi0_inb.root", false)
///                       analyse_eppi0("output_7/skim4_outb_pid_qa_andrey_selected", "output/analysed_eppi0_outb.root", true)                     
///
///
///                                               folder                  outputfile                       outb?
///
/// //////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h> 
#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TMath.h"
#include <vector>
#include <TAxis.h>
#include <TLorentzRotation.h>
using namespace std;

/// ////////////////////////////////////////////////////////////////////////////////////////////////////
/// settings:

double Ebeam = 10.6041;

/// particle momentum and theta limits:

double electron_p_min = 2.65;
double electron_theta_min = 0;
double electron_theta_max = 180;

double proton_p_min = 0;
double proton_p_max = 11.0;
double proton_theta_min = 0;
double proton_theta_max = 180;

double photon_p_min = 0.6;  
double photon_theta_min = 0;
double photon_theta_max = 180;

double pion_p_min = 0.6;  
double pion_theta_min = 0;
double pion_theta_max = 180;

// proton PID cuts:

double prot_chi2_pid_cut_low = -100;  //-2.64;
double prot_chi2_pid_cut_high = 100;  //2.64;

// pi0 ivariant mass cuts:

double pi0_min = 0.102485;
double pi0_max = 0.1644568;

//double pi0_min = 0.07;
//double pi0_max = 0.19;


bool correct = true;

/// //////////////////////////////

bool print = true;  // create output eventfile

/// //////////////////////////////
/// good run list

int run_in[] = {
5032, 5036, 5038, 5039, 5040, 5041, 5043, 5045, 5046, 5047, 5051, 5052, 5053, 5116, 5117, 5119, 5120, 5124, 5125, 5126, 
5127, 5128, 5129, 5130, 5137, 5138, 5139, 5153, 5158, 5159, 5160, 5162, 5163, 5164, 5165, 5166, 5167, 5168, 5169, 5180, 
5181, 5182, 5183, 5189, 5190, 5191, 5193, 5194, 5195, 5196, 5197, 5198, 5199, 5200, 5201, 5202, 5203, 5204, 5205, 5206, 
5208, 5211, 5212, 5215, 5216, 5219, 5220, 5221, 5222, 5223, 5225, 5229, 5230, 5231, 5232, 5233, 5234, 5235, 5237, 5238, 
5239, 5247, 5248, 5249, 5250, 5252, 5253, 5257, 5258, 5259, 5261, 5262, 5300, 5301, 5302, 5303, 5304, 5305, 5306, 5307, 
5310, 5311, 5315, 5316, 5317, 5318, 5319, 5320, 5323, 5324, 5325, 5333, 5334, 5335, 5336, 5339, 5340, 5341, 5342, 5343, 
5344, 5345, 5346, 5347, 5349, 5351, 5354, 5355, 5356, 5357, 5358, 5359, 5360, 5361, 5362, 5366, 5367, 5368, 5369, 5370, 
5371, 5372, 5373, 5374, 5375, 5376, 5377, 5378, 5379, 5380, 5381, 5382, 5383, 5386, 5390, 5391, 5392, 5393, 5394, 5398, 
5399, 5400, 5401, 5402, 5403, 5404, 5406, 5407
};  // 175 inbending runs


//int run_in[] = {5032, 5036, 5038, 5039, 5040, 5041, 5043, 5045, 5046, 5047};


int run_out[] = {
5422, 5423, 5424, 5425, 5426, 5428, 5429, 5430, 5431, 5432, 5434, 5435, 5436, 5437, 5438, 5439, 5440, 5441,
5442, 5443, 5444, 5445, 5447, 5448, 5449, 5450, 5451, 5452, 5453, 5454, 5455, 5456, 5457, 5460, 5462, 5464,
5465, 5466, 5467, 5468, 5469, 5470, 5471, 5472, 5473, 5474, 5475, 5476, 5478, 5479, 5480, 5481, 5482, 5483,
5485, 5486, 5487, 5495, 5496, 5497, 5498, 5499, 5500, 5504, 5505, 5507, 5516, 5517, 5518, 5519, 5520, 5521,
5522, 5523, 5524, 5525, 5526, 5527, 5528, 5530, 5532, 5533, 5534, 5535, 5536, 5537, 5538, 5540, 5541, 5542,
5543, 5544, 5545, 5546, 5547, 5548, 5549, 5550, 5551, 5552, 5554, 5555, 5556, 5557, 5558, 5559, 5561, 5562,
5564, 5565, 5566, 5567, 5569, 5570, 5571, 5572, 5573, 5574, 5577, 5578, 5581, 5584, 5586, 5589, 5590, 5591,
5592, 5594, 5595, 5597, 5598, 5600, 5601, 5602, 5603, 5604, 5606, 5607, 5609, 5610, 5611, 5612, 5613, 5614,
5615, 5616, 5617, 5618, 5619, 5620, 5621, 5623, 5624, 5625, 5626, 5627, 5628, 5629, 5630, 5631, 5632, 5633,
5634, 5635, 5637, 5638, 5639, 5641, 5643, 5644, 5645, 5646, 5647, 5648, 5649, 5650, 5651, 5652, 5654, 5655,
5656, 5662, 5663, 5664, 5665, 5666
};  // 186 outbending runs



/// /////////////////////////////

const static int BUFFER = 20;

const static int BUFFER_pi0   = 20;
const static int BUFFER_prot  = 20;
const static int BUFFER_gamma = 20;

/// ////////////////////////////////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////////////////////////////////

double kin_mismass_X2(TLorentzVector ele, TLorentzVector part1);
double kin_mismass2_X2(TLorentzVector ele, TLorentzVector part1);
double kin_mismass_X3(TLorentzVector ele, TLorentzVector part1, TLorentzVector part2);
double kin_mismass2_X3(TLorentzVector ele, TLorentzVector part1, TLorentzVector part2);

double kin_E_X3(TLorentzVector ele, TLorentzVector part1, TLorentzVector part2);
double kin_PT_X3(TLorentzVector ele, TLorentzVector part1, TLorentzVector part2);
double kin_Px_X3(TLorentzVector ele, TLorentzVector part1, TLorentzVector part2);
double kin_Py_X3(TLorentzVector ele, TLorentzVector part1, TLorentzVector part2);
double kin_Pz_X3(TLorentzVector ele, TLorentzVector part1, TLorentzVector part2);
double kin_cone_pi0(TLorentzVector ele, TLorentzVector part1, TLorentzVector pi0);

double kin_W(TLorentzVector ele);
double kin_Q2(TLorentzVector ele);
double kin_x(TLorentzVector ele);
double kin_y(TLorentzVector ele);
double kin_epsilon(TLorentzVector ele);

double kin_cmphi(TLorentzVector ele, TLorentzVector hadron);
double kin_cmcostheta(TLorentzVector ele, TLorentzVector hadron);
double kin_t(TLorentzVector ele, TLorentzVector hadron);

double alpha_p1p2(TLorentzVector particle1, TLorentzVector particle2);

double corr_inb_final(float px, float py, float pz, int sec, int ivec);
double corr_outb_final(float px, float py, float pz, int sec, int ivec);


/// ////////////////////////////////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////////////////////////////////


Int_t analyse_eppi0( Char_t *inFolder, Char_t *outputfile, bool outbending)
{

  double x0 = 0.15, y0 = 2.63;
  double x1 = 0.45, y1 = 4.5;
  double x01 = 0.34;
  double x4 = 0.45, y4 = 4.8;
  double x5 = 0.666, y5 = 6.2;

  int t_bins_bin1 = 3; 
  int t_bins_bin2 = 3; 
  int t_bins_bin3 = 3; 
  int t_bins_bin4 = 3; 
  int t_bins_bin5 = 3; 

  double t_bin1[] = {0.00, 0.42, 0.9, 2.0};
  double t_bin2[] = {0.00, 0.42, 0.9, 2.0};
  double t_bin3[] = {0.00, 0.46, 0.9, 2.0};
  double t_bin4[] = {0.05, 0.62, 0.98, 2.0};
  double t_bin5[] = {0.05, 0.82, 1.26, 2.0};

  double t_bin1_out[] = {0.00, 0.62, 1.09, 2.0};
  double t_bin2_out[] = {0.00, 0.60, 0.94, 2.0};
  double t_bin3_out[] = {0.00, 0.62, 1.02, 2.0};
  double t_bin4_out[] = {0.05, 0.70, 1.08, 2.0};
  double t_bin5_out[] = {0.05, 0.94, 1.38, 2.0};

  if(outbending == true){
    for(Int_t i = 0; i < 3; i++){ 
      t_bin1[i] = t_bin1_out[i];
      t_bin2[i] = t_bin2_out[i];
      t_bin3[i] = t_bin3_out[i];
      t_bin4[i] = t_bin4_out[i];
      t_bin5[i] = t_bin5_out[i];
    }
  }
 

/// ///////////////////////////////////////////////////////////////////////////
	
    Char_t tmpstr[80];
    Char_t file[200];
    Double_t fraction;
    ofstream myfile;
    
/// ///////////////////////////////////////////////////////////////////////////
///  initalize the input chain:
/// ///////////////////////////////////////////////////////////////////////////

    cout << "Initalize the input tree ... " << endl;

    Char_t *inTree="events_epgsX";
    TChain chain(inTree);

    int runs = sizeof(run_in)/sizeof(run_in[0]);
    if(outbending == true) runs = sizeof(run_out)/sizeof(run_out[0]);

    
    for(int i = 0; i < runs; i++){
      if(outbending == false) sprintf(file,"../../%s/skim4_selected_%i.root", inFolder, run_in[i]);
      if(outbending == true ) sprintf(file,"../../%s/skim4_selected_%i.root", inFolder, run_out[i]);
      chain.Add(file);
    }
   

    //chain.Add(inFolder);


    //chain.Add("../../output_7/skim4_inb_pid_qa_andrey_selected/skim4_selected_5051.root");  // inbending reference run
    //chain.Add("../../output_7/skim4_outb_pid_qa_andrey_selected/skim4_selected_5424.root");  // outbending reference run


    Double_t Pi = 3.14159265359;
    
    
/// /////////////////////////////////////////////////////////////////////////////    
///  create output-file for saving histograms:
/// /////////////////////////////////////////////////////////////////////////////

    TFile *out=new TFile(outputfile, "RECREATE");

    myfile.open ("eventfile.txt");
    myfile <<"eventfile for the exclsuive pi0 BSA analysis cross check" << endl;
    myfile <<"contact: Stefan Diehl (sdiehl@jlab.org)" << endl;
    myfile <<"---------------------------------------------------------------------------------------------------------------------" << endl;
    myfile <<"eventnumber	x	Q2	W	-t	phi[deg]	Mgg" << endl;
    myfile <<"---------------------------------------------------------------------------------------------------------------------" << endl;


/// ////////////////////////////////////////////////////////////////////////////
///  assign the tree branches and structs to variables:
/// ////////////////////////////////////////////////////////////////////////////

    int helicity;
    int eventnumber;
    double beam_charge;
    vector<double> *p4_ele_px = 0;
    vector<double> *p4_ele_py = 0;
    vector<double> *p4_ele_pz = 0;
    vector<double> *p4_ele_E = 0;
    vector<double> *p4_prot_px = 0;
    vector<double> *p4_prot_py = 0;
    vector<double> *p4_prot_pz = 0;
    vector<double> *p4_prot_E = 0;
    vector<double> *p4_phot_px = 0;
    vector<double> *p4_phot_py = 0;
    vector<double> *p4_phot_pz = 0;
    vector<double> *p4_phot_E = 0;
    vector<double> *p4_prot_chi2pid = 0;
    vector<int>    *p4_prot_det = 0;
    vector<int>    *p4_ele_sec = 0;
    vector<int>    *p4_phot_det = 0;
    vector<int>    *p4_phot_sec = 0;
    vector<double> *p4_prot_dcx1 = 0;
    vector<double> *p4_prot_dcy1 = 0;
    vector<double> *p4_prot_dcz1 = 0;

    if(print == true) chain.SetBranchAddress("eventnumber", &eventnumber);
    chain.SetBranchAddress("helicity", &helicity);
    chain.SetBranchAddress("beam_charge", &beam_charge);
    chain.SetBranchAddress("p4_ele_px", &p4_ele_px);
    chain.SetBranchAddress("p4_ele_py", &p4_ele_py);
    chain.SetBranchAddress("p4_ele_pz", &p4_ele_pz);
    chain.SetBranchAddress("p4_ele_E", &p4_ele_E);
    chain.SetBranchAddress("p4_prot_px", &p4_prot_px);
    chain.SetBranchAddress("p4_prot_py", &p4_prot_py);
    chain.SetBranchAddress("p4_prot_pz", &p4_prot_pz);
    chain.SetBranchAddress("p4_prot_E", &p4_prot_E);
    chain.SetBranchAddress("p4_phot_px", &p4_phot_px);
    chain.SetBranchAddress("p4_phot_py", &p4_phot_py);
    chain.SetBranchAddress("p4_phot_pz", &p4_phot_pz);
    chain.SetBranchAddress("p4_phot_E", &p4_phot_E);
    chain.SetBranchAddress("p4_prot_chi2pid", &p4_prot_chi2pid);
    chain.SetBranchAddress("phot_det", &p4_phot_det);
    chain.SetBranchAddress("prot_det", &p4_prot_det);
    chain.SetBranchAddress("ele_sec", &p4_ele_sec);
    chain.SetBranchAddress("phot_sec", &p4_phot_sec);
    chain.SetBranchAddress("p4_prot_dcx1", &p4_prot_dcx1);
    chain.SetBranchAddress("p4_prot_dcy1", &p4_prot_dcy1);
    chain.SetBranchAddress("p4_prot_dcz1", &p4_prot_dcz1);


/// ///////////////////////////////////////////////////////////////
///  create histograms
/// ///////////////////////////////////////////////////////////////

char name[400];
char title[400];

out->mkdir("event_information");				
out->cd ("event_information");

  TH1F *hist_helicity;
  TH1F *hist_beam_charge;
  TH1F *hist_ele_sec;
  TH1F *hist_phot_sec;
  TH1F *hist_phot_det;
  TH1F *hist_prot_det;
  TH1F *hist_prot_chi2pid;
  TH2F *hist_prot_chi2pid_vs_p;

  hist_helicity = new TH1F("helicity","helicity",5,-2.5,2.5);   
  hist_helicity->GetXaxis()->SetTitle("helicity");
  hist_helicity->GetYaxis()->SetTitle("counts");
  hist_beam_charge = new TH1F("beam_charge","beam charge", 10000, 0, 1000000);   
  hist_beam_charge->GetXaxis()->SetTitle("helicity");
  hist_beam_charge->GetYaxis()->SetTitle("counts");

  hist_ele_sec = new TH1F("ele_sec","sector of the electron", 6, 0.5, 6.5);   
  hist_ele_sec->GetXaxis()->SetTitle("sector");
  hist_ele_sec->GetYaxis()->SetTitle("counts");
  hist_phot_sec = new TH1F("phot_sec","sector of the photon", 6, 0.5, 6.5);   
  hist_phot_sec->GetXaxis()->SetTitle("sector");
  hist_phot_sec->GetYaxis()->SetTitle("counts");

  hist_phot_det = new TH1F("phot_det","detector of the  gamma", 4, 0.5, 4.5);   
  hist_phot_det->GetXaxis()->SetTitle("detector (2 = FD, 3 = CD)");
  hist_phot_det->GetYaxis()->SetTitle("counts");
  hist_prot_det = new TH1F("prot_det","detector of the proton", 4, 0.5, 4.5);   
  hist_prot_det->GetXaxis()->SetTitle("detector (2 = FD, 3 = CD)");
  hist_prot_det->GetYaxis()->SetTitle("counts");

  hist_prot_chi2pid = new TH1F("prot_chi2pid","chi2pid of protons", 300, -6, +6);   
  hist_prot_chi2pid->GetXaxis()->SetTitle("chi2pid");
  hist_prot_chi2pid->GetYaxis()->SetTitle("counts");
  hist_prot_chi2pid_vs_p = new TH2F("prot_chi2pid_vs_p","chi2pid vs p of protons",100, 0, 10, 300, -6, +6);   
  hist_prot_chi2pid_vs_p->GetXaxis()->SetTitle("p [GeV]");
  hist_prot_chi2pid_vs_p->GetYaxis()->SetTitle("chi2pid");



out->mkdir("particles_topology");
out->cd ("particles_topology");

  TH1F *hist_p_electron;
  TH1F *hist_phi_electron;
  TH1F *hist_theta_electron;
  TH2F *hist_theta_vs_p_electron;
  TH2F *hist_theta_vs_phi_electron;
  
  TH1F *hist_p_proton;
  TH1F *hist_phi_proton;
  TH1F *hist_theta_proton;
  TH2F *hist_theta_vs_p_proton;
  TH2F *hist_theta_vs_phi_proton;
  TH1F *hist_multiplicity_proton;
  
  TH1F *hist_p_phot;
  TH1F *hist_phi_phot;
  TH1F *hist_theta_phot;
  TH2F *hist_theta_vs_p_phot;
  TH2F *hist_theta_vs_phi_phot;
  
  TH1F *hist_p_pi0;
  TH1F *hist_phi_pi0;
  TH1F *hist_theta_pi0;
  TH2F *hist_theta_vs_p_pi0;
  TH2F *hist_theta_vs_phi_pi0;
  TH1F *hist_multiplicity_pi0;

  hist_p_electron = new TH1F("p_electron","p_{e^{-}}",220,0,11);   
  hist_p_electron->GetXaxis()->SetTitle("p_{e^{-}} [GeV]");
  hist_p_electron->GetYaxis()->SetTitle("counts");
  hist_phi_electron = new TH1F("phi_electron","#phi_{e^{-}}",360,-180,180);   
  hist_phi_electron->GetXaxis()->SetTitle("#phi_{e^{-}} [deg]");
  hist_phi_electron->GetYaxis()->SetTitle("counts");
  hist_theta_electron = new TH1F("theta_electron","#theta_{e^{-}}",180,0,45);   
  hist_theta_electron->GetXaxis()->SetTitle("#theta_{e^{-}} [deg]");
  hist_theta_electron->GetYaxis()->SetTitle("counts");
  hist_theta_vs_p_electron = new TH2F("electron_theta_vs_p","#theta_{e^{-}} vs p_{e^{-}}", 180, 1, 10, 90, 0, 45);   
  hist_theta_vs_p_electron->GetXaxis()->SetTitle("p_{e^{-}} [GeV]");
  hist_theta_vs_p_electron->GetYaxis()->SetTitle("#theta_{e^{-}} [deg]");
  hist_theta_vs_phi_electron = new TH2F("electron_theta_vs_phi","#theta_{e^{-}} vs #phi_{e^{-}}", 180, -180, 180, 140, 5, 40);   
  hist_theta_vs_phi_electron->GetXaxis()->SetTitle("#phi_{e^{-}} [deg]");
  hist_theta_vs_phi_electron->GetYaxis()->SetTitle("#theta_{e^{-}} [deg]");
  
  hist_p_proton = new TH1F("p_proton","p_{p}",200,0,8);   
  hist_p_proton->GetXaxis()->SetTitle("p_{p} [GeV]");
  hist_p_proton->GetYaxis()->SetTitle("counts");
  hist_phi_proton = new TH1F("phi_proton","#phi_{p}",180,-180,180);   
  hist_phi_proton->GetXaxis()->SetTitle("#phi_{p} [deg]");
  hist_phi_proton->GetYaxis()->SetTitle("counts");
  hist_theta_proton = new TH1F("theta_proton","#theta_{p}",270,0,135);   
  hist_theta_proton->GetXaxis()->SetTitle("#theta_{p} [deg]");
  hist_theta_proton->GetYaxis()->SetTitle("counts");
  hist_theta_vs_p_proton = new TH2F("proton_theta_vs_p","#theta_{p} vs p_{p}", 200, 0, 8, 270, 0, 135);   
  hist_theta_vs_p_proton->GetXaxis()->SetTitle("p_{p}");
  hist_theta_vs_p_proton->GetYaxis()->SetTitle("#theta_{p} [deg]");
  hist_theta_vs_phi_proton = new TH2F("proton_theta_vs_phi","#theta_{p} vs #phi_{p}", 180,-180, 180, 270, 0, 135);   
  hist_theta_vs_phi_proton->GetXaxis()->SetTitle("#phi_{p}");
  hist_theta_vs_phi_proton->GetYaxis()->SetTitle("#theta_{p} [deg]");
  
  hist_p_phot = new TH1F("p_phot","p_{#pi^{0}}",200,0,10);   
  hist_p_phot->GetXaxis()->SetTitle("p_{#pi^{0}} [GeV]");
  hist_p_phot->GetYaxis()->SetTitle("counts");
  hist_phi_phot = new TH1F("phi_phot","#phi_{#pi^{0}}",360,-180,180);   
  hist_phi_phot->GetXaxis()->SetTitle("#phi_{#pi^{0}} [deg]");
  hist_phi_phot->GetYaxis()->SetTitle("counts");
  hist_theta_phot = new TH1F("theta_phot","#theta_{#pi^{0}}",360,0,40);   
  hist_theta_phot->GetXaxis()->SetTitle("#theta_{#pi^{0}} [deg]");
  hist_theta_phot->GetYaxis()->SetTitle("counts");
  hist_theta_vs_p_phot = new TH2F("phot_theta_vs_p","#theta_{#pi^{0}} vs p_{#pi^{0}}", 200, 0, 10, 200, 0, 40);   
  hist_theta_vs_p_phot->GetXaxis()->SetTitle("p_{#pi^{0}} [GeV]");
  hist_theta_vs_p_phot->GetYaxis()->SetTitle("#theta_{#pi^{0}} [deg]");
  hist_theta_vs_phi_phot = new TH2F("phot_theta_vs_phi","#theta_{#pi^{0}} vs #phi_{#pi^{0}}", 180, -180, 180, 200, 0, 50);   
  hist_theta_vs_phi_phot->GetXaxis()->SetTitle("#phi_{#pi^{0}} [deg]");
  hist_theta_vs_phi_phot->GetYaxis()->SetTitle("#theta_{#pi^{0}} [deg]");
  
  hist_p_pi0 = new TH1F("p_pi0","p_{#pi^{0}}",200,0,10);   
  hist_p_pi0->GetXaxis()->SetTitle("p_{#pi^{0}} [GeV]");
  hist_p_pi0->GetYaxis()->SetTitle("counts");
  hist_phi_pi0 = new TH1F("phi_pi0","#phi_{#pi^{0}}",360,-180,180);   
  hist_phi_pi0->GetXaxis()->SetTitle("#phi_{#pi^{0}} [deg]");
  hist_phi_pi0->GetYaxis()->SetTitle("counts");
  hist_theta_pi0 = new TH1F("theta_pi0","#theta_{#pi^{0}}",360,0,40);   
  hist_theta_pi0->GetXaxis()->SetTitle("#theta_{#pi^{0}} [deg]");
  hist_theta_pi0->GetYaxis()->SetTitle("counts");
  hist_theta_vs_p_pi0 = new TH2F("pi0_theta_vs_p","#theta_{#pi^{0}} vs p_{#pi^{0}}", 200, 0, 10, 200, 0, 40);   
  hist_theta_vs_p_pi0->GetXaxis()->SetTitle("p_{#pi^{0}} [GeV]");
  hist_theta_vs_p_pi0->GetYaxis()->SetTitle("#theta_{#pi^{0}} [deg]");
  hist_theta_vs_phi_pi0 = new TH2F("pi0_theta_vs_phi","#theta_{#pi^{0}} vs #phi_{#pi^{0}}", 180, -180, 180, 200, 0, 50);   
  hist_theta_vs_phi_pi0->GetXaxis()->SetTitle("#phi_{#pi^{0}} [deg]");
  hist_theta_vs_phi_pi0->GetYaxis()->SetTitle("#theta_{#pi^{0}} [deg]");
  hist_multiplicity_pi0 = new TH1F("multiplicity_pi0","multiplicity pi0",9,-0.5,8.5);   
  hist_multiplicity_pi0->GetXaxis()->SetTitle("number of #pi^{0} per event");
  hist_multiplicity_pi0->GetYaxis()->SetTitle("counts");
  

out->mkdir("neutrals");				
out->cd ("neutrals");

  TH1F *hist_inv_mass;
  TH1F *hist_inv_mass2;
  
  hist_inv_mass = new TH1F("invariant_mass","M_{#gamma #gamma}",1200,0,1.2);   
  hist_inv_mass->GetXaxis()->SetTitle("M_{#gamma #gamma}");
  hist_inv_mass->GetYaxis()->SetTitle("counts");
  hist_inv_mass2 = new TH1F("invariant_mass2","M_{#gamma #gamma}^{2}",500,0.0,0.1);   
  hist_inv_mass2->GetXaxis()->SetTitle("M_{#gamma #gamma}^{2}");
  hist_inv_mass2->GetYaxis()->SetTitle("counts");


out->mkdir("missing_cuts_raw");
out->cd ("missing_cuts_raw");

  TH1F *hist_missing_mass2_eppi0;
  TH1F *hist_missing_mass2_epX;
  TH1F *hist_missing_mass_epi0X;

  TH1F *hist_missing_E_eppi0;
  TH1F *hist_missing_PT_eppi0;
  TH1F *hist_missing_Px_eppi0;
  TH1F *hist_missing_Py_eppi0;
  TH1F *hist_missing_Pz_eppi0;
  TH1F *hist_missing_cone_pi0;
  
  hist_missing_mass2_eppi0 = new TH1F("missing_mass2_eppi0","missing mass2 e p #pi^{0} raw", 2000, -0.1, 0.1);   
  hist_missing_mass2_eppi0->GetXaxis()->SetTitle("M^{2}_{miss} [GeV^{2}]");
  hist_missing_mass2_eppi0->GetYaxis()->SetTitle("counts");
  hist_missing_mass2_epX = new TH1F("missing_mass2_epX","missing mass2 e p X", 600, -1, 1);   
  hist_missing_mass2_epX->GetXaxis()->SetTitle("M^{2}_{miss} [GeV^{2}]");
  hist_missing_mass2_epX->GetYaxis()->SetTitle("counts");
  hist_missing_mass_epi0X = new TH1F("missing_mass_epi0X","missing mass e#pi^{0} X", 360, 0.6, 1.5);   
  hist_missing_mass_epi0X->GetXaxis()->SetTitle("M_{miss} [GeV]");
  hist_missing_mass_epi0X->GetYaxis()->SetTitle("counts");

  hist_missing_E_eppi0 = new TH1F("missing_E_eppi0","missing E  e p #pi^{0}", 400, -2, 4);   
  hist_missing_E_eppi0->GetXaxis()->SetTitle("E_{miss} [GeV]");
  hist_missing_E_eppi0->GetYaxis()->SetTitle("counts");
  hist_missing_PT_eppi0 = new TH1F("missing_PT_eppi0","missing P_{T}  e p #pi^{0}", 300, 0, 3);   
  hist_missing_PT_eppi0->GetXaxis()->SetTitle("P_{T} [GeV]");
  hist_missing_PT_eppi0->GetYaxis()->SetTitle("counts");
  hist_missing_Px_eppi0 = new TH1F("missing_Px_eppi0","missing P_{x}  e p #pi^{0}", 300, -3, 3);   
  hist_missing_Px_eppi0->GetXaxis()->SetTitle("P_{x} [GeV]");
  hist_missing_Px_eppi0->GetYaxis()->SetTitle("counts");
  hist_missing_Py_eppi0 = new TH1F("missing_Py_eppi0","missing P_{y}  e p #pi^{0}", 300, -3, 3);   
  hist_missing_Py_eppi0->GetXaxis()->SetTitle("P_{y} [GeV]");
  hist_missing_Py_eppi0->GetYaxis()->SetTitle("counts");
  hist_missing_Pz_eppi0 = new TH1F("missing_Pz_eppi0","missing P_{z}  e p #pi^{0}", 300, -3, 3);   
  hist_missing_Pz_eppi0->GetXaxis()->SetTitle("P_{z} [GeV]");
  hist_missing_Pz_eppi0->GetYaxis()->SetTitle("counts");
  hist_missing_cone_pi0 = new TH1F("missing_cone_pi0","missing cone pi0", 800, -20, 20);   
  hist_missing_cone_pi0->GetXaxis()->SetTitle("cone angle difference");
  hist_missing_cone_pi0->GetYaxis()->SetTitle("counts");


out->mkdir("missing_cuts_excl");
out->cd ("missing_cuts_excl");

  TH1F *hist_excl_missing_mass2_eppi0;
  TH1F *hist_excl_missing_mass2_epX;
  TH1F *hist_excl_missing_mass_epi0X;

  TH1F *hist_excl_missing_E_eppi0;
  TH1F *hist_excl_missing_PT_eppi0;
  TH1F *hist_excl_missing_Px_eppi0;
  TH1F *hist_excl_missing_Py_eppi0;
  TH1F *hist_excl_missing_Pz_eppi0;
  TH1F *hist_excl_missing_cone_pi0;
  
  hist_excl_missing_mass2_eppi0 = new TH1F("excl_missing_mass2_eppi0","exclusive missing mass2 e p #pi^{0} raw", 1000, -0.05, 0.05);   
  hist_excl_missing_mass2_eppi0->GetXaxis()->SetTitle("M^{2}_{miss} [GeV^{2}]");
  hist_excl_missing_mass2_eppi0->GetYaxis()->SetTitle("counts");
  hist_excl_missing_mass2_epX = new TH1F("excl_missing_mass2_epX","exclusive missing mass2 e p X", 500, -1, 1);   
  hist_excl_missing_mass2_epX->GetXaxis()->SetTitle("M^{2}_{miss} [GeV^{2}]");
  hist_excl_missing_mass2_epX->GetYaxis()->SetTitle("counts");
  hist_excl_missing_mass_epi0X = new TH1F("excl_missing_mass_epi0X","exclusive missing mass e#pi^{0} X", 360, 0.6, 1.5);   
  hist_excl_missing_mass_epi0X->GetXaxis()->SetTitle("M_{miss} [GeV]");
  hist_excl_missing_mass_epi0X->GetYaxis()->SetTitle("counts");

  hist_excl_missing_E_eppi0 = new TH1F("excl_missing_E_eppi0","exclusive missing E  e p #pi^{0}", 200, -1, 1);   
  hist_excl_missing_E_eppi0->GetXaxis()->SetTitle("E_{miss} [GeV]");
  hist_excl_missing_E_eppi0->GetYaxis()->SetTitle("counts");
  hist_excl_missing_PT_eppi0 = new TH1F("excl_missing_PT_eppi0","exclusive missing P_{T}  e p #pi^{0}", 250, 0, 0.5);   
  hist_excl_missing_PT_eppi0->GetXaxis()->SetTitle("P_{miss} [GeV]");
  hist_excl_missing_PT_eppi0->GetYaxis()->SetTitle("counts");
  hist_excl_missing_Px_eppi0 = new TH1F("excl_missing_Px_eppi0","exclusive missing P_{x}  e p #pi^{0}", 300, -3, 3);   
  hist_excl_missing_Px_eppi0->GetXaxis()->SetTitle("P_{miss} [GeV]");
  hist_excl_missing_Px_eppi0->GetYaxis()->SetTitle("counts");
  hist_excl_missing_Py_eppi0 = new TH1F("excl_missing_Py_eppi0","exclusive missing P_{y}  e p #pi^{0}", 300, -3, 3);   
  hist_excl_missing_Py_eppi0->GetXaxis()->SetTitle("P_{miss} [GeV]");
  hist_excl_missing_Py_eppi0->GetYaxis()->SetTitle("counts");
  hist_excl_missing_Pz_eppi0 = new TH1F("excl_missing_Pz_eppi0","exclusive missing P_{z}  e p #pi^{0}", 300, -3, 3);   
  hist_excl_missing_Pz_eppi0->GetXaxis()->SetTitle("P_{miss} [GeV]");
  hist_excl_missing_Pz_eppi0->GetYaxis()->SetTitle("counts");
  hist_excl_missing_cone_pi0 = new TH1F("excl_missing_cone_pi0","exclusive missing cone pi0", 800, -20, 20);   
  hist_excl_missing_cone_pi0->GetXaxis()->SetTitle("cone angle difference [deg]");
  hist_excl_missing_cone_pi0->GetYaxis()->SetTitle("counts");



out->mkdir("particles_exclusive_forward");
out->cd ("particles_exclusive_forward");

  TH1F *hist_forward_selected_p_electron;
  TH1F *hist_forward_selected_phi_electron;
  TH1F *hist_forward_selected_theta_electron;
  TH2F *hist_forward_selected_theta_vs_p_electron;
  TH2F *hist_forward_selected_theta_vs_phi_electron;
  
  TH1F *hist_forward_selected_p_proton;
  TH1F *hist_forward_selected_phi_proton;
  TH1F *hist_forward_selected_theta_proton;
  TH2F *hist_forward_selected_theta_vs_p_proton;
  TH2F *hist_forward_selected_theta_vs_phi_proton;
  
  TH1F *hist_forward_selected_p_pi0;
  TH1F *hist_forward_selected_phi_pi0;
  TH1F *hist_forward_selected_theta_pi0;
  TH2F *hist_forward_selected_theta_vs_p_pi0;
  TH2F *hist_forward_selected_theta_vs_phi_pi0;
  
  hist_forward_selected_p_electron = new TH1F("p_electron_forward","p_{e^{-}}",400,1,10);   
  hist_forward_selected_p_electron->GetXaxis()->SetTitle("p_{e^{-}} [GeV]");
  hist_forward_selected_p_electron->GetYaxis()->SetTitle("counts");
  hist_forward_selected_phi_electron = new TH1F("phi_electron_forward","#phi_{e^{-}}",360,-180,180);   
  hist_forward_selected_phi_electron->GetXaxis()->SetTitle("#phi_{e^{-}} [deg]");
  hist_forward_selected_phi_electron->GetYaxis()->SetTitle("counts");
  hist_forward_selected_theta_electron = new TH1F("theta_electron_forward","#theta_{e^{-}}",400,0,40);   
  hist_forward_selected_theta_electron->GetXaxis()->SetTitle("#theta_{e^{-}} [deg]");
  hist_forward_selected_theta_electron->GetYaxis()->SetTitle("counts");
  hist_forward_selected_theta_vs_p_electron = new TH2F("electron_theta_vs_p_forward","#theta_{e^{-}}vs p_{e^{-}}", 400, 1, 10, 400, 0, 40);   
  hist_forward_selected_theta_vs_p_electron->GetXaxis()->SetTitle("p_{e^{-}} [GeV]");
  hist_forward_selected_theta_vs_p_electron->GetYaxis()->SetTitle("#theta_{e^{-}} [deg]");
  hist_forward_selected_theta_vs_phi_electron = new TH2F("electron_theta_vs_phi_forward","#theta_{e^{-}}vs #phi_{e^{-}}", 360, -180, 180, 400, 0, 40);   
  hist_forward_selected_theta_vs_phi_electron->GetXaxis()->SetTitle("#phi_{e^{-}} [deg]");
  hist_forward_selected_theta_vs_phi_electron->GetYaxis()->SetTitle("#theta_{e^{-}} [deg]");

  hist_forward_selected_p_proton = new TH1F("p_proton_forward","p_{p}",400,0,4);   
  hist_forward_selected_p_proton->GetXaxis()->SetTitle("p_{p} [GeV]");
  hist_forward_selected_p_proton->GetYaxis()->SetTitle("counts");
  hist_forward_selected_phi_proton = new TH1F("phi_proton_forward","#phi_{p}",360,-180,180);   
  hist_forward_selected_phi_proton->GetXaxis()->SetTitle("#phi_{p} [deg]");
  hist_forward_selected_phi_proton->GetYaxis()->SetTitle("counts");
  hist_forward_selected_theta_proton = new TH1F("theta_proton_forward","#theta_{p}",200,0,100);   
  hist_forward_selected_theta_proton->GetXaxis()->SetTitle("#theta_{p} [deg]");
  hist_forward_selected_theta_proton->GetYaxis()->SetTitle("counts");
  hist_forward_selected_theta_vs_p_proton = new TH2F("proton_theta_vs_p_forward","#theta_{p} vs p_{p}", 240, 0, 10, 200, 0, 100);   
  hist_forward_selected_theta_vs_p_proton->GetXaxis()->SetTitle("p_{p} [GeV]");
  hist_forward_selected_theta_vs_p_proton->GetYaxis()->SetTitle("#theta_{p} [deg]");
  hist_forward_selected_theta_vs_phi_proton = new TH2F("proton_theta_vs_phi_forward","#theta_{p} vs #phi_{p}", 360, -180, 180, 200, 0, 100);   
  hist_forward_selected_theta_vs_phi_proton->GetXaxis()->SetTitle("#phi_{p} [deg]");
  hist_forward_selected_theta_vs_phi_proton->GetYaxis()->SetTitle("#theta_{p} [deg]");

  hist_forward_selected_p_pi0 = new TH1F("p_pi0_forward","p_{#pi^{0}}",400,0,9);   
  hist_forward_selected_p_pi0->GetXaxis()->SetTitle("p_{#pi^{0}} [GeV]");
  hist_forward_selected_p_pi0->GetYaxis()->SetTitle("counts");
  hist_forward_selected_phi_pi0 = new TH1F("phi_pi0_forward","#phi_{#pi^{0}}",360,-180,180);   
  hist_forward_selected_phi_pi0->GetXaxis()->SetTitle("#phi_{#pi^{0}} [deg]");
  hist_forward_selected_phi_pi0->GetYaxis()->SetTitle("counts");
  hist_forward_selected_theta_pi0 = new TH1F("theta_pi0_forward","#theta_{#pi^{0}}",400,0,40);   
  hist_forward_selected_theta_pi0->GetXaxis()->SetTitle("#theta_{#pi^{0}} [deg]");
  hist_forward_selected_theta_pi0->GetYaxis()->SetTitle("counts");
  hist_forward_selected_theta_vs_p_pi0 = new TH2F("pi0_theta_vs_p_forward","#theta_{#pi^{0}} vs p_{#pi^{0}}", 100, 0, 9, 200, 0, 40);   
  hist_forward_selected_theta_vs_p_pi0->GetXaxis()->SetTitle("p_{#pi^{0}} [GeV]");
  hist_forward_selected_theta_vs_p_pi0->GetYaxis()->SetTitle("#theta_{#pi^{0}} [deg]");
  hist_forward_selected_theta_vs_phi_pi0 = new TH2F("pi0_theta_vs_phi_forward","#theta_{#pi^{0}} vs #phi_{#pi^{0}}", 360, -180, 180, 200, 0, 40);   
  hist_forward_selected_theta_vs_phi_pi0->GetXaxis()->SetTitle("#phi_{#pi^{0}} [deg]");
  hist_forward_selected_theta_vs_phi_pi0->GetYaxis()->SetTitle("#theta_{#pi^{0}} [deg]");


out->mkdir("kinematics_exclusive_forward");				
out->cd ("kinematics_exclusive_forward");

  TH1F *hist_forward_selected_w;
  TH1F *hist_forward_selected_q2;
  TH1F *hist_forward_selected_x;
  TH1F *hist_forward_selected_t;
  TH1F *hist_forward_selected_cmphi;
  TH2F *hist_forward_selected_q2_vs_x;

  hist_forward_selected_w = new TH1F("W_forward","W",400,1.5,4.0);   
  hist_forward_selected_w->GetXaxis()->SetTitle("W [GeV]");
  hist_forward_selected_w->GetYaxis()->SetTitle("counts");
  hist_forward_selected_q2 = new TH1F("Q2_forward","Q^{2}",300,0,12);   
  hist_forward_selected_q2->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
  hist_forward_selected_q2->GetYaxis()->SetTitle("counts");
  hist_forward_selected_x = new TH1F("x_forward","x_{B}",200,0,1);   
  hist_forward_selected_x->GetXaxis()->SetTitle("x_{B}");
  hist_forward_selected_x->GetYaxis()->SetTitle("counts");
  hist_forward_selected_t = new TH1F("t_forward","t",400,0,4);   
  hist_forward_selected_t->GetXaxis()->SetTitle("-t [GeV^{2}]");
  hist_forward_selected_t->GetYaxis()->SetTitle("counts");
  hist_forward_selected_cmphi = new TH1F("phi_forward","#phi",180,0,360);   
  hist_forward_selected_cmphi->GetXaxis()->SetTitle("#phi /deg");
  hist_forward_selected_cmphi->GetYaxis()->SetTitle("counts");
  
  hist_forward_selected_q2_vs_x = new TH2F("Q2_vs_x_forward","Q^{2} vs x_{B}", 180, 0, 0.9, 600,0,12);
  hist_forward_selected_q2_vs_x->GetXaxis()->SetTitle("x_{B}");
  hist_forward_selected_q2_vs_x->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");


out->mkdir("kinematics_particles_correlations");
out->cd ("kinematics_particles_correlations");

  TH2F *hist_t_vs_p_electron;
  TH2F *hist_t_vs_theta_electron;
  TH2F *hist_t_vs_p_pion;
  TH2F *hist_t_vs_theta_pion;

  TH2F *hist_x_vs_p_electron;
  TH2F *hist_x_vs_theta_electron;
  TH2F *hist_q2_vs_p_electron;
  TH2F *hist_q2_vs_theta_electron;

  TH2F *hist_x_vs_p_electron_forward;
  TH2F *hist_x_vs_theta_electron_forward;
  TH2F *hist_q2_vs_p_electron_forward;
  TH2F *hist_q2_vs_theta_electron_forward;
  
  hist_t_vs_p_electron = new TH2F("t_vs_p_electron","t vs p electron", 200, 0, 10, 600,0,12);
  hist_t_vs_p_electron->GetXaxis()->SetTitle("p [GeV]");
  hist_t_vs_p_electron->GetYaxis()->SetTitle("-t [GeV^{2}]");
  hist_t_vs_theta_electron = new TH2F("t_vs_theta_electron","t vs theta electron", 200, 0, 40, 600,0,12);
  hist_t_vs_theta_electron->GetXaxis()->SetTitle("theta [deg]");
  hist_t_vs_theta_electron->GetYaxis()->SetTitle("-t [GeV^{2}]");
  hist_t_vs_p_pion = new TH2F("t_vs_p_pion","t vs p pion", 200, 0, 10, 600,0,12);
  hist_t_vs_p_pion->GetXaxis()->SetTitle("p [GeV]");
  hist_t_vs_p_pion->GetYaxis()->SetTitle("-t [GeV^{2}]");
  hist_t_vs_theta_pion = new TH2F("t_vs_theta_pion","t vs theta pion", 200, 0, 40, 600,0,12);
  hist_t_vs_theta_pion->GetXaxis()->SetTitle("theta [deg]");
  hist_t_vs_theta_pion->GetYaxis()->SetTitle("-t [GeV^{2}]");

  hist_x_vs_p_electron = new TH2F("x_vs_p_electron","x_{B} vs p electron", 200, 0, 10, 200,0,1);
  hist_x_vs_p_electron->GetXaxis()->SetTitle("p [GeV]");
  hist_x_vs_p_electron->GetYaxis()->SetTitle("x_{B}");
  hist_x_vs_theta_electron = new TH2F("x_vs_theta_electron","x_{B} vs theta electron", 200, 0, 40, 200,0,1);
  hist_x_vs_theta_electron->GetXaxis()->SetTitle("theta [deg]");
  hist_x_vs_theta_electron->GetYaxis()->SetTitle("x_{B}");
  hist_q2_vs_p_electron = new TH2F("q2_vs_p_electron","q2 vs p electron", 200, 0, 10, 600,0,12);
  hist_q2_vs_p_electron->GetXaxis()->SetTitle("p [GeV]");
  hist_q2_vs_p_electron->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
  hist_q2_vs_theta_electron = new TH2F("q2_vs_theta_electron","q2 vs theta electron", 200, 0, 40, 600,0,12);
  hist_q2_vs_theta_electron->GetXaxis()->SetTitle("theta [deg]");
  hist_q2_vs_theta_electron->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");

  hist_x_vs_p_electron_forward = new TH2F("x_vs_p_electron_forward","x_{B} vs p electron forward", 200, 0, 10, 200,0,1);
  hist_x_vs_p_electron_forward->GetXaxis()->SetTitle("p [GeV]");
  hist_x_vs_p_electron_forward->GetYaxis()->SetTitle("x_{B}");
  hist_x_vs_theta_electron_forward = new TH2F("x_vs_theta_electron_forward","x_{B} vs theta electron forward", 200, 0, 40, 200,0,1);
  hist_x_vs_theta_electron_forward->GetXaxis()->SetTitle("theta [deg]");
  hist_x_vs_theta_electron_forward->GetYaxis()->SetTitle("x_{B}");
  hist_q2_vs_p_electron_forward = new TH2F("q2_vs_p_electron_forward","q2 vs p electron forward", 200, 0, 10, 600,0,12);
  hist_q2_vs_p_electron_forward->GetXaxis()->SetTitle("p [GeV]");
  hist_q2_vs_p_electron_forward->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");
  hist_q2_vs_theta_electron_forward = new TH2F("q2_vs_theta_electron_forward","q2 vs theta electron forward", 200, 0, 40, 600,0,12);
  hist_q2_vs_theta_electron_forward->GetXaxis()->SetTitle("theta [deg]");
  hist_q2_vs_theta_electron_forward->GetYaxis()->SetTitle("Q^{2} [GeV^{2}]");



out->mkdir("beam_asymmetry_bins");
out->cd ("beam_asymmetry_bins");

  TH1F *histbin_t_bin1;
  TH1F *histbin_t_bin2;
  TH1F *histbin_t_bin3;
  TH1F *histbin_t_bin4;
  TH1F *histbin_t_bin5;

  histbin_t_bin1 = new TH1F("t_bins_bin1", "t bins bin1", t_bins_bin1, t_bin1);   
  histbin_t_bin1->GetXaxis()->SetTitle("-t");
  histbin_t_bin1->GetYaxis()->SetTitle("counts");
  histbin_t_bin2 = new TH1F("t_bins_bin2", "t bins bin2", t_bins_bin2, t_bin2);   
  histbin_t_bin2->GetXaxis()->SetTitle("-t");
  histbin_t_bin2->GetYaxis()->SetTitle("counts");
  histbin_t_bin3 = new TH1F("t_bins_bin3", "t bins bin3", t_bins_bin3, t_bin3);   
  histbin_t_bin3->GetXaxis()->SetTitle("-t");
  histbin_t_bin3->GetYaxis()->SetTitle("counts");
  histbin_t_bin4 = new TH1F("t_bins_bin4", "t bins bin4", t_bins_bin4, t_bin4);   
  histbin_t_bin4->GetXaxis()->SetTitle("-t");
  histbin_t_bin4->GetYaxis()->SetTitle("counts");
  histbin_t_bin5 = new TH1F("t_bins_bin5", "t bins bin5", t_bins_bin5, t_bin5);   
  histbin_t_bin5->GetXaxis()->SetTitle("-t");
  histbin_t_bin5->GetYaxis()->SetTitle("counts");
  
  
out->mkdir("kinematics_bins");
out->cd ("kinematics_bins");

  TH1F *hist_t_bin1;
  TH1F *hist_t_bin2;
  TH1F *hist_t_bin3;
  TH1F *hist_t_bin4;
  TH1F *hist_t_bin5;

  TH1F *hist_q2_bin1;
  TH1F *hist_q2_bin2;
  TH1F *hist_q2_bin3;
  TH1F *hist_q2_bin4;
  TH1F *hist_q2_bin5;
  
  TH1F *hist_x_bin1;
  TH1F *hist_x_bin2;
  TH1F *hist_x_bin3;
  TH1F *hist_x_bin4;
  TH1F *hist_x_bin5;
  
  TH1F *hist_epsilon_bin1;
  TH1F *hist_epsilon_bin2;
  TH1F *hist_epsilon_bin3;
  TH1F *hist_epsilon_bin4;
  TH1F *hist_epsilon_bin5;
  
  TH2F *hist_q2_vs_x_bin1;
  TH2F *hist_q2_vs_x_bin2;
  TH2F *hist_q2_vs_x_bin3;
  TH2F *hist_q2_vs_x_bin4;
  TH2F *hist_q2_vs_x_bin5;

  TH1F *hist_q2_t_bin1[t_bins_bin1];
  TH1F *hist_x_t_bin1[t_bins_bin1];
  TH1F *hist_t_t_bin1[t_bins_bin1];
  TH1F *hist_epsilon_t_bin1[t_bins_bin1];

  TH1F *hist_q2_t_bin2[t_bins_bin2];
  TH1F *hist_x_t_bin2[t_bins_bin2];
  TH1F *hist_t_t_bin2[t_bins_bin2];
  TH1F *hist_epsilon_t_bin2[t_bins_bin2];
  
  TH1F *hist_q2_t_bin3[t_bins_bin3];
  TH1F *hist_x_t_bin3[t_bins_bin3];
  TH1F *hist_t_t_bin3[t_bins_bin3];
  TH1F *hist_epsilon_t_bin3[t_bins_bin3];

  TH1F *hist_q2_t_bin4[t_bins_bin4];
  TH1F *hist_x_t_bin4[t_bins_bin4];
  TH1F *hist_t_t_bin4[t_bins_bin4];
  TH1F *hist_epsilon_t_bin4[t_bins_bin4];

  TH1F *hist_q2_t_bin5[t_bins_bin5];
  TH1F *hist_x_t_bin5[t_bins_bin5];
  TH1F *hist_t_t_bin5[t_bins_bin5];
  TH1F *hist_epsilon_t_bin5[t_bins_bin5];


  hist_t_bin1 = new TH1F("t_bin1", "-t bin 1", 500, 0, 8);   
  hist_t_bin2 = new TH1F("t_bin2", "-t bin 2", 500, 0, 8);
  hist_t_bin3 = new TH1F("t_bin3", "-t bin 3", 500, 0, 8);
  hist_t_bin4 = new TH1F("t_bin4", "-t bin 4", 500, 0, 8);
  hist_t_bin5 = new TH1F("t_bin5", "-t bin 5", 500, 0, 8);
  
  hist_q2_bin1 = new TH1F("q2_bin1", "q2 bin 1", 300, 0, 12);   
  hist_q2_bin2 = new TH1F("q2_bin2", "q2 bin 2", 300, 0, 12);
  hist_q2_bin3 = new TH1F("q2_bin3", "q2 bin 3", 300, 0, 12);
  hist_q2_bin4 = new TH1F("q2_bin4", "q2 bin 4", 300, 0, 12);
  hist_q2_bin5 = new TH1F("q2_bin5", "q2 bin 5", 300, 0, 12);
  
  hist_x_bin1 = new TH1F("x_bin1", "x bin 1", 200, 0, 1);   
  hist_x_bin2 = new TH1F("x_bin2", "x bin 2", 200, 0, 1);
  hist_x_bin3 = new TH1F("x_bin3", "x bin 3", 200, 0, 1);
  hist_x_bin4 = new TH1F("x_bin4", "x bin 4", 200, 0, 1);
  hist_x_bin5 = new TH1F("x_bin5", "x bin 5", 200, 0, 1);
  
  hist_epsilon_bin1 = new TH1F("epsilon_bin1", "epsilon bin 1", 200, 0, 1);   
  hist_epsilon_bin2 = new TH1F("epsilon_bin2", "epsilon bin 2", 200, 0, 1);
  hist_epsilon_bin3 = new TH1F("epsilon_bin3", "epsilon bin 3", 200, 0, 1);
  hist_epsilon_bin4 = new TH1F("epsilon_bin4", "epsilon bin 4", 200, 0, 1);
  hist_epsilon_bin5 = new TH1F("epsilon_bin5", "epsilon bin 5", 200, 0, 1);
  
  hist_q2_vs_x_bin1 = new TH2F("q2_vs_x_bin1", "q2 vs x bin 1", 200, 0, 1, 300, 0, 12); 
  hist_q2_vs_x_bin2 = new TH2F("q2_vs_x_bin2", "q2 vs x bin 2", 200, 0, 1, 300, 0, 12); 
  hist_q2_vs_x_bin3 = new TH2F("q2_vs_x_bin3", "q2 vs x bin 3", 200, 0, 1, 300, 0, 12); 
  hist_q2_vs_x_bin4 = new TH2F("q2_vs_x_bin4", "q2 vs x bin 4", 200, 0, 1, 300, 0, 12); 
  hist_q2_vs_x_bin5 = new TH2F("q2_vs_x_bin5", "q2 vs x bin 5", 200, 0, 1, 300, 0, 12); 


  for(int i = 0; i < t_bins_bin1; i++){
    sprintf(name,"q2_t_bin1_%02d", i);
    hist_q2_t_bin1[i] = new TH1F(name,name,650,0,12);   
    sprintf(name,"x_t_bin1_%02d", i);
    hist_x_t_bin1[i] = new TH1F(name,name,100,0,1);    
    sprintf(name,"t_t_bin1_%02d", i);
    hist_t_t_bin1[i] = new TH1F(name,name,500,0,8.0);      
    sprintf(name,"epsilon_t_bin1_%02d", i);
    hist_epsilon_t_bin1[i] = new TH1F(name,name,200,0,1.0);   
  }

  for(int i = 0; i < t_bins_bin2; i++){
    sprintf(name,"q2_t_bin2_%02d", i);
    hist_q2_t_bin2[i] = new TH1F(name,name,650,0,12);   
    sprintf(name,"x_t_bin2_%02d", i);
    hist_x_t_bin2[i] = new TH1F(name,name,100,0,1);    
    sprintf(name,"t_t_bin2_%02d", i);
    hist_t_t_bin2[i] = new TH1F(name,name,500,0,8.0); 
    sprintf(name,"epsilon_t_bin2_%02d", i);
    hist_epsilon_t_bin2[i] = new TH1F(name,name,200,0,1.0);   
  }

  for(int i = 0; i < t_bins_bin3; i++){
    sprintf(name,"q2_t_bin3_%02d", i);
    hist_q2_t_bin3[i] = new TH1F(name,name,650,0,12);   
    sprintf(name,"x_t_bin3_%02d", i);
    hist_x_t_bin3[i] = new TH1F(name,name,100,0,1);    
    sprintf(name,"t_t_bin3_%02d", i);
    hist_t_t_bin3[i] = new TH1F(name,name,500,0,8.0); 
    sprintf(name,"epsilon_t_bin3_%02d", i);
    hist_epsilon_t_bin3[i] = new TH1F(name,name,200,0,1.0);   
  }

  for(int i = 0; i < t_bins_bin4; i++){
    sprintf(name,"q2_t_bin4_%02d", i);
    hist_q2_t_bin4[i] = new TH1F(name,name,650,0,12);   
    sprintf(name,"x_t_bin4_%02d", i);
    hist_x_t_bin4[i] = new TH1F(name,name,100,0,1);    
    sprintf(name,"t_t_bin4_%02d", i);
    hist_t_t_bin4[i] = new TH1F(name,name,500,0,8.0);   
    sprintf(name,"epsilon_t_bin4_%02d", i);
    hist_epsilon_t_bin4[i] = new TH1F(name,name,200,0,1.0);   
  }

  for(int i = 0; i < t_bins_bin5; i++){
    sprintf(name,"q2_t_bin5_%02d", i);
    hist_q2_t_bin5[i] = new TH1F(name,name,650,0,12);   
    sprintf(name,"x_t_bin5_%02d", i);
    hist_x_t_bin5[i] = new TH1F(name,name,100,0,1);    
    sprintf(name,"t_t_bin5_%02d", i);
    hist_t_t_bin5[i] = new TH1F(name,name,500,0,8.0);
    sprintf(name,"epsilon_t_bin5_%02d", i);
    hist_epsilon_t_bin5[i] = new TH1F(name,name,200,0,1.0);   
  }
 
 

out->mkdir("beam_asymmetry");
out->cd ("beam_asymmetry");

  TH1F *hist_CMPhi_t_bin1[t_bins_bin1][2];
  TH1F *hist_CMPhi_t_bin2[t_bins_bin2][2];     
  TH1F *hist_CMPhi_t_bin3[t_bins_bin3][2];     
  TH1F *hist_CMPhi_t_bin4[t_bins_bin4][2];  
  TH1F *hist_CMPhi_t_bin5[t_bins_bin5][2];
  
  for(int i = 0; i < t_bins_bin1; i++){
    sprintf(name,"CMPhi_t_bin1_m1_%02d", i);
    hist_CMPhi_t_bin1[i][0] = new TH1F(name,name, 9, 0, 360);     
    hist_CMPhi_t_bin1[i][0]->GetXaxis()->SetTitle("#phi_{CM} [deg]");
    hist_CMPhi_t_bin1[i][0]->GetYaxis()->SetTitle("counts");
    sprintf(name,"CMPhi_t_bin1_p1_%02d", i);
    hist_CMPhi_t_bin1[i][1] = new TH1F(name,name, 9, 0, 360);    
    hist_CMPhi_t_bin1[i][1]->GetXaxis()->SetTitle("#phi_{CM} [deg]");
    hist_CMPhi_t_bin1[i][1]->GetYaxis()->SetTitle("counts");
  }

  for(int i = 0; i < t_bins_bin2; i++){
    sprintf(name,"CMPhi_t_bin2_m1_%02d", i);
    hist_CMPhi_t_bin2[i][0] = new TH1F(name,name, 9, 0, 360);     
    hist_CMPhi_t_bin2[i][0]->GetXaxis()->SetTitle("#phi_{CM} [deg]");
    hist_CMPhi_t_bin2[i][0]->GetYaxis()->SetTitle("counts");
    sprintf(name,"CMPhi_t_bin2_p1_%02d", i);
    hist_CMPhi_t_bin2[i][1] = new TH1F(name,name, 9, 0, 360);    
    hist_CMPhi_t_bin2[i][1]->GetXaxis()->SetTitle("#phi_{CM} [deg]");
    hist_CMPhi_t_bin2[i][1]->GetYaxis()->SetTitle("counts");
  }

  for(int i = 0; i < t_bins_bin3; i++){
    sprintf(name,"CMPhi_t_bin3_m1_%02d", i);
    hist_CMPhi_t_bin3[i][0] = new TH1F(name,name, 9, 0, 360);     
    hist_CMPhi_t_bin3[i][0]->GetXaxis()->SetTitle("#phi_{CM} [deg]");
    hist_CMPhi_t_bin3[i][0]->GetYaxis()->SetTitle("counts");
    sprintf(name,"CMPhi_t_bin3_p1_%02d", i);
    hist_CMPhi_t_bin3[i][1] = new TH1F(name,name, 9, 0, 360);    
    hist_CMPhi_t_bin3[i][1]->GetXaxis()->SetTitle("#phi_{CM} [deg]");
    hist_CMPhi_t_bin3[i][1]->GetYaxis()->SetTitle("counts");
  }

  for(int i = 0; i < t_bins_bin4; i++){
    sprintf(name,"CMPhi_t_bin4_m1_%02d", i);
    hist_CMPhi_t_bin4[i][0] = new TH1F(name,name, 9, 0, 360);     
    hist_CMPhi_t_bin4[i][0]->GetXaxis()->SetTitle("#phi_{CM} [deg]");
    hist_CMPhi_t_bin4[i][0]->GetYaxis()->SetTitle("counts");
    sprintf(name,"CMPhi_t_bin4_p1_%02d", i);
    hist_CMPhi_t_bin4[i][1] = new TH1F(name,name, 9, 0, 360);    
    hist_CMPhi_t_bin4[i][1]->GetXaxis()->SetTitle("#phi_{CM} [deg]");
    hist_CMPhi_t_bin4[i][1]->GetYaxis()->SetTitle("counts");
  }

  for(int i = 0; i < t_bins_bin5; i++){
    sprintf(name,"CMPhi_t_bin5_m1_%02d", i);
    hist_CMPhi_t_bin5[i][0] = new TH1F(name,name, 9, 0, 360);     
    hist_CMPhi_t_bin5[i][0]->GetXaxis()->SetTitle("#phi_{CM} [deg]");
    hist_CMPhi_t_bin5[i][0]->GetYaxis()->SetTitle("counts");
    sprintf(name,"CMPhi_t_bin5_p1_%02d", i);
    hist_CMPhi_t_bin5[i][1] = new TH1F(name,name, 9, 0, 360);    
    hist_CMPhi_t_bin5[i][1]->GetXaxis()->SetTitle("#phi_{CM} [deg]");
    hist_CMPhi_t_bin5[i][1]->GetYaxis()->SetTitle("counts");
  }


out->mkdir("invmass_bins");
out->cd ("invmass_bins");

  TH1F *hist_invmass_t_bin1[t_bins_bin1];
  TH1F *hist_invmass_t_bin2[t_bins_bin2];     
  TH1F *hist_invmass_t_bin3[t_bins_bin3];     
  TH1F *hist_invmass_t_bin4[t_bins_bin4];  
  TH1F *hist_invmass_t_bin5[t_bins_bin5];
  
  for(int i = 0; i < t_bins_bin1; i++){
    sprintf(name,"invmass_t_bin1_%02d", i);
    hist_invmass_t_bin1[i] = new TH1F(name,name, 200, 0.07, 0.19);     
    hist_invmass_t_bin1[i]->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV]");
    hist_invmass_t_bin1[i]->GetYaxis()->SetTitle("counts");
  }

  for(int i = 0; i < t_bins_bin2; i++){
    sprintf(name,"invmass_t_bin2_%02d", i);
    hist_invmass_t_bin2[i] = new TH1F(name,name, 200, 0.07, 0.19);     
    hist_invmass_t_bin2[i]->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV]");
    hist_invmass_t_bin2[i]->GetYaxis()->SetTitle("counts");
  }

  for(int i = 0; i < t_bins_bin3; i++){
    sprintf(name,"invmass_t_bin3_%02d", i);
    hist_invmass_t_bin3[i] = new TH1F(name,name, 200, 0.07, 0.19);     
    hist_invmass_t_bin3[i]->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV]");
    hist_invmass_t_bin3[i]->GetYaxis()->SetTitle("counts");
  }

  for(int i = 0; i < t_bins_bin4; i++){
    sprintf(name,"invmass_t_bin4_%02d", i);
    hist_invmass_t_bin4[i] = new TH1F(name,name, 200, 0.07, 0.19);     
    hist_invmass_t_bin4[i]->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV]");
    hist_invmass_t_bin4[i]->GetYaxis()->SetTitle("counts");
  }

  for(int i = 0; i < t_bins_bin5; i++){
    sprintf(name,"invmass_t_bin5_%02d", i);
    hist_invmass_t_bin5[i] = new TH1F(name,name, 200, 0.07, 0.19);     
    hist_invmass_t_bin5[i]->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV]");
    hist_invmass_t_bin5[i]->GetYaxis()->SetTitle("counts");
  }




/// ///////////////////////////////////////////////////////////////////////////////
///  define variables and vectors for the calculations:
/// /////////////////////////////////////////////////////////////  


TLorentzVector target; 
TLorentzVector beam; 
TLorentzVector ele; 
TLorentzVector proton; 
TLorentzVector pi0;

TLorentzVector proton_cand[BUFFER_prot];
TLorentzVector proton_excl[BUFFER_prot];
TLorentzVector gamma[BUFFER_gamma];
TLorentzVector pi0_cand[BUFFER_pi0];
TLorentzVector pi0_excl[BUFFER_pi0];

double prot_dcx1[BUFFER_prot];
double prot_dcy1[BUFFER_prot];
double prot_dcz1[BUFFER_prot];

double w, q2, x, y, s, tmin;
double t, cmphi, cmphi_shift, epsilon;


/// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///  start of the event loop     ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

cout << "Analysing Tree: " << inTree << endl;
cout << "Event Loop starting ... " << endl;
cout << endl;


for(Int_t k=0; k<chain.GetEntriesFast();k++){    

  chain.GetEntry(k);

  if(k % 100000 == 0){

      double events = chain.GetEntries();
      double percent = k/(events/100);
      
      printf("Analysing event number %i of %.00f (%.01f percent)\n", k, events, percent);
  }


/// ///////////////////////////////////////////////////////////////////
/// initialize variables and vectors:

  target.SetPxPyPzE(0,0,0,0.9384);
  beam.SetPxPyPzE(0,0,Ebeam,Ebeam);
  ele.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
  proton.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
  pi0.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

  t = 0; epsilon = 0; cmphi = 0; cmphi_shift = 0; 
  

  for(Int_t i = 0; i < BUFFER_pi0; i++){ 
    pi0_excl[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    pi0_cand[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
  }
  
  for(Int_t i = 0; i < BUFFER_gamma; i++){ 
    gamma[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
  }
  
  for(Int_t i = 0; i < BUFFER_prot; i++){ 
    proton_excl[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    proton_cand[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
  }
  

/// asign the 4 vectors of the particles:

  double ele_sec = -1;

  if(p4_ele_E->size() > 0){ ele.SetPxPyPzE(p4_ele_px->at(0),p4_ele_py->at(0),p4_ele_pz->at(0),p4_ele_E->at(0)); ele_sec = p4_ele_sec->at(0);}
 
  for(Int_t i = 0; i < BUFFER_gamma; i++){
    int phot_size = p4_phot_E->size();
    if(phot_size > i){
      if(p4_phot_det->at(i) == 2 && p4_phot_sec->at(i) != ele_sec){
        gamma[i].SetPxPyPzE(p4_phot_px->at(i),p4_phot_py->at(i),p4_phot_pz->at(i),p4_phot_E->at(i));
      }
    }
  }
  
  for(Int_t i = 0; i < BUFFER_prot; i++){
    int prot_size = p4_prot_E->size();
    if(prot_size > i){
      if(p4_prot_chi2pid->at(i) > prot_chi2_pid_cut_low && p4_prot_chi2pid->at(i) < prot_chi2_pid_cut_high && p4_prot_det->at(i) == 2){
        proton_cand[i].SetPxPyPzE(p4_prot_px->at(i),p4_prot_py->at(i),p4_prot_pz->at(i),p4_prot_E->at(i));
        prot_dcx1[i] = p4_prot_dcx1->at(i);
        prot_dcy1[i] = p4_prot_dcy1->at(i);
        prot_dcz1[i] = p4_prot_dcz1->at(i);
      }
    } 
  }


/// ///////////////////////////////////////////////////////////////////////////////////
/// determine proton sector

  double sec_prot[BUFFER];

  for(Int_t i = 0; i < BUFFER; i++){ 
     sec_prot[i] = 0;
  }

  for(Int_t i = 0; i < BUFFER; i++){ 
      if(proton_cand[i].Phi()*180/Pi > -45  && proton_cand[i].Phi()*180/Pi < 15)   sec_prot[i] = 1;
      if(proton_cand[i].Phi()*180/Pi > 15   && proton_cand[i].Phi()*180/Pi < 75)   sec_prot[i] = 2;
      if(proton_cand[i].Phi()*180/Pi > 75   && proton_cand[i].Phi()*180/Pi < 135)  sec_prot[i] = 3;
      if(proton_cand[i].Phi()*180/Pi > 135  && proton_cand[i].Phi()*180/Pi < 180)  sec_prot[i] = 4;
      if(proton_cand[i].Phi()*180/Pi > -180 && proton_cand[i].Phi()*180/Pi < -165) sec_prot[i] = 4;
      if(proton_cand[i].Phi()*180/Pi > -165 && proton_cand[i].Phi()*180/Pi < -115) sec_prot[i] = 5;
      if(proton_cand[i].Phi()*180/Pi > -115 && proton_cand[i].Phi()*180/Pi < -45)  sec_prot[i] = 6;
  }


/// ////////////////////////////////////////////////////////////////////////////////////
/// electron momentum correction:

  double fe = 1;
  if(outbending == false) fe = corr_inb_final(ele.Px(), ele.Py(), ele.Pz(), ele_sec, 0) + 1;
  if(outbending == true)  fe = corr_outb_final(ele.Px(), ele.Py(), ele.Pz(), ele_sec, 0) + 1;
  double dp = fe;
  
  if(correct == true) ele.SetPxPyPzE(dp * ele.Px(), dp * ele.Py(), dp * ele.Pz(), sqrt(pow(dp * ele.P(),2) + pow(0.000511,2)));


/// //////////////////////////////////////////////////////////////////////////////////////
/// proton corrections:


  for(Int_t i = 0; i < BUFFER_prot; i++){
  
    double pp = proton_cand[i].P();
    double psec = sec_prot[i];
    double dc1th = atan2(sqrt(prot_dcx1[i]*prot_dcx1[i] + prot_dcy1[i]*prot_dcy1[i]), prot_dcz1[i])*TMath::RadToDeg();
    bool lowband = dc1th < (-53.14680163254601 + 79.61307254040804*pow(pp-0.3, 0.05739232362022314));
    double fp = 1;
    
    if(outbending == false){
    double dploss = exp(-2.739 - 3.932*pp) + 0.002907;
    if(!lowband) dploss = exp(-1.2 - 4.228*pp) + 0.007502;
      double pfac = 0.5;
      if(psec==1) pfac = 0.8;
      else if(psec==2) pfac = 0.7;
      else if(psec==6) pfac = 1;
      fp = (pp + dploss + dploss*pfac)/pp;
    } 
    else{
      double dploss = exp(-2.739 - 3.932*pp) + 0.002907;
      if(!lowband) dploss = exp(-1.871 - 3.063*pp) + 0.007517;
      double pfac = 0;
      if(psec==2) pfac = 0.3;
      else if(psec==5) pfac = -0.5;
      fp = (pp + dploss + dploss*pfac)/pp;
    }
     
    if(correct == true) proton_cand[i].SetPxPyPzE(fp * proton_cand[i].Px(), fp * proton_cand[i].Py(), fp * proton_cand[i].Pz(), 
                                                  sqrt(pow(fp * proton_cand[i].P(),2) + pow(0.938272,2)));

  }
 


/// /////////////////////////////////////////////////////////////////////////////////////
/// reject particles which are outside the kienmatic limits: 

   
  if(ele.P() < electron_p_min || ele.Theta()*180/Pi < electron_theta_min || ele.Theta()*180/Pi > electron_theta_max || kin_W(ele) < 2 || kin_Q2(ele) < 2){ ele.SetPxPyPzE(0,0,0,0);}
  
  for(Int_t i = 0; i < BUFFER_gamma; i++){  
    if(gamma[i].P() < photon_p_min || gamma[i].Theta()*180/Pi < photon_theta_min || gamma[i].Theta()*180/Pi > photon_theta_max ){  
      pi0_cand[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    }
  }

  for(Int_t i = 0; i < BUFFER_prot; i++){  
    if(proton_cand[i].P() < proton_p_min || proton_cand[i].P() > proton_p_max || proton_cand[i].Theta()*180/Pi < proton_theta_min || proton_cand[i].Theta()*180/Pi > proton_theta_max ){  
      proton_cand[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    }
  }



/// //////////////////////////////////////////////////////////////////
/// find the correct two photon pair for eth pi0:

  for(Int_t i = 0; i < BUFFER; i++){ 
    if(gamma[i].P() > 0) hist_p_phot->Fill(gamma[i].P());
    if(gamma[i].P() > 0) hist_phi_phot->Fill(gamma[i].Phi()*180/Pi);
    if(gamma[i].P() > 0) hist_theta_phot->Fill(gamma[i].Theta()*180/Pi);
    if(gamma[i].P() > 0) hist_theta_vs_p_phot->Fill(gamma[i].P(), gamma[i].Theta()*180/Pi);
    if(gamma[i].P() > 0) hist_theta_vs_phi_phot->Fill(gamma[i].Phi()*180/Pi, gamma[i].Theta()*180/Pi);
  }

  TLorentzVector neutral_iter[28];  
 
  for(Int_t i = 0; i <= 27; i++){ 
    neutral_iter[i].SetPxPyPzE(0,0,0,0);
    pi0_cand[i].SetPxPyPzE(0,0,0,0);
  }

  /// //////////////////////////////

  double g_thr1 = 0.6;
  double g_thr2 = 0.6;

  /// /////////////////////////////

  /*
  for(Int_t i = 0; i < BUFFER; i++){ 
    if(alpha_p1p2(ele, gamma[i])*180/Pi < 8 || gamma[i].Theta()*180/Pi < 5){     // reject radiative photons
      gamma[i].SetPxPyPzE(0,0,0,0);
    }
  }
  */

  /// select all 28 pairs of the 8 detected gammas with the highest energy and calculate the opening angle:

  if(gamma[0].E() > g_thr1 && gamma[1].E() > g_thr2 && p4_phot_sec->at(0) == p4_phot_sec->at(1)){ neutral_iter[0] = gamma[0] + gamma[1];}
  if(gamma[0].E() > g_thr1 && gamma[2].E() > g_thr2 && p4_phot_sec->at(0) == p4_phot_sec->at(2)){ neutral_iter[1] = gamma[0] + gamma[2];}
  if(gamma[1].E() > g_thr1 && gamma[2].E() > g_thr2 && p4_phot_sec->at(1) == p4_phot_sec->at(2)){ neutral_iter[2] = gamma[1] + gamma[2];}
  if(gamma[0].E() > g_thr1 && gamma[3].E() > g_thr2 && p4_phot_sec->at(0) == p4_phot_sec->at(3)){ neutral_iter[3] = gamma[0] + gamma[3];}
  if(gamma[1].E() > g_thr1 && gamma[3].E() > g_thr2 && p4_phot_sec->at(1) == p4_phot_sec->at(3)){ neutral_iter[4] = gamma[1] + gamma[3];}
  if(gamma[2].E() > g_thr1 && gamma[3].E() > g_thr2 && p4_phot_sec->at(2) == p4_phot_sec->at(3)){ neutral_iter[5] = gamma[2] + gamma[3];}
  if(gamma[0].E() > g_thr1 && gamma[4].E() > g_thr2 && p4_phot_sec->at(0) == p4_phot_sec->at(4)){ neutral_iter[6] = gamma[0] + gamma[4];}
  if(gamma[1].E() > g_thr1 && gamma[4].E() > g_thr2 && p4_phot_sec->at(1) == p4_phot_sec->at(4)){ neutral_iter[7] = gamma[1] + gamma[4];}
  if(gamma[2].E() > g_thr1 && gamma[4].E() > g_thr2 && p4_phot_sec->at(2) == p4_phot_sec->at(4)){ neutral_iter[8] = gamma[2] + gamma[4];}
  if(gamma[3].E() > g_thr1 && gamma[4].E() > g_thr2 && p4_phot_sec->at(3) == p4_phot_sec->at(4)){ neutral_iter[9] = gamma[3] + gamma[4];}
  if(gamma[0].E() > g_thr1 && gamma[5].E() > g_thr2 && p4_phot_sec->at(0) == p4_phot_sec->at(5)){ neutral_iter[10] = gamma[0] + gamma[5];}
  if(gamma[1].E() > g_thr1 && gamma[5].E() > g_thr2 && p4_phot_sec->at(1) == p4_phot_sec->at(5)){ neutral_iter[11] = gamma[1] + gamma[5];}
  if(gamma[2].E() > g_thr1 && gamma[5].E() > g_thr2 && p4_phot_sec->at(2) == p4_phot_sec->at(5)){ neutral_iter[12] = gamma[2] + gamma[5];}
  if(gamma[3].E() > g_thr1 && gamma[5].E() > g_thr2 && p4_phot_sec->at(3) == p4_phot_sec->at(5)){ neutral_iter[13] = gamma[3] + gamma[5];}
  if(gamma[4].E() > g_thr1 && gamma[5].E() > g_thr2 && p4_phot_sec->at(4) == p4_phot_sec->at(5)){ neutral_iter[14] = gamma[4] + gamma[5];}
  if(gamma[0].E() > g_thr1 && gamma[6].E() > g_thr2 && p4_phot_sec->at(0) == p4_phot_sec->at(6)){ neutral_iter[15] = gamma[0] + gamma[6];}
  if(gamma[1].E() > g_thr1 && gamma[6].E() > g_thr2 && p4_phot_sec->at(1) == p4_phot_sec->at(6)){ neutral_iter[16] = gamma[1] + gamma[6];}
  if(gamma[2].E() > g_thr1 && gamma[6].E() > g_thr2 && p4_phot_sec->at(2) == p4_phot_sec->at(6)){ neutral_iter[17] = gamma[2] + gamma[6];}
  if(gamma[3].E() > g_thr1 && gamma[6].E() > g_thr2 && p4_phot_sec->at(3) == p4_phot_sec->at(6)){ neutral_iter[18] = gamma[3] + gamma[6];}
  if(gamma[4].E() > g_thr1 && gamma[6].E() > g_thr2 && p4_phot_sec->at(4) == p4_phot_sec->at(6)){ neutral_iter[19] = gamma[4] + gamma[6];}
  if(gamma[5].E() > g_thr1 && gamma[6].E() > g_thr2 && p4_phot_sec->at(5) == p4_phot_sec->at(6)){ neutral_iter[20] = gamma[5] + gamma[6];}
  if(gamma[0].E() > g_thr1 && gamma[7].E() > g_thr2 && p4_phot_sec->at(0) == p4_phot_sec->at(7)){ neutral_iter[21] = gamma[0] + gamma[7];}
  if(gamma[1].E() > g_thr1 && gamma[7].E() > g_thr2 && p4_phot_sec->at(1) == p4_phot_sec->at(7)){ neutral_iter[22] = gamma[1] + gamma[7];}
  if(gamma[2].E() > g_thr1 && gamma[7].E() > g_thr2 && p4_phot_sec->at(2) == p4_phot_sec->at(7)){ neutral_iter[23] = gamma[2] + gamma[7];}
  if(gamma[3].E() > g_thr1 && gamma[7].E() > g_thr2 && p4_phot_sec->at(3) == p4_phot_sec->at(7)){ neutral_iter[24] = gamma[3] + gamma[7];}
  if(gamma[4].E() > g_thr1 && gamma[7].E() > g_thr2 && p4_phot_sec->at(4) == p4_phot_sec->at(7)){ neutral_iter[25] = gamma[4] + gamma[7];}
  if(gamma[5].E() > g_thr1 && gamma[7].E() > g_thr2 && p4_phot_sec->at(5) == p4_phot_sec->at(7)){ neutral_iter[26] = gamma[5] + gamma[7];}
  if(gamma[6].E() > g_thr1 && gamma[7].E() > g_thr2 && p4_phot_sec->at(6) == p4_phot_sec->at(7)){ neutral_iter[27] = gamma[6] + gamma[7];}


  for(Int_t i = 0; i <= 27; i++){ 
    if(neutral_iter[i].M() > 0) hist_inv_mass->Fill(neutral_iter[i].M());
    if(neutral_iter[i].M2() != 0) hist_inv_mass2->Fill(neutral_iter[i].M2());
  }


  /// /////////////////////////////////////////////////////////////////////////////////////////////////
  /// select and sort the pi0s which will be used for the analyis

  for(Int_t i = 0; i <= 27; i++){ 
    if(neutral_iter[i].M() > pi0_min && neutral_iter[i].M() < pi0_max){
      pi0_cand[i].SetPxPyPzE(neutral_iter[i].Px(), neutral_iter[i].Py(), neutral_iter[i].Pz(), neutral_iter[i].E());
    }
    else{ pi0_cand[i].SetPxPyPzE(0.0, 0.0, 0.0, 0.0); }
  }


/// //////////////////////////////////////////////////////////////////////////////
/// fill event information:

    hist_helicity->Fill(helicity);
    hist_beam_charge->Fill(beam_charge);
    hist_ele_sec->Fill(ele_sec);


    for(Int_t i = 0; i < BUFFER_prot; i++){ 
      int prot_size = p4_prot_E->size();
      if(prot_size > i){
        if(p4_prot_det->at(i) > 0) hist_prot_det->Fill(p4_prot_det->at(i));
        if(p4_prot_chi2pid->at(i) != 0) hist_prot_chi2pid->Fill(p4_prot_chi2pid->at(i));
        if(proton_cand[i].P() > 0 && p4_prot_chi2pid->at(i) != 0) hist_prot_chi2pid_vs_p->Fill(proton_cand[i].P(), p4_prot_chi2pid->at(i));
      }
    } 


    for(Int_t i = 0; i < BUFFER_gamma; i++){
      int phot_size = p4_phot_E->size();
      if(phot_size > i){
        if(p4_phot_det->at(i) > 0) hist_phot_det->Fill(p4_phot_det->at(i));
        if(p4_phot_sec->at(i) > 0) hist_phot_sec->Fill( p4_phot_sec->at(i));
      }
    }


/// ///////////////////////////////////////////////////////////////////////////////
/// consider topology particles first:

  // count the particles

  int pi0_count = 0;
  for(Int_t i = 0; i < BUFFER_pi0; i++){ if(pi0_cand[i].P() > 0.0){pi0_count += 1;}}
  if(pi0_count > 0) hist_multiplicity_pi0->Fill(pi0_count);

  int prot_count = 0;
  for(Int_t i = 0; i < BUFFER_prot; i++){ if(proton_cand[i].P() > 0){ prot_count += 1;}}
  if(prot_count > 0) hist_multiplicity_proton->Fill(prot_count, 1);

  // fill raw particle properties
  
  if(ele.P() > 0){
    hist_p_electron->Fill(ele.P());
    hist_phi_electron->Fill(ele.Phi()*180/Pi);
    hist_theta_electron->Fill(ele.Theta()*180/Pi);
    hist_theta_vs_p_electron->Fill(ele.P(), ele.Theta()*180/Pi);
    hist_theta_vs_phi_electron->Fill(ele.Phi()*180/Pi, ele.Theta()*180/Pi);
  }
  for(Int_t i = 0; i < BUFFER_prot; i++){ 
    if(proton_cand[i].P() > 0){
      hist_p_proton->Fill(proton_cand[i].P());
      hist_phi_proton->Fill(proton_cand[i].Phi()*180/Pi);
      hist_theta_proton->Fill(proton_cand[i].Theta()*180/Pi);
      hist_theta_vs_p_proton->Fill(proton_cand[i].P(), proton_cand[i].Theta()*180/Pi);
      hist_theta_vs_phi_proton->Fill(proton_cand[i].Phi()*180/Pi, proton_cand[i].Theta()*180/Pi);
    }
  }
  for(Int_t i = 0; i < BUFFER_pi0; i++){ 
    if(pi0_cand[i].P() > 0){
      hist_p_pi0->Fill(pi0_cand[i].P());
      hist_phi_pi0->Fill(pi0_cand[i].Phi()*180/Pi);
      hist_theta_pi0->Fill(pi0_cand[i].Theta()*180/Pi);
      hist_theta_vs_p_pi0->Fill(pi0_cand[i].P(), pi0_cand[i].Theta()*180/Pi);
      hist_theta_vs_phi_pi0->Fill(pi0_cand[i].Phi()*180/Pi, pi0_cand[i].Theta()*180/Pi);
    }
  }

  // raw kinematics:

  w = 0; q2 = 0; x = 0; y = 0; s = 0; tmin = 0;

  if(ele.P() > 0){
    w  = kin_W(ele);
    q2 = kin_Q2(ele);
    x  = kin_x(ele);
    y  = kin_y(ele);
    epsilon = kin_epsilon(ele);
    
    s = (q2/x) - q2 + pow(0.938272, 2.0);
    tmin = -1. / (2 * s) * ( pow(s, 2.0) - s * (-q2 + 2 * pow(0.938272, 2.0))
                           + ( -q2 - pow(0.938272, 2.0)) * -1.0 * pow(0.938272, 2.0) 
                           - sqrt( pow(s, 2.0) + pow(0.938272, 4.0) 
                                   + pow(q2, 2.0) - 2.0 * s * (pow(0.938272, 2.0) - q2) + 2 * q2 * pow(0.938272, 2.0)
                                 ) * sqrt(pow(s, 2.0) + pow(0.938272, 4.0) - 2.0 * s * pow(0.938272, 2.0))
                         ); 
  }


/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// selection of exclusive events:
  
  int prot_count_excl = 0;
  int gamma_count_excl = 0;
  int pi0_count_excl   = 0;
  
  double t_raw[BUFFER_pi0];

  if(ele.P() > 0){  // check if an electron within the kienmatic limits exists

    for(Int_t i = 0; i < BUFFER_prot; i++){ 
      for(Int_t j = 0; j < BUFFER_pi0; j++){ 

        if(proton_cand[i].P() > 0 && pi0_cand[j].P() > 0){    // if all particles exist in the combination

           hist_missing_mass2_eppi0->Fill(kin_mismass2_X3(ele, proton_cand[i], pi0_cand[j]));
           hist_missing_mass2_epX->Fill(kin_mismass2_X2(ele, proton_cand[i]));
           hist_missing_mass_epi0X->Fill(kin_mismass_X2(ele, pi0_cand[j]));
           hist_missing_E_eppi0->Fill(kin_E_X3(ele, proton_cand[i], pi0_cand[j]));
           hist_missing_PT_eppi0->Fill(kin_PT_X3(ele, proton_cand[i], pi0_cand[j]));
           hist_missing_Px_eppi0->Fill(kin_Px_X3(ele, proton_cand[i], pi0_cand[j]));
           hist_missing_Py_eppi0->Fill(kin_Py_X3(ele, proton_cand[i], pi0_cand[j]));
           hist_missing_Pz_eppi0->Fill(kin_Pz_X3(ele, proton_cand[i], pi0_cand[j]));
           hist_missing_cone_pi0->Fill(kin_cone_pi0(ele, proton_cand[i], pi0_cand[j])*180/Pi);
           
           
           if(   abs(kin_Px_X3(ele, proton_cand[i], pi0_cand[j])) < 0.3
              && abs(kin_Py_X3(ele, proton_cand[i], pi0_cand[j])) < 0.3 
              && abs(kin_cone_pi0(ele, proton_cand[i], pi0_cand[j])*180/Pi) < 4
              && kin_Pz_X3(ele, proton_cand[i], pi0_cand[j]) > -0.5 && kin_Pz_X3(ele, proton_cand[i], pi0_cand[j]) < 0.9 
              && kin_mismass2_X2(ele, proton_cand[i]) > -0.3 && kin_mismass2_X2(ele, proton_cand[i]) < 0.4
             ){
              proton_excl[i].SetPxPyPzE(proton_cand[i].Px(), proton_cand[i].Py(), proton_cand[i].Pz(), proton_cand[i].E()); 
              pi0_excl[j].SetPxPyPzE(pi0_cand[j].Px(), pi0_cand[j].Py(), pi0_cand[j].Pz(), pi0_cand[j].E()); 
           }
         }
       }
     }
  }


  for(Int_t i = 0; i < BUFFER_prot; i++){ if(proton_excl[i].P() > 0){prot_count_excl += 1;}}
  for(Int_t i = 0; i < BUFFER_pi0; i++){ if(pi0_excl[i].P() > 0){pi0_count_excl += 1;}}



  if(ele.P() > 0 && prot_count_excl == 1){

    for(Int_t i = 0; i < BUFFER_prot; i++){ 
      if(proton_excl[i].P() > 0) proton.SetPxPyPzE(proton_excl[i].Px(), proton_excl[i].Py(), proton_excl[i].Pz(), proton_excl[i].E()); 
    }      
    for(Int_t i = 0; i < BUFFER_pi0; i++){ 
      if(pi0_excl[i].P() > 0) pi0.SetPxPyPzE(pi0_excl[i].Px(), pi0_excl[i].Py(), pi0_excl[i].Pz(), pi0_excl[i].E()); 
    }             
  }


/// ///////////////////////////////////////////////////////////////////////////////////////////
/// From here on we have exclsuive events:



  if(ele.P() > 0 && proton.P() > 0 && pi0.P() > 0 && (proton.Theta()*180/Pi < (44.11 - 6.625 * proton.P() + 1.438 * pow(proton.P(),2)))){  

     hist_excl_missing_mass2_eppi0->Fill(kin_mismass2_X3(ele, proton, pi0));
     hist_excl_missing_mass2_epX->Fill(kin_mismass2_X2(ele, proton));
     hist_excl_missing_mass_epi0X->Fill(kin_mismass_X2(ele, pi0));
     hist_excl_missing_E_eppi0->Fill(kin_E_X3(ele, proton, pi0));
     hist_excl_missing_PT_eppi0->Fill(kin_PT_X3(ele, proton, pi0));
     hist_excl_missing_Px_eppi0->Fill(kin_Px_X3(ele, proton, pi0));
     hist_excl_missing_Py_eppi0->Fill(kin_Py_X3(ele, proton, pi0));
     hist_excl_missing_Pz_eppi0->Fill(kin_Py_X3(ele, proton, pi0));
     hist_excl_missing_cone_pi0->Fill(kin_cone_pi0(ele, proton, pi0)*180/Pi);


    /// calculate kienmatics:

    t = kin_t(ele, proton);
    cmphi = kin_cmphi(ele, proton);
    
    //if(cmphi < 0){cmphi_shift= cmphi + 2*3.141592654;}
    //else{cmphi_shift = cmphi;}
    //cmphi = cmphi_shift;


  /// ///////////////////////////////////////////////////////////////////////////////////////////////////
  /// invariant mass distributions in the single bins:
  
    if(x < x1 && q2 > 2 && q2 < (y0 + (x-x0)/(x1-x0)*(y1-y0)) && x < x01){
      for(int i = 0; i < t_bins_bin1; i++){
        if(-t > t_bin1[i] && -t < t_bin1[i+1]){hist_invmass_t_bin1[i]->Fill(pi0.M());}
      }
    }

    if(x < x1 && q2 > 2 && q2 < (y0 + (x-x0)/(x1-x0)*(y1-y0)) && x > x01){
      for(int i = 0; i < t_bins_bin2; i++){
        if(-t > t_bin2[i] && -t < t_bin2[i+1]){hist_invmass_t_bin2[i]->Fill(pi0.M());}
      }
    }

    if(x < x1 && q2 > 2 && q2 > (y0 + (x-x0)/(x1-x0)*(y1-y0))){
      for(int i = 0; i < t_bins_bin3; i++){
        if(-t > t_bin3[i] && -t < t_bin3[i+1]){hist_invmass_t_bin3[i]->Fill(pi0.M());}
      }
    }

    if(x > x1 && q2 > 2 && q2 < (y4 + (x-x4)/(x5-x4)*(y5-y4))){
      for(int i = 0; i < t_bins_bin4; i++){
        if(-t > t_bin4[i] && -t < t_bin4[i+1]){hist_invmass_t_bin4[i]->Fill(pi0.M());}
      }
    }

    if(x > x1 && q2 > 2 && q2 > (y4 + (x-x4)/(x5-x4)*(y5-y4))){
      for(int i = 0; i < t_bins_bin5; i++){
        if(-t > t_bin5[i] && -t < t_bin5[i+1]){hist_invmass_t_bin5[i]->Fill(pi0.M());}
      }
    }


/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(abs(t) < 2.0){  
      hist_forward_selected_p_electron->Fill(ele.P());
      hist_forward_selected_phi_electron->Fill(ele.Phi()*180/Pi);
      hist_forward_selected_theta_electron->Fill(ele.Theta()*180/Pi);
      hist_forward_selected_theta_vs_p_electron->Fill(ele.P(), ele.Theta()*180/Pi);
      hist_forward_selected_theta_vs_phi_electron->Fill(ele.Phi()*180/Pi, ele.Theta()*180/Pi);
       
      hist_forward_selected_p_proton->Fill(proton.P());
      hist_forward_selected_phi_proton->Fill(proton.Phi()*180/Pi);
      hist_forward_selected_theta_proton->Fill(proton.Theta()*180/Pi);
      hist_forward_selected_theta_vs_p_proton->Fill(proton.P(), proton.Theta()*180/Pi);
      hist_forward_selected_theta_vs_phi_proton->Fill(proton.Phi()*180/Pi, proton.Theta()*180/Pi);

      hist_forward_selected_p_pi0->Fill(pi0.P());
      hist_forward_selected_phi_pi0->Fill(pi0.Phi()*180/Pi);
      hist_forward_selected_theta_pi0->Fill(pi0.Theta()*180/Pi);
      hist_forward_selected_theta_vs_p_pi0->Fill(pi0.P(), pi0.Theta()*180/Pi);
      hist_forward_selected_theta_vs_phi_pi0->Fill(pi0.Phi()*180/Pi, pi0.Theta()*180/Pi);

      if(w  > 0) hist_forward_selected_w->Fill(w);
      if(q2 > 0) hist_forward_selected_q2->Fill(q2);
      if(x  > 0) hist_forward_selected_x->Fill(x);
      if(t != 0)  hist_forward_selected_t->Fill(-t);
      if(cmphi != 0) hist_forward_selected_cmphi->Fill(cmphi*180/Pi);
      if(x > 0 && q2 > 0) hist_forward_selected_q2_vs_x->Fill(x, q2);

      if(x  > 0) hist_x_vs_p_electron_forward->Fill(ele.P(), x);
      if(x  > 0) hist_x_vs_theta_electron_forward->Fill(ele.Theta()*180/Pi, x);
      if(q2 > 0) hist_q2_vs_p_electron_forward->Fill(ele.P(), q2);
      if(q2 > 0) hist_q2_vs_theta_electron_forward->Fill(ele.Theta()*180/Pi, q2);
    }


   //// write output eventfile:
   if(print==true){ myfile <<eventnumber << "	" << x << "	" << q2 << "	" << w << "	" << -t <<  "	" << cmphi*180/Pi <<  "	" << pi0.M() <<  endl;}
   //////////////////////  

    /// ///////////////////////////////////////////////////////////////////////////////////////////////////
    /// beam asymmetry and kinematics in x, q2 and t bins
    ///

    if(x < x1 && q2 > 2 && q2 < (y0 + (x-x0)/(x1-x0)*(y1-y0)) && x < x01){
      hist_t_bin1->Fill(-t); 
      histbin_t_bin1->Fill(-t);
      hist_q2_bin1->Fill(q2);
      hist_x_bin1->Fill(x);
      hist_epsilon_bin1->Fill(epsilon);
      hist_q2_vs_x_bin1->Fill(x, q2);
      for(int i = 0; i < t_bins_bin1; i++){
        if(-t > t_bin1[i] && -t < t_bin1[i+1]){
          hist_q2_t_bin1[i]->Fill(q2);  
          hist_x_t_bin1[i]->Fill(x);  
          hist_t_t_bin1[i]->Fill(-t);  
          hist_epsilon_t_bin1[i]->Fill(epsilon); 
          if(helicity == +1) hist_CMPhi_t_bin1[i][1]->Fill(cmphi*180/Pi);
          if(helicity == -1) hist_CMPhi_t_bin1[i][0]->Fill(cmphi*180/Pi);
        }
      }
    }

    if(x < x1 && q2 > 2 && q2 < (y0 + (x-x0)/(x1-x0)*(y1-y0)) && x > x01){
      hist_t_bin2->Fill(-t); 
      histbin_t_bin2->Fill(-t);
      hist_q2_bin2->Fill(q2);
      hist_x_bin2->Fill(x);
      hist_epsilon_bin2->Fill(epsilon);
      hist_q2_vs_x_bin2->Fill(x, q2);
      for(int i = 0; i < t_bins_bin2; i++){
        if(-t > t_bin2[i] && -t < t_bin2[i+1]){
          hist_q2_t_bin2[i]->Fill(q2);  
          hist_x_t_bin2[i]->Fill(x);  
          hist_t_t_bin2[i]->Fill(-t);  
          hist_epsilon_t_bin2[i]->Fill(epsilon);
          if(helicity == +1) hist_CMPhi_t_bin2[i][1]->Fill(cmphi*180/Pi);
          if(helicity == -1) hist_CMPhi_t_bin2[i][0]->Fill(cmphi*180/Pi);
        }
      }
    }

    if(x < x1 && q2 > 2 && q2 > (y0 + (x-x0)/(x1-x0)*(y1-y0))){
      hist_t_bin3->Fill(-t);
      histbin_t_bin3->Fill(-t);
      hist_q2_bin3->Fill(q2);
      hist_x_bin3->Fill(x);
      hist_epsilon_bin3->Fill(epsilon);
      hist_q2_vs_x_bin3->Fill(x, q2);
      for(int i = 0; i < t_bins_bin3; i++){
        if(-t > t_bin3[i] && -t < t_bin3[i+1]){
          hist_q2_t_bin3[i]->Fill(q2);  
          hist_x_t_bin3[i]->Fill(x);  
          hist_t_t_bin3[i]->Fill(-t);  
          hist_epsilon_t_bin3[i]->Fill(epsilon);
          if(helicity == +1) hist_CMPhi_t_bin3[i][1]->Fill(cmphi*180/Pi);
          if(helicity == -1) hist_CMPhi_t_bin3[i][0]->Fill(cmphi*180/Pi);
        }
      }
    }

    if(x > x1 && q2 > 2 && q2 < (y4 + (x-x4)/(x5-x4)*(y5-y4))){
      hist_t_bin4->Fill(-t);
      histbin_t_bin4->Fill(-t);
      hist_q2_bin4->Fill(q2);
      hist_x_bin4->Fill(x);
      hist_epsilon_bin4->Fill(epsilon);
      hist_q2_vs_x_bin4->Fill(x, q2);
      for(int i = 0; i < t_bins_bin4; i++){
        if(-t > t_bin4[i] && -t < t_bin4[i+1]){
          hist_q2_t_bin4[i]->Fill(q2);  
          hist_x_t_bin4[i]->Fill(x);  
          hist_t_t_bin4[i]->Fill(-t);
          hist_epsilon_t_bin4[i]->Fill(epsilon);
          if(helicity == +1) hist_CMPhi_t_bin4[i][1]->Fill(cmphi*180/Pi);
          if(helicity == -1) hist_CMPhi_t_bin4[i][0]->Fill(cmphi*180/Pi);
        }
      }
    }

    if(x > x1 && q2 > 2 && q2 > (y4 + (x-x4)/(x5-x4)*(y5-y4))){
      hist_t_bin5->Fill(-t);
      histbin_t_bin5->Fill(-t);
      hist_q2_bin5->Fill(q2);
      hist_x_bin5->Fill(x);
      hist_epsilon_bin5->Fill(epsilon);
      hist_q2_vs_x_bin5->Fill(x, q2);
      for(int i = 0; i < t_bins_bin5; i++){
        if(-t > t_bin5[i] && -t < t_bin5[i+1]){
          hist_q2_t_bin5[i]->Fill(q2);  
          hist_x_t_bin5[i]->Fill(x);  
          hist_t_t_bin5[i]->Fill(-t);
          hist_epsilon_t_bin5[i]->Fill(epsilon);
          if(helicity == +1) hist_CMPhi_t_bin5[i][1]->Fill(cmphi*180/Pi);
          if(helicity == -1) hist_CMPhi_t_bin5[i][0]->Fill(cmphi*180/Pi);
        }
      }
    }


 // } // end of tmin cut

  }  // end of exclusivity condition
  

} // end of event loop
   
  
  /// /////////////////////////////////////////////////////////////

    cout << endl;
    cout << "Tree successfully analysed!" << endl;
    cout << "Writing the output file ... " << endl;
    out->Write(); // Saving Histograms
    cout << "Histograms saved in File: " << outputfile << endl;
    out->Close(); // Closing Output File
    cout << "... Completed!" << endl;

    return 1;
    
/// ///////////////////////////////////////////////////////////////////////////////
}   /// end of main
/// ///////////////////////////////////////////////////////////////////////////////



/// ///////////////////////////////////////////////////////////////////////////////


double kin_cmphi(TLorentzVector ele, TLorentzVector hadron){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  
  auto vqq = beam - ele;
  auto lnorm = vqq.Vect().Cross(beam.Vect());
  auto hnorm = hadron.Vect().Cross(vqq.Vect());
  double phistar = lnorm.Angle(hnorm);
  if( lnorm.Dot(hadron.Vect()) > 0 ) phistar = 2*3.141592654 - phistar;
  return phistar;
}


double kin_cmcostheta(TLorentzVector ele, TLorentzVector hadron){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);

  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;
  TVector3 boost = -fCM.BoostVector();
  TLorentzVector hadronClone(hadron); 
  hadronClone.Boost(boost);

  return TMath::Cos(hadronClone.Theta()); 
}


double kin_W(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;
  return fCM.M();
}

double kin_Q2(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  return -fGamma.M2();
}

double kin_x(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  return -fGamma.M2()/(2*0.938272*fGamma.E());
}

double kin_y(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  return fGamma.E()/Ebeam;
}

double kin_t(TLorentzVector ele, TLorentzVector hadron){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector transfer = hadron - target;
  return transfer.M2();
}

double kin_epsilon(TLorentzVector ele){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector fGamma = beam - ele;
  TLorentzVector fCM = fGamma + target;
  double W = fCM.M();
  double Q2 = -fGamma.M2();
  double x = -fGamma.M2()/(2*0.938272*fGamma.E());
  double gamma = 2*x*0.93827/sqrt(Q2);
  double eta = (W*W-0.93827*0.93827+Q2)/(2*0.93827*0.93827);
  double y = eta/Ebeam;
  double epsilon = (1-y-0.25*y*y*gamma*gamma)/(1-y+0.5*y*y+0.25*y*y*gamma*gamma);
  return epsilon; 
}

double kin_mismass_X2(TLorentzVector ele, TLorentzVector part1){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - part1;
  return sqrt(miss.E()*miss.E() - miss.Px()*miss.Px() - miss.Py()*miss.Py() - miss.Pz()*miss.Pz());
}
double kin_mismass2_X2(TLorentzVector ele, TLorentzVector part1){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - part1;
  return (miss.E()*miss.E() - miss.Px()*miss.Px() - miss.Py()*miss.Py() - miss.Pz()*miss.Pz());
}
double kin_mismass_X3(TLorentzVector ele, TLorentzVector part1, TLorentzVector part2){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - part1 - part2;
  return sqrt(miss.E()*miss.E() - miss.Px()*miss.Px() - miss.Py()*miss.Py() - miss.Pz()*miss.Pz());
}
double kin_mismass2_X3(TLorentzVector ele, TLorentzVector part1, TLorentzVector part2){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - part1 - part2;
  return (miss.E()*miss.E() - miss.Px()*miss.Px() - miss.Py()*miss.Py() - miss.Pz()*miss.Pz());
}

double kin_E_X3(TLorentzVector ele, TLorentzVector part1, TLorentzVector part2){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - part1 - part2;
  return miss.E();
}

double kin_PT_X3(TLorentzVector ele, TLorentzVector part1, TLorentzVector part2){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - part1 - part2;
  return miss.Pt();
}

double kin_Px_X3(TLorentzVector ele, TLorentzVector part1, TLorentzVector part2){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - part1 - part2;
  return miss.Px();
}

double kin_Py_X3(TLorentzVector ele, TLorentzVector part1, TLorentzVector part2){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - part1 - part2;
  return miss.Py();
}

double kin_Pz_X3(TLorentzVector ele, TLorentzVector part1, TLorentzVector part2){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - part1 - part2;
  return miss.Pz();
}

double kin_cone_pi0(TLorentzVector ele, TLorentzVector part1, TLorentzVector pi0){
  TLorentzVector beam(0,0,Ebeam,Ebeam);
  TLorentzVector target(0,0,0,0.93827);
  TLorentzVector miss   = beam + target - ele - part1;
  return pi0.Phi() - miss.Phi();
}

double alpha_p1p2(TLorentzVector particle1, TLorentzVector particle2){
  return particle1.Vect().Angle(particle2.Vect());
}


/// ///////////////////////////////////
/// momentum correction

double corr_inb_final(float px, float py, float pz, int sec, int ivec){

  int phibins = 3;

  double xx[] = { 0.0263375, 0.0158871, 0.0130852, 0.0252757, 0.0156601, 0.00984872, 0.0171495, 0.00359637, -0.0046115, 0.0189465, 0.0131816, 0.0262004, 
                  0.0116485, 0.0105681, 0.0149848, 0.0213057, 0.0112999, 0.0100216, -0.00366006, 0.00694866, 0.0197195, 0.00244435, 0.00681414, 0.0294068, 
                  0.00314739, 0.0136338, 0.0768753, 0.00375165, 0.00907457, 0.0486894, 0.000318094, -0.00480124, 0.0395545, 0.000653685, 0.0093174, 0.0822385,
                  0.0290127, 0.019353, 0.00619702, 0.0250068, 0.0127755, 0.00356361, 0.0258017, 0.011097, -0.00706104, 0.0314799, 0.017144, 0.00642617, 0.0322597, 
                  0.0214737, 0.0123113, 0.0294037, 0.0235748, 0.0127779};

  double pars[3][6][3];
  int ipar=0;

  for(int ivec=0;ivec<3;ivec++)
  for(int isec=0;isec<6;isec++)
  {
    if(ivec!=2) {
        double dp1=xx[ipar++], dp5=xx[ipar++], dp9=xx[ipar++];
        pars[ivec][isec][0] = (dp1 - 2*dp5 + dp9)/32.;
        pars[ivec][isec][1] = (-7*dp1)/16. + (5*dp5)/8. - (3*dp9)/16.;
        pars[ivec][isec][2] = (45*dp1)/32. - (9*dp5)/16. + (5*dp9)/32.;
    } else {
        double dp1=xx[ipar++], dp2=xx[ipar++], dp4=xx[ipar++];
        double a = (dp4 - 3*dp2 + 2*dp1) / 6.0;
        double b = dp2 - dp1 - 3*a;
        double c = dp1 - a -b;
        pars[ivec][isec][0] = a;
        pars[ivec][isec][1] = b;
        pars[ivec][isec][2] = c;
    }
  }

  double pp = sqrt(px*px + py*py + pz*pz);
  double a=pars[ivec][sec-1][0],
         b=pars[ivec][sec-1][1],
         c=pars[ivec][sec-1][2];
  double dp = a*pp*pp + b*pp + c; //pol2 corr func

  double fie = TMath::RadToDeg()*atan2(py,px);
  double phi = fie + (fie<0 && sec>1)*360 - (sec-1)*60;
  phi = phi - 30/pp;

  if(ivec == 0){
    if(sec == 1){
      dp = 0.45*b*(pp-9)+0.1*c;     
      if(phibins == 3){
        if(phi < -5) dp = -0.01*b*(pp-9)+1.1*c; //phi<-5
        if(phi < 5 && phi > -5) dp = 0.45*b*(pp-9)+0.1*c; //-5<phi<5
        if(phi > 5) dp = 1.4*b*(pp-9)-1.15*c; //phi>5
      }
    }
    if(sec == 2){
      dp = -0.15*b*(pp-8.0)+0.1*c;
      if(phibins == 3){
        if(phi < -5) dp = -0.7*b*(pp-8.0)+0.5*c; //phi<-5
        if(phi < 5 && phi > -5) dp = -0.15*b*(pp-8.0)+0.1*c; //-5<phi<5
        if(phi > 5) dp = 1.7*b*(pp-8.0)-0.55*c; //phi>5
      }
    }
    if(sec == 3){
      dp = 2.*b*(pp-5.4)-0.6*c;
      if(phibins == 3){
        if(phi < -5) dp = 2.75*b*(pp-5.4)-1.0*c; //phi<-5
        if(phi < 5 && phi > -5) dp = 2.*b*(pp-5.4)-0.6*c; //-5<phi<5
        if(phi > 5) dp = 1.25*b*(pp-5.4)+0.1*c; //phi>5
      }
    }
    if(sec == 4){
      dp = 0.25*b*(pp-9.25)+0.5*c;
      if(phibins == 3){
        if(phi < -5) dp = 0.25*b*(pp-9.25)+0.01*c; //phi<-5
        if(phi < 5 && phi > -5) dp = 0.25*b*(pp-9.25)+0.5*c; //-5<phi<5
        if(phi > 5) dp = 0.1*b*(pp-9.25)+1.1*c; //phi>5
      }
    }
    if(sec == 5){
      dp = 2.2*b*(pp-7.5);
      if(phibins == 3){
        if(phi < -5) dp = 2.2*b*(pp-7.5); //phi<-5
        if(phi < 5 && phi > -5) dp = 2.2*b*(pp-7.5); //-5<phi<5
        if(phi > 5) dp = 2.2*b*(pp-7.5); //phi>5
      }
    }
    if(sec == 6){
      dp = 0.5*b*(pp-7)-0.1*c;
      if(phibins == 3){
        if(phi < -5) dp = 0.95*b*(pp-7)+0.25*c; //phi<-5
        if(phi < 5 && phi > -5) dp = 0.5*b*(pp-7)-0.1*c; //-5<phi<5
        if(phi > 5) dp = 1.25*b*(pp-7)-0.7*c; //phi>5  
      }    
    }
  }

  return dp/pp;
};


double corr_outb_final(float px, float py, float pz, int sec, int ipart){    // 20.09.2021

  int phibins = 3;

  double xx[] = { 0.0219879, 0.00406117, 0.000287491, 0.0244179, 0.0169383, 0.000121762,
                  0.0209204, -0.000675913, -0.00874854, 0.025209, 0.0113607, -0.0104661,
                  0.0211029, 0.00524283, 0.0116993, 0.0242328, 0.00706621, -0.0185997,
                  0.0240847, 0.0054933, 0.00358604, 0.0264154, 0.0111607, -0.00691424,
                  0.0243936, -1.30348e-06, -0.0157793, 0.0222698, 0.0123583, -0.00728148,
                  0.0224922, 0.0200913, 0.0337443, 0.0262862, 0.0170036, -0.00152548};

  double pars[6][2][3];
  int ipar=0;
  for(int isec=0;isec<6;isec++)
  for(int ivec=0;ivec<2;ivec++) {
      if (ivec==0){
      double dp1=xx[ipar++], dp5=xx[ipar++], dp9=xx[ipar++];
      pars[isec][ivec][0] = (dp1 - 2*dp5 + dp9)/32.;
      pars[isec][ivec][1] = (-7*dp1)/16. + (5*dp5)/8. - (3*dp9)/16.;
      pars[isec][ivec][2] = (45*dp1)/32. - (9*dp5)/16. + (5*dp9)/32.;
    } else{
      double dp1=xx[ipar++], dp2=xx[ipar++], dp4=xx[ipar++];
      double a = (dp4 - 3*dp2 + 2*dp1) / 6.0;
      double b = dp2 - dp1 - 3*a;
      double c = dp1 - a -b;
      pars[isec][ivec][0] = a;
      pars[isec][ivec][1] = b;
      pars[isec][ivec][2] = c;
    }
  }

  double pp = sqrt(px*px + py*py + pz*pz);
  int ivec = ipart==1 ? 1 : 0;
  double a=pars[sec-1][ivec][0],
         b=pars[sec-1][ivec][1],
         c=pars[sec-1][ivec][2];
  double dp = a*pp*pp + b*pp + c;

  double fie = TMath::RadToDeg()*atan2(py,px);
  double phi = fie + (fie<0 && sec>1)*360 - (sec-1)*60;
  phi = phi - 30/pp;


  if(ivec == 0){
    if(sec == 1){
      dp = 0.5*b*pp + 1.7*c;
      if(phibins == 3){
        if(phi < -5) dp = 1.0*b*pp + 1.55*c; //phi<-5
        if(phi < 5 && phi > -5) dp = 1.0*b*pp + 2.4*c; //5<phi<-5
        if(phi > 5) dp = 0.5*b*pp + 2.5*c; //phi>5
      }
    }
    if(sec == 2){
      dp = b*pp + 2.5*c;
      if(phibins == 3){
        if(phi < -5) dp = 1.2*b*pp + 2.3*c; //phi<-5
        if(phi < 5 && phi > -5) dp = 0.9*b*pp + 2.25*c; //5<phi<-5
        if(phi > 5) dp = 0.8*b*pp + 2.8*c; //phi>5
      }
    }
    if(sec == 3){
      dp = 0.5*b*pp + 2.*c;
      if(phibins == 3){
        if(phi < -5) dp = 0.5*b*pp + 2.25*c; //phi<-5
        if(phi < 5 && phi > -5) dp = 0.5*b*pp + 2.*c; //5<phi<-5
        if(phi > 5) dp = 0.5*b*pp + 1.75*c; //phi>5
      }
    }
    if(sec == 4){
      dp = b*pp + 2.45*c;
      if(phibins == 3){
        if(phi < -5) dp = 0.8*b*pp + 2.3*c; //phi<-5
        if(phi < 5 && phi > -5) dp = b*pp + 2.4*c; //5<phi<-5
        if(phi > 5) dp = b*pp + 2.35*c; //phi>5
      }
    }
    if(sec == 5){
      dp = b*pp + 1.85*c;
      if(phibins == 3){
        if(phi < -5) dp = 1.2*b*pp + 2.0*c; //phi<-5
        if(phi < 5 && phi > -5) dp = b*pp + 1.65*c; //5<phi<-5
        if(phi > 5) dp = b*pp + 2.0*c; //phi>5
      }
    }
    if(sec == 6){
      dp = 0.01*b*pp + 1.75*c;
      if(phibins == 3){
        if(phi < -5) dp = 0.5*b*pp + 2.3*c; //phi<-5
        if(phi < 5 && phi > -5) dp = 0.9*b*pp + 2.5*c; //5<phi<-5
        if(phi > 5) dp = 0.01*b*pp + 2.0*c; //phi>5
      }
    }
  }

  return dp/pp;
};


