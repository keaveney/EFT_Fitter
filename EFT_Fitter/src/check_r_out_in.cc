#include <iostream>
#include <sstream>
#include <string.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"

#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TRandom3.h"

using namespace std;

int main(int argc, const char * argv[]){

  double tt_xsec = 831.76;
  double dy_xsec = (6025.2 + 22635.1);
  double int_lumi = 35900.0;

  double nominal_Ntt_incl = (tt_xsec * int_lumi);
  double nominal_Ndy_incl = (dy_xsec * int_lumi);
    
  double n_tt_out_ee     = 41775.3;
  double n_tt_in_ee      = 35829.0;
  double n_tt_out_mumu   = 54508.8;
  double n_tt_emu        = 148556;
  double n_tt_in_emu     = 39046;
  double n_dy_out_ee     = 51413.0;
  double n_dy_out_mumu   = 7584.1;
  double n_dy_emu        = 1228.2;
  double n_dy_in_ee      = 406151;
  double n_dy_in_mumu    = 890443.0;

  double eff_tt_out_ee      =   (n_tt_out_ee  / nominal_Ntt_incl );
  double eff_tt_in_ee       =   (n_tt_in_ee  / nominal_Ntt_incl );
  double eff_tt_out_mumu    =   (n_tt_out_mumu  / nominal_Ntt_incl );
  double eff_tt_emu         =   (n_tt_emu / nominal_Ntt_incl );
  double eff_tt_in_emu      =   (n_tt_in_emu / nominal_Ntt_incl );
  double eff_dy_out_ee      =   (n_dy_out_ee  / nominal_Ndy_incl );
  double eff_dy_out_mumu    =   (n_dy_out_mumu  / nominal_Ndy_incl );
  double eff_dy_emu         =   (n_dy_emu / nominal_Ndy_incl );
  double eff_dy_in_ee       =   (n_dy_in_ee  / nominal_Ndy_incl );
  double eff_dy_in_mumu     =   (n_dy_in_mumu  / nominal_Ndy_incl );
    
  double k_mumu = 1.44;
  double k_ee =  0.70;
  double R_out_in = 1.3;
    
  cout <<" eff_tt_ee =   "   << eff_tt_out_ee    << endl;
  cout <<" eff_tt_mumu = "   << eff_tt_out_mumu  << endl;
  cout <<" eff_tt_emu =  "   << eff_tt_emu   << endl;

  cout <<"******* RUNNING TOYS *********"<< endl;
    
  TRandom3 * gRand = new TRandom3();
  TH2F * h_trivial =  new TH2F("trivial","trivial", 20, 700.0, 1000.0, 20, 700.0, 1000.0);
  TH2F * h_real =  new TH2F("real","real", 40, 825.0, 840.0, 40, 825.0, 840.0);

  TH1F * h_real_pull =  new TH1F("real_pull","real_pull", 150, -1.0, 1.0);

  TH1F * h_dy_in_ee =  new TH1F("h_dy_in_ee","h_dy_in_ee", 30, 300000.0, 600000.0);
  TH1F * h_dy_out_ee =  new TH1F("h_dy_out_ee","h_dy_out_ee", 20, 0.0, 5000.0);
  TH1F * h_tt_out_ee =  new TH1F("h_tt_out_ee","h_tt_out_ee", 20, 5000.0, 50000.0);
  TH1F * h_tt_emu =  new TH1F("h_tt_emu","h_tt_emu", 20, 50000.0, 200000.0);

  TH2F * h_dy_in_ee_dy_out_ee =  new TH2F("h_dy_in_ee_dy_out_ee","h_dy_in_ee_dy_out_ee",30, 300000.0, 600000.0, 20, 0.0, 5000.0);
  TH2F * h_tt_emu_tt_out_ee =  new TH2F("h_tt_emu_tt_out_ee","h_tt_emu_tt_out_ee",20, 50000.0, 200000.0, 20, 5000.0, 50000.0);

  int n_toys = 10000;
    
  TGraph *g_real = new TGraph(n_toys);
  TGraph *g_trivial = new TGraph(n_toys);

  for (int toy = 0 ; toy < n_toys; toy++){
      
    double toy_ntt_incl = gRand->Gaus(tt_xsec*int_lumi,pow(tt_xsec*int_lumi,  0.5));
    double toy_ndy_incl = gRand->Gaus(dy_xsec*int_lumi,pow(dy_xsec*int_lumi,  0.5));

    double toy_n_tt_out_ee          =  toy_ntt_incl * eff_tt_out_ee;
    double toy_n_tt_in_ee           =  toy_ntt_incl * eff_tt_in_ee;
    double toy_n_tt_out_mumu        =  toy_ntt_incl * eff_tt_out_mumu;
    double toy_n_tt_emu             =  toy_ntt_incl * eff_tt_emu;
    double toy_n_tt_in_emu          =  toy_ntt_incl * eff_tt_in_emu;
    double toy_n_dy_out_ee          =  toy_ndy_incl * eff_dy_out_ee;
    double toy_n_dy_out_mumu        =  toy_ndy_incl * eff_dy_out_mumu;
    double toy_n_dy_emu             =  toy_ndy_incl * eff_dy_emu;
    double toy_n_dy_in_ee           =  toy_ndy_incl * eff_dy_in_ee;
    double toy_n_dy_in_mumu         =  toy_ndy_incl * eff_dy_in_mumu;
      

    h_dy_in_ee->Fill(toy_n_dy_in_ee);
    h_dy_out_ee->Fill(toy_n_dy_out_ee);
    h_tt_out_ee->Fill(toy_n_tt_out_ee);
    h_tt_emu->Fill(toy_n_tt_emu);

    h_dy_in_ee_dy_out_ee->Fill(toy_n_dy_in_ee,toy_n_dy_out_ee);
    h_tt_emu_tt_out_ee->Fill(toy_n_tt_emu,toy_n_tt_out_ee);

    // 1. Completely trivial check
    double toy_xsec_trivial_ee = (toy_n_tt_out_ee )/(eff_tt_out_ee);
    double toy_xsec_trivial_mumu = (toy_n_tt_out_mumu )/(eff_tt_out_mumu);
    double toy_xsec_trivial_emu = (toy_n_tt_emu )/(eff_tt_emu);
    double toy_xsec_trivial =  ( toy_xsec_trivial_ee  +  toy_xsec_trivial_emu ) / (2.0 * int_lumi);

    h_trivial->Fill(toy_ntt_incl,toy_xsec_trivial);
      
    cout <<"  "<< endl;
    cout <<" toy xsec = "<<  toy_ntt_incl/int_lumi   << endl;
    cout <<" trivial xsec = "<<  toy_xsec_trivial   << endl;

    //2. Simulate real analysis with R_out_in
    
    double r_out_in_ee = 0.119;
    double r_out_in_mumu = 0.124;
      
    double N_out_ee =  (r_out_in_ee) * ( toy_n_dy_in_ee + toy_n_tt_in_ee - (0.5 * toy_n_tt_in_emu * k_ee));
      
    cout <<" toy n out ee = "<<  toy_n_dy_out_ee   << "  estimated n out ee "<< N_out_ee  <<endl;

    double N_incl_ee = ((toy_n_tt_out_ee + toy_n_dy_out_ee)  - N_out_ee )  / (eff_tt_out_ee) ;
    double N_incl_emu = (toy_n_tt_emu )  / (eff_tt_emu) ;

    double toy_xsec_real =  ( N_incl_ee + N_incl_emu ) / (2.0 * int_lumi);
    double toy_xsec =  ( toy_ntt_incl  ) / (int_lumi);

    cout <<" toy xsec, DD xsec = "<<  toy_xsec  <<"  " <<toy_xsec_real   << endl;

    h_real_pull->Fill( 100*(toy_xsec_real - toy_xsec)/  (toy_xsec)  );
    h_real->Fill(toy_xsec,toy_xsec_real);
    g_real->SetPoint(toy,toy_xsec,toy_xsec_real);
    g_trivial->SetPoint(toy,toy_xsec,toy_xsec_trivial);

    }
    
    TFile * f = new TFile("toy_study_Routin.root", "RECREATE");
    
    h_dy_in_ee->Write();
    h_dy_out_ee->Write();
    h_tt_out_ee->Write();
    h_tt_emu->Write();
    h_dy_in_ee_dy_out_ee->Write();
    h_tt_emu_tt_out_ee->Write();
    g_real->Write();
    g_trivial->Write();
    h_real_pull->Write();

    TCanvas * c = new TCanvas();
    h_real->SetMarkerColor(kRed);
    h_real->Draw("AP");
    h_trivial->Draw("PSAME");
    c->Write();
    h_trivial->Write();
    h_real->Write();
    
    f->Close();
    
    return 0;
    
}



