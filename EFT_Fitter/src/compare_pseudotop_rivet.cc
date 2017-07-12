#include <iostream>
#include <sstream>

#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"


using namespace std;

void makeComparisonPlots(string, string, string, string );


int main(int argc, const char * argv[]) {

    std::cout <<"comparing histos "<< std::endl;
    
    makeComparisonPlots("files/emu_ttbarsignalplustau_1M_noALLfilters.root","VisPseudoToppT", "files/pwhg_p8_1m.root", "CMS_dilepton_diff/t_pT");
    makeComparisonPlots("files/emu_ttbarsignalplustau_1M_noALLfilters.root","VisPseudoTTBarMass", "files/pwhg_p8_1m.root", "CMS_dilepton_diff/ttbar_mass");
    makeComparisonPlots("files/emu_ttbarsignalplustau_1M_noALLfilters.root","VisPseudoTTBarRapidity", "files/pwhg_p8_1m.root", "CMS_dilepton_diff/ttbar_y");
    makeComparisonPlots("files/emu_ttbarsignalplustau_1M_noALLfilters.root","VisPseudoTTBarpT", "files/pwhg_p8_1m.root", "CMS_dilepton_diff/ttbar_pT");
    makeComparisonPlots("files/emu_ttbarsignalplustau_1M_noALLfilters.root","VisPseudoTopRapidity", "files/pwhg_p8_1m.root", "CMS_dilepton_diff/t_y");

    
    makeComparisonPlots("files/emu_ttbarsignalplustau_1M_noALLfilters.root","VisGenToppT", "files/pwhg_p8_1m.root", "CMS_dilepton_diff/t_pT_parton");
    makeComparisonPlots("files/emu_ttbarsignalplustau_1M_noALLfilters.root","VisGenTTBarMass", "files/pwhg_p8_1m.root", "CMS_dilepton_diff/ttbar_mass_parton");
    makeComparisonPlots("files/emu_ttbarsignalplustau_1M_noALLfilters.root","VisGenTTBarRapidity", "files/pwhg_p8_1m.root", "CMS_dilepton_diff/ttbar_y_parton");
    makeComparisonPlots("files/emu_ttbarsignalplustau_1M_noALLfilters.root","VisGenTTBarpT", "files/pwhg_p8_1m.root", "CMS_dilepton_diff/ttbar_pT_parton");
    makeComparisonPlots("files/emu_ttbarsignalplustau_1M_noALLfilters.root","VisGenTopRapidity", "files/pwhg_p8_1m.root", "CMS_dilepton_diff/t_y_parton");
    
    

    makeComparisonPlots("files/emu_ttbarsignalplustau_1M_noALLfilters.root","VisPseudoLLBarDPhi", "files/pwhg_p8_1m.root", "CMS_dilepton_diff/ll_delphi_abs");
    makeComparisonPlots("files/emu_ttbarsignalplustau_1M_noALLfilters.root","VisPseudoLLBarMass", "files/pwhg_p8_1m.root", "CMS_dilepton_diff/ll_mass");
    makeComparisonPlots("files/emu_ttbarsignalplustau_1M_noALLfilters.root","VisPseudoLeptonpT", "files/pwhg_p8_1m.root", "CMS_dilepton_diff/l_pT");
  //  makeComparisonPlots("files/emu_ttbarsignalplustau_1M_noALLfilters.root","HypjetMulti", "files/pwhg_p8_1m.root", "CMS_dilepton_diff/Njets");
    makeComparisonPlots("files/emu_ttbarsignalplustau_1M_noALLfilters.root","VisPseudoTTBarDeltaPhi", "files/pwhg_p8_1m.root", "CMS_dilepton_diff/ttbar_delphi_abs");

        makeComparisonPlots("files/emu_ttbarsignalplustau_1M_noALLfilters.root","VisGenLeptonBjetMass", "files/pwhg_p8_1m.root", "CMS_dilepton_diff/lb_mass");
    
    
}



void makeComparisonPlots(string pt_filename, string pt_plotname, string rivet_filename, string rivet_plotname){

    
    std::cout <<" makeComparisonPlots "<< pt_filename << "  "  <<  pt_plotname  << "  "  << rivet_filename << "  " << rivet_plotname <<std::endl;

    TFile * f_pseudotop = new TFile(pt_filename.c_str());
    TFile * f_rivet = new TFile(rivet_filename.c_str());
    
    
    TH1F * pseudotop_histo = (TH1F*)f_pseudotop->Get(pt_plotname.c_str());
    TH1F * rivet_histo = (TH1F*)f_rivet->Get(rivet_plotname.c_str());
    
    TH1F *hnew;
    
    
    if (!f_pseudotop)     std::cout <<"cant find pseudotop file "<< std::endl;
    if (!f_rivet  )       std::cout <<"cant find rivet file "<< std::endl;
    if (!rivet_histo)     std::cout <<"cant find pseudotop file "<< std::endl;
    if (!pseudotop_histo) std::cout <<"cant find pseudotop histo "<< std::endl;
    
    if (pt_plotname == "VisPseudoToppT" ||  pt_plotname == "VisGenToppT"   ){
    
        std::cout <<"Making comparison of  "<< pt_plotname <<std::endl;
        double edges_pt[7] = { 0.0, 65.0, 124.5, 200.0, 290.0, 400.0, 550.00 };
        hnew = (TH1F*)pseudotop_histo->Rebin(6,"hnew",edges_pt);
        hnew->SetTitle("");
        hnew->SetXTitle("PT_{t}");
        hnew->SetYTitle("A.U.");


    }else if (pt_plotname == "VisPseudoTTBarMass" ||  pt_plotname == "VisGenTTBarMass" ){
        std::cout <<"Making comparison of  "<< pt_plotname <<std::endl;

        double edges_mass[8] = { 300.00, 380.0,  470.0, 620.0,  820.0,  1100.0,  1550.0, 2500.0};
        hnew = (TH1F*)pseudotop_histo->Rebin(7,"hnew",edges_mass);
        hnew->SetTitle("");
        hnew->SetXTitle("M_{tt}");
        hnew->SetYTitle("A.U.");

    
    }else if (pt_plotname == "VisPseudoTTBarRapidity"  ||  pt_plotname == "VisGenTTBarRapidity"  ){
        std::cout <<"Making comparison of  "<< pt_plotname <<std::endl;

        double edges_tty[9] = { -2.5, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.5};
        hnew = (TH1F*)pseudotop_histo->Rebin(8,"hnew",edges_tty);
        hnew->SetTitle("");
        hnew->SetXTitle("Y_{tt}");
        hnew->SetYTitle("A.U.");
        
    }else if (pt_plotname == "VisPseudoTTBarpT"  ||  pt_plotname == "VisGenTTBarpT"  ){
        std::cout <<"Making comparison of  "<< pt_plotname <<std::endl;
        double edges_ttpt[6] = { 0.0, 30.0, 80.0, 170.0, 300.0, 500.00 };
        hnew = (TH1F*)pseudotop_histo->Rebin(5,"hnew",edges_ttpt);
        hnew->SetTitle("");
        hnew->SetXTitle("PT_{tt}");
        hnew->SetYTitle("A.U.");

        
    }else if (pt_plotname == "VisPseudoTopRapidity"  ||  pt_plotname == "VisGenTopRapidity"){
        std::cout <<"Making comparison of  "<< pt_plotname <<std::endl;
        double edges_ty[9] = { -2.5, -1.6, -1.0, -0.5, 0.0, 0.5, 1.0, 1.6, 2.5};
        hnew = (TH1F*)pseudotop_histo->Rebin(8,"hnew",edges_ty);
        hnew->SetTitle("");
        hnew->SetXTitle("Y_{t}");
        hnew->SetYTitle("A.U.");

    }else if (pt_plotname == "VisPseudoLLBarDPhi"){
        std::cout <<"Making comparison of  "<< pt_plotname <<std::endl;
//double edges_lldphi[16] = { -3.14, -2.722, -2.304, -1.886,-1.486,-1.05, -0.632, -0.214, 0.204,0.622, 1.04, 1.458, 1.876, 2.294, 2.712 , 3.14};
        
        //double edges_lldphi[21] = {0.0, 0.16, 0.32, 0.48, 0.64, 0.80, 0.96, 1.12, 1.28, 1.44, 1.60, 1.76, 1.92, 2.08, 2.24, 2.40, 2.56, 2.72, 2.88, 3.04, 3.20};
        //double edges_lldphi[14] = {0.0, 0.24,0.48, 0.73, 0.97, 1.21, 1.45, 1.7, 1.94, 2.18, 2.42, 2.67, 2.91, 3.15};

        double edges_lldphi[7] = {0.0, 0.475, 0.975, 1.45, 1.94, 2.425, 2.91};

        hnew = (TH1F*)pseudotop_histo->Rebin(6,"hnew",edges_lldphi);
        hnew->SetTitle("");
        hnew->SetXTitle("DelPhi_{ll}");
        hnew->SetYTitle("A.U.");
        
        rivet_histo->Rebin(2);
        
    }else if (pt_plotname == "HypjetMulti"){
        std::cout <<"Making comparison of  "<< pt_plotname <<std::endl;
        double edges_njets[11] = { 0,1,2,3,4,5,6,7,8,9,10};
       // hnew = (TH1F*)pseudotop_histo->Rebin(15,"hnew",edges_njets);
        hnew->SetTitle("");
        hnew->SetXTitle("N_Jets");
        hnew->SetYTitle("A.U.");
        
    }
    
    else if (pt_plotname == "VisPseudoLLBarMass"){
        std::cout <<"Making comparison of  "<< pt_plotname <<std::endl;
        double edges_llmass[26] = { 0.0, 8.0, 16.0, 24.0, 32.0, 40.0, 48.0, 56.0, 64.0, 72.0, 80.0, 88.0, 96.0, 104.0, 112.0, 120.0, 128.0, 136.0, 144.0, 152.0, 160.0, 168.0, 176.0, 184.0, 192.0, 200.0};
        hnew = (TH1F*)pseudotop_histo->Rebin(25,"hnew",edges_llmass);
        hnew->SetTitle("");
        hnew->SetXTitle("M_ll");
        hnew->SetYTitle("A.U.");
        
    }  else if (pt_plotname == "VisPseudoLeptonpT"){
        std::cout <<"Making comparison of  "<< pt_plotname <<std::endl;
        double edges_lpt[31] = {0.0,10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0,160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0};
        hnew = (TH1F*)pseudotop_histo->Rebin(15,"hnew",edges_lpt);
        hnew->SetTitle("");
        hnew->SetXTitle("pT_l");
        hnew->SetYTitle("A.U.");
        
    }else if (pt_plotname == "VisPseudoTTBarDeltaPhi"){
        std::cout <<"Making comparison of  "<< pt_plotname <<std::endl;
        //double edges_ttdelphi[9] = {0.0, 0.3925 , 0.785, 1.1775, 1.57, 1.9625, 2.355, 2.7475, 3.14  };
        double edges_ttdelphi[5] = {0.0, 1.89,2.77, 3.05, 3.15 };
        hnew = (TH1F*)pseudotop_histo->Rebin(4,"hnew",edges_ttdelphi);
        hnew->SetTitle("");
        hnew->SetXTitle("delphi_tt");
        hnew->SetYTitle("A.U.");
        
    }
    
    else if (pt_plotname == "VisGenLeptonBjetMass"){
        std::cout <<"Making comparison of  "<< pt_plotname <<std::endl;
        double edges_mlb[31] = {0.0,10.0, 20.0, 30.0, 40.0, 50.0,  60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0,160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 260.0,  270.0, 280.0, 290.0, 300.0 };
        hnew = (TH1F*)pseudotop_histo->Rebin(30,"hnew",edges_mlb);
        hnew->SetTitle("");
        hnew->SetXTitle("mlb");
        hnew->SetYTitle("A.U.");
        
    }
    
    
    
    
     std::cout <<"making canvas"<< std::endl;

    
    TCanvas * c1  = new TCanvas();
    gStyle->SetOptStat(0);
    hnew->SetLineColor(2);
    hnew->SetLineWidth(2.5);
    hnew->SetMarkerStyle(22);
    hnew->SetMarkerColor(2);

    rivet_histo->SetLineWidth(2.5);


    
    hnew->DrawNormalized("E1p");
    rivet_histo->DrawNormalized("SAME");
    
    TLegend * leg = new TLegend(0.2,0.6,0.4,0.8);
    leg->AddEntry(hnew,"PseudoTop Producer","Epl");
    leg->AddEntry(rivet_histo,"RIVET","Epl");
    leg->Draw();
    
    
    string canvas_name = pt_plotname + ".png";
    
    c1->SaveAs(canvas_name.c_str());


}
