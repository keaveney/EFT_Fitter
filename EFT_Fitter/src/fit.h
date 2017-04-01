#include <iostream>
#include <sstream>
#include <string.h>
#include "TH1F.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLegendEntry.h"
#include <utility>

using namespace std;

bool debug = false;


class Fitter {
    
    public:
        Fitter(string, string);
        void run_extraction(int, float bins[], std::string, std::string, std::string, std::string ,bool );
        std::pair <TH1F*, vector<TH1F *>> initialise(std::string, std::string, std::string, int, float bins[], std::string, bool);
        void scan_couplings(std::string, std::pair <TH1F*, vector<TH1F *>>,  std::string mode  );
        void make_summary_plot(vector <TGraphErrors*>);
        void create_dummy_fiducial_measurement(double, double);

    private:
        double calculate_test_statistic(TH1F*, TH1F*);
};


Fitter* f_EFT;

vector <TGraphErrors*> scans;
vector <std::string> obs_names;

TH1F * dummy_fiducial_data;
vector<TH1F *> mc_histos_fiducial;


float bins_mtt[] = { 0.0, 380.0, 470.0, 620.0, 820.0, 1100.0, 1600.0 };
float bins_ptt[] = { 0.0, 65.0, 125.0, 200.0, 290.0, 400.0, 545.0 };
float bins_pttt[] = { 0.0, 30.0, 80.0, 170.0, 300.0, 500.0 };
float bins_ytt[] = { -2.5, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.5};
float bins_yt[] = { -2.5, -1.6, -1.0, -0.5, 0.0, 0.5, 1.0, 1.6, 2.5};



//scale basis histos according to cross sections (re-check these numbers)
//NLO SM Sigma
double sigma_sm_nlo = 674.0;
double sigma_CtG_neg2_nlo = 373.6;
double sigma_CtG_pos2_nlo = 975.0;
double k_factor = 831.76/674.1;
double BR = 0.045;
//acceptances (to be re-checked with higher stats)
double acc_sm_nlo = 0.3785;
double acc_CtG_neg2_nlo = 0.3899;
double acc_CtG_pos2_nlo = 0.3854;
double sigma_sm_fid = sigma_sm_nlo * k_factor * BR *acc_sm_nlo;
double sigma_CtG_neg2_fid = sigma_CtG_neg2_nlo * k_factor * BR * acc_CtG_neg2_nlo;
double sigma_CtG_pos2_fid = sigma_CtG_pos2_nlo * k_factor * BR * acc_CtG_pos2_nlo;



void Fitter::create_dummy_fiducial_measurement(double result, double rel_uncertainty){
    dummy_fiducial_data = new TH1F("","",1,0.0,1.0);
    dummy_fiducial_data->SetBinContent(1, result);
    dummy_fiducial_data->SetBinError(1, result*rel_uncertainty);
    TFile * f_dummy_fiducial =  new TFile("dummy_fiducial.root", "RECREATE");
    dummy_fiducial_data->Write();
    f_dummy_fiducial->Close();
}


std::pair <TH1F*, vector<TH1F *>> Fitter::initialise(std::string graphname_data, std::string filename_data, std::string histoname_pred, int nbins, float bins[], string mode, bool closure_test){
    
    
    TH1F* data_histo;
    vector<TH1F *> mc_histos;

    double running_total = 0.0;
    TH1F * h_CtG_pred;
    TH1F * h_CtG_pred_fiducial;
    double CtG_vals[11] = {-2.0,-1.6,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.6,2.0};
    double scaling =1.0, CtG;
    
    if (closure_test){
        
        TFile* f_data = new TFile("files/EFT_tt_rivet_small4.root");
        data_histo = (TH1F*)f_data->Get("CMS_dilepton_diff/ttbar_mass");
        if(debug)    cout << "Closure test\n";
        
    }else{
        if(debug)    cout << "data mode\n";
        if(debug)    cout << "Extracting data graph from " <<  filename_data.c_str()    <<"\n";
        
        TFile *f_data = new TFile(filename_data.c_str());
        TGraphAsymmErrors* g_data = (TGraphAsymmErrors*) f_data->Get(graphname_data.c_str() );
        
        data_histo = new TH1F("n","t", nbins, bins);
        double bin_centre, bin_height;
        
        for (int point = 0; point <= nbins; point ++){
            g_data->GetPoint(point, bin_centre, bin_height);
             cout << "looping on graph points, "<< point <<  "  "<<  bin_height  << " running_total " << running_total  <<"\n";

            data_histo->SetBinContent(point+1, bin_height);
            double bin_error = (g_data->GetErrorYhigh(point) + g_data->GetErrorYlow(point))/2.0 ; //hmmmmm need to set asymmetric errors in the histo or use TGraph throughout
            data_histo->SetBinError(point+1, bin_error);
            running_total = running_total + bin_height;
        }
    }
    
    if(debug)    cout << "Running total for data histo  =" <<  running_total  <<"\n";

    string filename_neg2 = "files/CTG_-2_nominal_2M.root";
    string filename_pos2 = "files/CTG_2_nominal_2M.root";
    string filename_0 = "files/CTG_0_nominal_2M.root";
    
    TFile * f_neg2 = new TFile(filename_neg2.c_str());
    TFile * f_0 = new TFile(filename_0.c_str());
    TFile * f_pos2 = new TFile(filename_pos2.c_str());
    
    TH1F * mc_histo_neg2 = (TH1F*)f_neg2->Get(histoname_pred.c_str());
    TH1F * mc_histo_0 = (TH1F*)f_0->Get(histoname_pred.c_str());
    TH1F * mc_histo_pos2 = (TH1F*)f_pos2->Get(histoname_pred.c_str());
    
    if (!mc_histo_neg2) cout << "mc histo: "<<  filename_neg2  <<" not found" << endl;
    
    
    double running_fid_xs_neg2 = 0.0;
    double running_fid_xs_0 = 0.0;
    double running_fid_xs_pos2 = 0.0;

    
    //divide bin heights by bin widths
    for (int bin = 1 ; bin <= nbins; bin++){
        double bin_width = bins[bin] - bins[bin-1];
        double bin_xsec_0  = mc_histo_0->GetBinContent(bin) / bin_width;
        double bin_xsec_neg2  = mc_histo_neg2->GetBinContent(bin) / bin_width;
        double bin_xsec_pos2  = mc_histo_pos2->GetBinContent(bin) / bin_width;
        running_fid_xs_neg2 = running_fid_xs_neg2 + mc_histo_neg2->GetBinContent(bin)  ;
        running_fid_xs_0 = running_fid_xs_0  + mc_histo_0->GetBinContent(bin) ;
        running_fid_xs_pos2 = running_fid_xs_pos2 + mc_histo_pos2->GetBinContent(bin);
        mc_histo_neg2->SetBinContent(bin, bin_xsec_neg2);
        mc_histo_0->SetBinContent(bin, bin_xsec_0);
        mc_histo_pos2->SetBinContent(bin, bin_xsec_pos2);

    }
    
    

       //now scale the basis histos to the right XS
        mc_histo_0->Scale(sigma_sm_fid/running_fid_xs_0);
        mc_histo_neg2->Scale(sigma_CtG_neg2_fid/running_fid_xs_neg2);
        mc_histo_pos2->Scale(sigma_CtG_pos2_fid/running_fid_xs_pos2);
    
    
    
       //make histo of pure CtG contribution
        TH1F * h_pure_Ctg = (TH1F*)mc_histo_pos2->Clone();
        h_pure_Ctg->Add(mc_histo_neg2, -1);
        h_pure_Ctg->Scale(0.25); //factor to get pure interference corresponding to CtG==1

       TFile * f_pred = new TFile("preds.root", "RECREATE");

        //loop over chosen CtG values and make specific prediction histos
        for (int scale = 0 ; scale < 11; scale++){
            
            CtG = CtG_vals[scale];
            h_CtG_pred = (TH1F*)mc_histo_0->Clone();
            h_CtG_pred->Add(h_pure_Ctg, CtG);
            
            //create histo containing prediction for fiducial cross section
            h_CtG_pred_fiducial = new TH1F("","", 1, 0.0,1.0);
            double fid_xsec_pred =0.0;
            
            for (int bin = 1 ; bin <= nbins; bin++){
                double bin_width = bins[bin] - bins[bin-1];
                fid_xsec_pred  = fid_xsec_pred  + ( h_CtG_pred->GetBinContent(bin) * bin_width );
            }
            
            h_CtG_pred_fiducial->SetBinContent(1,fid_xsec_pred);
            h_CtG_pred_fiducial->SetBinError(1,(fid_xsec_pred * 0.05));
            mc_histos_fiducial.push_back(h_CtG_pred_fiducial);
            
           // cout <<"VAR: "<< histoname_pred << " data integral =" <<  running_total  <<" Pred integral = "<<  running_fid_xs  <<"\n";

        
            //if running in norm or norm_fid, rescale the diff. predictions to match the measured data.
            if(mode == "norm_only" || mode == "norm_fid") {
                //if doing norm only  or norm + fid analysis, just set pred integral to data integral
                scaling = running_total/h_CtG_pred->Integral();
                h_CtG_pred->Scale(scaling);
            }
        
            mc_histos.push_back(h_CtG_pred);
        
            std::string pred_name_s = "CtG_" + to_string(CtG_vals[scale]);
            std::string pred_name_fid_s = "fid_CtG_" + to_string(CtG_vals[scale]);
            h_CtG_pred->SetName(pred_name_s.c_str());
            h_CtG_pred->Write();
            h_CtG_pred_fiducial->SetName(pred_name_fid_s.c_str());
            h_CtG_pred_fiducial->Write();
            }
    
    //now write basis histos to a file for debugging/validation
    std::stringstream ss;
    ss << histoname_pred;
    
    std::string segment;
    std::vector<std::string> seglist;
    
    while(std::getline(ss, segment, '/'))
    {
        seglist.push_back(segment);
    }
    
    std::string basis_name = "basis_histos/" + seglist[1] + "_basis.root";
    cout <<"basis name "<<   basis_name <<endl;
    
    TFile * f_basis = new TFile(basis_name.c_str(), "RECREATE");
    mc_histo_neg2->SetName("Neg2");
    mc_histo_neg2->Write();
    mc_histo_0->SetName("0");
    mc_histo_0->Write();
    mc_histo_pos2->SetName("Pos2");
    mc_histo_pos2->Write();
    h_pure_Ctg->SetName("pureInterference");
    h_pure_Ctg->Write();
    data_histo->Write();
    
    
    //f_basis->Close();
    // f_neg2->Close();
    // f_pos2->Close();
    // f_0->Close();
    
    std::pair <TH1F*, vector<TH1F *>> histos;
    histos = std::make_pair (data_histo, mc_histos);
    
    return histos;
}





