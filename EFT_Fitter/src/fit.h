#include <iostream>
#include <sstream>
#include <string.h>
#include "TH1F.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TSpline.h"
#include "TLine.h"
#include "TLegendEntry.h"
#include <utility>
#include <tuple>


using namespace std;

bool debug = false;

const int n_preds = 51;

//double CtG_vals[n_preds] = {-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0};
double CtG_vals[n_preds];


class Fitter {
    
    public:
        Fitter(string, string);
        void run_extraction(int, float bins[], std::string, std::string, std::string, std::string ,bool, bool );
        std::tuple <TH1F*, vector<TH1F *>, vector<TGraphAsymmErrors *> > initialise(std::string, std::string, std::string, int, float bins[], std::string, bool, bool);
        void scan_couplings(std::string, std::tuple <TH1F*, vector<TH1F *> , vector<TGraphAsymmErrors *> >,  std::string mode , bool add_pwhg );
        void make_summary_plot(vector <TGraphErrors*>);
        void create_dummy_fiducial_measurement(double, double);

    private:
        double calculate_test_statistic(TH1F*, TH1F*, TGraphAsymmErrors*);
};


Fitter* f_EFT;

vector <TGraphErrors*> scans;
vector <std::string> obs_names;

TH1F * dummy_fiducial_data;
vector<TH1F *> mc_histos_fiducial;


float bins_mtt[] = { 345.0, 380.0, 470.0, 620.0, 820.0, 1100.0, 1600.0 };
float bins_ptt[] = { 0.0, 65.0, 125.0, 200.0, 290.0, 400.0, 545.0 };
float bins_pttt[] = { 0.0, 30.0, 80.0, 170.0, 300.0, 500.0 };
float bins_ytt[] = { -2.5, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.5};
float bins_yt[] = { -2.5, -1.6, -1.0, -0.5, 0.0, 0.5, 1.0, 1.6, 2.5};
//float bins_delphill[] = { -3.14, -2.72, -2.30, -1.88, -1.47, -1.05, -0.63, -0.21, 0.21, 0.63, 1.05, 1.47, 1.88, 2.30, 2.72, 3.14};
//float bins_delphill[] = { -3.14, -2.72, -2.30, -1.88, -1.47, -1.05, -0.63, -0.21, 0.21, 0.63, 1.05, 1.47, 1.88, 2.30, 2.72, 3.14 };
//float bins_delphill[] = {0.0, 0.16, 0.32, 0.48, 0.64, 0.80, 0.96, 1.12, 1.28, 1.44, 1.60, 1.76, 1.92, 2.08, 2.24, 2.40, 2.56, 2.72, 2.88, 3.04, 3.20};

float bins_delphill[] ={0.0, 0.24,0.48, 0.73, 0.97, 1.21, 1.45, 1.7, 1.94, 2.18, 2.42, 2.67, 2.91, 3.15};

float bins_delphibb[] = {0.0, 0.16, 0.32, 0.48, 0.64, 0.80, 0.96, 1.12, 1.28, 1.44, 1.60, 1.76, 1.92, 2.08, 2.24, 2.40, 2.56, 2.72, 2.88, 3.04, 3.20};

float bins_delphitt[] = {0.0, 1.89,2.77, 3.05, 3.15};
float bins_dephilb[] = { 0.0, 0.314, 0.628, 0.942, 1.256, 1.57 , 1.884, 2.198, 2.512, 2.826, 3.14};
float bins_costhetastar[] = {-1.0, -0.834, -0.668, -0.502, -0.336, -0.17, -0.004 , 0.162, 0.328 , 0.494, 0.66,  0.826, 1.0};

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


std::tuple <TH1F*, vector<TH1F *>, vector<TGraphAsymmErrors *> > Fitter::initialise(std::string graphname_data, std::string filename_data, std::string histoname_pred, int nbins, float bins[], string mode, bool closure_test, bool add_pwhg){
    
    TH1F* data_histo;
    vector<TH1F *> mc_histos;
    vector<TGraphAsymmErrors *> error_graphs;

    double running_total = 0.0;
    TH1F * h_CtG_pred;
    TH1F * h_CtG_pred_scaledown;
    TH1F * h_CtG_pred_scaleup;
    TH1F * h_CtG_pred_fiducial;
    


    double scaling =1.0, CtG;
    
    if (closure_test){
        
        TFile* f_data = new TFile("files/CtG_0_nominal_v2.root");
        data_histo = (TH1F*)f_data->Get(histoname_pred.c_str());
        
        if(debug)    cout << "Closure test\n";
        
        for (int bin = 1 ; bin <= nbins; bin++){
            double bin_width = bins[bin] - bins[bin-1];
            double bin_xsec_0  = data_histo->GetBinContent(bin) / bin_width;
            data_histo->SetBinContent(bin, bin_xsec_0);
            data_histo->SetBinError(bin, 0.5);

        }
        
        double data_int = data_histo->Integral();
        running_total = 1.0;
        data_histo->Scale(running_total/data_int);
    
        
    }else{
        if(debug)    cout << "data mode\n";
        if(debug)    cout << "Extracting data graph from " <<  filename_data.c_str()    <<"\n";
        
        TFile *f_data = new TFile(filename_data.c_str());
        TGraphAsymmErrors* g_data = (TGraphAsymmErrors*) f_data->Get(graphname_data.c_str() );
        
        data_histo = new TH1F("n","t", nbins, bins);
        double bin_centre, bin_height;
        
        for (int point = 0; point < nbins; point ++){
            g_data->GetPoint(point, bin_centre, bin_height);
            data_histo->SetBinContent(point+1, bin_height);
            double bin_error = (g_data->GetErrorYhigh(point) + g_data->GetErrorYlow(point))/2.0 ; //hmmmmm need to set asymmetric errors in the histo or use TGraph throughout
            data_histo->SetBinError(point+1, bin_error);
              cout << "looping on graph points, "<< point <<  "  "<<  bin_height  <<"  bin error "<<bin_error   << " running_total " << running_total  <<"\n";

            
            running_total = running_total + bin_height;
        }
//        data_histo->Scale(running_total/data_histo->Integral());

    }
    
    if(debug)    cout << "Running total for data histo  =" <<  running_total  <<"\n";
//CTG_pos2_1M.root
    
    
    //string filename_neg2 = "files/CTG_-2_nominal_v2.root";
    //string filename_pos2 = "files/CtG_2_nominal_v2.root";
    //string filename_0 = "files/CtG_0_nominal_v2.root";
    
    //nominal files
  //  string filename_neg2 = "files/CtG_-2_scaledown_v2.root";
   // string filename_pos2 = "files/CtG_2_scaledown_v2.root";
   // string filename_0 = "files/CtG_0_scaledown_v2.root";
    
    
    //nominal files
    string filename_neg2 = "files/CtG_-2_scaleup_v2.root";
    string filename_pos2 = "files/CtG_2_scaleup_v2.root";
    string filename_0 = "files/CtG_0_scaleup_v2.root";
    
    
    
    //scale down files
    string filename_neg2_scaledown = "files/CtG_-2_scaledown_v2.root";
    string filename_pos2_scaledown = "files/CtG_2_scaledown_v2.root";
    string filename_0_scaledown = "files/CtG_0_scaledown_v2.root";
    
    //scale up files
    string filename_neg2_scaleup = "files/CtG_-2_scaleup_v2.root";
    string filename_pos2_scaleup = "files/CtG_2_scaleup_v2.root";
    string filename_0_scaleup = "files/CtG_0_scaleup_v2.root";
    
    TFile * f_neg2 = new TFile(filename_neg2.c_str());
    TFile * f_0 = new TFile(filename_0.c_str());
    TFile * f_pos2 = new TFile(filename_pos2.c_str());

    TFile * f_neg2_scaledown = new TFile(filename_neg2_scaledown.c_str());
    TFile * f_0_scaledown = new TFile(filename_0_scaledown.c_str());
    TFile * f_pos2_scaledown = new TFile(filename_pos2_scaledown.c_str());
    
    TFile * f_neg2_scaleup = new TFile(filename_neg2_scaleup.c_str());
    TFile * f_0_scaleup = new TFile(filename_0_scaleup.c_str());
    TFile * f_pos2_scaleup = new TFile(filename_pos2_scaleup.c_str());
    
    TH1F * mc_histo_neg2 = (TH1F*)f_neg2->Get(histoname_pred.c_str());
    TH1F * mc_histo_0 = (TH1F*)f_0->Get(histoname_pred.c_str());
    TH1F * mc_histo_pos2 = (TH1F*)f_pos2->Get(histoname_pred.c_str());
    
    TH1F * mc_histo_neg2_scaledown = (TH1F*)f_neg2_scaledown->Get(histoname_pred.c_str());
    TH1F * mc_histo_0_scaledown = (TH1F*)f_0_scaledown->Get(histoname_pred.c_str());
    TH1F * mc_histo_pos2_scaledown = (TH1F*)f_pos2_scaledown->Get(histoname_pred.c_str());
    
    TH1F * mc_histo_neg2_scaleup = (TH1F*)f_neg2_scaleup->Get(histoname_pred.c_str());
    TH1F * mc_histo_0_scaleup = (TH1F*)f_0_scaleup->Get(histoname_pred.c_str());
    TH1F * mc_histo_pos2_scaleup = (TH1F*)f_pos2_scaleup->Get(histoname_pred.c_str());
    
    TH1F* mc_histo_pwhg;
    TFile * f_pwhg = new TFile("files/test_2M.root");
    if (add_pwhg) mc_histo_pwhg = (TH1F*)f_pwhg->Get(histoname_pred.c_str());

    
    if (!mc_histo_neg2) cout << "mc histo: "<<  filename_neg2  <<" not found" << endl;

    
    double running_fid_xs_neg2 = 0.0;
    double running_fid_xs_0 = 0.0;
    double running_fid_xs_pos2 = 0.0;
    
    double running_fid_xs_neg2_scaledown = 0.0;
    double running_fid_xs_0_scaledown = 0.0;
    double running_fid_xs_pos2_scaledown = 0.0;
    
    double running_fid_xs_neg2_scaleup = 0.0;
    double running_fid_xs_0_scaleup = 0.0;
    double running_fid_xs_pos2_scaleup = 0.0;
    
    double running_fid_xs_pwhg = 0.0;

    
    //divide bin heights by bin widths
    for (int bin = 1 ; bin <= nbins; bin++){
        double bin_width = bins[bin] - bins[bin-1];
       if (debug) cout << "bin width  "<<  bin_width  << endl;
        
        
        double bin_xsec_0  = mc_histo_0->GetBinContent(bin) / bin_width;
        double bin_xsec_neg2  = mc_histo_neg2->GetBinContent(bin) / bin_width;
        double bin_xsec_pos2  = mc_histo_pos2->GetBinContent(bin) / bin_width;
        
        double bin_xsec_0_scaledown  = mc_histo_0_scaledown->GetBinContent(bin) / bin_width;
        double bin_xsec_neg2_scaledown  = mc_histo_neg2_scaledown->GetBinContent(bin) / bin_width;
        double bin_xsec_pos2_scaledown  = mc_histo_pos2_scaledown->GetBinContent(bin) / bin_width;
        
        double bin_xsec_0_scaleup  = mc_histo_0_scaleup->GetBinContent(bin) / bin_width;
        double bin_xsec_neg2_scaleup  = mc_histo_neg2_scaleup->GetBinContent(bin) / bin_width;
        double bin_xsec_pos2_scaleup  = mc_histo_pos2_scaleup->GetBinContent(bin) / bin_width;
        double bin_xsec_pwhg;
        
        if (add_pwhg)  bin_xsec_pwhg  = mc_histo_pwhg->GetBinContent(bin) / bin_width;

        running_fid_xs_neg2 = running_fid_xs_neg2 + mc_histo_neg2->GetBinContent(bin)  ;
        running_fid_xs_0 = running_fid_xs_0  + mc_histo_0->GetBinContent(bin) ;
        running_fid_xs_pos2 = running_fid_xs_pos2 + mc_histo_pos2->GetBinContent(bin);
        
        running_fid_xs_neg2_scaledown = running_fid_xs_neg2_scaledown + mc_histo_neg2_scaledown->GetBinContent(bin)  ;
        running_fid_xs_0_scaledown = running_fid_xs_0_scaledown  + mc_histo_0_scaledown->GetBinContent(bin) ;
        running_fid_xs_pos2_scaledown = running_fid_xs_pos2_scaledown + mc_histo_pos2_scaledown->GetBinContent(bin);
        
        running_fid_xs_neg2_scaleup = running_fid_xs_neg2_scaleup + mc_histo_neg2_scaleup->GetBinContent(bin)  ;
        running_fid_xs_0_scaleup = running_fid_xs_0_scaleup  + mc_histo_0_scaleup->GetBinContent(bin) ;
        running_fid_xs_pos2_scaleup = running_fid_xs_pos2_scaleup + mc_histo_pos2_scaleup->GetBinContent(bin);
        
        
        if (add_pwhg)   running_fid_xs_pwhg = running_fid_xs_pwhg + mc_histo_pwhg->GetBinContent(bin);

        mc_histo_neg2->SetBinContent(bin, bin_xsec_neg2);
        mc_histo_0->SetBinContent(bin, bin_xsec_0);
        mc_histo_pos2->SetBinContent(bin, bin_xsec_pos2);
        
        mc_histo_neg2_scaledown->SetBinContent(bin, bin_xsec_neg2_scaledown);
        mc_histo_0_scaledown->SetBinContent(bin, bin_xsec_0_scaledown);
        mc_histo_pos2_scaledown->SetBinContent(bin, bin_xsec_pos2_scaledown);
        
        mc_histo_neg2_scaleup->SetBinContent(bin, bin_xsec_neg2_scaleup);
        mc_histo_0_scaleup->SetBinContent(bin, bin_xsec_0_scaleup);
        mc_histo_pos2_scaleup->SetBinContent(bin, bin_xsec_pos2_scaleup);
        

        
      if (add_pwhg)  mc_histo_pwhg->SetBinContent(bin, bin_xsec_pwhg);
    }
    
       //now scale the basis histos to the right XS
       //need to check the XS for scale variations, right now assuming the
       //dependence of XS on CtG doesnt depend on scale
        mc_histo_0->Scale(sigma_sm_fid/running_fid_xs_0);
        mc_histo_neg2->Scale(sigma_CtG_neg2_fid/running_fid_xs_neg2);
        mc_histo_pos2->Scale(sigma_CtG_pos2_fid/running_fid_xs_pos2);
    
    
        mc_histo_0_scaledown->Scale(sigma_sm_fid/running_fid_xs_0_scaledown);
        mc_histo_neg2_scaledown->Scale(sigma_CtG_neg2_fid/running_fid_xs_neg2_scaledown);
        mc_histo_pos2_scaledown->Scale(sigma_CtG_pos2_fid/running_fid_xs_pos2_scaledown);
    
        mc_histo_0_scaleup->Scale(sigma_sm_fid/running_fid_xs_0_scaleup);
        mc_histo_neg2_scaleup->Scale(sigma_CtG_neg2_fid/running_fid_xs_neg2_scaleup);
        mc_histo_pos2_scaleup->Scale(sigma_CtG_pos2_fid/running_fid_xs_pos2_scaleup);
    
    
     if (add_pwhg)    mc_histo_pwhg->Scale(sigma_sm_fid/running_fid_xs_pwhg);
    
    
       //make histo of pure CtG contribution
        TH1F * h_pure_Ctg = (TH1F*)mc_histo_pos2->Clone();
        TH1F * h_pure_Ctg_scaledown = (TH1F*)mc_histo_pos2_scaledown->Clone();
        TH1F * h_pure_Ctg_scaleup = (TH1F*)mc_histo_pos2_scaleup->Clone();
    
        h_pure_Ctg->Add(mc_histo_neg2, -1);
        h_pure_Ctg->Scale(0.25); //factor to get pure interference corresponding to CtG==1
        h_pure_Ctg->SetLineColor(kBlue); // the 0.25 comes from the fact that CtG = +/- 2 in the samples
        h_pure_Ctg->SetLineStyle(2);
    
        h_pure_Ctg_scaledown->Add(mc_histo_neg2_scaledown, -1);
        h_pure_Ctg_scaledown->Scale(0.25);
        h_pure_Ctg_scaledown->SetLineColor(kBlue);
        h_pure_Ctg_scaledown->SetLineStyle(2);
    
        h_pure_Ctg_scaleup->Add(mc_histo_neg2_scaleup, -1);
        h_pure_Ctg_scaleup->Scale(0.25);
        h_pure_Ctg_scaleup->SetLineColor(kBlue);
        h_pure_Ctg_scaleup->SetLineStyle(2);

        mc_histo_neg2->SetLineColor(kBlack);
        mc_histo_0->SetLineColor(kBlack);
        mc_histo_pos2->SetLineColor(kBlack);
    
        mc_histo_neg2_scaledown->SetLineColor(kGreen);
        mc_histo_0_scaledown->SetLineColor(kGreen);
        mc_histo_pos2_scaledown->SetLineColor(kGreen);
    
        mc_histo_neg2_scaleup->SetLineColor(kRed);
        mc_histo_0_scaleup->SetLineColor(kRed);
        mc_histo_pos2_scaleup->SetLineColor(kRed);

    
        TFile * f_pred = new TFile("preds.root", "RECREATE");
        TCanvas * c_discrimination =  new TCanvas();
       // h_pure_Ctg->DrawNormalized();
       // h_pure_Ctg_scaledown->DrawNormalized();
       // h_pure_Ctg_scaleup->DrawNormalized();

        mc_histo_0->DrawNormalized();
        mc_histo_neg2->DrawNormalized("SAME");
        mc_histo_pos2->DrawNormalized("SAME");
    
        mc_histo_0_scaledown->DrawNormalized("SAME");
        mc_histo_neg2_scaledown->DrawNormalized("SAME");
        mc_histo_pos2_scaledown->DrawNormalized("SAME");
    
        mc_histo_0_scaleup->DrawNormalized("SAME");
        mc_histo_neg2_scaleup->DrawNormalized("SAME");
        mc_histo_pos2_scaleup->DrawNormalized("SAME");


        //loop over chosen CtG values and make specific prediction histos
        for (int scale = 0 ; scale < n_preds; scale++){
            
            CtG = CtG_vals[scale];
            h_CtG_pred = (TH1F*)mc_histo_0->Clone();
            h_CtG_pred_scaledown = (TH1F*)mc_histo_0_scaledown->Clone();
            h_CtG_pred_scaleup   = (TH1F*)mc_histo_0_scaleup->Clone();

            h_CtG_pred->Add(h_pure_Ctg, CtG);
            h_CtG_pred_scaledown->Add(h_pure_Ctg_scaledown, CtG);
            h_CtG_pred_scaleup->Add(h_pure_Ctg_scaleup, CtG);
            
            
            //create histo containing prediction for fiducial cross section
            h_CtG_pred_fiducial = new TH1F("","", 1, 0.0,1.0);
            double fid_xsec_pred =0.0;
            
            TGraphAsymmErrors * gr_errors = new TGraphAsymmErrors(nbins);
            
            for (int bin = 1 ; bin <= nbins; bin++){
                double bin_width = bins[bin] - bins[bin-1];
                fid_xsec_pred  = fid_xsec_pred  + ( h_CtG_pred->GetBinContent(bin) * bin_width );
                h_CtG_pred->SetBinError(bin, h_CtG_pred->GetBinContent(bin) * 0.02);

            }
            
            h_CtG_pred_fiducial->SetBinContent(1,fid_xsec_pred);
            h_CtG_pred_fiducial->SetBinError(1,(fid_xsec_pred * 0.10));
            mc_histos_fiducial.push_back(h_CtG_pred_fiducial);
            
            
           if (debug) cout <<"BEFORE SCALING: "<< histoname_pred << " data integral =" <<  running_total  <<" Pred integral = "<<  h_CtG_pred->Integral()  <<"\n";

        
            //if running in norm or norm_fid, rescale the diff. predictions to match the measured data.
            if(mode == "norm_only" || mode == "norm_fid") {
                //if doing norm only  or norm + fid analysis, just set pred integral to data integral
                scaling = running_total/h_CtG_pred->Integral();
                h_CtG_pred->Scale(scaling);
                scaling = running_total/h_CtG_pred_scaledown->Integral();
                h_CtG_pred_scaledown->Scale(scaling);
                scaling = running_total/h_CtG_pred_scaleup->Integral();
                h_CtG_pred_scaleup->Scale(scaling);
            }
            
           if (debug) cout <<"AFTER SCALING: "<< histoname_pred << " data integral =" <<  running_total  <<" Pred integral = "<<  h_CtG_pred->Integral()  <<"\n";

            
            //fill error graph
            
            //NEXT: get errors from scale up/down plots instead of arbitary values...
            
            
            for ( int bin = 1 ; bin <= nbins; bin++){
                 double bin_width = bins[bin] - bins[bin-1];
                 double bin_centre = bins[bin-1] + (( bins[bin] - bins[bin-1] ) /2.0);
                gr_errors->SetPoint(bin-1,    bin_centre ,  h_CtG_pred->GetBinContent(bin) );
                
                //unflipped
                if (h_CtG_pred_scaleup->GetBinContent(bin) > h_CtG_pred->GetBinContent(bin) &&  h_CtG_pred->GetBinContent(bin) > h_CtG_pred_scaledown->GetBinContent(bin) ){
                    gr_errors->SetPointEYhigh(bin-1,3.0* ( h_CtG_pred_scaleup->GetBinContent(bin) - h_CtG_pred->GetBinContent(bin) ) );
                    gr_errors->SetPointEYlow(bin-1,3.0* (h_CtG_pred->GetBinContent(bin) - h_CtG_pred_scaledown->GetBinContent(bin) ) );
                }
                //flipped
                else if(h_CtG_pred_scaleup->GetBinContent(bin) < h_CtG_pred->GetBinContent(bin) &&  h_CtG_pred->GetBinContent(bin) < h_CtG_pred_scaledown->GetBinContent(bin)){
                    gr_errors->SetPointEYhigh(bin-1, 3.0*( h_CtG_pred_scaledown->GetBinContent(bin) - h_CtG_pred->GetBinContent(bin) ) );
                    gr_errors->SetPointEYlow(bin-1, 3.0*(h_CtG_pred->GetBinContent(bin) - h_CtG_pred_scaleup->GetBinContent(bin) ) );
                
                }
                //both higher
                else if(h_CtG_pred_scaleup->GetBinContent(bin) > h_CtG_pred->GetBinContent(bin) && h_CtG_pred_scaledown->GetBinContent(bin) > h_CtG_pred->GetBinContent(bin) ){
                   
                    if (h_CtG_pred_scaleup->GetBinContent(bin) > h_CtG_pred_scaledown->GetBinContent(bin) ){
                        gr_errors->SetPointEYhigh(bin-1, 3.0*( h_CtG_pred_scaleup->GetBinContent(bin) - h_CtG_pred->GetBinContent(bin) ) );
                    }else{
                        gr_errors->SetPointEYhigh(bin-1,3.0* ( h_CtG_pred_scaledown->GetBinContent(bin) - h_CtG_pred->GetBinContent(bin) ) );
                    }
                    
                }            //both lower
                else if(h_CtG_pred_scaleup->GetBinContent(bin) < h_CtG_pred->GetBinContent(bin) && h_CtG_pred_scaledown->GetBinContent(bin) < h_CtG_pred->GetBinContent(bin) ){
                    
                    if (h_CtG_pred_scaleup->GetBinContent(bin) > h_CtG_pred_scaledown->GetBinContent(bin) ){
                        gr_errors->SetPointEYlow(bin-1, 3.0*( h_CtG_pred->GetBinContent(bin) - h_CtG_pred_scaledown->GetBinContent(bin) ) );
                    }else{
                        gr_errors->SetPointEYlow(bin-1,3.0* ( h_CtG_pred->GetBinContent(bin) - h_CtG_pred_scaleup->GetBinContent(bin) ) );
                    }
                    
                }
                
                gr_errors->SetPointEXlow(bin-1, bin_width/2.0  );
                gr_errors->SetPointEXhigh(bin-1, bin_width/2.0 );
            }
            
            
            //TODO: update to new rivet files,
            // check why scale uncertainty is smaller for Ctg preds that nominal (if still the case)
            // make global variable for N CtG vals to avoid tedious changing
            // add dummy covariance matrices
            // keep data in graph rather than converting to TH1
            // start on slides
            
        
            mc_histos.push_back(h_CtG_pred);
            error_graphs.push_back(gr_errors);

        
            std::string pred_name_s = "CtG_" + to_string(CtG_vals[scale]);
            std::string pred_name_fid_s = "fid_CtG_" + to_string(CtG_vals[scale]);
            h_CtG_pred->SetName(pred_name_s.c_str());
            h_CtG_pred->Write();
            h_CtG_pred_fiducial->SetName(pred_name_fid_s.c_str());
            h_CtG_pred_fiducial->Write();
            }
    
    if(mode == "norm_only" || mode == "norm_fid") {
        //if doing norm only  or norm + fid analysis, just set pred integral to data integral
        if (add_pwhg) {scaling = running_total/mc_histo_pwhg->Integral();
        mc_histo_pwhg->Scale(scaling);
        }
    }
    
    if (add_pwhg) mc_histos.push_back(mc_histo_pwhg);

    
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
    std::string basis_canvas_name = "basis_histos/" + seglist[1] + "_canvas.pdf";
    std::string discrimination_canvas_name = "basis_histos/" + seglist[1] + "discrim_canvas.pdf";
    c_discrimination->SaveAs(discrimination_canvas_name.c_str());
    
    
    cout <<"basis name "<<   basis_name <<endl;
    
    TCanvas * c_basis =  new TCanvas();
    TPad *pad1 = new TPad("pad1","pad1",0,0.45,1,1);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    
    mc_histo_neg2->SetLineColor(kGreen);
    mc_histo_neg2->SetMarkerColor(kGreen);
  //  mc_histo_neg2->SetFillColor(kGreen);

    mc_histo_pos2->SetLineColor(kRed);
    mc_histo_pos2->SetMarkerColor(kRed);
//    mc_histo_pos2->SetFillColor(kRed);

    mc_histo_0->SetMarkerColor(kBlack);
    mc_histo_0->SetLineColor(kBlack);
 //   mc_histo_0->SetFillColor(kBlack);
    
    mc_histo_0->GetYaxis()->SetRangeUser(0.0, 5.0);


    mc_histo_0->Draw("E1p");
    mc_histo_pos2->Draw("E1pSAME");
    mc_histo_neg2->Draw("E1pSAME");
    
    
  // pad1->SetLogy();
    c_basis->cd();
    
    TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.45);
    pad2->SetTopMargin(0);
    pad2->Draw();
    pad2->cd();
  //  pad2->SetBottomMargin(0.3);
   // pad2->SetTopMargin(0.0);
    
    TH1F* mc_temp_neg = (TH1F*)  mc_histo_neg2->Clone();
    mc_temp_neg->Sumw2();
    mc_temp_neg->SetStats(0);
    mc_temp_neg->Divide(mc_histo_0);
    mc_temp_neg->SetMarkerStyle(21);
    mc_temp_neg->SetMarkerColor(kGreen);
    mc_temp_neg->SetLineColor(kGreen);
    mc_temp_neg->GetYaxis()->SetRangeUser(0.3, 2.5);
    mc_temp_neg->GetYaxis()->SetLabelSize(0.1);
    mc_temp_neg->GetXaxis()->SetLabelSize(0.1);
    mc_temp_neg->GetXaxis()->SetTitleSize(0.16);
    mc_temp_neg->GetXaxis()->SetTitleOffset(0.8);
    mc_temp_neg->GetYaxis()->SetTitleOffset(0.2);
    mc_temp_neg->GetYaxis()->SetTitleSize(0.12);
    mc_temp_neg->SetYTitle("THEORY/DATA");
    
    
    TH1F* mc_temp_pos = (TH1F*)  mc_histo_pos2->Clone();
    mc_temp_pos->Sumw2();
    mc_temp_pos->SetStats(0);
    mc_temp_pos->Divide(mc_histo_0);
    mc_temp_pos->SetMarkerStyle(21);
    mc_temp_pos->SetMarkerColor(kRed);
    mc_temp_pos->SetLineColor(kRed);

    
    mc_temp_neg->Draw();
    mc_temp_pos->Draw("SAME");
    pad2->SetGridy();

    
    
    c_basis->SaveAs(basis_canvas_name.c_str());

    
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
    
   // std::pair <TH1F*, vector<TH1F *>> histos;
   // histos = std::make_pair (data_histo, mc_histos);
    
    std::tuple<TH1F*, vector<TH1F*>, vector<TGraphAsymmErrors * >>  histos (  data_histo, mc_histos, error_graphs    );

    
    return histos;
}


void Fitter::scan_couplings(std::string var_name, std::tuple <TH1F*, vector<TH1F *>, vector<TGraphAsymmErrors *>> histos, std::string mode, bool add_pwhg){
    
    if (debug) cout << "Fitter::scan_couplings" << endl;
    
    if (!std::get<0>(histos)) cout << "data histo not found" << endl;
    
    if (debug) cout << "N histos = " <<  std::get<1>(histos).size()  <<endl;
    
    
    TGraphErrors * g = new TGraphErrors();
    TCanvas * c_compare_dists = new TCanvas();
    TPad *pad1 = new TPad("pad1","pad1",0,0.45,1,1);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    //pad1->SetLogy();
 
    
    int n_hists;
    if (add_pwhg){
        n_hists = std::get<1>(histos).size() -1;
    }
    else{
        n_hists = std::get<1>(histos).size();
    }
    
    
    if (debug) cout << "Fitter::scan_couplings:: data integral " << std::get<0>(histos)->Integral()   <<endl;

    for (int weight = 0 ; weight < n_hists ; weight++){
        if (  std::get<0>(histos) && std::get<1>(histos)[weight] ){
            
            if (debug) cout << "Fitter::scan_couplings:: pred integral " << std::get<1>(histos)[weight]->Integral()   <<endl;

            std::get<1>(histos)[weight]->SetLineColor(weight+1);
            std::get<1>(histos)[weight]->SetLineWidth(1.0);
            std::get<2>(histos)[weight]->SetLineColor(weight+1);
            std::get<2>(histos)[weight]->SetFillColor(weight+1);
            std::get<2>(histos)[weight]->SetFillStyle(3001);
            
            //             std::get<2>(histos)[weight]->SetFillStyle(3001);
            //histos.second[weight]->SetFillColor(weight+1);
            //histos.second[weight]->SetFillStyle(3004);
            
            if (weight == 0) {
                std::get<1>(histos)[weight]->SetStats(kFALSE);
                std::get<1>(histos)[weight]->GetYaxis()->SetTitleSize(0.09);
                std::get<1>(histos)[weight]->GetYaxis()->SetTitleOffset(0.35);
                std::get<1>(histos)[weight]->SetYTitle("#frac{1}{#sigma} #frac{#delta(#sigma_tt)}{#delta(x)}");
                std::get<1>(histos)[weight]->Draw("HIST");
                //                  std::get<1>(histos)[weight]->GetYaxis()->SetRangeUser(-0.008, 0.6);
                //std::get<1>(histos)[weight]->GetYaxis()->SetRangeUser(0.16, 0.5);
                std::get<1>(histos)[weight]->GetYaxis()->SetRangeUser(0.06, 0.1);

                std::get<2>(histos)[weight]->Draw("E2SAME");
                
            }else {
                std::get<1>(histos)[weight]->Draw("HISTSAME");
                std::get<2>(histos)[weight]->Draw("E2SAME");
                
               if (debug) cout << "Drawing histo # " << weight  <<endl;
                
            }
            double chi2 = -1.0;
            
            if (mode=="norm_only" || mode=="abs_only") {
                chi2 = this->calculate_test_statistic( std::get<0>(histos), std::get<1>(histos)[weight], std::get<2>(histos)[weight]);
            }
            else if(mode=="norm_fid"){
                chi2 = (this->calculate_test_statistic( std::get<0>(histos),  std::get<1>(histos)[weight], std::get<2>(histos)[weight]) + this->calculate_test_statistic(dummy_fiducial_data, mc_histos_fiducial[weight] , std::get<2>(histos)[weight]  ));
            }
            g->SetPoint(weight, CtG_vals[weight], chi2);
        }
    }
    
    if (add_pwhg)std::get<1>(histos)[std::get<1>(histos).size()-1]->SetLineStyle(2);
    //std::get<1>(histos)[std::get<1>(histos).size()-1]->Draw("HISTSAME");
    
    
    if (  std::get<0>(histos)) {
        std::get<0>(histos)->SetMarkerSize(0.5);
        std::get<0>(histos)->SetMarkerStyle(20);
        std::get<0>(histos)->SetFillColor(kRed);
        std::get<0>(histos)->SetFillStyle(3010);
        std::get<0>(histos)->Draw("E1psame");
    }
    
    TLegend *leg = new TLegend(0.7,0.1,0.85,0.45);
    //leg->SetTextFont(65);
    leg->AddEntry(   std::get<0>(histos) ,"DATA","E1p");
    for (int ctg = 0; ctg < n_preds;  ctg++){
        float ctg_val =CtG_vals[ctg];
        std::string pred_name_leg = "CtG = " +  std::to_string(ctg_val).substr(0,4);
        //    leg->AddEntry( std::get<1>(histos)[ctg], pred_name_leg.c_str(),"l");
        
        leg->AddEntry( std::get<2>(histos)[ctg], pred_name_leg.c_str(),"lf");
        
    }
    
   if (add_pwhg) leg->AddEntry( std::get<1>(histos)[std::get<1>(histos).size()-1], "PWHG+P8","fl");
    leg->Draw();
    
    c_compare_dists->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.45);
    pad2->SetTopMargin(0);
    pad2->Draw();
    pad2->cd();
    pad2->SetBottomMargin(0.3);
    pad2->SetTopMargin(0.0);
    
    if (debug) cout << "Fitter::scan_couplings::making ratio plot" << endl;
    
    for (int histo = 0; histo < std::get<1>(histos).size(); histo++){
        if (debug) cout << "Fitter::scan_couplings::making ratio plot, looping" << endl;
        
        TH1F* mc_temp = (TH1F*)  std::get<1>(histos)[histo]->Clone();
        mc_temp->Sumw2();
        mc_temp->SetStats(0);
        mc_temp->Divide(std::get<0>(histos));
        if (debug) cout << "Fitter::scan_couplings::making ratio plot" <<  ", Nbins data = "   << std::get<0>(histos)->GetNbinsX() << ", Nbins pred = "   <<  mc_temp->GetNbinsX() <<endl;
        
        //mc_temp->SetMarkerStyle(21);
        //  mc_temp->SetMarkerColor(histo+1);
        //  mc_temp->SetLineColor(histo+1);
//        mc_temp->GetYaxis()->SetRangeUser(0.75, 1.4);
        mc_temp->GetYaxis()->SetRangeUser(0.9, 1.2);
        mc_temp->GetYaxis()->SetLabelSize(0.05);
        mc_temp->GetXaxis()->SetLabelSize(0.1);
        mc_temp->GetXaxis()->SetTitleSize(0.16);
        mc_temp->GetXaxis()->SetTitleOffset(0.8);
        mc_temp->GetYaxis()->SetTitleOffset(0.35);
        mc_temp->GetYaxis()->SetTitleSize(0.09);
        mc_temp->SetYTitle("THEORY/DATA");
//        mc_temp->SetXTitle("#Delta#Phi(l^{+}l^{-})");
        mc_temp->SetXTitle("Cos #theta *");
      
        
        std::stringstream ss;
        ss << var_name;
        
        std::string segment;
        std::vector<std::string> seglist;
        
        while(std::getline(ss, segment, '/'))
        {
            seglist.push_back(segment);
        }
        
        std::string xtitle =seglist[1] + "  (GeV)";
        // mc_temp->SetXTitle(xtitle.c_str());
        
        if (add_pwhg) {
        if(histo == std::get<1>(histos).size() -2) {
           mc_temp->SetLineStyle(2);
        }
        }
        
        if(histo == 0){
            mc_temp->Draw("HIST");
        }else{
            mc_temp->Draw("HISTSAME");
        }
    }
    
    c_compare_dists->cd();
    //  pad1->SetLogy();
    
    std::stringstream ss;
    ss << var_name;
    
    std::string segment;
    std::vector<std::string> seglist;
    
    while(std::getline(ss, segment, '/'))
    {
        seglist.push_back(segment);
    }
    
    std::string compare_canvas_name = "compare_" +  seglist[1] + "_.pdf";
    
    //TFile * f_compare = new TFile("compare.root", "RECREATE");
    //c_compare_dists->Write();
    c_compare_dists->SaveAs(compare_canvas_name.c_str());
    
    TCanvas * c1 = new TCanvas();
    // g->SetMinimum(0.002);
    // g->SetMaximum(0.0026);
    g->GetHistogram()->GetXaxis()->SetTitle("CtG");
    g->GetHistogram()->GetYaxis()->SetTitle("#chi^{2} (DATA-THEORY)");
    g->GetHistogram()->GetXaxis()->SetRangeUser(-5.0, 5.0);
    
    g->Draw("AL");
    scans.push_back(g);
    obs_names.push_back(seglist[1]);
    
    c1->SetLogy();
std:string scan_canvas_name = "scans/scan_" +  seglist[1]  + "_.pdf";
    c1->SaveAs(scan_canvas_name.c_str());
    
    std::string scan_rootfile_name = "scan_" +  seglist[1]  + ".root";
    TFile * f_out = new TFile(scan_rootfile_name.c_str(), "RECREATE");
    std::get<0>(histos)->Write();
    for ( int histo = 0; histo <   std::get<1>(histos).size(); histo++){
        std::get<1>(histos)[histo]->Write();
    }
    c1->Write();
    c_compare_dists->Write();
    f_out->Close();
    
    return;
}






