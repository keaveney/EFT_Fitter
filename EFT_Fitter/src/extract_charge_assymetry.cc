#include <iostream>
#include <sstream>
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"

#include "TFile.h"
#include "TFrame.h"

#include "TLine.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TF1.h"
#include <iostream>
#include <string>

#include <armadillo>

#include <boost/algorithm/string.hpp>

using namespace std;

using namespace arma;

bool debug = true;

void calculate_CA(TGraph *, std::string);
void run_extraction(std::string);
void summary_plot();


int main(int argc, const char * argv[]) {
    
    gStyle->SetOptStat(0);

    
  //  run_extraction("files/July3/DiffXS_HypTTBarDeltaRapidity_source_norm.root");
  //  run_extraction("files/July3/DiffXS_HypTTBarDeltaRapidity_source_abs.root");
  //  run_extraction("files/July3/DiffXS_HypTTBarDeltaRapidity_source_normparticle.root");
  //  run_extraction("files/July3/DiffXS_HypTTBarDeltaRapidity_source_absparticle.root");
    
    run_extraction("files/Oct11/DiffXS_HypLLBarDEta_source_normparticlell.root");

    summary_plot();
}



void run_extraction(std::string filename){

    std::cout <<" "<< std::endl;
    std::cout <<" "<< std::endl;
    std::cout <<" "<< std::endl;
    std::cout <<"*********************************************"<< std::endl;
    std::cout <<"*********************************************"<< std::endl;
    std::cout <<"*********************************************"<< std::endl;
    std::cout <<"********  Extracting charge asymmetry  ******"<< std::endl;
    std::cout <<"*********************************************"<< std::endl;
    std::cout <<"*********************************************"<< std::endl;
    std::cout <<"*********************************************"<< std::endl;
    std::cout <<" filename  = "<< filename  << std::endl;

    std::string mc = "none";
    
    double x, y;
    double bin_pred;
    TFile * f_CA = new TFile(filename.c_str());
    TGraphAsymmErrors * toppt_delphi = (TGraphAsymmErrors*)f_CA->Get("data");
    TH1F * toppt_delphi_mc = (TH1F*)f_CA->Get("mc");
    TCanvas * c_result = (TCanvas*)f_CA->Get("canvas;1");
    TH1F * h_amCNLO = (TH1F*)c_result->GetPrimitive("POWHEGplot");
    
    
    if (c_result) std::cout <<" got canvas  "<< std::endl;
    if (h_amCNLO) std::cout <<" got primitive  "<< std::endl;
    
        
        
        
    std::cout <<" waMC pred  "<< bin_pred <<std::endl;

    
    if (mc == "pwg"){
        for (int bin = 1; bin <= toppt_delphi_mc->GetNbinsX(); bin++ ){
            toppt_delphi->GetPoint(bin-1, x, y);
            toppt_delphi->SetPoint(bin-1, x, toppt_delphi_mc->GetBinContent(bin));
        }
    }
    else if (mc == "aMC"){
        for (int bin = 1; bin <= toppt_delphi_mc->GetNbinsX(); bin++ ){
            toppt_delphi->GetPoint(bin-1, x, y);
            bin_pred = h_amCNLO->GetBinContent(bin);
            toppt_delphi->SetPoint(bin-1,x, bin_pred);
            std::cout <<" waMC pred  "<< bin_pred <<std::endl;

        }
    }

    std::vector<std::string> results;
    boost::split(results, filename, [](char c){return c == '_';});
    calculate_CA(toppt_delphi, results[3] );
}



void calculate_CA(TGraph * toppt_delphi, std::string write) {
    
    std::cout <<" write string "<< write <<std::endl;


    int n = toppt_delphi->GetN();
    int mid_point = (n/2)-1;
    TFile * f_CA_result;
    
    double x, xsec_bin,  xsec_bin_errup ,xsec_bin_errdown, xsec_neg, xsec_neg_errup, xsec_neg_errdown, xsec_pos, xsec_pos_errup, xsec_pos_errdown, CA_errup, CA_errdown;
    
    TFile * f_cov = new TFile("Unfolding_combined_TtBar_Rapidity_HypTTBarDeltaRapidity.root");
    TH2F * h_cov = (TH2F*)f_cov->Get("SVD_combined_TtBar_Rapidity_HypTTBarDeltaRapidity_STATCORR;1");

    
    for(int p = 0; p < n; p++){
 
        toppt_delphi->GetPoint(p, x, xsec_bin);
        xsec_bin_errup = toppt_delphi->GetErrorYhigh(p);
        xsec_bin_errdown = toppt_delphi->GetErrorYlow(p);
 
        if(p <= mid_point){
            xsec_neg = xsec_neg + xsec_bin;
            xsec_neg_errup = xsec_neg_errup + (xsec_bin + xsec_bin_errup);
            xsec_neg_errdown = xsec_neg_errdown + (xsec_bin - xsec_bin_errdown);
        } else {
            xsec_pos = xsec_pos + xsec_bin;
            xsec_pos_errup = xsec_pos_errup + (xsec_bin + xsec_bin_errup);
            xsec_pos_errdown = xsec_pos_errdown + (xsec_bin - xsec_bin_errdown);
        }
    }
 
    
    double dCA_dbinx = 0.0;
    double dCA_dbinj = 0.0;
    double CA_var = 0.0;
    double bin_var = 0.0;
    double cov_ij;
    
    //propagate uncertainties
    for (int x = 0; x < n; x++){
        for (int j = 0; j < n; j++){
            
            if(x <= mid_point){
                dCA_dbinx = ( - 2*xsec_pos ) / ( pow(xsec_pos + xsec_neg, 2.0) );
            } else if (x > mid_point){
                dCA_dbinx = ( 2*xsec_neg ) / ( pow(xsec_pos + xsec_neg, 2.0) );
            }
        
            if(j <= mid_point){
                dCA_dbinj = ( - 2*xsec_pos ) / ( pow(xsec_pos + xsec_neg, 2.0) );
            } else if (j > mid_point){
                dCA_dbinj = ( 2*xsec_neg ) / ( pow(xsec_pos + xsec_neg, 2.0) );
            }
            
            if (x == j){
            
            cov_ij = 1.0 * (toppt_delphi->GetErrorYhigh(x)) * (toppt_delphi->GetErrorYhigh(j)) ;
                      } else{
                
            cov_ij = (h_cov->GetBinContent(x+2, j+2)/100.0) * (toppt_delphi->GetErrorYhigh(x)) * (toppt_delphi->GetErrorYhigh(j)) ;

                      }

          //  cout <<"i = "<< x <<" j = "<< j <<" COR ij  = "<< h_cov->GetBinContent(x+2, j+2)  << endl;
            
            bin_var = (dCA_dbinx)*(dCA_dbinj) * cov_ij;
            
            CA_var += bin_var;
            
        }
    }
    

    double CA_sd = sqrt(CA_var);
    double CA =  ( xsec_pos - xsec_neg   ) / ( xsec_pos + xsec_neg );
 
 
 
        //'Coherent' errors
        // CA_errup =  ( xsec_pos_errup - xsec_neg_errup   ) / ( xsec_pos_errup + xsec_neg_errup );
        // CA_errdown =( xsec_pos_errdown - xsec_neg_errdown   ) / ( xsec_pos_errdown + xsec_neg_errdown );
    
        //'Antogonistic' errors (max shape deviation, probably overconservative)
        CA_errup =  ( xsec_pos_errup - xsec_neg_errdown   ) / ( xsec_pos_errup + xsec_neg_errdown );
        CA_errdown = ( xsec_pos_errdown - xsec_neg_errup   ) / ( xsec_pos_errdown + xsec_neg_errup );
    
        TGraphAsymmErrors * g_CA_68 = new TGraphAsymmErrors();
        TGraphAsymmErrors * g_CA_95 = new TGraphAsymmErrors();
        TGraphAsymmErrors * g_CA_central = new TGraphAsymmErrors();
   
    
        std::cout <<"********   CA =  "  <<  CA         <<  std::endl;
        std::cout <<"********   CA+ =  " <<  CA_errup   <<  std::endl;
        std::cout <<"********   CA- =  " <<  CA_errdown <<  std::endl;
        std::cout <<"********   CA SD =  " <<  CA_sd <<  std::endl;

    /*
        if(CA_errup > CA && CA_errdown < CA ){
        
            std::cout <<"********   CA =  "<<  CA << " + "<< fabs(CA_errup - CA) << " - "<< fabs(CA - CA_errdown) <<std::endl;
            g_CA->SetPoint(0, CA, 0.0);
            g_CA->SetPointEXhigh(0, fabs(CA_errup - CA));
     g_CA->SetPointEXhigh(0, fabs(CA_errup - CA));
     
        } else if (CA_errup < CA && CA_errdown > CA){
        
            std::cout <<"********   CA =  "<<  CA << " + "<< fabs(CA_errdown - CA) << " - "<< fabs(CA - CA_errup) <<std::endl;
            g_CA->SetPoint(0, CA, 0.0);
            g_CA->SetPointEXhigh(0, fabs(CA_errdown - CA));
            g_CA->SetPointEXlow(0, fabs(CA - CA_errup));
        }
     
     */
        
    double CA_sd_95 = 2.0 * CA_sd;

    
    if(write == "norm.root"){
        g_CA_68->SetPoint(0, CA, 0.0);
        g_CA_95->SetPoint(0, CA, 0.0);
        g_CA_central->SetPoint(0, CA, 0.0);
        g_CA_68->SetPointEXhigh(0, CA_sd);
        g_CA_68->SetPointEXlow(0, CA_sd);
        g_CA_95->SetPointEXhigh(0, CA_sd_95);
        g_CA_95->SetPointEXlow(0, CA_sd_95);

    }else if(write == "abs.root"){
        g_CA_68->SetPoint(0, CA, 1.0);
        g_CA_95->SetPoint(0, CA, 1.0);
        g_CA_central->SetPoint(0, CA, 1.0);
        g_CA_68->SetPointEXhigh(0, CA_sd);
        g_CA_68->SetPointEXlow(0, CA_sd);
        g_CA_95->SetPointEXhigh(0, CA_sd_95);
        g_CA_95->SetPointEXlow(0, CA_sd_95);
    }
    
    else if(write == "normparticle.root"){
        g_CA_68->SetPoint(0, CA, 0.5);
        g_CA_95->SetPoint(0, CA, 0.5);
        g_CA_central->SetPoint(0, CA, 0.5);
        g_CA_68->SetPointEXhigh(0, CA_sd);
        g_CA_68->SetPointEXlow(0, CA_sd);
        g_CA_95->SetPointEXhigh(0, CA_sd_95);
        g_CA_95->SetPointEXlow(0, CA_sd_95);
    }
    else if(write == "absparticle.root"){
        g_CA_68->SetPoint(0, CA, 3.0);
        g_CA_95->SetPoint(0, CA, 3.0);
        g_CA_central->SetPoint(0, CA, 3.0);
        g_CA_68->SetPointEXhigh(0, CA_sd);
        g_CA_68->SetPointEXlow(0, CA_sd);
        g_CA_95->SetPointEXhigh(0, CA_sd_95);
        g_CA_95->SetPointEXlow(0, CA_sd_95);
    }
    else if(write == "normparticlell.root"){
        g_CA_68->SetPoint(0, CA, 1.0);
        g_CA_95->SetPoint(0, CA, 1.0);
        g_CA_central->SetPoint(0, CA, 1.0);
        g_CA_68->SetPointEXhigh(0, CA_sd);
        g_CA_68->SetPointEXlow(0, CA_sd);
        g_CA_95->SetPointEXhigh(0, CA_sd_95);
        g_CA_95->SetPointEXlow(0, CA_sd_95);
    }
    
        g_CA_central->SetPointEXhigh(0, 0.0);
        g_CA_central->SetPointEXlow(0, 0.0);
        g_CA_central->SetPointEYhigh(0, 0.2);
        g_CA_central->SetPointEYlow(0, 0.2);
    
        g_CA_68->SetPointEYhigh(0, 0.2);
        g_CA_68->SetPointEYlow(0, 0.2);
        g_CA_95->SetPointEYhigh(0, 0.2);
        g_CA_95->SetPointEYlow(0, 0.2);
    
        std::string filename = "f_CA_result_" + write;
        f_CA_result = new TFile(filename.c_str(), "RECREATE");
    
        g_CA_68->Write("68");
        g_CA_95->Write("95");
        g_CA_central->Write("central");
        f_CA_result->Write();
    
        return ;

    
}

    
    void summary_plot(){
    
        
        std::cout <<"********   Making summry plot  "   <<  std::endl;

      // vector<std::string> runs = {"norm.root", "abs.root", "normparticle.root","absparticle.root" };
        
        vector<std::string> runs = {"norm.root", "normparticle.root", "normparticlell.root" };

       // vector<std::string> runs = {"norm.root"};

        
        TCanvas * c_results = new TCanvas();
        
        TGraphAsymmErrors * gr_68;
        TGraphAsymmErrors * gr_95;
        TGraphAsymmErrors * gr_central;
        
        TH1F * h_base = new TH1F("","", 1, -0.03, 0.05);
        
        h_base->SetLineWidth(0.0);
        h_base->GetXaxis()->SetRangeUser(-0.07, 0.07);
        h_base->GetYaxis()->SetRangeUser(-0.5, 1.5);
        h_base->GetYaxis()->SetTickLength(0.);
        h_base->GetYaxis()->SetLabelSize(0.);
        h_base->GetXaxis()->SetTitle("A_{C}^{X}");
        h_base->GetYaxis()->SetLabelOffset(999);
        h_base->GetYaxis()->SetAxisColor(0,100.0);
        h_base->Draw();
        
        for(int r =0 ;  r < runs.size(); r++){
            
        string filename = "f_CA_result_" + runs[r];
            
        TFile * f_results = new TFile(filename.c_str());
            
        gr_68 = (TGraphAsymmErrors*)f_results->Get("68;1");
        gr_95 = (TGraphAsymmErrors*)f_results->Get("95;1");
        gr_central = (TGraphAsymmErrors*)f_results->Get("central;1");
        f_results->Close();
            
        gr_68->SetFillColor(kGreen);
        gr_95->SetFillColor(kYellow);
            
            if (r == 0){
        gr_95->Draw("SAME2");
        gr_68->Draw("SAME2");
        gr_central->Draw("E1Z");

            } else{
                gr_95->Draw("SAME2");
                gr_68->Draw("SAME2");
                gr_central->Draw("SAMEE1Z");
            
            }
    
        }
        
        h_base->GetXaxis()->SetRangeUser(-0.03, 0.03);

        
        TLine *line_mg5h6_parton = new TLine(0.0045,-0.2,0.0045,0.2);
        line_mg5h6_parton->SetLineColor(kRed);
        line_mg5h6_parton->Draw();
        
        TLine *line_pwhgp8_parton = new TLine(0.003442,-0.2,0.003442,0.2);
        line_pwhgp8_parton->SetLineColor(kBlue);
        line_pwhgp8_parton->SetLineStyle(2);
        line_pwhgp8_parton->Draw();
        
        TLine *line_pwhgp8_particle = new TLine(0.001551,0.3,0.001551,0.7);
        line_pwhgp8_particle->SetLineColor(kBlue);
        line_pwhgp8_particle->SetLineStyle(2);
        line_pwhgp8_particle->Draw();
        
        TLine *line_pwhgp8_particle_ll = new TLine(-0.0018,0.8,-0.0018,1.2);
        line_pwhgp8_particle_ll->SetLineColor(kBlue);
        line_pwhgp8_particle_ll->SetLineStyle(2);
        line_pwhgp8_particle_ll->Draw();
    
        
        TLine *line_aMC_particle_ll = new TLine(-0.00148,0.8,-0.00148351,1.2);
        line_aMC_particle_ll->SetLineColor(kViolet);
        line_aMC_particle_ll->SetLineStyle(2);
        line_aMC_particle_ll->Draw();
        
        
        TLatex Tl;
        Tl.SetTextAlign(12);
        Tl.SetTextSize(0.03);
        Tl.DrawLatex(-0.027, 0.0,"A_{c}^{tt} parton level");
        Tl.DrawLatex(-0.027, 0.5,"A_{c}^{tt} particle level");
        Tl.DrawLatex(-0.027, 1.0,"A_{c}^{ll}  particle level");
        
        auto legend = new TLegend(0.69,0.45,0.99,0.75);
        legend->AddEntry(line_mg5h6_parton,"MG5_aMC@NLO+HERWIG6","E1");
        legend->AddEntry(line_pwhgp8_parton,"POWHEGv2+PYTHIA8","E1");
        legend->AddEntry(line_aMC_particle_ll,"MG5_aMC@NLO+PYTHIA8","E1");

        legend->AddEntry(gr_68,"Data","E1");
        legend->SetTextSize(0.028);
        legend->Draw();
        
        gStyle->SetOptStat(00000);
        
        TFile * f_results_summary = new TFile("results_summary.root", "RECREATE");
       // c_results->SetFillColor(0);
       // c_results->SetFillStyle(4000);
        //c_results->SetBorderSize(2);
       // c_results->SetFrameFillStyle(0);
        //c_results->SetFrameLineColor(0);
       // c_results->SetFrameFillStyle(0);
        
      //  c_results->SetLeftMargin(0.07);
      //  c_results->SetRightMargin(0.02);
      //  c_results->SetTopMargin(0.06);
      //  c_results->SetBottomMargin(0.16);
        gPad->RedrawAxis();

        c_results->SaveAs("CA_results.pdf");
        c_results->Write();

        
    
    }


