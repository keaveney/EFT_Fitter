#include <iostream>
#include <sstream>
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TMatrixD.h"
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

void calculate_CA(TGraph *, TH2F*, std::string);
void run_extraction(std::string, std::string);
void summary_plot();


int main(int argc, const char * argv[]) {
    
gStyle->SetOptStat(0);
    
run_extraction("files/Nov1/particle/normalised/results/DiffXS_HypTTBarDeltaRapidity_source.root", "files/Nov1/particle/normalised/covariance/HypTTBarDeltaRapidity_totCovEnvXSMtrxFile.root");
    
run_extraction("files/Nov1/parton/normalised/results/DiffXS_HypTTBarDeltaRapidity_source.root", "files/Nov1/parton/normalised/covariance/HypTTBarDeltaRapidity_totCovEnvXSMtrxFile.root");
    
run_extraction("files/Nov1/particle/normalised/results/DiffXS_HypLLBarDEta_source.root", "files/Nov1/particle/normalised/covariance/HypLLBarDEta_totCovEnvXSMtrxFile.root");

//run_extraction("files/Nov1/particle/absolute/results/DiffXS_HypLLBarDEta_source.root", "files/Nov1/particle/absolute/covariance/HypLLBarDEta_totCovEnvXSMtrxFile.root");
    
//    run_extraction("files/Nov1/parton/absolute/results/DiffXS_HypTTBarDeltaRapidity_source.root", "files/Nov1/parton/absolute/covariance/HypTTBarDeltaRapidity_totCovEnvXSMtrxFile.root");
    
 //   run_extraction("files/Nov1/particle/absolute/results/DiffXS_HypTTBarDeltaRapidity_source.root", "files/Nov1/particle/absolute/covariance/HypTTBarDeltaRapidity_totCovEnvXSMtrxFile.root");
    
    
    summary_plot();
}



void run_extraction(std::string filename, std::string cov_filename){

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
    TGraphAsymmErrors * data = (TGraphAsymmErrors*)f_CA->Get("data");
    TH1F * pwhg_mc = (TH1F*)f_CA->Get("mc");
    TCanvas * c_result = (TCanvas*)f_CA->Get("canvas;1");
    TH1F * aMCNLO_mc = (TH1F*)c_result->GetPrimitive("POWHEGplot");
    
    TFile * f_cov = new TFile(cov_filename.c_str());
    TH2F * h_cov = (TH2F*)f_cov->Get("cov");
    
    if (c_result) std::cout <<" got canvas  "<< std::endl;
    if (aMCNLO_mc) std::cout <<" got primitive  "<< std::endl;
    
    //std::cout <<" waMC pred  "<< bin_pred <<std::endl;
    
    if (mc == "pwg"){
        for (int bin = 1; bin <= pwhg_mc->GetNbinsX(); bin++ ){
            data->GetPoint(bin-1, x, y);
            data->SetPoint(bin-1, x, pwhg_mc->GetBinContent(bin));
        }
    }
    else if (mc == "aMC"){
        for (int bin = 1; bin <= aMCNLO_mc->GetNbinsX(); bin++ ){
            data->GetPoint(bin-1, x, y);
            bin_pred = aMCNLO_mc->GetBinContent(bin);
            data->SetPoint(bin-1,x, bin_pred);
            std::cout <<" waMC pred  "<< bin_pred <<std::endl;
        }
    }

    std::vector<std::string> results1, results2;
    boost::split(results1, filename, [](char c){return c == '/';});
    boost::split(results2, results1[5], [](char c){return c == '_';});

    std::string tag;
    
    if (results1[3] == "normalised" && results1[2] == "particle" && results2[1] == "HypTTBarDeltaRapidity" ){
    tag = "normparticle.root";
    } else if (results1[3] == "normalised" && results1[2] == "parton" && results2[1] == "HypTTBarDeltaRapidity"){
    tag = "normparton.root";
    }else if (results1[3] == "normalised" && results1[2] == "particle" && results2[1] == "HypLLBarDEta"){
        tag = "normparticlell.root";
    }else if (results1[3] == "absolute" && results1[2] == "particle" && results2[1] == "HypTTBarDeltaRapidity"){
        tag = "absparticle.root";
    }
    else if (results1[3] == "absolute" && results1[2] == "parton" && results2[1] == "HypTTBarDeltaRapidity"){
        tag = "absparton.root";
    }
    else if (results1[3] == "absolute" && results1[2] == "particle" && results2[1] == "HypLLBarDEta"){
        tag = "absparticlell.root";
    }
    
    
    std::cout <<"mode is  "<< tag << std::endl;
    
    calculate_CA(data, h_cov, tag );
}



void calculate_CA(TGraph * data, TH2F * h_cov, std::string write) {
    
    int n = data->GetN();
    int mid_point = (n/2);
    TFile * f_CA_result;
    
    std::cout <<" N points "<< n  <<  "  mid point = "<< mid_point  <<std::endl;
    
    double x = 0.0; double xsec_bin = 0.0; double xsec_bin_errup = 0.0; double xsec_bin_errdown = 0.0; double xsec_neg = 0.0; double xsec_neg_errup = 0.0; double xsec_neg_errdown = 0.0; double xsec_pos = 0.0; double xsec_pos_errup = 0.0; double xsec_pos_errdown = 0.0; double CA_errup = 0.0; double CA_errdown = 0.0; double rel_unc = 0.0; double rel_unc_corr = 0.0;
    
    for(int p = 0; p < n; p++){
 
        data->GetPoint(p, x, xsec_bin);
        xsec_bin_errup = data->GetErrorYhigh(p);
        xsec_bin_errdown = data->GetErrorYlow(p);
        
        rel_unc = rel_unc + pow((xsec_bin_errup / xsec_bin ), 2.0);
        
        if(p < mid_point){
            xsec_neg = xsec_neg + xsec_bin;
            xsec_neg_errup = xsec_neg_errup + (xsec_bin + xsec_bin_errup);
            xsec_neg_errdown = xsec_neg_errdown + (xsec_bin - xsec_bin_errdown);
            std::cout <<" neg side "<< p  <<  "  xsec=  = "<<  xsec_bin <<" xsec neg "<< xsec_neg <<std::endl;

        } else {
            xsec_pos = xsec_pos + xsec_bin;
            xsec_pos_errup = xsec_pos_errup + (xsec_bin + xsec_bin_errup);
            xsec_pos_errdown = xsec_pos_errdown + (xsec_bin - xsec_bin_errdown);
            std::cout <<" pos side "<< p  <<  "  xsec=  = "<<  xsec_bin <<" xsec pos "<< xsec_pos <<std::endl;
        }
    }
    
    double rel_unc_sqr = pow(rel_unc, 0.5);
 
    double dCA_dbinx = 0.0;
    double dCA_dbinj = 0.0;
    double CA_var = 0.0;
    double bin_var = 0.0;
    double cov_ij;
    
    TMatrixD mat_cov = TMatrixD(n,n);
    TMatrixD mat_derivatives(n, 1);
    TMatrixD mat_derivatives_tp(1,n);
    TMatrixD mat_sd(1,1);

    for (int x = 0; x < n; x++){

        if(x < mid_point){
            mat_derivatives[x][0]  = ( -2*xsec_pos ) / ( pow(xsec_pos + xsec_neg, 2.0) );
            mat_derivatives_tp[0][x]  = ( -2*xsec_pos ) / ( pow(xsec_pos + xsec_neg, 2.0) );
        } else if (x >= mid_point){
            mat_derivatives[x][0] = ( 2*xsec_neg ) / ( pow(xsec_pos + xsec_neg, 2.0) );
            mat_derivatives_tp[0][x] = ( 2*xsec_neg ) / ( pow(xsec_pos + xsec_neg, 2.0) );
        }
    }
    
    //propagate uncertainties
    for (int x = 0; x < n; x++){
        for (int j = 0; j < n; j++){
            
            if(x < mid_point){
                dCA_dbinx = ( -2*xsec_pos ) / ( pow(xsec_pos + xsec_neg, 2.0) );
            } else if (x >= mid_point){
                dCA_dbinx = ( 2*xsec_neg ) / ( pow(xsec_pos + xsec_neg, 2.0) );
            }
        
            if(j < mid_point){
                dCA_dbinj = ( -2*xsec_pos ) / ( pow(xsec_pos + xsec_neg, 2.0) );
            } else if (j >= mid_point){
                dCA_dbinj = ( 2*xsec_neg ) / ( pow(xsec_pos + xsec_neg, 2.0) );
            }

            cov_ij = h_cov->GetBinContent(x+1, j+1) ;
            
            mat_cov[x][j] = cov_ij;

          //  cout <<"i = "<< x <<" j = "<< j <<" COR ij  = "<< h_cov->GetBinContent(x+2, j+2)  << endl;
            rel_unc_corr = rel_unc_corr + pow((xsec_bin_errup / xsec_bin ), 2.0);
            
            bin_var = (dCA_dbinx)*(dCA_dbinj) * cov_ij;
            CA_var += bin_var;
        }
    }

        mat_derivatives.Print();
        mat_derivatives_tp.Print();

        mat_sd = mat_derivatives_tp*mat_cov*mat_derivatives;
        mat_sd.Print();

        double CA_matrix_var = mat_sd[0][0];
        double CA_matrix_sd = pow(CA_matrix_var,0.5);
    
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
        std::cout <<"********   CA SD (uncorrelated) =  " <<  fabs(rel_unc_sqr*CA) <<  std::endl;
        std::cout <<"********   CA SD (matrix) =  " <<  CA_matrix_sd  <<  std::endl;


    
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

    
    if(write == "normparton.root"){
        g_CA_68->SetPoint(0, CA, 0.0);
        g_CA_95->SetPoint(0, CA, 0.0);
        g_CA_central->SetPoint(0, CA, 0.0);
        g_CA_68->SetPointEXhigh(0, CA_sd);
        g_CA_68->SetPointEXlow(0, CA_sd);
        g_CA_95->SetPointEXhigh(0, CA_sd_95);
        g_CA_95->SetPointEXlow(0, CA_sd_95);

    }else if(write == "normparticle.root"){
        g_CA_68->SetPoint(0, CA, 0.5);
        g_CA_95->SetPoint(0, CA, 0.5);
        g_CA_central->SetPoint(0, CA, 0.5);
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
    
    
    else if(write == "absparton.root"){
        g_CA_68->SetPoint(0, CA, 3.0);
        g_CA_95->SetPoint(0, CA, 3.0);
        g_CA_central->SetPoint(0, CA, 3.0);
        g_CA_68->SetPointEXhigh(0, CA_sd);
        g_CA_68->SetPointEXlow(0, CA_sd);
        g_CA_95->SetPointEXhigh(0, CA_sd_95);
        g_CA_95->SetPointEXlow(0, CA_sd_95);
    }
    
    else if(write == "absparticle.root"){
        g_CA_68->SetPoint(0, CA, 4.0);
        g_CA_95->SetPoint(0, CA, 4.0);
        g_CA_central->SetPoint(0, CA, 4.0);
        g_CA_68->SetPointEXhigh(0, CA_sd);
        g_CA_68->SetPointEXlow(0, CA_sd);
        g_CA_95->SetPointEXhigh(0, CA_sd_95);
        g_CA_95->SetPointEXlow(0, CA_sd_95);
    }
    
    else if(write == "absparticlell.root"){
        g_CA_68->SetPoint(0, CA, 5.0);
        g_CA_95->SetPoint(0, CA, 5.0);
        g_CA_central->SetPoint(0, CA, 5.0);
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

      //vector<std::string> runs = {"normparton.root", "absparton.root", "normparticle.root","absparticle.root", "normparticlell.root","absparticlell.root" };
        
        vector<std::string> runs = {"normparton.root", "normparticle.root", "normparticlell.root" };

       // vector<std::string> runs = {"norm.root"};

        TCanvas * c_results = new TCanvas();
        c_results->SetBottomMargin(0.15);
        
        TGraphAsymmErrors * gr_68;
        TGraphAsymmErrors * gr_95;
        TGraphAsymmErrors * gr_central;
        
        TH1F * h_base = new TH1F("","", 1, -0.035, 0.04);
        
        h_base->SetLineWidth(0.0);
        h_base->GetXaxis()->SetRangeUser(-0.07, 0.07);
        h_base->GetYaxis()->SetRangeUser(-0.3, 1.7);
        h_base->GetYaxis()->SetTickLength(0.);
        h_base->GetYaxis()->SetLabelSize(0.);
        h_base->GetXaxis()->SetTitle("A_{C}^{X}");
        h_base->GetXaxis()->SetTitleSize(0.05);
        h_base->GetXaxis()->SetTitleOffset(1.0);

        h_base->GetYaxis()->SetLabelOffset(999);
        //h_base->GetYaxis()->SetAxisColor(0,100.0);
        h_base->Draw();
        
        for(int r =0 ;  r < runs.size(); r++){
            
        string filename = "f_CA_result_" + runs[r];
            
        TFile * f_results = new TFile(filename.c_str());
            
        gr_68 = (TGraphAsymmErrors*)f_results->Get("68;1");
        gr_95 = (TGraphAsymmErrors*)f_results->Get("95;1");
        gr_central = (TGraphAsymmErrors*)f_results->Get("central;1");
        f_results->Close();
            
            
        gr_central->SetLineWidth(3);

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

       // TLine *line_mg5h6_parton = new TLine(0.0045,-0.2,0.0045,0.2);
       // line_mg5h6_parton->SetLineColor(kRed);
       // line_mg5h6_parton->Draw();
        
        TLine *line_pwhgp8_parton = new TLine(0.003442,-0.2,0.003442,0.2);
        line_pwhgp8_parton->SetLineColor(kBlue);
        line_pwhgp8_parton->SetLineStyle(2);
        line_pwhgp8_parton->SetLineWidth(3);
        line_pwhgp8_parton->Draw();
        
        TLine *line_aMC_parton = new TLine(0.003007,-0.2,0.003007,0.2);
        line_aMC_parton->SetLineColor(kRed);
        line_aMC_parton->SetLineStyle(3);
        line_aMC_parton->SetLineWidth(3);
        line_aMC_parton->Draw();
        
        TLine *line_bernreuther_parton = new TLine(0.0075,-0.2,0.0075,0.2);
        line_bernreuther_parton->SetLineColor(kViolet);
        line_bernreuther_parton->SetLineStyle(8);
        line_bernreuther_parton->SetLineWidth(3);
        line_bernreuther_parton->Draw();
        
        TLine *line_pwhgp8_particle = new TLine(0.001551,0.3,0.001551,0.7);
        line_pwhgp8_particle->SetLineColor(kBlue);
        line_pwhgp8_particle->SetLineStyle(2);
        line_pwhgp8_particle->SetLineWidth(3);
        line_pwhgp8_particle->Draw();
        
        TLine *line_aMC_particle = new TLine(0.00221,0.3,0.00221,0.7);
        line_aMC_particle->SetLineColor(kRed);
        line_aMC_particle->SetLineStyle(3);
        line_aMC_particle->SetLineWidth(3);
        line_aMC_particle->Draw();
        
        TLine *line_pwhgp8_particle_ll = new TLine(-0.0018,0.8,-0.0018,1.2);
        line_pwhgp8_particle_ll->SetLineColor(kBlue);
        line_pwhgp8_particle_ll->SetLineStyle(7);
        line_pwhgp8_particle_ll->SetLineWidth(3);
        line_pwhgp8_particle_ll->Draw();
    
        TLine *line_aMC_particle_ll = new TLine(-0.00148351,0.8,-0.00148351,1.2);
        line_aMC_particle_ll->SetLineColor(kRed);
        line_aMC_particle_ll->SetLineStyle(3);
        line_aMC_particle_ll->SetLineWidth(3);
        line_aMC_particle_ll->Draw();

        
        TLatex Tl;
        Tl.SetTextAlign(12);
        Tl.SetTextSize(0.045);
        Tl.DrawLatex(-0.034, 0.0,"A_{c}^{tt} parton level");
        Tl.DrawLatex(-0.034, 0.5,"A_{c}^{tt} particle level");
        Tl.DrawLatex(-0.034, 1.0,"A_{c}^{ll}  particle level");
        
       // Tl.DrawLatex(-0.03, 3.0,"A_{c}^{tt} parton level");
       // Tl.DrawLatex(-0.03, 4.0,"A_{c}^{tt} particle level");
       // Tl.DrawLatex(-0.03, 5.0,"A_{c}^{ll}  particle level");
        
        float H = c_results->GetWh();
        float W = c_results->GetWw();
        float l = c_results->GetLeftMargin();
        float t = c_results->GetTopMargin();
        float r = c_results->GetRightMargin();
        float b = c_results->GetBottomMargin();
        float extraOverCmsTextSize  = 0.8;
        
        TString cmsText, extraText, lumiText;
        cmsText += "CMS";
        extraText += "";
        lumiText += "35.9 fb^{-1} (13 TeV)";
        TLatex latex;
        latex.SetNDC();
        latex.SetTextAngle(0);
        latex.SetTextSize(0.46*t);
        latex.SetTextColor(kBlack);
        latex.SetTextFont(61);
        latex.SetTextAlign(31);
        latex.DrawLatex(0.17,0.916,cmsText);
        
        latex.SetTextFont(52);
        latex.SetTextSize(0.71*t*extraOverCmsTextSize);
        latex.DrawLatex(0.3,0.912,extraText);
        
        latex.SetTextFont(42);
        latex.SetTextSize(0.46*t);
        latex.DrawLatex(0.9,0.916,lumiText);
        
        gr_68->SetLineWidth(0);
        gr_95->SetLineWidth(0);
        
        auto legend = new TLegend(0.49,0.55,0.79,0.85);
        legend->AddEntry(gr_central,"Data","E1Z");
        legend->AddEntry(gr_68,"68% CI","f");
        legend->AddEntry(gr_95,"95% CI","f");
        legend->AddEntry(line_pwhgp8_parton,"POWHEGV2 + PYTHIA8","E1");
        legend->AddEntry(line_aMC_particle_ll,"MG5_aMC@NLO + PYTHIA8 [FxFx]","E1");
        legend->AddEntry(line_bernreuther_parton,"NLO+EW","E1");
        legend->SetTextSize(0.032);
        legend->SetBorderSize(0);
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


