#include <iostream>
#include <sstream>
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"

#include "TFile.h"
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

double calculate_CA(TGraph *, std::string);
void run_extraction(std::string);
void summary_plot();


int main(int argc, const char * argv[]) {
    
    run_extraction("files/July3/DiffXS_HypTTBarDeltaRapidity_source_norm.root");
    run_extraction("files/July3/DiffXS_HypTTBarDeltaRapidity_source_abs.root");
    summary_plot();
}



void run_extraction(std::string filename){

    std::cout <<"*********************************************"<< std::endl;
    std::cout <<"*********************************************"<< std::endl;
    std::cout <<"*********************************************"<< std::endl;
    std::cout <<"********  Extracting charge asymmetry  ******"<< std::endl;
    std::cout <<"*********************************************"<< std::endl;
    std::cout <<"*********************************************"<< std::endl;
    std::cout <<"*********************************************"<< std::endl;
    
    TFile * f_CA = new TFile(filename.c_str());
    TGraphAsymmErrors * toppt_delphi = (TGraphAsymmErrors*)f_CA->Get("data");
    double x, x_err,  xsec_bin, xsec_bin_errup, xsec_bin_errdown,  xsec_neg = 0.0, xsec_pos = 0.0, xsec_neg_errup = 0.0,xsec_neg_errdown = 0.0, xsec_pos_errup = 0.0, xsec_pos_errdown = 0.0,  xsec_pos_err = 0.0, CA, CA_errup, CA_errdown, CA_err,CA_err_abs;
    
    int n = toppt_delphi->GetN();
    int mid_point = (n)/2.0;
    double x1_ind, x2_ind, x1_unc, x2_unc;
    
    TRandom3 * gRand = new TRandom3();
    TH1F * h_CA_DATA = new TH1F("","", 30, -0.02, 0.02);
    
    
    std::vector<std::string> results;
    boost::split(results, filename, [](char c){return c == '_';});
    
    cout <<"RESULT TYPE =  "<< results[3]  << endl;

    
    //1 extract CA from DATA
    double CA_DATA =  calculate_CA(toppt_delphi, "write");
    //h_CA_DATA->Fill(CA_DATA);
    
    cout <<" CA DATA "<< " = "<< CA_DATA << endl;

   //2 check coverage
    
    
    /*
    mat A(8,8);
    mat C(8,8);


    //read real cov matrix
    double cov1, cov2, corr;
    TFile * f_cov = new TFile("Unfolding_combined_TtBar_Rapidity_HypTTBarDeltaRapidity.root");
    TH2F * h_cov = (TH2F*)f_cov->Get("SVD_combined_TtBar_Rapidity_HypTTBarDeltaRapidity_STATCORR;1");
    
    //Generate correlated toys (not fully working yet)
    for (double x1 = 0; x1 <= 7; x1++ ){
        
        toppt_delphi->GetPoint(x1, x, xsec_bin);
        xsec_bin_errup = toppt_delphi->GetErrorYhigh(x1+1);
        xsec_bin_errdown = toppt_delphi->GetErrorYlow(x1+1);
        x1_unc = (xsec_bin_errup + xsec_bin_errdown)/(2.0);
        
        
        for (double x2 = 0; x2 <= 7; x2++ ){
            
            toppt_delphi->GetPoint(x2, x, xsec_bin);
            xsec_bin_errup = toppt_delphi->GetErrorYhigh(x2+1);
            xsec_bin_errdown = toppt_delphi->GetErrorYlow(x2+1);
            x2_unc = (xsec_bin_errup + xsec_bin_errdown)/(2.0);
            
            if (fabs (x1 - x2) == 0   ){
                C(x1,x2) = 1.0;
                x1_ind = (1.0)*(x1_unc * x2_unc);
                
            } else{
                corr = (h_cov->GetBinContent(x1+2,x2+2)/100.0);
                C(x1,x2) = corr;
                x1_ind = (corr)*(x1_unc * x2_unc);
            }
            
            A(x1,x2) = x1_ind;

        }
    }
    
    
    mat R1 = chol(A, "lower");
    C.print("COR:");
    A.print("COV:");
    R1.print("chol(COV):");
    
    mat T1(8,1);
    mat T2(8,1);
    
    
    T2 = R1 * T1;
    T1.print("Toy:");
    T2.print("Corr. Toy:");

    */

  /*
    TGraph * toy_gr = new TGraph();
    TH1F * h_CA_toys = new TH1F("","", 30, -0.15, 0.15);

    for (int toy = 0; toy < 10000; toy++){
        
        toy_gr->Set(0);
        
        for(int p = 0; p < n; p++){
            
            toppt_delphi->GetPoint(p, x, xsec_bin);
            xsec_bin_errup = toppt_delphi->GetErrorYhigh(p);
            xsec_bin_errdown = toppt_delphi->GetErrorYlow(p);
            
            double toy_x1 = gRand->Gaus(xsec_bin,(xsec_bin_errup + xsec_bin_errdown)/2.0);
            toy_gr->SetPoint(p, x,toy_x1 );
            
        }
        
        double CA_toy = calculate_CA(toy_gr, "calc");
        h_CA_toys->Fill(CA_toy);
        cout <<" CA toy # "<< toy << " = "<< CA_toy << endl;
        cout <<" "<< endl;

    }
    
    
    h_CA_toys->Fit("gaus");
    TF1 *myfunc = (TF1*)h_CA_toys->GetFunction("gaus");
    
    Double_t p1 = myfunc->GetParameter(1);
    Double_t e1 = myfunc->GetParError(1);
    Double_t p2 = myfunc->GetParameter(2);
    Double_t e2 = myfunc->GetParError(2);
    
    
    
    std::string results_filename = "results_" + std::string(results[3]);
    
    double y_pos = 0.0;
    
    if (results[3] == "norm.root"){
        y_pos = 0.0;
        
    }else if (results[3] == "abs.root") {
        y_pos = 1.0;
    }
    
    TFile * f_results = new TFile(results_filename.c_str(), "RECREATE");
    TGraphAsymmErrors * gr_68 = new TGraphAsymmErrors();
    TGraphAsymmErrors * gr_95 = new TGraphAsymmErrors();
    TGraphAsymmErrors * gr_central = new TGraphAsymmErrors();
    

    gr_central->SetPoint(0, CA_DATA, y_pos);
    gr_central->SetPointEYlow(0, 0.2);
    gr_central->SetPointEYhigh(0, 0.2);


    gr_68->SetPoint(0, CA_DATA, y_pos);
    gr_68->SetPointEXlow(0, p2);
    gr_68->SetPointEXhigh(0, p2);
    gr_68->SetPointEYlow(0, 0.2);
    gr_68->SetPointEYhigh(0, 0.2);

    
    
    gr_95->SetPoint(0, CA_DATA, y_pos);
    gr_95->SetPointEXlow(0, p2*2.0);
    gr_95->SetPointEXhigh(0, p2*2.0);
    gr_95->SetPointEYlow(0, 0.2);
    gr_95->SetPointEYhigh(0, 0.2);

    gr_68->Write("68");
    gr_95->Write("95");
    gr_central->Write("central");

    
    h_CA_toys->Write();
    h_CA_DATA->Write();
    f_results->Write();
    f_results->Close();
   
   */
    
    
    /*
    //more playing with correlated random numbers
    TH2F * h_uncorr = new TH2F("","", 50, -2.0, 2.0, 50, -2.0, 2.0);
    TH2F * h_corr = new TH2F("","", 50, -2.0, 2.0, 50, -2.0, 2.0);
    
    double mean1 = 0.0 , width1 = 0.2, mean2 = 0.0 , width2 = 1.0;
    
    for (int bin = 1; bin <= 8; bin++){
        
        double toy_x1 = gRand->Gaus(mean1, width1);
        double toy_x2 = gRand->Gaus(mean2, width2);
        
        h_uncorr->Fill(toy_x1, toy_x2);
        
        double rho = 0.5;
        
        
        //x_rotated = ((x - x_origin) * cos(angle)) - ((y_origin - y) * sin(angle)) + x_origin
        //y_rotated = ((y_origin - y) * cos(angle)) - ((x - x_origin) * sin(angle)) + y_origin
        
        double x_rotated = ((toy_x1 - mean1) * cos(rho)) - ((mean2 - toy_x2) * sin(rho)) + mean1;
        double y_rotated = ((mean2 - toy_x2) * cos(rho)) - ((toy_x1 - mean1) * sin(rho)) + mean2;
        
        h_corr->Fill(x_rotated,y_rotated);
        
    }
    
    TCanvas * c_corrtoys = new TCanvas();
    h_corr->SetMarkerColor(kRed);
    h_uncorr->Draw();
    h_corr->Draw("SAME");

    c_corrtoys->SaveAs("corrtoys.pdf");
    */

    
}




double calculate_CA(TGraph * toppt_delphi, std::string write) {

    int n = toppt_delphi->GetN();
    int mid_point = (n/2)-1;
    
    
    double x, xsec_bin,  xsec_bin_errup ,xsec_bin_errdown, xsec_neg, xsec_neg_errup, xsec_neg_errdown, xsec_pos, xsec_pos_errup, xsec_pos_errdown, CA_errup, CA_errdown;
    
    TH2F * h_cov = (TH2F*)f_cov->Get("SVD_combined_TtBar_Rapidity_HypTTBarDeltaRapidity_STATCORR;1");
    double cov[8][8];

    
    for(int p = 0; p < n; p++){
 
        toppt_delphi->GetPoint(p, x, xsec_bin);
        xsec_bin_errup = toppt_delphi->GetErrorYhigh(p);
        xsec_bin_errdown = toppt_delphi->GetErrorYlow(p);
        
 
        if(p <= mid_point){
            xsec_neg = xsec_neg + xsec_bin;
            xsec_neg_errup = xsec_neg_errup + (xsec_bin + xsec_bin_errup);
            xsec_neg_errdown = xsec_neg_errdown + (xsec_bin - xsec_bin_errdown);
            cout <<"Neg side  ="<<  xsec_bin << endl;
 
        } else {
            xsec_pos = xsec_pos + xsec_bin;
            xsec_pos_errup = xsec_pos_errup + (xsec_bin + xsec_bin_errup);
            xsec_pos_errdown = xsec_pos_errdown + (xsec_bin - xsec_bin_errdown);
            cout <<"Pos side  ="<<  xsec_bin  << endl;

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
                
            cov_ij = (h_cov->GetBinContent(x+2, j+2)) * (toppt_delphi->GetErrorYhigh(x)) * (toppt_delphi->GetErrorYhigh(j)) ;

                      }
                      
                      
            bin_var = (dCA_dbinx)*(dCA_dbinj) * cov_ij;
            
            
                      CA_var += bin_var;
            
        }
    }
    

    double CA_sd = sqrt(CA_var);
    double CA =  ( xsec_pos - xsec_neg   ) / ( xsec_pos + xsec_neg );
 
 
    if(write == "write"){
 
        //'Coherent' errors
        // CA_errup =  ( xsec_pos_errup - xsec_neg_errup   ) / ( xsec_pos_errup + xsec_neg_errup );
        // CA_errdown =( xsec_pos_errdown - xsec_neg_errdown   ) / ( xsec_pos_errdown + xsec_neg_errdown );
    
        //'Antogonistic' errors
        CA_errup =  ( xsec_pos_errup - xsec_neg_errdown   ) / ( xsec_pos_errup + xsec_neg_errdown );
        CA_errdown = ( xsec_pos_errdown - xsec_neg_errup   ) / ( xsec_pos_errdown + xsec_neg_errup );
    
        TGraphAsymmErrors * g_CA = new TGraphAsymmErrors();
    
        std::cout <<"********   CA =  "  <<  CA         <<  std::endl;
        std::cout <<"********   CA+ =  " <<  CA_errup   <<  std::endl;
        std::cout <<"********   CA- =  " <<  CA_errdown <<  std::endl;
        std::cout <<"********   CA SD =  " <<  CA_sd <<  std::endl;

    
        if(CA_errup > CA && CA_errdown < CA ){
        
            std::cout <<"********   CA =  "<<  CA << " + "<< fabs(CA_errup - CA) << " - "<< fabs(CA - CA_errdown) <<std::endl;
            g_CA->SetPoint(0, CA, 0.0);
            g_CA->SetPointEXhigh(0, fabs(CA_errup - CA));
            g_CA->SetPointEXlow(0, fabs(CA - CA_errdown));
        
        } else if (CA_errup < CA && CA_errdown > CA){
        
            std::cout <<"********   CA =  "<<  CA << " + "<< fabs(CA_errdown - CA) << " - "<< fabs(CA - CA_errup) <<std::endl;
            g_CA->SetPoint(0, CA, 0.0);
            g_CA->SetPointEXhigh(0, fabs(CA_errdown - CA));
            g_CA->SetPointEXlow(0, fabs(CA - CA_errup));
        }
    
        TCanvas * c_CA = new TCanvas();
        TFile * f_CA_result = new TFile("f_CA_result.root", "RECREATE");
    
        g_CA->GetHistogram()->GetXaxis()->SetTitle("A_{C}");

    
        g_CA->SetMarkerStyle(20);
        g_CA->Draw("AE1p");
    
        g_CA->GetHistogram()->GetXaxis()->SetRangeUser(-0.1, 0.1);
        g_CA->GetHistogram()->GetYaxis()->SetRangeUser(-1.0, 1.0);
    
        TLine *line = new TLine(0.0045,-0.2,0.0045,1.2);
        line->SetLineColor(kRed);
        line->Draw();
    
        auto legend = new TLegend(0.2,0.6,0.53,0.85);
        legend->AddEntry(line,"MG5_aMC@NLO+HERWIG6","l");
        legend->AddEntry(g_CA,"Data","E1p");
        legend->Draw();
    
        g_CA->Write();
        c_CA->Write();
    
        c_CA->SaveAs("CA.pdf");
        }
    
    return CA;

    
}

    
    void summary_plot(){
    
        vector<std::string> runs = {"norm.root", "abs.root"};
        
        TCanvas * c_results = new TCanvas();
        
        TGraphAsymmErrors * gr_68;
        TGraphAsymmErrors * gr_95;
        TGraphAsymmErrors * gr_central;
        
        TH1F * h_base = new TH1F("","", 1, -0.15, 0.15);
        h_base->SetLineWidth(0.0);
        h_base->GetXaxis()->SetRangeUser(-0.07, 0.07);
        h_base->GetYaxis()->SetRangeUser(-1.0, 2.0);
        h_base->GetYaxis()->SetTickLength(0.);
        h_base->GetYaxis()->SetLabelSize(0.);
        h_base->Draw();
        
        for(int r;  r < runs.size(); r++){
            
            
        string filename = "results_" + runs[r];
            
        TFile * f_results = new TFile(filename.c_str());
        gr_68 = (TGraphAsymmErrors*)f_results->Get("68");
        gr_95 = (TGraphAsymmErrors*)f_results->Get("95");
        gr_central = (TGraphAsymmErrors*)f_results->Get("central");
        f_results->Close();
            
            if (r == 0){
        gr_95->Draw("SAME2");
        gr_68->Draw("SAME2");
        gr_central->Draw("SAMEP");
            } else{
            
                gr_95->Draw("SAME2");
                gr_68->Draw("SAME2");
                gr_central->Draw("SAMEP");
            
            }
        gr_68->SetFillColor(kGreen);
        gr_95->SetFillColor(kYellow);

        gr_95->GetHistogram()->GetXaxis()->SetTitle("A_{C}");
       // gr_95->GetHistogram()->GetYaxis()->SetTickLength(0.);
       // gr_95->GetHistogram()->GetYaxis()->SetLabelSize(0.);

    
        }
        
        
        TLine *line = new TLine(0.0045,-0.5,0.0045,0.5);
        line->SetLineColor(kRed);
        line->Draw();
        
        TLatex Tl;
        Tl.SetTextAlign(12);
        Tl.SetTextSize(0.03);
        Tl.DrawLatex(-0.14, 0.0,"normalised parton level");
        Tl.DrawLatex(-0.14, 1.0,"absolute parton level");
        
        auto legend = new TLegend(0.2,0.6,0.53,0.85);
        legend->AddEntry(line,"MG5_aMC@NLO+HERWIG6","l");
        legend->AddEntry(gr_68,"Data","E1p");
        legend->Draw();
        
        TFile * f_results_summary = new TFile("results_summary.root", "RECREATE");
        c_results->Write();

        
    
    }


