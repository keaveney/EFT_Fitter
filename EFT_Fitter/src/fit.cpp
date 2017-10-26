//
//  fit.cpp
//  EFT_Fitter
//
//  Created by James Keaveney on 19/01/2017.
//  Copyright © 2017 James Keaveney. All rights reserved.
//

#include "fit.h"
using namespace std;

int main(int argc, const char * argv[]){

    gStyle->SetOptStat(0);
    
    f_EFT->create_dummy_fiducial_measurement(14.1649, 0.02);

    std::string mode = "abs_only";
    
    for(int pred = 0; pred < n_preds; pred++){
    CtG_vals[pred] = (2.0/(n_preds-1))*pred - 1.0;
    }

    if (mode == "abs_only"){
       make_covariance_matrix("HypLLBarDPhi_totCovMtrxFile.txt", "files/Oct4/DiffXS_HypLLBarDPhi_source.root");
       f_EFT->run_extraction(10, bins_delphill, "data", "files/Oct4/DiffXS_HypLLBarDPhi_source.root", "CMS_dilepton_diff/ll_delphi_abs", mode, false , false);
    }
    else if (mode == "norm_fid" || mode == "norm_only" || mode == "fid_only"){
        f_EFT->run_extraction( 13, bins_delphill,"data", "files/April20/Norm/DiffXS_HypLLBarDPhi_source.root", "CMS_dilepton_diff/ll_delphi_abs", mode, false , false);
        f_EFT->run_extraction( 6, bins_ptt, "data_staterror_only","files/April20/Norm/DiffXS_HypToppT_source.root", "CMS_dilepton_diff/t_pT", mode, false , true );
    }

    f_EFT->make_summary_plot(scans);
    if(debug) std::cout <<" make_summary_plot done "<< "\n";

    return 0;
}

void Fitter::run_extraction(int nbins, float bins[], std::string graphname_data, std::string filename_data,std::string histoname_pred,  std::string mode, bool closure_test, bool add_pwhg){
    
     tuple<TH1F*, vector<TH1F *>, vector<TGraphAsymmErrors *>>  histos  = f_EFT->initialise(graphname_data, filename_data, histoname_pred, nbins, bins, mode, closure_test, add_pwhg, "nom");
        // f_EFT->toy_study( histos, 1);
    f_EFT->scan_couplings("data", histoname_pred, histos ,mode, add_pwhg) ;
}





void Fitter::make_summary_plot(vector<TGraphErrors*> scans){
    
   if(debug) std::cout <<" make_summary_plot "<< "\n";
    TFile * all_scans = new TFile("all_scans.root", "RECREATE");

    double x, y, min_y;
    min_y = 99999999.9;
    vector <double> min_vals;
    min_vals.clear();
    
    vector<TGraphErrors*> rel_scans;

    //first find min chi2 of each var for sensitivity comparison
    for (int scan = 0; scan< scans.size(); scan ++){
        
        //int n = scans[scan]->GetN(); //get ploted array dimension
        min_y = 99999999.9;

        for(Int_t i=0; i< scans[scan]->GetN(); i++) {
            scans[scan]->GetPoint(i,x,y);
            if (y < min_y) min_y = y;
           if (debug) std::cout <<" points   = = "<< y <<"\n";
        }
          if (debug) std::cout <<" "<<"\n";
          if (debug) std::cout <<" found min chi2,  var  "<< scan <<"  " <<min_y<<"\n";
          if (debug) std::cout <<" "<<"\n";

      //  if (min_y < 0.000001) min_y = 0.0;
        min_vals.push_back(min_y);
    }
    
    TGraphErrors * gr_pval = new TGraphErrors();
    double min_relchi =999999.9, best_fit, cl_68_lo, cl_95_lo , cl_68_hi, cl_95_hi;
    bool  hi_68 =true , lo_68= true , lo_95 = true, hi_95 =true;

    for (int scan = 0; scan< scans.size(); scan ++){
        TGraphErrors * gr_rel = new TGraphErrors();

        for(Int_t i=0; i < scans[scan]->GetN(); i++) {
            scans[scan]->GetPoint(i,x,y);
            double pval = TMath::Prob(y, 9);
            
            double rel_y  = y - min_vals[scan];
            gr_rel->SetPoint(i,x,rel_y);
            gr_pval->SetPoint(i,x,pval);
            
            if (rel_y < min_relchi){
                min_relchi = rel_y;
                best_fit = x;
            }
            
            if (rel_y <= 3.84 && lo_95) {
                cl_95_lo = x;
                lo_95 = false;
                cout <<"found lo 95"<< endl;

            }
        if (rel_y <= 1.0 && lo_68 && !(lo_95)) {
                cl_68_lo = x;
                lo_68 = false;
                cout <<"found lo 68"<< endl;
            }
            if (rel_y >= 1.0 && hi_68 && !lo_95 && !lo_68) {
                cl_68_hi = x;
                hi_68 = false;
                cout <<"found hi 68"<< endl;


            }
            if (rel_y >= 3.84 && hi_95 && !lo_95 && !lo_68 && !hi_68) {
                cl_95_hi = x;
                hi_95 = false;
                cout <<"found hi 95"<< endl;

            }
            
            

          //  std::cout <<" point " << i  <<" min val =   " << min_vals[scan] <<"\n";
          //  std::cout <<" point " << i  <<", x = "<< x  <<" rel chi =   " <<rel_y<<"\n";
        }
        rel_scans.push_back(gr_rel);
        if (debug) std::cout <<" added graph "<<"\n";

        std::cout <<" BEST fit  =  " << best_fit <<" with rel chi =   " << min_relchi  <<"\n";
        std::cout << "CLs = "<<cl_95_lo<< "  "  <<  cl_68_lo <<"  "<< cl_68_hi << "  "  <<  cl_95_hi<< endl;
        std::cout << "Unc = "<<fabs(best_fit - cl_95_lo) << "  "  <<  fabs( best_fit - cl_68_lo) <<"  "<< fabs( best_fit - cl_68_hi) << "  "  << fabs( best_fit -  cl_95_hi)<< endl;

    }
    
    TCanvas * all_pvalscans_c = new TCanvas();
    gr_pval->Draw("AC");
    all_pvalscans_c->SaveAs("pvals_scan.pdf");

    
    TCanvas * all_relscans_c = new TCanvas("all_relscans","",800,600);
    
    

    TH1F * base_histo = new TH1F("","", 1000, -10.0, 10.0);
    base_histo->GetYaxis()->SetRangeUser(0.0, 5.0);
    base_histo->GetXaxis()->SetRangeUser(-0.5, 0.5);
    
    base_histo->GetYaxis()->SetTitle("#Delta #chi^{2}");
    base_histo->GetXaxis()->SetTitle("CtG/#Lambda^{2}");
    
    base_histo->GetYaxis()->SetLabelSize(0.04);
    base_histo->GetXaxis()->SetLabelSize(0.04);
    base_histo->GetXaxis()->SetTitleSize(0.04);
    base_histo->GetYaxis()->SetTitleSize(0.04);
    base_histo->GetXaxis()->SetTitleOffset(0.95);
    base_histo->GetYaxis()->SetTitleOffset(0.95);
    
    base_histo->Draw();
    

    all_relscans_c->SetTopMargin(0.11);
    all_relscans_c->SetBottomMargin(0.15);

    float H = all_relscans_c->GetWh();
    float W = all_relscans_c->GetWw();
    float l = all_relscans_c->GetLeftMargin();
    float t = all_relscans_c->GetTopMargin();
    float r = all_relscans_c->GetRightMargin();
    float b = all_relscans_c->GetBottomMargin();
    float extraOverCmsTextSize  = 0.76;
    
    TString cmsText, extraText, lumiText;
    cmsText += "CMS";
    extraText += "Preliminary";
    lumiText += "35.9 fb^{-1} (13 TeV)";
    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextSize(0.4*t);
    latex.SetTextColor(kBlack);
    latex.SetTextFont(61);
    latex.SetTextAlign(31);
    latex.DrawLatex(0.17,0.9,cmsText);
    
    latex.SetTextFont(52);
    latex.SetTextSize(0.4*t*extraOverCmsTextSize);
    latex.DrawLatex(0.3,0.9,extraText);
    
    latex.SetTextFont(42);
    latex.SetTextSize(0.3*t);
    latex.DrawLatex(0.9,0.9,lumiText);
    
    
    TGraphAsymmErrors * g_68 = new TGraphAsymmErrors();
    TGraphAsymmErrors * g_95 = new TGraphAsymmErrors();
    TGraphAsymmErrors * g_central = new TGraphAsymmErrors();
    
    
    g_68->SetPoint(0, 0.216, 2.42);
    g_95->SetPoint(0, 0.216, 2.42);
    g_central->SetPoint(0, 0.216, 2.42);
    g_68->SetPointEXhigh(0, 0.116);
    g_68->SetPointEXlow(0, 0.112);
    g_95->SetPointEXhigh(0, 0.224);
    g_95->SetPointEXlow(0, 0.22);
    
    g_central->SetPointEXhigh(0, 0.0);
    g_central->SetPointEXlow(0, 0.0);
    g_central->SetPointEYhigh(0, 1.42);
    g_central->SetPointEYlow(0, 2.42);
    
    g_68->SetPointEYhigh(0, 1.42);
    g_68->SetPointEYlow(0, 2.42);
    g_95->SetPointEYhigh(0, 1.42);
    g_95->SetPointEYlow(0, 2.42);
    
    g_68->SetFillColor(kGreen);
    g_95->SetFillColor(kYellow);
    
    g_68->SetLineWidth(0);
    g_95->SetLineWidth(0);
    
   // g_68->SetFillColorAlpha(kGreen, 0.35);
   // g_95->SetFillColorAlpha(kYellow, 0.35);

    g_95->Draw("SAME2");
    g_68->Draw("SAME2");
    g_central->Draw("SAMEE1Z");
    
    TLine *line_68 = new TLine(-0.5,1.0,0.332,1.0);
    line_68->SetLineColor(kGreen);
    line_68->SetLineStyle(4);
    line_68->SetLineWidth(3);
    
    TLine *line_95 = new TLine(-0.5,3.84,0.44,3.84);
    line_95->SetLineColor(kYellow);
    line_95->SetLineStyle(4);
    line_95->SetLineWidth(3);
    
    line_68->Draw();
    line_95->Draw();
    

    TLegend *leg_rel = new TLegend(0.15,0.33,0.4,0.63);
    
    for (int scan = 0; scan< rel_scans.size(); scan ++){
        if (debug) std::cout <<" looping relscans "<<"\n";

        rel_scans[scan]->SetLineColor(scan+1);
        rel_scans[scan]->SetLineWidth(2);
        rel_scans[scan]->SetLineStyle(2);

        //rel_scans[0]->GetHistogram()->GetYaxis()->SetRangeUser(0.0, 5.0);
        //rel_scans[0]->GetHistogram()->GetXaxis()->SetRangeUser(-0.2, 0.5);


        rel_scans[scan]->SetMarkerSize(1);
        rel_scans[scan]->SetMarkerColor(scan+1);
        rel_scans[scan]->SetMarkerStyle(22);
        
        
        std::string gr_name_rel = obs_names[scan] + "_relscan";
        rel_scans[scan]->Write(gr_name_rel.c_str());
        if (debug) std::cout <<" relscan written "<<"\n";

        
        TSpline3 *s = new TSpline3("grs",rel_scans[scan]);
        s->SetLineColor(scan+1);
        
        if (scan ==0){
            rel_scans[scan]->Draw("CSAME");
            // use a cubic spline to smooth the graph
           // s->Draw("same");
        }else{
            rel_scans[scan]->Draw("CSAME");
           // s->Draw("same");

        }
        leg_rel->AddEntry( rel_scans[scan], "#Delta#chi^{2} vs. CtG/#Lambda^{2} ","l");
    }

    
    leg_rel->AddEntry( g_central, "Best-fit value","E1");
    leg_rel->AddEntry( g_68, "68% CI","f");
    leg_rel->AddEntry( g_95, "95% CI","f");
    leg_rel->SetBorderSize(0);
    
    
    //leg_rel->AddEntry( line_lower , "ar#Chiiv:1503.08841" ,"l");

    
    //all_relscans_c->SetLogy();
    leg_rel->Draw();
    
    
    all_relscans_c->RedrawAxis();
    
    all_relscans_c->SaveAs("all_relscans.pdf");

    
    TCanvas * allscans_c = new TCanvas();
    TLegend *leg = new TLegend(0.5,0.6,0.8,0.8);

    scans[0]->GetHistogram()->GetYaxis()->SetRangeUser(0, 120.0);

    for (int scan = 0; scan< scans.size(); scan ++){
        if (debug) std::cout <<" drawing scans "<<"\n";

        scans[scan]->SetLineColor(scan+1);
        scans[scan]->SetLineWidth(2);
        
       std:string gr_name = obs_names[scan] + "_scan";
        scans[scan]->Write(gr_name.c_str());
        

        if (scan ==0){
            scans[scan]->Draw("AL");
            }else{
                scans[scan]->Draw("L");
            }
            leg->AddEntry( scans[scan], obs_names[scan].c_str(),"l");
            }
    
    leg->Draw();
    //allscans_c->SetLogy();
    allscans_c->SaveAs("all_scans.pdf");
    allscans_c->Write("Scans");
    all_relscans_c->Write("RelScans");
    all_scans->Close();
    
    
    
    TCanvas * limit_c = new TCanvas();
    TFile * f_scan_nominal = new TFile("all_scans_nominal.root");
    TFile * f_scan_scaledown = new TFile("all_scans_scaledown.root");
    TFile * f_scan_scaleup = new TFile("all_scans_scaleup.root");

    if (debug) std::cout <<" got all scan variations "<<"\n";

    
  //  TGraphErrors * gr_nom;
  //  TGraphErrors * gr_down;
   // TGraphErrors * gr_up;
    
    if (  f_scan_nominal && f_scan_scaledown && f_scan_scaleup  ){
        for (int scan = 0; scan< scans.size(); scan ++){
            TCanvas * limit_c = new TCanvas();
            std::string gr_name = obs_names[scan] + "_relscan";
            std::string c_name = obs_names[scan] + "_limit.pdf";

        TGraphErrors * gr_nom =  (TGraphErrors*)f_scan_nominal->Get(gr_name.c_str());
        TGraphErrors * gr_down =  (TGraphErrors*)f_scan_scaledown->Get(gr_name.c_str());
        TGraphErrors * gr_up =    (TGraphErrors*)f_scan_scaleup->Get(gr_name.c_str());
            
            if (debug) std::cout <<" got graph variations "<<"\n";

      //   gr_down->SetLineStyle(2);
      //   gr_up->SetLineStyle(2);
      //   gr_nom->Draw("AL");
      //   gr_down->Draw("LSAME");
      //   gr_up->Draw("LSAME");
         limit_c->SaveAs(c_name.c_str());
            if (debug) std::cout <<" graph variations drawn "<<"\n";
        }
    }
    
    ////////////////////////////////////////////////
    ////////// Make coverage plot //////////////////
    ////////////////////////////////////////////////
    
    TFile * f_toys = new TFile("toy.root");
    TH1F * h_coverage = (TH1F*)f_toys->Get("coverage");
    h_coverage->SetXTitle("best fit CtG/#Lambda^{2}");
    h_coverage->SetYTitle("N_{toys}");

    
    TCanvas * c_coverage = new TCanvas();
    h_coverage->Fit("gaus");
    TF1 *myfunc = (TF1*)h_coverage ->GetFunction("gaus");

    Double_t p1 = myfunc->GetParameter(1);
    Double_t e1 = myfunc->GetParError(1);
    Double_t p2 = myfunc->GetParameter(2);
    Double_t e2 = myfunc->GetParError(2);

    
    TLatex Tl;
    Tl.SetTextAlign(12);
    Tl.SetTextSize(0.03);
  
    TFile * f_cov = new TFile("coverage.root", "RECREATE");
    
    std::string fit_result_1 = "Mean = " + std::to_string(p1) + "+/-" + std::to_string(e1);
    std::string fit_result_2 = "#sigma = " + std::to_string(p2) + "+/-" + std::to_string(e2);
    

    Tl.DrawLatex(0.04,1.0,fit_result_1.c_str());
    Tl.DrawLatex(0.04,0.9,fit_result_2.c_str());

    Float_t ymax = h_coverage->GetMaximum();

    TLine *line1 = new TLine(( p1 - p2 ),0.0,( p1 - p2 ),ymax);
    line1->SetLineColor(kGreen);
    line1->Draw();
    
    TLine *line2 = new TLine(( p1 + p2 ),0.0,( p1 + p2 ),ymax);
    line2->SetLineColor(kGreen);
    line2->Draw();
    
    TLine *line3 = new TLine(( p1 - 2*p2 ),0.0,( p1 - 2*p2 ),ymax);
    line3->SetLineColor(kYellow);
    line3->Draw();
    
    TLine *line4 = new TLine(( p1 + 2*p2 ),0.0,( p1 + 2*p2 ),ymax);
    line4->SetLineColor(kYellow);
    line4->Draw();
    
    Tl.SetTextColor(kGreen);
    Tl.DrawLatex( (p1 - p2), (ymax +1.0), "-1 #sigma");
    Tl.DrawLatex( (p1 + p2), (ymax +1.0), "+1 #sigma");

    Tl.SetTextColor(kYellow);
    Tl.DrawLatex( (p1 - 2*p2), (ymax +1.0), "-2 #sigma");
    Tl.DrawLatex( (p1 + 2*p2), (ymax +1.0), "+2 #sigma");

    
    
    c_coverage->Write();
    
    
    cout <<"MEAN = "<<  p1  <<"  SIGMA = "<< p2<< endl;

}




