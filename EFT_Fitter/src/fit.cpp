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
           //make_covariance_matrix("files/Oct11/HypLLBarDPhi_totCovMtrxFile.txt", "files/Nov1/particle/absolute/results/DiffXS_HypLLBarDPhi_source.root");
       make_covariance_matrix("files/Nov1/particle/absolute/covariance/HypLLBarDPhi_totCovEnvXSMtrxFile.txt", "files/Nov1/particle/absolute/results/DiffXS_HypLLBarDPhi_source.root", "#Delta #Phi l #bar{l}");
       f_EFT->run_extraction(10, bins_delphill, "data", "files/Nov1/particle/absolute/results/DiffXS_HypLLBarDPhi_source.root", "CMS_dilepton_diff/ll_delphi_abs", mode, false , false);
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

    double x, y, min_y, chi2_at_bestfit;
    vector <double> min_vals;
    min_vals.clear();
    TGraphErrors * gr_rel;
    vector<TGraphErrors*> rel_scans;
    vector<double> boundaries;
    TCanvas * all_relscans_c;
    
    for (int scan = 0; scan< scans.size(); scan ++){
        int n = scans[scan]->GetN();
        double* ny = scans[scan]->GetY();
        int locmax = TMath::LocMin(n,ny);
        double min_y = ny[locmax];
        
        min_vals.push_back(min_y);
        gr_rel = new TGraphErrors();

        for(Int_t i=0; i < scans[scan]->GetN(); i++){
            scans[scan]->GetPoint(i,x,y);
            double rel_y  = y - min_y;
            gr_rel->SetPoint(i,x,rel_y);
            cout <<"rel chi2 vals  = "<< rel_y  << endl;
        }
        rel_scans.push_back(gr_rel);
        boundaries = get_bestfit_CL_boundaries(gr_rel);
        std::cout<<"BEST fit =  "<<boundaries[2]<<" with chi2 = "<<  min_y <<"\n";
        std::cout<<"CL boundaries = "<< boundaries[0]<<" "<<boundaries[1] <<" "<<boundaries[3]<<" "<<boundaries[4]<< endl;
        std::cout<<"Unc = "<<fabs(boundaries[2] - boundaries[0])<<"  "<<fabs(boundaries[2]-boundaries[1])<<"  "<<fabs(boundaries[2]- boundaries[3])<<"  "<<fabs(boundaries[2]-boundaries[4])<<endl;
    }
    
    //now make final plot
    //1. define base histo and make canvas with CMS text
    std::cout <<"1. define base histo and make canvas with CMS text"<< std::endl;

    TH1F * base_histo = new TH1F("","", 1000, -10.0, 10.0);
    all_relscans_c =  make_results_canvas(base_histo, boundaries);
    
    // 2. Add shaded CL bands, best-fit point and 'guide lines'
    std::cout <<"2. Add shaded CL bands, best-fit point and 'guide lines'"<< std::endl;

    TGraphAsymmErrors * g_68 = new TGraphAsymmErrors();
    TGraphAsymmErrors * g_95 = new TGraphAsymmErrors();
    TGraphAsymmErrors * g_central = new TGraphAsymmErrors();
    g_central->SetLineStyle(7);
    g_central->SetLineWidth(2);

    
    g_68->SetPoint(0, boundaries[2], 2.42);
    g_95->SetPoint(0, boundaries[2], 2.42);
    g_central->SetPoint(0, boundaries[2], 2.42);
    
    g_68->SetPointEXhigh(0, fabs(boundaries[2] - boundaries[3]));
    g_68->SetPointEXlow(0, fabs(boundaries[2] - boundaries[1]));
    g_95->SetPointEXhigh(0, fabs(boundaries[2] - boundaries[4]));
    g_95->SetPointEXlow(0, fabs(boundaries[2] - boundaries[0]));
    
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
    
    g_95->Draw("SAME2");
    g_68->Draw("SAME2");
    g_central->Draw("SAMEE1Z");
    
    TLine *line_68 = new TLine(-0.5,1.0,boundaries[3],1.0);
    line_68->SetLineColor(kGreen);
    line_68->SetLineStyle(4);
    line_68->SetLineWidth(4);
    
    TLine *line_95 = new TLine(-0.5,3.84,boundaries[4],3.84);
    line_95->SetLineColor(kYellow);
    line_95->SetLineStyle(4);
    line_95->SetLineWidth(4);
    line_68->Draw();
    line_95->Draw();
    
    //3. Add legend and nominal fit parabola
    TLegend *leg_rel = new TLegend(0.15,0.33,0.38,0.63);
    TFile * f_scan_nom = new TFile("all_scans.root", "RECREATE");
    std::cout <<"3. Add legend and nominal fit parabola"<< std::endl;
    for (int scan = 0; scan< rel_scans.size(); scan ++){
        if (debug) std::cout <<" looping relscans "<<"\n";
        rel_scans[scan]->SetLineColor(scan+1);
        rel_scans[scan]->SetLineWidth(3);
//        rel_scans[scan]->SetLineStyle(2);
        rel_scans[scan]->SetMarkerSize(1);
        rel_scans[scan]->SetMarkerColor(scan+1);
        rel_scans[scan]->SetMarkerStyle(22);
        
        std::string gr_name_rel = obs_names[scan] + "_relscan";
        
        if (scan ==0){
            rel_scans[scan]->Draw("CSAME");
        }else{
            rel_scans[scan]->Draw("CSAME");
        }
        
        //leg_rel->AddEntry( rel_scans[scan], "#Delta#chi^{2} vs. CtG/#Lambda^{2} ","l");
        leg_rel->AddEntry( rel_scans[scan], "nominal fit","l");

        //5. Write parabola graphs (need to manually chnage this file name after running sys. variation)
        rel_scans[scan]->Write();
    }
    
    f_scan_nom->Close();

    //4. Add systematic fit parabolas
    std::cout <<"4. Add systematic fit parabolas"<< std::endl;

    TFile * f_scan_scaledown = new TFile("all_scans_scaledown.root");
    TFile * f_scan_scaleup = new TFile("all_scans_scaleup.root");
    TGraphErrors * gr_down;
    TGraphErrors * gr_up;
    
    if ( f_scan_scaledown && f_scan_scaleup  ){
        for (int scan = 0; scan< scans.size(); scan ++){

            std::string gr_name = obs_names[scan] + "_relscan";
            
            gr_down =  (TGraphErrors*)f_scan_scaledown->Get(";1");
            gr_up =    (TGraphErrors*)f_scan_scaleup->Get(";1");
            gr_down->SetLineStyle(1);
            gr_up->SetLineStyle(1);
            gr_down->SetLineWidth(2);
            gr_up->SetLineWidth(2);
            gr_down->SetLineColorAlpha(kBlue, 0.55);
            gr_up->SetLineColorAlpha(kRed, 0.55);
            gr_down->Draw("CSAME");
            gr_up->Draw("CSAME");
        }
    }
    
    f_scan_scaledown->Close();
    f_scan_scaleup->Close();

    leg_rel->AddEntry( g_central, "best-fit value","E1");
    leg_rel->AddEntry( gr_down, "+ theory syst","l");
    leg_rel->AddEntry( gr_up, "- theory syst","l");
    leg_rel->AddEntry( g_68, "68% CI","f");
    leg_rel->AddEntry( g_95, "95% CI","f");
    leg_rel->SetBorderSize(0);
    leg_rel->SetTextSize(0.035);
    //leg_rel->AddEntry( line_lower , "ar#Chiiv:1503.08841" ,"l");
    leg_rel->Draw();
    
    
    all_relscans_c->RedrawAxis();
    all_relscans_c->SaveAs("all_relscans.pdf");

    /*
    TCanvas * allscans_c = new TCanvas();
    TLegend *leg = new TLegend(0.5,0.6,0.8,0.8);

    scans[0]->GetHistogram()->GetYaxis()->SetRangeUser(0, 120.0);

    for (int scan = 0; scan< scans.size(); scan ++){
        if (debug) std::cout <<" drawing scans "<<"\n";
        scans[scan]->SetLineColor(scan+1);
        scans[scan]->SetLineWidth(2);
        
        std:string gr_name = obs_names[scan] + "_scan";

        if (scan ==0){
            scans[scan]->Draw("AL");
            }else{
                scans[scan]->Draw("L");
            }
            leg->AddEntry( scans[scan], obs_names[scan].c_str(),"l");
            }
    leg->Draw();

    allscans_c->SaveAs("all_scans.pdf");
    
    allscans_c->Write("Scans");
     */

    
    ////////////////////////////////////////////////
    ////////// Make coverage plot //////////////////
    ////////////////////////////////////////////////
    
    /*
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
     */

}




