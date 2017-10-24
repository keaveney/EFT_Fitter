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


double Fitter::calculate_test_statistic(TH1F* data, TH1F* pred, TGraphAsymmErrors * gr_errors){
    
    double test_statistic = 0.0;
    int nbins  = data->GetSize() - 2;
    bool schmitt_fit = true;
    vector<double> deltas;
    vector<double> errors;

    double corr_coff = 1.0;
    bool data_errors_only = true;
    double chisq = 0.0;
    double chi2_bin = 0.0;
    double theory_error =0.0;
    double total_error =0.0;
    
    //Retrieve covariance matrix from text file
    TFile *cov_file = new TFile("HypLLBarDPhi_totCovMtrxFile.root");
    TH2D * cov = (TH2D*)cov_file->Get("inv_cov");
   /// cov = Fitter::make_covariance_matrix("HypTTBarDeltaRapidity_totCovMtrxFile.txt");

    deltas.clear();
    errors.clear();
    
    for (int i=1;i <=  data->GetSize() -2; i++){
        
        double data_bin = data->GetBinContent(i);
        double pred_bin = pred->GetBinContent(i);
        double delta = data_bin - pred_bin;
    
        deltas.push_back(delta);
        
        if (data_errors_only){
       //std::cout << "bin  "<< i<<", data_bin "<< data_bin <<", pred_bin = "<< pred_bin <<  " delta =  " << delta << " delta_sq  " << delta*delta <<  ", error  = "<< data->GetBinError(i)  <<  " chi2 in bin =  " <<  ((delta)*(delta)) /  data->GetBinError(i)  <<"\n";
            errors.push_back(data->GetBinError(i));
           // errors.push_back(pred_bin);

        }else{
        //implmement quadrature addition of data and pred errors
            
            if (data_bin > pred_bin ){
            theory_error = gr_errors->GetErrorYhigh(i-1);//accroding to aMC@NLO manual, the scale uncertainty of the NLO predicitons is 10%
            }else if (data_bin < pred_bin ){
                theory_error = gr_errors->GetErrorYlow(i-1);
            }
            total_error =  pow(  pow( (theory_error +  data->GetBinError(i) ) ,2.0), 0.5);
            errors.push_back( total_error);
        }
    }
    
    for (int i=1;i<=nbins; i++) {
        for (int j=1;j<=nbins; j++){
         //   std::cout <<"  "<< std::endl;
            double corr_coff = cov->GetBinContent(i,j);
            if (schmitt_fit){
                chisq += deltas[i-1]*deltas[j-1]*corr_coff;
            //std::cout <<"delta i "<< deltas[i-1]<< " delta j " <<  deltas[j-1] <<"  corr coff "<< corr_coff <<"chi2 "<< chisq << "\n";

                
            }
            else {
                if(i == j){
       //std::cout<<"Calc chi2 " << "bin  "<< i << " delta =  " << deltas[i] <<  ", error  = "<< errors[i]  <<"\n";
                    double delta_sq = deltas[i-1]*deltas[i-1];
                    //double chi2_bin = delta_sq/errors[i-1];
//                    double chi2_bin = delta_sq/(errors[i-1]*errors[i-1]);
                    chi2_bin = delta_sq/(errors[i-1]*errors[i-1]);
                  //  std::cout <<"chi2 "<< chisq << "\n";

                    chisq += chi2_bin;
                }
            }
        }
    }
// std::cout <<"chi2 final "<< chisq << "\n";
    
    cov_file->Close();

    return chisq;
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
    

    for (int scan = 0; scan< scans.size(); scan ++){
        TGraphErrors * gr_rel = new TGraphErrors();
        for(Int_t i=0; i < scans[scan]->GetN(); i++) {
            scans[scan]->GetPoint(i,x,y);
            double rel_y  = y - min_vals[scan];
            gr_rel->SetPoint(i,x,rel_y);
          //  std::cout <<" point " << i  <<" min val =   " << min_vals[scan] <<"\n";
          //  std::cout <<" point " << i  <<" rel chi =   " <<rel_y<<"\n";
        }
        rel_scans.push_back(gr_rel);
        if (debug) std::cout <<" added graph "<<"\n";

    }
    
    TCanvas * all_relscans_c = new TCanvas();
    TLegend *leg_rel = new TLegend(0.5,0.6,0.8,0.8);
    
    // TLine *line_lower = new TLine(-0.42,0,-0.42,2.7);
    // TLine *line_upper = new TLine(0.3,0,0.3,2.7);
    TLine *line_lower = new TLine(-0.32,0,-0.32,2.7);
    TLine *line_upper = new TLine(0.3,0,0.3,2.7);
    
    line_lower->SetLineColor(kGreen);
    line_upper->SetLineColor(kGreen);
    line_lower->SetLineStyle(4);
    line_upper->SetLineStyle(4);
    line_lower->SetLineWidth(3);
    line_upper->SetLineWidth(3);

    
    for (int scan = 0; scan< rel_scans.size(); scan ++){
        if (debug) std::cout <<" looping relscans "<<"\n";

        rel_scans[scan]->SetLineColor(scan+1);
        rel_scans[scan]->SetLineWidth(2);
        rel_scans[0]->GetHistogram()->GetYaxis()->SetRangeUser(-1.0, 10.0);
        rel_scans[0]->GetHistogram()->GetXaxis()->SetRangeUser(-0.5, 0.5);
        rel_scans[scan]->SetMarkerSize(1);
        rel_scans[scan]->SetMarkerColor(scan+1);
        rel_scans[scan]->SetMarkerStyle(22);
        
        
        std::string gr_name_rel = obs_names[scan] + "_relscan";
        rel_scans[scan]->Write(gr_name_rel.c_str());
        if (debug) std::cout <<" relscan written "<<"\n";

        
        TSpline3 *s = new TSpline3("grs",rel_scans[scan]);
        s->SetLineColor(scan+1);
        
        if (scan ==0){
            rel_scans[scan]->Draw("AC");
            // use a cubic spline to smooth the graph
           // s->Draw("same");
        }else{
            rel_scans[scan]->Draw("CSAME");
           // s->Draw("same");

        }
        leg_rel->AddEntry( rel_scans[scan], obs_names[scan].c_str(),"l");
    }
    
    leg_rel->AddEntry( line_lower , "ar#Chiiv:1503.08841" ,"l");

    
    //all_relscans_c->SetLogy();
    leg_rel->Draw();
    
    line_lower->Draw();
    line_upper->Draw();
    
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




