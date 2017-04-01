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
    
f_EFT->create_dummy_fiducial_measurement(14.1649, 0.02);

    
std::string mode = "abs_only";
    
    
    if (mode == "abs_only"){
f_EFT->run_extraction(6,bins_mtt,"data","files/Abs/DiffXS_HypTTBarMass_source.root","CMS_dilepton_diff/ttbar_mass",mode,false);
f_EFT->run_extraction(6,bins_ptt,"data","files/Abs/DiffXS_HypToppT_source.root","CMS_dilepton_diff/t_pT",mode,false );
f_EFT->run_extraction(5,bins_pttt,"data","files/Abs/DiffXS_HypTTBarpT_source.root","CMS_dilepton_diff/ttbar_pT",mode,false);
f_EFT->run_extraction(8,bins_ytt,"data","files/Abs/DiffXS_HypTTBarRapidity_source.root","CMS_dilepton_diff/ttbar_y",mode,false);
f_EFT->run_extraction(8,bins_yt,"data","files/Abs/DiffXS_HypTopRapidity_source.root","CMS_dilepton_diff/t_y",mode,false);
    }
    else if (mode == "norm_fid" || mode == "norm_only" || mode == "fid_only"){
//f_EFT->run_extraction( 6, bins_mtt,"data", "files/Norm/DiffXS_HypTTBarMass_source.root", "CMS_dilepton_diff/ttbar_delphi", mode, false );
//f_EFT->run_extraction( 6, bins_mtt,"data", "files/Norm/DiffXS_HypTTBarMass_source.root", "CMS_dilepton_diff/ll_delphi", mode, false );
f_EFT->run_extraction( 6, bins_mtt,"data", "files/Norm/DiffXS_HypTTBarMass_source.root", "CMS_dilepton_diff/ttbar_mass", mode, false );
f_EFT->run_extraction( 6, bins_ptt, "data","files/Norm/DiffXS_HypToppT_source.root", "CMS_dilepton_diff/t_pT", mode, false  );
f_EFT->run_extraction( 5, bins_pttt,"data", "files/Norm/DiffXS_HypTTBarpT_source.root ", "CMS_dilepton_diff/ttbar_pT", mode, false );
f_EFT->run_extraction( 8, bins_ytt,"data", "files/Norm/DiffXS_HypTTBarRapidity_source.root", "CMS_dilepton_diff/ttbar_y",mode, false );
f_EFT->run_extraction( 8, bins_yt, "data","files/Norm/DiffXS_HypTopRapidity_source.root", "CMS_dilepton_diff/t_y", mode, false );
    }
    
    
    
f_EFT->make_summary_plot(scans);
    
    return 0;
}

void Fitter::run_extraction(int nbins, float bins[], std::string graphname_data, std::string filename_data,std::string histoname_pred, std::string mode, bool closure_test){
    
     pair<TH1F*, vector<TH1F *>>  histos  = f_EFT->initialise(graphname_data, filename_data, histoname_pred, nbins, bins, mode, closure_test);
    f_EFT->scan_couplings(histoname_pred, histos ,mode) ;
}


void Fitter::scan_couplings(std::string var_name,std::pair <TH1F*, vector<TH1F *>> histos, std::string mode){

    if (debug) cout << "Fitter::scan_couplings" << endl;

    if (!  histos.first) cout << "data histo not found" << endl;

    TGraphErrors * g = new TGraphErrors();
    TCanvas * c_compare_dists = new TCanvas();
    TPad *pad1 = new TPad("pad1","pad1",0,0.45,1,1);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    
    double CtG_vals[11] = {-2.0,-1.6,-1.2,-0.8,-0.4,0.0,0.4,0.8,1.2,1.6,2.0};

    for (int weight = 0 ; weight < histos.second.size(); weight++){
        if (  histos.first && histos.second[weight] ){
              histos.second[weight]->SetLineColor(weight);
            if (weight ==0) {
                  histos.second[weight]->SetStats(kFALSE);
                  histos.second[weight]->SetYTitle("#frac{1}{#sigma} \; #frac{#delta(#sigma_tt)}{#delta(x)}");
                  histos.second[weight]->Draw("HIST");
                  histos.second[weight]->GetYaxis()->SetRangeUser(0.0001, 50.0);
            }else {
                  histos.second[weight]->Draw("HISTSAME");
            }
            double chi2 = -1.0;
            
            if (mode=="norm_only" || mode=="abs_only") {
                chi2 = this->calculate_test_statistic(  histos.first,  histos.second[weight]);
            }
            else if(mode=="norm_fid"){
              chi2 = (this->calculate_test_statistic(  histos.first,  histos.second[weight]) + this->calculate_test_statistic(dummy_fiducial_data,mc_histos_fiducial[weight])   );
            }
           g->SetPoint(weight, CtG_vals[weight], chi2);
        }
    }
    
    if (  histos.first) {
          histos.first->SetMarkerSize(1.0);
          histos.first->SetMarkerStyle(22);
          histos.first->SetFillColor(kRed);
          histos.first->SetFillStyle(3005);
          histos.first->Draw("E2psame");
    }
    
    TLegend *leg = new TLegend(0.5,0.6,0.8,0.8);
    //leg->SetTextFont(65);
    leg->AddEntry(   histos.first ,"DATA","E2p");
    leg->AddEntry(   histos.second[1] ,"EFT predictions","fl");
    leg->Draw();
    
    c_compare_dists->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.45);
    pad2->SetTopMargin(0);
    pad2->Draw();
    pad2->cd();
    pad2->SetBottomMargin(0.3);
    pad2->SetTopMargin(0.0);
    
    if (debug) cout << "Fitter::scan_couplings::making ratio plot" << endl;

    for (int histo = 0; histo <   histos.second.size(); histo++){
        if (debug) cout << "Fitter::scan_couplings::making ratio plot, looping" << endl;

    TH1F* mc_temp = (TH1F*)  histos.second[histo]->Clone();
    mc_temp->Sumw2();
    mc_temp->SetStats(0);
    mc_temp->Divide(histos.first);
    mc_temp->SetMarkerStyle(21);
    mc_temp->SetMarkerColor(histo);
    mc_temp->SetLineColor(histo);
    mc_temp->GetYaxis()->SetRangeUser(0.3, 2.5);
    mc_temp->GetYaxis()->SetLabelSize(0.1);
    mc_temp->GetXaxis()->SetLabelSize(0.1);
    mc_temp->GetXaxis()->SetTitleSize(0.16);
    mc_temp->GetXaxis()->SetTitleOffset(0.8);
    mc_temp->GetYaxis()->SetTitleOffset(0.2);
    mc_temp->GetYaxis()->SetTitleSize(0.12);
    mc_temp->SetYTitle("THEORY/DATA");
        
        std::stringstream ss;
        ss << var_name;
        
        std::string segment;
        std::vector<std::string> seglist;
        
        while(std::getline(ss, segment, '/'))
        {
            seglist.push_back(segment);
        }
        
        std::string xtitle =seglist[1] + "  (GeV)";
        mc_temp->SetXTitle(xtitle.c_str());

        if(histo == 0){
            mc_temp->Draw("E1p");
        }else{
            mc_temp->Draw("E1pSAME");
        }
    }
    
    c_compare_dists->cd();
    pad1->SetLogy();
    
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

    g->Draw("PAC");
    scans.push_back(g);
    obs_names.push_back(seglist[1]);
    
    c1->SetLogy();
    std:string scan_canvas_name = "scans/scan_" +  seglist[1]  + "_.pdf";
    c1->SaveAs(scan_canvas_name.c_str());
    
    std::string scan_rootfile_name = "scan_" +  seglist[1]  + ".root";
    TFile * f_out = new TFile(scan_rootfile_name.c_str(), "RECREATE");
      histos.first->Write();
    for ( int histo = 0; histo <   histos.second.size(); histo++){
          histos.second[histo]->Write();
    }
    c1->Write();
    c_compare_dists->Write();
    f_out->Close();
    
    return;
}


double Fitter::calculate_test_statistic(TH1F* data, TH1F* pred){
    
    double test_statistic = 0.0;
    int nbins  = data->GetSize() - 2;
    bool schmitt_fit = false;
    vector<double> deltas;
    vector<double> errors;

    double corr_coff = 1.0;
    bool data_errors_only = true;
    double chisq = 0.0;
    
    deltas.clear();
    errors.clear();
    
    for (int i=1;i<=  data->GetSize() -2   ; i++){
        double data_bin = data->GetBinContent(i);
        double pred_bin = pred->GetBinContent(i);
        double delta = data_bin - pred_bin;
        deltas.push_back(delta);
        
        if (data_errors_only){
      //  std::cout << "bin  "<< i<<", data_bin "<< data_bin <<", pred_bin = "<< pred_bin <<  " delta =  " << delta << " delta_sq  " << delta*delta <<  ", error  = "<< data->GetBinError(i)  <<  "  running chi2 =  " <<  ((delta)*(delta)) /  data->GetBinError(i)  <<"\n";
            errors.push_back(data->GetBinError(i));
        }else{
        //implmement quadrature addition of data and pred errors
            errors.push_back(data->GetBinError(i));
        }
    }
    
    for (int i=1;i<=nbins; i++) {
        for (int j=1;j<=nbins; j++){
            
           // double corr_coff = cov->GetBinContent(i+1,j+1);
            if (schmitt_fit) {
                chisq += deltas[i]*deltas[j]*corr_coff;
            }
            else {
                if(i == j){
       //std::cout<<"Calc chi2 " << "bin  "<< i << " delta =  " << deltas[i] <<  ", error  = "<< errors[i]  <<"\n";
                    double delta_sq = deltas[i-1]*deltas[i-1];
                    //double chi2_bin = delta_sq/errors[i-1];
//                    double chi2_bin = delta_sq/(errors[i-1]*errors[i-1]);
                    double chi2_bin = delta_sq/(0.1* errors[i-1]*errors[i-1]);

                    chisq += chi2_bin;
                }
            }
        }
    }
   if(debug) std::cout <<" chi2 = "<< chisq << "\n";
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
           // std::cout <<" points   = = "<< y <<"\n";
        }
        //std::cout <<" "<<"\n";
        //std::cout <<" found min chi2,  var  "<< scan <<"  " <<min_y<<"\n";
        //std::cout <<" "<<"\n";

        min_vals.push_back(min_y);
    }
    

    for (int scan = 0; scan< scans.size(); scan ++){
        TGraphErrors * gr_rel = new TGraphErrors();
        for(Int_t i=0; i < scans[scan]->GetN(); i++) {
            scans[scan]->GetPoint(i,x,y);
            double rel_y  = y - min_vals[scan];
            gr_rel->SetPoint(i,x,rel_y);
        
    //        std::cout <<" point " << i  <<" rel chi =   " <<rel_y<<"\n";
        }
        rel_scans.push_back(gr_rel);
    }
    

    TCanvas * all_relscans_c = new TCanvas();
    TLegend *leg_rel = new TLegend(0.5,0.6,0.8,0.8);
    
    for (int scan = 0; scan< rel_scans.size(); scan ++){
        rel_scans[scan]->SetLineColor(scan+1);
        rel_scans[scan]->SetLineWidth(2);
        rel_scans[0]->GetHistogram()->GetYaxis()->SetRangeUser(0.0, 80.0);

        if (scan ==0){
            rel_scans[scan]->Draw("AC");
        }else{
            rel_scans[scan]->Draw("C");
        }
        leg_rel->AddEntry( rel_scans[scan], obs_names[scan].c_str(),"l");
    }
    leg_rel->Draw();
    
    //all_relscans_c->SetLogy();

    TLine *line_lower = new TLine(-0.42,0,-0.42,0.10);
    TLine *line_upper = new TLine(0.3,0,0.3,0.10);
    line_lower->SetLineColor(kGreen);
    line_upper->SetLineColor(kGreen);
    line_lower->SetLineStyle(4);
    line_upper->SetLineStyle(4);
    line_lower->SetLineWidth(3);
    line_upper->SetLineWidth(3);
    line_lower->Draw();
    line_upper->Draw();
    
    all_relscans_c->SaveAs("all_relscans.pdf");
    TCanvas * allscans_c = new TCanvas();
    TLegend *leg = new TLegend(0.5,0.6,0.8,0.8);

    scans[0]->GetHistogram()->GetYaxis()->SetRangeUser(100, 20000.0);

    for (int scan = 0; scan< scans.size(); scan ++){
        scans[scan]->SetLineColor(scan+1);
        scans[scan]->SetLineWidth(2);

        if (scan ==0){
            scans[scan]->Draw("AC");
            }else{
                scans[scan]->Draw("C");
            }
            leg->AddEntry( scans[scan], obs_names[scan].c_str(),"l");
            }
    
    leg->Draw();
    //allscans_c->SetLogy();
    allscans_c->SaveAs("all_scans.pdf");
    allscans_c->Write();
    all_relscans_c->Write();
    all_scans->Close();
    
}
