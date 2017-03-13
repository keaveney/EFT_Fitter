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

    float bins_mtt[] = { 0.0, 380.0, 470.0, 620.0, 820.0, 1100.0, 1600.0 };
    float bins_ptt[] = { 0.0, 65.0, 125.0, 200.0, 290.0, 400.0, 545.0 };
    float bins_pttt[] = { 0.0, 30.0, 80.0, 170.0, 300.0, 500.0 };
    float bins_ytt[] = { -2.5, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.5};
    float bins_yt[] = { -2.5, -1.6, -1.0, -0.5, 0.0, 0.5, 1.0, 1.6, 2.5};

    f_EFT->run_extraction( 6, bins_mtt, "files/DiffXS_HypTTBarMass_source.root", "CMS_dilepton_diff/ttbar_mass", true, "data" );
    f_EFT->run_extraction( 6, bins_ptt, "files/DiffXS_HypToppT_source.root", "CMS_dilepton_diff/t_pT", true, "data" );
    f_EFT->run_extraction( 5, bins_pttt, "files/DiffXS_HypTTBarpT_source.root ", "CMS_dilepton_diff/ttbar_pT", true, "data" );
    f_EFT->run_extraction( 8, bins_ytt, "files/DiffXS_HypTTBarRapidity_source.root", "CMS_dilepton_diff/ttbar_y", true, "data" );
    f_EFT->run_extraction( 8, bins_yt, "files/DiffXS_HypTopRapidity_source.root", "CMS_dilepton_diff/t_y", true, "data" );

    f_EFT->make_summary_plot(scans);
    
    return 0;
}

void Fitter::run_extraction(int nbins, float bins[], std::string filename_data,std::string histoname_pred, bool norm, std::string mode){
    
    //First prepare histos
    pair<TH1F*, vector<TH1F *>>  histos  = f_EFT->initialise(filename_data, histoname_pred, nbins, bins, norm, mode);
    
    if(debug)    cout << "Histos prepared\n";

    //Second run scan over predictions
    f_EFT->scan_couplings(histoname_pred, histos) ;
}


std::pair <TH1F*, vector<TH1F *>> Fitter::initialise(std::string filename_data, std::string histoname_pred, int nbins, float bins[], bool norm, string mode){
   
    TH1F* data_histo;
    double running_total = 0.0;

    if (mode == "closuretest"){
        
        TFile* f_data = new TFile("files/EFT_tt_rivet_small4.root");
        data_histo = (TH1F*)f_data->Get("CMS_dilepton_diff/ttbar_mass");
        if(debug)    cout << "Closure test\n";
        
    }else if (mode == "data"){
        if(debug)    cout << "data mode\n";
        if(debug)    cout << "Extracting data graph from " <<  filename_data.c_str()    <<"\n";

        TFile *f_data = new TFile(filename_data.c_str());
        TGraphAsymmErrors* g_data = (TGraphAsymmErrors*) f_data->Get("data");
        
        data_histo = new TH1F("n","t", nbins, bins);
        double bin_centre, bin_height;

        for (int point = 0; point <= nbins; point ++){
            if(debug)    cout << "looping on graph points, "<<  " running_total " << running_total  <<"\n";
            g_data->GetPoint(point, bin_centre, bin_height);
            data_histo->SetBinContent(point+1, bin_height);
            double bin_error = g_data->GetErrorYhigh(point);
            data_histo->SetBinError(point+1, bin_error);
            running_total = running_total + bin_height;
            }
        }
    
        if(debug)    cout << "Running total for data histo  =" <<  running_total  <<"\n";
    
        //first extract raw histos to make basis histos
        string filename_neg2 = "files/EFTNLO_CtGNeg2.root";
        string filename_pos2 = "files/EFTNLO_CtGPos2.root";
        string filename_0 = "files/EFTNLO_CtG0.root";
        
        TFile * f_neg2 = new TFile(filename_neg2.c_str());
        TFile * f_0 = new TFile(filename_0.c_str());
        TFile * f_pos2 = new TFile(filename_pos2.c_str());
        
        TH1F * mc_histo_neg2 = (TH1F*)f_neg2->Get(histoname_pred.c_str());
        TH1F * mc_histo_0 = (TH1F*)f_0->Get(histoname_pred.c_str());
        TH1F * mc_histo_pos2 = (TH1F*)f_pos2->Get(histoname_pred.c_str());

        if (!mc_histo_neg2) cout << "mc histo: "<<  filename_neg2  <<" not found" << endl;
    
        //make histo of pure CtG contribution
        TH1F * h_pure_Ctg = (TH1F*)mc_histo_pos2->Clone();
        h_pure_Ctg->Add(mc_histo_neg2, -1);
 
    
      //  std::string basis_name = "basis_histos/" + histoname_pred + "_basis.root";
       // TFile * f_basis = new TFile(basis_name.c_str(), "RECREATE");
       // mc_histo_neg2->SetName("Neg2");
      //  mc_histo_neg2->Write();
      //  mc_histo_0->SetName("0");
      //  mc_histo_0->Write();
      //  mc_histo_pos2->SetName("Pos2");
      //  mc_histo_pos2->Write();
      //  h_pure_Ctg->Write();
      //  f_basis->Close();
      //  f_neg2->Close();
      //  f_pos2->Close();
      //  f_0->Close();

    
        double CtG_vals[11] = {-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0};
        double scaling =1.0, CtG;
        vector<TH1F *> mc_histos;
        TH1F * h_CtG_pred;
    
        //loop over chosen CtG values and make prediction histos
        for (int scale = 0 ; scale < 11; scale++){
            CtG = CtG_vals[scale];
            h_CtG_pred = (TH1F*)mc_histo_0->Clone();
            h_CtG_pred->Add(h_pure_Ctg, CtG/2.0);
            
            for (int bin = 1 ; bin <= nbins; bin++){
                double bin_width = bins[bin] - bins[bin-1];
                double bin_xsec  = h_CtG_pred->GetBinContent(bin) / bin_width;
                h_CtG_pred->SetBinContent(bin, bin_xsec);
            }
            
            scaling = running_total/h_CtG_pred->Integral();
            if(norm) h_CtG_pred->Scale(scaling);
            mc_histos.push_back(h_CtG_pred);
       }

        std::pair <TH1F*, vector<TH1F *>> histos;
        histos = std::make_pair (data_histo, mc_histos);

        return histos;
}



void Fitter::scan_couplings(std::string var_name,  std::pair <TH1F*, vector<TH1F *>> histos){

    if (debug) cout << "Fitter::scan_couplings" << endl;

    if (!histos.first) cout << "data histo not found" << endl;

    TGraphErrors * g = new TGraphErrors();
    TCanvas * c_compare_dists = new TCanvas();
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    
double CtG_vals[11] = {-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0};

    for (int weight = 0 ; weight < histos.second.size() ; weight++){

        if (!histos.second[weight]) cout << "mc histo not found" << endl;
        if (!histos.first) cout << "data histo not found" << endl;

        if (histos.first && histos.second[weight] ){
           // histos.second[weight]->Write();
            histos.second[weight]->SetLineColor(weight);
            
            if (weight ==0) {

                histos.second[weight]->SetStats(kFALSE);
                histos.second[weight]->SetYTitle("#frac{1}{#sigma} \; #frac{#delta(#sigma_tt)}{#delta(x)}");

                histos.second[weight]->Draw("HIST");
                histos.second[weight]->GetYaxis()->SetRangeUser(0.00001, 1.0);

            }else {
            
                histos.second[weight]->Draw("HISTSAME");
            }
        
            //Standard chi2 calculation from ROOT
            //double chi2 = histos.first->Chi2Test(histos.second[weight],"WWCHI2/NDOF");
            double chi2 = this->calculate_test_statistic(histos.first,histos.second[weight]);
            double weight_d = (weight/5.0) * 2.5  ;
            g->SetPoint(weight, CtG_vals[weight], chi2);
        }
    }
    
    if (histos.first) {
        histos.first->SetMarkerSize(1.0);
        histos.first->SetMarkerStyle(22);
        histos.first->SetFillColor(kBlue);
        histos.first->SetFillStyle(3005);
        histos.first->Draw("E2psame");
    }
    
    TLegend *leg = new TLegend(0.5,0.6,0.8,0.8);
    //leg->SetTextFont(65);
    leg->AddEntry( histos.first ,"DATA","E1p");
    leg->AddEntry( histos.second[1] ,"EFT predictions","fl");
    leg->Draw();
    
    c_compare_dists->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
    pad2->SetTopMargin(0);
    pad2->Draw();
    pad2->cd();
    pad2->SetBottomMargin(0.3);
    pad2->SetTopMargin(0.0);
    
    for (int histo = 0; histo < histos.second.size(); histo++){

    TH1F* mc_temp = (TH1F*)histos.second[histo]->Clone("");
    mc_temp->Sumw2();
    mc_temp->SetStats(0);
    mc_temp->Divide(histos.first);
    mc_temp->SetMarkerStyle(21);
    mc_temp->SetMarkerColor(histo);
    mc_temp->SetLineColor(histo);
    mc_temp->GetYaxis()->SetRangeUser(0.3, 1.7);
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
            mc_temp->Draw("HIST");
        }else{
            mc_temp->Draw("HISTSAME");
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
    for ( int histo = 0; histo < histos.second.size(); histo++){
        histos.second[histo]->Write();
    }
    c1->Write();
    c_compare_dists->Write();
    f_out->Close();
    
    
}


double Fitter::calculate_test_statistic(TH1F* data, TH1F* pred){

    std::cout << "\n";
    std::cout << "\n";
    std::cout << "Fitter::calculate_test_statistic \n";
    
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
    
    for (int i=1;i<=nbins; i++){
        double data_bin = data->GetBinContent(i);
        double pred_bin = pred->GetBinContent(i);
        double delta = data_bin - pred_bin;
        deltas.push_back(delta);
        
        if (data_errors_only){
std::cout << "bin  "<< i<<", data_bin "<< data_bin <<", pred_bin = "<< pred_bin <<  " delta =  " << delta << " delta_sq  " << delta*delta <<  ", error  = "<< data->GetBinError(i)  <<  "  running chi2 =  " <<  ((delta)*(delta)) /  data->GetBinError(i)  <<"\n";

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
       // std::cout << "bin  "<< i << " delta =  " << deltas[i] <<  ", error  = "<< errors[i]  <<"\n";
                    double delta_sq = deltas[i-1]*deltas[i-1];
                    double chi2_bin = delta_sq/errors[i-1];
                    chisq += chi2_bin;
                }
            }
        }
    }
    
    std::cout <<" chi2 = "<< chisq << "\n";
    return chisq;
    
    
}




void Fitter::make_summary_plot(vector<TGraphErrors*> scans){
    
    std::cout <<" make_summary_plot "<< "\n";
    
    TFile * all_scans = new TFile("all_scans.root", "RECREATE");

    double x, y, min_y;
    min_y = 99999999.9;
    
    vector <double> min_vals;
    min_vals.clear();
    
    vector<TGraphErrors*> rel_scans;

    //first find min chi2 of each var for sensitivity comparison
    for (int scan = 0; scan< scans.size(); scan ++){
        
        //int n = scans[scan]->GetN(); //get ploted array dimention
        min_y = 99999999.9;

        for(Int_t i=0; i< scans[scan]->GetN(); i++) {
            scans[scan]->GetPoint(i,x,y);
            if (y < min_y) min_y = y;
            std::cout <<" points   = = "<< y <<"\n";
        }
        std::cout <<" "<<"\n";
        std::cout <<" found min chi2,  var  "<< scan <<"  " <<min_y<<"\n";
        std::cout <<" "<<"\n";

        
        min_vals.push_back(min_y);
    }

    
    std::cout <<" found min chi2 "<< "\n";

    
    for (int scan = 0; scan< scans.size(); scan ++){
        TGraphErrors * gr_rel = new TGraphErrors();
        for(Int_t i=0; i < scans[scan]->GetN(); i++) {
            scans[scan]->GetPoint(i,x,y);
            double rel_y  = y / min_vals[scan];
            gr_rel->SetPoint(i,x,rel_y);
        }
        rel_scans.push_back(gr_rel);

    }
    
    
    std::cout <<" make vec of min chi2 "<< "\n";

        TCanvas * all_relscans_c = new TCanvas();
    TLegend *leg_rel = new TLegend(0.5,0.6,0.8,0.8);
    
    for (int scan = 0; scan< rel_scans.size(); scan ++){
        rel_scans[scan]->SetLineColor(scan+1);
        rel_scans[scan]->SetLineWidth(2);
    
        if (scan ==0){
            rel_scans[scan]->Draw("AC");
        }else{
            rel_scans[scan]->Draw("C");
        }
        leg_rel->AddEntry( rel_scans[scan], obs_names[scan].c_str(),"l");

    
    }
    leg_rel->Draw();
    
    std::cout <<" all relscans drawn "<< "\n";
    all_relscans_c->SaveAs("all_relscans.pdf");
    

    
    
    TCanvas * allscans_c = new TCanvas();
    TLegend *leg = new TLegend(0.5,0.6,0.8,0.8);


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
    
    scans[0]->GetHistogram()->GetYaxis()->SetRangeUser(0.0001, 1.0);

    
    leg->Draw();
    allscans_c->SetLogy();

    allscans_c->SaveAs("all_scans.pdf");
    
    allscans_c->Write();
    all_relscans_c->Write();

    all_scans->Close();
    
}


