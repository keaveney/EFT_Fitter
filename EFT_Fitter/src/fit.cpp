//
//  fit.cpp
//  EFT_Fitter
//
//  Created by James Keaveney on 19/01/2017.
//  Copyright © 2017 James Keaveney. All rights reserved.
//

#include "fit.h"


using namespace std;


int main(int argc, const char * argv[]) {
    
    cout << "Running fit\n";
    
    Fitter *f_EFT;
    
    vector<TH1F *> mc_histos;
    TFile * f_data;
    TGraphAsymmErrors * g_data;
    TH1F* data_histo;
    
    double bin_centre, bin_height, running_total =0.0;

    string mode = "data";
    
    if (mode == "closuretest"){
    
        f_data = new TFile("files/EFT_tt_rivet_small4.root");
        data_histo = (TH1F*)f_data->Get("CMS_dilepton_diff/ttbar_mass");
        
        
    }else if (mode == "data") {
    
     f_data = new TFile("files/DiffXS_HypTTBarMass_source.root");
   // f_data = new TFile("files/DiffXS_HypToppT_source.root");

     g_data = (TGraphAsymmErrors*) f_data->Get("data");
        
    Float_t bins[] = { 0.0, 380.0, 470.0, 620.0, 820.0, 1100.0, 1600.0 };
    Int_t  binnum = sizeof(bins)/sizeof(Float_t) - 1;
    data_histo = new TH1F("n","t", binnum, bins);
    
    
    if (!data_histo) cout << "data histo not found" << endl;
    
    for (int point = 0; point < 7; point ++){
        g_data->GetPoint(point, bin_centre, bin_height);
        data_histo->SetBinContent(point+1, bin_height);
        running_total = running_total + bin_height;
    }
    }
    
    std::cout << "Data histo has integral " << running_total;


    mc_histos = f_EFT->initialise(running_total);

    f_EFT->scan_couplings( mc_histos, data_histo) ;

    return 0;
}

vector <TH1F*> Fitter::initialise(double integral){
    vector<TH1F *> mc_histos;

    std::cout << "Initialise... \n";

   string var_name = "CMS_dilepton_diff/ttbar_mass";
    // string var_name = "CMS_dilepton_diff/t_pT";

    stringstream ss;
    
    double CtG_vals[11] = {  -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0  };

    
    //first extract raw histos
        string filename_neg2 = "files/EFTNLO_CtGNeg2.root";
        string filename_pos2 = "files/EFTNLO_CtGPos2.root";
        string filename_0 = "files/EFTNLO_CtG0.root";
        
        TFile * f_neg2 = new TFile(filename_neg2.c_str());
        TFile * f_0 = new TFile(filename_0.c_str());
        TFile * f_pos2 = new TFile(filename_pos2.c_str());
        
        TH1F * mc_histo_neg2 = (TH1F*)f_neg2->Get(var_name.c_str());
        TH1F * mc_histo_0 = (TH1F*)f_0->Get(var_name.c_str());
        TH1F * mc_histo_pos2 = (TH1F*)f_pos2->Get(var_name.c_str());

        if (!mc_histo_neg2) cout << "mc histo: "<<  filename_neg2  <<" not found" << endl;
    
   //make histo of pure CtG contribution
       TH1F * h_pure_Ctg = (TH1F*)mc_histo_pos2->Clone();
       h_pure_Ctg->Add(mc_histo_neg2, -1);
 
    
    TFile * f_basis = new TFile("basis_histos.root", "RECREATE");
    mc_histo_neg2->SetName("Neg2");
    mc_histo_neg2->Write();
    mc_histo_0->SetName("0");
    mc_histo_0->Write();
    mc_histo_pos2->SetName("Pos2");
    mc_histo_pos2->Write();
    h_pure_Ctg->Write();
    f_basis->Close();
    
    
    //loop over chosen CtG values and make prediction
       for (int scale = 0 ; scale < 11; scale++){
          double CtG = CtG_vals[scale];
           TH1F * h_CtG_pred = (TH1F*)mc_histo_0->Clone();
           h_CtG_pred->Add(h_pure_Ctg, CtG/2.0);
           
           double scaling = integral/h_CtG_pred->Integral();
           h_CtG_pred->Scale(scaling);
           mc_histos.push_back(h_CtG_pred);

       }
    
    
    

    
    return mc_histos;

}



void Fitter::scan_couplings(vector<TH1F *> mc_histos,TH1F * data_histo ){
    
    
    if (!data_histo) cout << "data histo not found" << endl;


    TGraphErrors * g = new TGraphErrors();

    TCanvas * c_compare_dists = new TCanvas();
  
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
    pad1->SetBottomMargin(0);
    pad1->Draw();
    pad1->cd();
    
    std::cout << "Scanning couplings with "<<   mc_histos.size() <<" weights \n";
    
    double CtG_vals[11] = {  -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0  };

    TFile * f_preds = new TFile("pred_histos.root", "RECREATE");

    for (int weight = 0 ; weight < mc_histos.size() ; weight++){
        
        std::cout << "Looping histos...1 \n";

        if (!mc_histos[weight]) cout << "mc histo not found" << endl;
        if (!data_histo) cout << "data histo not found" << endl;



        if (data_histo && mc_histos[weight] ){
            
            std::cout << "got histos \n";
            mc_histos[weight]->Write();
            

            mc_histos[weight]->SetLineColor(weight);
            
            if (weight ==0) {
                std::cout << "first hsto 1 \n";

                mc_histos[weight]->SetStats(kFALSE);
                mc_histos[weight]->SetYTitle("A.U.");
                mc_histos[weight]->Draw("HIST");
                mc_histos[weight]->GetYaxis()->SetRangeUser(0.00001, 0.01);

                std::cout << "first hsto N\n";

            }else {
            
                mc_histos[weight]->Draw("HISTSAME");
            }
            std::cout << "gcalculating chi2t ....\n";

            
            double chi2 = data_histo->Chi2Test(mc_histos[weight],"WWCHI2/NDOF");
            std::cout << "Weight = "<< weight <<" CtG "<<  CtG_vals[weight]  <<"  chi2 = "<< chi2 << "\n";
            std::cout << "Weight = "<< "\n";
            double weight_d = (weight/5.0) * 2.5  ;
            g->SetPoint(weight, CtG_vals[weight], chi2);
            std::cout << "gsetting point ....\n";

            
        }
        
        
    }
    
    f_preds->Close();
    
     cout << "DRAWING" << endl;

    
    if (data_histo) {
        data_histo->SetMarkerSize(1.0);
        data_histo->SetMarkerStyle(22);
        data_histo->Draw("psame");
    }

    
    
    
    TLegend *leg = new TLegend(0.5,0.6,0.8,0.8);
    //leg->SetTextFont(65);
    leg->AddEntry( data_histo ,"DATA","E1p");
    leg->AddEntry( mc_histos[9] ,"EFT tt predictions","l");

   // TLegendEntry* l2 = leg->AddEntry("TOP_16_008_NNLO_DATA.yoda","lepton + jets (arXiv:1610.04191, sub. to PRD)","E1p");
   // l2->SetMarkerColor(4);
   // l2->SetMarkerStyle(23);
   // l2->SetLineColor(4);
    //TLegendEntry* l3 = leg->AddEntry("fit","fit","lf");
   // leg->SetBorderSize(0.0);
    leg->Draw();
    
    
    
    c_compare_dists->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
    pad2->SetTopMargin(0);
    pad2->Draw();
    pad2->cd();
    pad2->SetBottomMargin(0.3);
    pad2->SetTopMargin(0.0);
    
    for (int histo = 0; histo < mc_histos.size(); histo++){

    TH1F* mc_temp = (TH1F*)mc_histos[histo]->Clone("");
    mc_temp->Sumw2();
    mc_temp->SetStats(0);
    mc_temp->Divide(data_histo);
    mc_temp->SetMarkerStyle(21);
    mc_temp->SetMarkerColor(histo);
    mc_temp->SetLineColor(histo);
    mc_temp->GetYaxis()->SetRangeUser(-0.0, 5.0);
    mc_temp->GetYaxis()->SetLabelSize(0.1);
    mc_temp->GetXaxis()->SetLabelSize(0.1);
    mc_temp->GetXaxis()->SetTitleSize(0.16);
    mc_temp->GetXaxis()->SetTitleOffset(0.8);
    mc_temp->GetYaxis()->SetTitleOffset(0.2);
    mc_temp->GetYaxis()->SetTitleSize(0.12);

    mc_temp->SetYTitle("THEORY/DATA");
    mc_temp->SetXTitle("M_{tt} GeV");

        if(histo == 0){
            mc_temp->Draw("HIST");
        }else{
            mc_temp->Draw("HISTSAME");
        }
    }
    
    
    c_compare_dists->cd();
    
    
    
   pad1->SetLogy();
   c_compare_dists->SaveAs("compare_dists.pdf");
    
    TCanvas * c1 = new TCanvas();
   // g->SetMinimum(0.002);
   // g->SetMaximum(0.0026);
    g->GetHistogram()->GetXaxis()->SetTitle("CtG");
    g->GetHistogram()->GetYaxis()->SetTitle("#chi^{2} (DATA-THEORY)");
    g->GetHistogram()->GetXaxis()->SetRangeUser(-5.0, 5.0);


    g->Draw("PAC");
    c1->SetLogy();
    c1->SaveAs("chi2scan.png");
    
    TFile * f_out = new TFile("scan_result.root", "RECREATE");
    data_histo->Write();
    
    for ( int histo = 0; histo < mc_histos.size(); histo++){
        mc_histos[histo]->Write();
        std::cout << "writing mc histo\n";
        
    }
    c1->Write();
    c_compare_dists->Write();
    
    f_out->Close();
    
    
    
}




double Fitter::calculate_test_statistic(TH1F* h_data, TH1F* h_mc){
    
    std::cout << "Fitter::calculate_test_statistic \n";
    
    double test_statistic = 0;
    
    return test_statistic;
}
