//
//  TopMassExtractor.cpp
//  EFT_Fitter
//
//  Created by James Keaveney on 08/03/2017.
//  Copyright © 2017 James Keaveney. All rights reserved.
//

#include "TopMassExtractor.hpp"


int main(int argc, const char * argv[]) {

    std::cout << "Top Mass Extractor... "<< std::endl;


    
    TFile * f_data = new TFile("files/DiffXS_HypToppT_source_parton.root");
    TGraphAsymmErrors * g_data = (TGraphAsymmErrors*)f_data->Get("data");
    f_data->Close();
    
    TFile * f_difftop_scan = new TFile("Difftop_scan.root", "RECREATE");
    
    TH1F * h;
    vector<std::string> mass_points = {"170", "172", "173.3", "174", "176"};
    vector<TH1F*> histos;

    TGraphAsymmErrors * g_chi2_scan = new TGraphAsymmErrors(mass_points.size());

    
    for (int pred = 0; pred < mass_points.size(); pred++){
    
        std::string pred_name = "difftop_predictions/" + mass_points[pred] + ".txt";
        h = make_histo(pred_name);
        histos.push_back(h);
    
        std::cout << "Writing histo... "<< std::endl;

        h->SetName(pred_name.c_str());
        h->Write();
    }
    f_difftop_scan->Write();


    for( int pred = 0; pred < mass_points.size(); pred++){
    
        
        double mass_val = atof(mass_points[pred].c_str());
        double running_chi2 = calc_chi2(histos[pred], g_data) / 5.0;
    
        g_chi2_scan->SetPoint(pred, mass_val, running_chi2);

         std::cout << "Running chi2... "<< running_chi2 <<  std::endl;
        
    }
    g_chi2_scan->GetXaxis()->SetTitle("M^{pole}_{t} (GeV)");
    g_chi2_scan->GetYaxis()->SetTitle("#chi^{2} DATA-THEORY");

    TCanvas * c1 = new TCanvas();
    g_chi2_scan->Draw("AC*");
    g_chi2_scan->Write();
    c1->SaveAs("Chi2Scan.pdf");

}








