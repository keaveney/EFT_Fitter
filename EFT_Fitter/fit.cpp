//
//  main.cpp
//  EFT_Fitter
//
//  Created by James Keaveney on 19/01/2017.
//  Copyright © 2017 James Keaveney. All rights reserved.
//

#include "fit.h"

using namespace std;



int main(int argc, const char * argv[]) {
    // insert code here...

    cout << "Running fit\n";
    
    Fitter *f ;
    
    vector<TH1F *> mc_histos;
    TFile * f_data = new TFile("files/EFT_tt_rivet_small0.root");
    TH1F * data_histo;
    data_histo =  (TH1F*)f_data->Get("CMS_LesHouches2015/ttbar_mass");
    
    mc_histos = f->initialise();
    f->scan_couplings( mc_histos, data_histo) ;

    return 0;
}

vector <TH1F*> Fitter::initialise(){
    vector<TH1F *> mc_histos;

    std::cout << "Initialising fitter\n";

    string var_name = "CMS_LesHouches2015/ttbar_mass";
    stringstream ss;
    
    for (int file = 0 ; file < 11; file++){
        
        ss.str("");
        ss << file;
        string file_s = ss.str();
        
        string filename = "files/EFT_tt_rivet_small" + file_s + ".root";
        cout << filename<< "\n";
        
        TFile * f = new TFile(filename.c_str());
        
        TH1F * mc_histo = (TH1F*)f->Get(var_name.c_str());
        
        mc_histos.push_back(mc_histo);
}
    
    return mc_histos;

}



void Fitter::scan_couplings(vector<TH1F *> mc_histos,TH1F * data_histo ){

    TGraphErrors * g = new TGraphErrors();

    
    std::cout << "Scanning couplings with "<<   mc_histos.size() <<" weights \n";
    
    for (int weight = 0 ; weight < mc_histos.size() ; weight++){
        
        if (data_histo && mc_histos[weight] ){
            double chi2 = data_histo->Chi2Test(mc_histos[weight],"WWCHI2/NDOF");
            std::cout << "Weight = "<< weight <<"  chi2 = "<< chi2 << "\n";
            std::cout << "Weight = "<< "\n";
            double weight_d = weight;
            g->SetPoint(weight, weight_d, chi2);
            
        }
        
    }
    
    TCanvas * c1 = new TCanvas();
    g->Draw("AC");
    c1->SaveAs("chi2scan.png");
    
}




double Fitter::calculate_test_statistic(TH1F* h_data, TH1F* h_mc){
    
    std::cout << "Fitter::calculate_test_statistic \n";
    
    double test_statistic = 0;
    
    return test_statistic;
}
