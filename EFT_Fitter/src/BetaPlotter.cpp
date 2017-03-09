//
//  BetaPlotter.cpp
//  EFT_Fitter
//
//  Created by James Keaveney on 07/03/2017.
//  Copyright © 2017 James Keaveney. All rights reserved.
//

#include "BetaPlotter.hpp"


int main(int argc, const char * argv[]) {
    
    std::cout << "Beta plotter... "<< std::endl;
    
    double lambda = 1000.00;
    double beta_1 = 239.5;
    double beta_2 = 43.2;
    double sm = 839.0;

    TGraph * g_sm = new TGraph();
    TGraph * g_beta1 = new TGraph();
    TGraph * g_beta2 = new TGraph();


    for (double i = -100; i < 100; i = i + 1 ){
    
     double CtG = i*(10000);
        
     double cross_section_b1 =  ((CtG / (lambda*lambda)) * beta_1);
     double cross_section_b2 =  (        pow (  (CtG / (lambda*lambda)), 2 )     )*(beta_2);
        
        std::cout << "cross_section_b1 =  "<< cross_section_b1 << std::endl;
        std::cout << "cross_section_b2 =  "<< cross_section_b2 << std::endl;

        
     g_sm->SetPoint(i, CtG, sm );
     g_beta1->SetPoint(i, (CtG/lambda*lambda), cross_section_b1 );
     g_beta2->SetPoint(i, (CtG/lambda*lambda), cross_section_b2 );
    
    }
    
    TFile * f = new TFile("Beta.root", "RECREATE");

    TCanvas * c1 = new TCanvas();
    
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(g_sm);
    mg->Add(g_beta1);
    mg->Add(g_beta2);
    mg->Draw("AL");
    mg->GetXaxis()->SetTitle("CtG / #Lambda^{2}");
    mg->GetYaxis()->SetTitle("#sigma_{tt}");


    c1->SetLogy();

    c1->SaveAs("BetaPlot.pdf");
    c1->Write();
    f->Write();


}
