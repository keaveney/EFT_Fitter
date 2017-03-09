//
//  TopMassExtractor.hpp
//  EFT_Fitter
//
//  Created by James Keaveney on 08/03/2017.
//  Copyright © 2017 James Keaveney. All rights reserved.
//

#ifndef TopMassExtractor_hpp
#define TopMassExtractor_hpp

#include <stdio.h>

#include <iostream>
#include <sstream>
#include <stdio.h>
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMultiGraph.h"
#include "TVectorF.h"
#include <fstream>
#include <stdlib.h>     /* atof */
#include <math.h>       /* exp */




using namespace std;

TH1F * h;

TH1F * make_histo(std::string filename){
    
    cout <<"Making histo from prediction" << filename <<endl;
    
    int nbins = 6;

    
    TVectorF vx(nbins), vy(nbins), vexl(nbins), vexh(nbins), veyl(nbins), veyh(nbins), clhi(nbins), cllo(nbins), line_y(nbins);
    Double_t bin_edges[nbins+1];
    
    vector<string> tokens;
    string line;
    ifstream myfile (filename);
    int row = 0;
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
            tokens.clear();
            istringstream iss(line);
            copy(istream_iterator<string>(iss),
                 istream_iterator<string>(),
                 back_inserter(tokens));
            
            cout <<"N tokens = " << tokens.size() << endl;
            if (tokens.size() == 13){
                vexl[row] = stof(tokens[0]);
                vexh[row] = stof(tokens[1]);
                vx[row] =   (vexl[row] + vexh[row])/2.0 ;

                bin_edges[row] = stof(tokens[0]);
                
cout << "Bin " << row <<" , bin centre  "<< vx[row]  <<"bin low edge "<<  vexl[row] <<" bin high edge " << vexh[row] << endl;
                
                //  vexl[row] = 0.0;
                //  vexh[row] = 0.0;
                
                vy[row] = stof(tokens[7]);
                
                //veyl[row] = stof(tokens[3])*0.02;
                //veyh[row] = stof(tokens[3])*0.02;
                
                veyl[row] = stof(tokens[9]);
                veyh[row] = stof(tokens[8]);
                
                // cout <<"VX = "<<  vx[row]  <<" VY "<<  vy[row] <<  " VEYL  " << veyl[row]  << " VEYH  " << veyh[row]<<endl;
                
                row++;
            }
        }
        myfile.close();
        bin_edges[nbins] = (vexh[nbins-1]);
        
        cout <<"closed file"<<endl;
        
    }
    else         cout <<"unable to open file"<<endl;
    
    
    h = new TH1F("","",6,bin_edges);
    
    bool correct_pred = true;
    
    for (int i =0; i < nbins; i++){
        
        double bin_content = vy[i];
        double bin_error   = veyl[i];
        
        double correction = exp ( 0.0416   - (0.0003*vx[i])   );
       if (correct_pred) bin_content = correction*bin_content;
        
        cout <<"Correction "<<  correction << endl;
        h->SetBinContent(i+1, bin_content );
        h->SetBinError(i+1, bin_error);
    }
    return h;

    
}

double calc_chi2(TH1F* pred, TGraphAsymmErrors * data){
    
    std::cout << "Calculating chi2... "<< std::endl;
    
    double x_data, y_data, x_data_unc, y_data_unc;
    double running_chi2= 0.0;
    double bin_unc, final_unc;
    
    for (int bin = 0; bin < 6 ; bin++){
    
        double bin_pred = pred->GetBinContent(bin+1);
        double bin_error = pred->GetBinError(bin+1);

        data->GetPoint(bin, x_data, y_data);
        
        double sq_diff = pow( ( bin_pred - y_data   )   ,   2);
        
        //calculating appropriate data uncertainty
        if (bin_pred > y_data){
        
            bin_unc = data->GetErrorYhigh(bin);
        }else{
        bin_unc = data->GetErrorYlow(bin);
        
        }
        
        //DATA & THEORY in quadrature
        //final_unc =  pow (pow ( bin_unc,2)  + pow ( bin_error,2) ,0.5)  ;

        //DATA only
        final_unc =   bin_unc;
                
        running_chi2 = running_chi2 + (    sq_diff / (final_unc*final_unc));
    
        cout <<"Bin = "<< bin  <<" bin_pred  "<< bin_pred << " bin_data "  << y_data << " sq_diff " <<  sq_diff <<"  bin_unc"<< bin_unc<<endl;

    }
    running_chi2 = running_chi2/5.0;
    
    return running_chi2;

    
}


#endif /* TopMassExtractor_hpp */
