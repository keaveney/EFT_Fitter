//This header file contains a few simple tools that were needed for the EFT analysis
//but ended up being useful for the charge asymmetry analysis and creation
//of pval tables and plots for the TOP-17-014 paper.
//Therefore they are now decoupled from those class and live here
// as standalone methods that can be intgreated into a similar analysis.

#include <iostream>
#include <sstream>
#include <fstream>

#include "TMatrixD.h"

#include "TRandom.h"
#include "TRandom3.h"

TH2D* make_covariance_matrix(std::string, std::string);
TH1* make_poisson_toy(TH1*,TH1*, int, double);
double calculate_test_statistic(TH1F*, TH1F*, TGraphAsymmErrors*,std::string);


double calculate_test_statistic(TH1F* data, TH1F* pred, TGraphAsymmErrors * gr_errors, std::string cov_file){
    
    double test_statistic = 0.0;
    int nbins  = data->GetSize() - 2;
    bool schmitt_fit = true;
    std::vector<double> deltas;
    std::vector<double> errors;
    
    double corr_coff = 1.0;
    bool data_errors_only = true;
    double chisq = 0.0;
    double chi2_bin = 0.0;
    double theory_error =0.0;
    double total_error =0.0;
    TH2D * cov;
    
    //Retrieve covariance matrix from text file
    //TFile *cov_file = new TFile("HypLLBarDPhi_totCovMtrxFile.root");
    //TH2D * cov = (TH2D*)cov_file->Get("inv_cov");
    cov = make_covariance_matrix(cov_file.c_str());
    
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



TH2D* make_covariance_matrix(std::string filename, std::string data_file){
    
    std::ifstream readFile(filename.c_str());
    
    std::vector<std::vector<std::string>> matrix;
    std::vector<std::string> tokens;
    
    std::string line;
    std::string delimiter = " ";
    
    while(getline(readFile,line)){
        tokens.clear();
        std::istringstream iss(line);
        
        //copy(istream_iterator<string>(iss),istream_iterator<string>(), back_inserter(tokens));
        size_t pos = 0;
        std::string token;
        
        //std::cout <<" " << std::endl;
        //std::cout <<"line: " << std::endl;
        
        while ((pos = line.find(delimiter)) != std::string::npos){
            token = line.substr(0, pos);
            //   std::cout <<" mini-loop token "<< token << std::endl;
            line.erase(0, pos + delimiter.length());
            tokens.push_back(token);
        }
        //std::cout <<" mini-loop token "<< line << std::endl;
        
        tokens.push_back(line);
        
        matrix.push_back(tokens);
        // cout <<"line-- = "<< line << " " <<  matrix[0][0] <<endl;
    }
    
    readFile.close();
    
    std::cout <<" rows "<< matrix.size() <<"  columns  "<< matrix[0].size() <<std::endl;
    
    TMatrixD mat =  TMatrixD(matrix.size(), matrix.size());
    
    TFile * f_data_file = new TFile(data_file.c_str());
    TGraphAsymmErrors * g_data = (TGraphAsymmErrors*)f_data_file->Get("data");
    

    //make parameters of mth2d auto from graph
    TH2D* cov = new TH2D("cov","cov", matrix.size() , 0.0, 3.142, matrix.size() , 0.0, 3.142);
    TH2D* inv_cov = new TH2D("inv_cov","inv_cov", matrix.size() , 0.0, 3.142, matrix.size() , 0.0, 3.142);
    
    double point_i_x,point_i_y, point_j_x, point_j_y;
    
    
    for (int x = 0; x < matrix.size(); x++){
        for (int y = 0; y < matrix[x].size(); y++){
            std::cout <<" xm "<< matrix[x][y] <<std::endl;

            std::string::size_type sz;
            std::string cov_elem_str = matrix[x][y];
            double cov_elem = std::stod(cov_elem_str, &sz);
            
            g_data->GetPoint(x, point_i_x, point_i_y );
            g_data->GetPoint(y, point_j_x, point_j_y );
            cov_elem = (cov_elem)*(point_i_y)*(point_j_y);
            
           std::cout <<" x = "<< x <<"  y "<< y <<"  "<<  cov_elem_str <<std::endl;
            cov->SetBinContent(x+1, y+1, cov_elem);
            mat[x][y] = cov_elem;
        }
    }
    
    std::string rootfile_name = filename.substr(0, filename.length() - 3) + "root";
    
    //invert convariance matrix and write
    Double_t det2;
    mat.Invert(&det2);
    
    for (int x = 0; x < matrix.size(); x++){
        for (int y = 0; y < matrix[x].size(); y++){
            inv_cov->SetBinContent(x+1, y+1, mat[x][y]);
        }
    }
    
    
    TFile * file = new TFile(rootfile_name.c_str(), "RECREATE");
    cov->Write();
    inv_cov->Write();
    mat.Write();
    
    file->Close();
    
    return cov;
    
}


TH1* make_poisson_toy(TH1* hh, TH1* data, int nevents, double data_integral){
    
    //cout <<"making poisson toy " << endl;
    
    std::string opt = "manual";
    
    //double n_toys = 523850.625;
    
    //This function make a Poisson toy for a given prediction at particle level
    // corresponding to the same number of events as observed in data.
    
    //Histos to contain raw distributions for prediction and data.
    
    // TH1F* hh_t = new TH1F("","", hh->GetNbinsX() , 0.0, 3.15);
    // TH1F* hh_data = new TH1F("","", hh->GetNbinsX() , 0.0, 3.15);
    
    //TH1F* hh_t = (TH1F*)hh->Clone();
    //TH1F* hh_data = (TH1F*)hh->Clone();
    
    TH1 * hh_f;
    hh_f =(TH1F*)hh->Clone();
    hh_f->Reset();
    
    
    //hh_t->Reset();
    //hh_data->Reset();
    double bin_height_pred;
    double bin_height_data;
    TRandom3 * gRand = new TRandom3();

    
    
    
    if (opt == "manual"){
        for (int bin = 1; bin <= hh->GetNbinsX(); bin++){
            double mean = hh->GetBinContent(bin);
            double width = data->GetBinError(bin);
            double toy_bin_height = gRand->Gaus(mean,width);
            double toy_bin_error = width;
            //  cout <<"TOY =  "<< "mean " << mean <<" width  "<<width<< " toy height " <<  toy_bin_height <<endl;
            hh_f->SetBinContent(bin, toy_bin_height);
            hh_f->SetBinError(bin, toy_bin_error);
            //hh_t->SetBinContent(bin, bin_height_pred);
            //hh_data->SetBinContent(bin, bin_height_data);
        }
    }
    else if (opt=="roofit"){
        //Multiply back by bin widths to get cross sections as bin heights
        for (int bin = 1; bin <= hh->GetNbinsX(); bin++){
            if (bin == hh->GetNbinsX() ){
                bin_height_pred = (0.24)*hh->GetBinContent(bin);
                bin_height_data = (0.24)*data->GetBinContent(bin);
            }else{
                bin_height_pred = (0.25)*hh->GetBinContent(bin);
                bin_height_data = (0.25)*data->GetBinContent(bin);
            }
            //hh_t->SetBinContent(bin, bin_height_pred);
            //     hh_data->SetBinContent(bin, bin_height_data);
        }
        
        
        //now scale by integrated lumi
        //double scaling = hh->Integral() / hh_t->Integral();
        
        //    double int_lumi = 37500000.0;
        //double int_lumi = 3750.0;
        
        double int_lumi = 6000.0;
        
        //hh_t->Scale(int_lumi);
        
        //TH1F* h_forPDF = new TH1F("","", hh->GetNbinsX() , -1.0,1.0);
        TH1F* h_forPDF = (TH1F*)hh->Clone();
        h_forPDF->Reset();
        
        
        //fill histo with raw event counts for this prediction
        // for (int bin = 1; bin <= hh->GetNbinsX(); bin++){
        //     h_forPDF->SetBinContent(bin, int_lumi *  hh_t->GetBinContent(bin));
        // }
        
        
        //some neccessary Roofit declarations
        RooRealVar x("x","x",0,3.15);
        x.setBins(hh->GetNbinsX());
        
        
        //convert raw prediction into pdf
        RooDataHist dh("hh_t","hh_t",x,RooFit::Import(*h_forPDF)) ;
        RooHistPdf histpdf1("histpdf1","histpdf1",x,dh,0) ;
        
        //generate toy event count
        TRandom * r = new TRandom();
        double n_expected = r->Poisson(h_forPDF->Integral());
        
        //generate binned toys according the PDF with Poisson errors
        RooDataHist* toy_data = histpdf1.generateBinned(x,n_expected);
        
        //convert it back to a TH1
        hh_f = toy_data->createHistogram("toy", x);
        
        
        double bbin;
        double bbin_error;
        
        //convert back to 'cross section result' histo
        //divide by bin width
        for (int bin = 1; bin <= hh->GetNbinsX(); bin++){
            
            if (bin == hh->GetNbinsX() ){
                bbin = hh_f->GetBinContent(bin)/(0.24);
                // bbin_error = hh_f->GetBinError(bin)/(0.24);
                bbin_error = bbin *(0.01);
                
            }else{
                bbin = hh_f->GetBinContent(bin)/(0.25);
                // bbin_error = hh_f->GetBinError(bin)/(0.25);
                bbin_error = bbin *(0.01);
                
            }
            
            hh_f->SetBinContent(bin,bbin);
            hh_f->SetBinError(bin,bbin_error);
        }
        
        //divide by lumi
        hh_f->Scale(1.0/int_lumi);
        
        
        
        // double scaling_2 = data_integral/hh_f->Integral();
        double scaling_2 = hh->Integral()/hh_f->Integral();
        
        //hh_f->Scale(scaling_2);
        
    }
    //cout <<" toy integral  = "<<  hh_f->Integral()<< endl;
    
    return hh_f;
}


