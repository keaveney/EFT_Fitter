//make pvals
// script to write tables of p-values, chi2/ndof between// differential
// cross section results and theory/MC predictions.

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TMath.h"

#include "TLatex.h"

#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TF1.h"
#include <utility>
#include <tuple>

using namespace std;

void write_latex(string, string, vector<string>, float[] , int[], float[]);
std::tuple <float, int, float > process_result(string,string);

TH1F * h_model;

ofstream myfile;


int main(int argc, const char * argv[]){
    
    string mode = "norm_parton";
    vector<string> filenames;

    myfile.precision(2);
    string filename = mode + ".tex";
    myfile.open (filename);
    myfile << "\\begin{table}"<<endl;
    myfile << "\\centering"<<endl;
    myfile << "\\makebox[0pt][c]{\\parbox{1.0\\textwidth}{%"<<endl;
    
    vector<string> modelnames = {
        "\\Powheg+\\Pythia",
       // "\\MGaMCatNLO+\\Pythia"
        "\\Powheg+\\Herwigpp"

    };

    if (mode == "norm_parton"){
    filenames = {
        "files/July3/norm_parton/DiffXS_HypToppT_source.root",
        "files/July3/norm_parton/DiffXS_HypAntiToppT_source.root",
        "files/July3/norm_parton/DiffXS_HypToppTLead_source.root",
        "files/July3/norm_parton/DiffXS_HypToppTNLead_source.root",
        "files/July3/norm_parton/DiffXS_HypToppTTTRestFrame_source.root",
        "files/July3/norm_parton/DiffXS_HypAntiToppTTTRestFrame_source.root",
        "files/July3/norm_parton/DiffXS_HypTopRapidity_source.root",
        "files/July3/norm_parton/DiffXS_HypAntiTopRapidity_source.root",
        "files/July3/norm_parton/DiffXS_HypTopRapidityLead_source.root",
        "files/July3/norm_parton/DiffXS_HypTopRapidityNLead_source.root",
        "files/July3/norm_parton/DiffXS_HypTTBarpT_source.root",
        "files/July3/norm_parton/DiffXS_HypTTBarRapidity_source.root",
        "files/July3/norm_parton/DiffXS_HypTTBarMass_source.root",
        "files/July3/norm_parton/DiffXS_HypTTBarDeltaRapidity_source.root",
        "files/July3/norm_parton/DiffXS_HypTTBarDeltaPhi_source.root"
    };
    } else if (mode == "abs_parton"){
        filenames = {

        "files/July3/abs_parton/DiffXS_HypToppT_source.root",
        "files/July3/abs_parton/DiffXS_HypAntiToppT_source.root",
        "files/July3/abs_parton/DiffXS_HypToppTLead_source.root",
        "files/July3/abs_parton/DiffXS_HypToppTNLead_source.root",
        "files/July3/abs_parton/DiffXS_HypToppTTTRestFrame_source.root",
        "files/July3/abs_parton/DiffXS_HypAntiToppTTTRestFrame_source.root",
        "files/July3/abs_parton/DiffXS_HypTopRapidity_source.root",
        "files/July3/abs_parton/DiffXS_HypAntiTopRapidity_source.root",
        "files/July3/abs_parton/DiffXS_HypTopRapidityLead_source.root",
        "files/July3/abs_parton/DiffXS_HypTopRapidityNLead_source.root",
        "files/July3/abs_parton/DiffXS_HypTTBarpT_source.root",
        "files/July3/abs_parton/DiffXS_HypTTBarRapidity_source.root",
        "files/July3/abs_parton/DiffXS_HypTTBarMass_source.root",
        "files/July3/abs_parton/DiffXS_HypTTBarDeltaRapidity_source.root",
        "files/July3/abs_parton/DiffXS_HypTTBarDeltaPhi_source.root"
        };
    }else if (mode == "norm_particle"){
        filenames = {
            "files/July3/norm_particle/DiffXS_HypToppT_source.root",
            "files/July3/norm_particle/DiffXS_HypAntiToppT_source.root",
            "files/July3/norm_particle/DiffXS_HypToppTLead_source.root",
            "files/July3/norm_particle/DiffXS_HypToppTNLead_source.root",
            "files/July3/norm_particle/DiffXS_HypToppTTTRestFrame_source.root",
            "files/July3/norm_particle/DiffXS_HypAntiToppTTTRestFrame_source.root",
            "files/July3/norm_particle/DiffXS_HypTopRapidity_source.root",
            "files/July3/norm_particle/DiffXS_HypAntiTopRapidity_source.root",
            "files/July3/norm_particle/DiffXS_HypTopRapidityLead_source.root",
            "files/July3/norm_particle/DiffXS_HypTopRapidityNLead_source.root",
            "files/July3/norm_particle/DiffXS_HypTTBarpT_source.root",
            "files/July3/norm_particle/DiffXS_HypTTBarRapidity_source.root",
            "files/July3/norm_particle/DiffXS_HypTTBarMass_source.root",
            "files/July3/norm_particle/DiffXS_HypTTBarDeltaRapidity_source.root",
            "files/July3/norm_particle/DiffXS_HypTTBarDeltaPhi_source.root"
        };
    }
    
        else if (mode == "abs_particle"){
            filenames = {
                "files/July3/abs_particle/DiffXS_HypToppT_source.root",
                "files/July3/abs_particle/DiffXS_HypAntiToppT_source.root",
                "files/July3/abs_particle/DiffXS_HypToppTLead_source.root",
                "files/July3/abs_particle/DiffXS_HypToppTNLead_source.root",
                "files/July3/abs_particle/DiffXS_HypToppTTTRestFrame_source.root",
                "files/July3/abs_particle/DiffXS_HypAntiToppTTTRestFrame_source.root",
                "files/July3/abs_particle/DiffXS_HypTopRapidity_source.root",
                "files/July3/abs_particle/DiffXS_HypAntiTopRapidity_source.root",
                "files/July3/abs_particle/DiffXS_HypTopRapidityLead_source.root",
                "files/July3/abs_particle/DiffXS_HypTopRapidityNLead_source.root",
                "files/July3/abs_particle/DiffXS_HypTTBarpT_source.root",
                "files/July3/abs_particle/DiffXS_HypTTBarRapidity_source.root",
                "files/July3/abs_particle/DiffXS_HypTTBarMass_source.root",
                "files/July3/abs_particle/DiffXS_HypTTBarDeltaRapidity_source.root",
                "files/July3/abs_particle/DiffXS_HypTTBarDeltaPhi_source.root"
            };
        
        }
    

    vector<string> vars = {
        "\\pt top",
        "\\pt anti-top",
        "\\pt top  (leading)",
        "\\pt top  (sub-leading)",
        "\\pt top (\\ttbar\\ rest frame)",
        "\\pt anti-top (\\ttbar\\ rest frame)",
        "Y top",
        "Y anti-top",
        "Y top (leading)",
        "Y top (sub-leading)",
        "\\pt \\ttbar\\",
        "Y \\ttbar\\",
        "M_{\\ttbar}",
        "\\Delta\\ Y \\ttbar\\",
        "\\Delta\\ \\Phi\\ \\ttbar\\"
    };
    
   // vector<  vector<float>  > chisq;
   // vector<  vector<int>    > ndof;
   // vector<  vector<float>  > pval;
    float   chisq[modelnames.size()][vars.size()];
    int     ndof[modelnames.size()][vars.size()];
    float   pval[modelnames.size()][vars.size()];

    
    std::tuple <float, int, float > gof;
    
    for (int m = 0; m < modelnames.size(); m++){
        for (int f = 0; f < filenames.size(); f++){
            gof = process_result(modelnames[m], filenames[f]);
         //   chisq[m].push_back(std::get<0>(gof));
         //   ndof[m].push_back(std::get<1>(gof));
         //   pval[m].push_back(std::get<2>(gof));
            
            chisq[m][f] = std::get<0>(gof);
            ndof[m][f] = std::get<1>(gof);
            pval[m][f] = std::get<2>(gof);

           // cout <<"RETURNED chisq , ndof, pval  =  "<< chisq[m][f] <<"  "<< ndof[m][f] <<"  "<<  pval[m][f]  << endl;
        }

        write_latex(mode, modelnames[m], vars, chisq[m], ndof[m], pval[m]);

    }
    
    myfile << "}}"<<endl;
    myfile << "\\end{table}"<<endl;
    myfile.close();
    
    return 1;

}


std::tuple <float, int, float > process_result(string modelname, string filename){
    
    string modelhistoname, dataname;
    double chisq_running = 0.0, pval = 0.0;
    double bin_x, bin_data, e_data;
    int ndof;
    dataname = "Graph";
    TFile * f_results = new TFile(filename.c_str());
    TCanvas * c = (TCanvas*)f_results->Get("canvas;1");
    
    if (f_results) cout <<" got file"<<" "<< filename <<endl;

    if (c) cout <<" got canvas"<<endl;
    
    if (modelname == "\\Powheg+\\Pythia"){
        modelhistoname = "Nominalplot";
    } else if (modelname == "\\Powheg+\\Herwigpp"){
        modelhistoname = "MCATNLOplot";
    }
    else if (modelname == "\\MGaMCatNLO+\\Pythia"){
        modelhistoname = "MCATNLOplot";
    }
    
    cout <<"modelhist name is "<<modelhistoname<<endl;
    
    //TGraphAsymmErrors * g_data = (TGraphAsymmErrors*)c->GetPrimitive(dataname.c_str());
    TGraphAsymmErrors * g_data = (TGraphAsymmErrors*)f_results->Get("data");

    h_model = (TH1F*)c->GetPrimitive(modelhistoname.c_str());
    
    if (g_data) cout <<" got data graph"<<endl;
    if (h_model) cout <<" got model histo"<<endl;

    
    cout <<"  "<<endl;
    
    
    for (int bin = 0; bin < h_model->GetNbinsX(); bin++){
        
        double bin_model = h_model->GetBinContent(bin+1);
        g_data->GetPoint(bin, bin_x, bin_data);
        e_data = g_data->GetErrorY(bin);
        
        cout <<" data  =  "<<  bin_data << "  data error = "<< e_data  <<"  model " <<  bin_model <<endl;

        //Extact term in covariance matrix here (need extra loop also)
        
        chisq_running += pow(bin_data - bin_model, 2.0) / (  pow(e_data, 2.0 ) ) ;
        
    }
    
    ndof = (h_model->GetNbinsX()); //need to reduce this by 1 for normalised results
    pval = TMath::Prob(chisq_running, ndof);
    
    cout <<"chisq , ndof, pval  =  "<< chisq_running <<"  "<< ndof <<"  "<<  pval << endl;
    
    std::tuple<float, int, float>  results ( chisq_running, ndof, pval);
    
    
    cout <<"returning results" << endl;

    return results;
    
}





void write_latex(string mode, string model, vector<string> vars, float chisq[], int ndof[], float pval[]){

    string mode_string;
    
    if (mode == "norm_parton"){
        mode_string = "normalised, parton-level";
    }else if (mode == "abs_parton"){
        mode_string = "absolute, parton-level";
    }
    else if (mode == "norm_particle"){
        mode_string = "normalised, particle-level";
    }
    else if (mode == "abs_particle"){
        mode_string = "absolute, particle-level";
    }
    
    myfile << "\\begin{minipage}[b]{0.42\\hsize}\\centering \\small"<<endl;
    myfile << "\\begin{tabular}{lcc}"<<endl;
    myfile << "\\hline"<<endl;
    myfile << " &  	$\\chi^{2}$ / ndof      & p-value  \\\\"<<endl;
    for (int i  = 0; i< vars.size() ; i++ ) myfile <<  vars[i] << "   &  " << chisq[i] <<"/"<<  ndof[i] <<"	&  " <<pval[i] << " \\\\"<<endl;
    myfile << "\\hline"<<endl;
    myfile << "\\end{tabular}"<<endl;
    myfile << "\\caption{\\chi^{2}/ndof$ and p-values between data and " << model << "\\ for "<< mode_string <<" measurements}"<< endl;
    myfile << "\\label{tab:"<< model <<mode<< "}" << endl;
    myfile << "\\end{minipage}"<<endl;
    myfile << "\\hspace{1cm}"<<endl;

    //myfile << "\\hfill"<<endl;
   // myfile << "{\\small"<<endl;
   // myfile << "\\begin{tabular}{lcc}"<<endl;
   // myfile << "\\hline"<<endl;
   // myfile << " &  	$\\chi^{2}$ / ndof      & p-value  \\\\"<<endl;
   // for (int i  = 0; i< vars.size() ; i++ ) myfile <<  vars[i] << "   &  " << chisq[i] <<"/"<<  ndof[i] <<"	&  " <<pval[i] << " \\\\"<<endl;
   // myfile << "\\hline"<<endl;
   // myfile << "\\end{tabular}"<<endl;
   // myfile << "\\end{minipage}"<<endl;


}
