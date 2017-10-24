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
#include "TFrame.h"
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

void write_latex(string, vector<string>, vector<string>, float[2][32] , int[2][32], float[2][32]);
void summary_plot(string, string);
std::tuple <float, int, float > process_result(string,string);

TH1F * h_model;

ofstream myfile;


int main(int argc, const char * argv[]){
    
    string mode = "norm_particle";
    
    vector<string> filenames;

    myfile.precision(2);
    string filename = mode + ".tex";
    myfile.open (filename);
    myfile << "\\begin{table}"<<endl;
    myfile << "\\centering"<<endl;
    myfile << "\\begin{tabular}{| l | c | c | c | c |}"<<endl;
    myfile << "\\hline"<<endl;
    
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
            "files/July3/norm_particle/DiffXS_HypLeptonpT_source.root",
            "files/July3/norm_particle/DiffXS_HypAntiLeptonpT_source.root",
            "files/July3/norm_particle/DiffXS_HypLeptonpTLead_source.root",
            "files/July3/norm_particle/DiffXS_HypLeptonpTNLead_source.root",
            "files/July3/norm_particle/DiffXS_HypLeptonEta_source.root",
            "files/July3/norm_particle/DiffXS_HypAntiLeptonEta_source.root",
            "files/July3/norm_particle/DiffXS_HypLeptonEtaLead_source.root",
            "files/July3/norm_particle/DiffXS_HypLeptonEtaNLead_source.root",
            "files/July3/norm_particle/DiffXS_HypLLBarpT_source.root",
            "files/July3/norm_particle/DiffXS_HypLLBarMass_source.root",
            "files/July3/norm_particle/DiffXS_HypLLBarDPhi_source.root",
            "files/July3/norm_particle/DiffXS_HypJetMultpt30_source.root",
            "files/July3/norm_particle/DiffXS_HypBJetpTLead_source.root",
            "files/July3/norm_particle/DiffXS_HypBJetpTNLead_source.root",
            "files/July3/norm_particle/DiffXS_HypBJetEtaLead_source.root",
            "files/July3/norm_particle/DiffXS_HypBJetEtaNLead_source.root",
            "files/July3/norm_particle/DiffXS_HypBBBarpT_source.root",
            "files/July3/norm_particle/DiffXS_HypBBBarMass_source.root"
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
                
                "files/July3/abs_particle/DiffXS_HypLeptonpT_source.root",
                "files/July3/abs_particle/DiffXS_HypAntiLeptonpT_source.root",
                "files/July3/abs_particle/DiffXS_HypLeptonpTLead_source.root",
                "files/July3/abs_particle/DiffXS_HypLeptonpTNLead_source.root",
                "files/July3/abs_particle/DiffXS_HypLeptonEta_source.root",
                "files/July3/abs_particle/DiffXS_HypAntiLeptonEta_source.root",
                "files/July3/abs_particle/DiffXS_HypLeptonEtaLead_source.root",
                "files/July3/abs_particle/DiffXS_HypLeptonEtaNLead_source.root",
                "files/July3/abs_particle/DiffXS_HypLLBarpT_source.root",
                "files/July3/abs_particle/DiffXS_HypLLBarMass_source.root",
                "files/July3/abs_particle/DiffXS_HypLLBarDPhi_source.root",
                "files/July3/abs_particle/DiffXS_HypJetMultpt30_source.root",
                "files/July3/abs_particle/DiffXS_HypBJetpTLead_source.root",
                "files/July3/abs_particle/DiffXS_HypBJetpTNLead_source.root",
                "files/July3/abs_particle/DiffXS_HypBJetEtaLead_source.root",
                "files/July3/abs_particle/DiffXS_HypBJetEtaNLead_source.root",
                "files/July3/abs_particle/DiffXS_HypBBBarpT_source.root",
                "files/July3/abs_particle/DiffXS_HypBBBarMass_source.root"
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
    //   "\\Delta\\ \\Phi\\ \\ttbar\\",
        "\\pt^{l}",
        "\\pt^{\\lbar}",
        "\\pt^{l (leading)}",
        "\\pt^{l (sub-leading)}",
        "\\eta^{l}",
        "\\eta^{\\lbar}",
        "\\eta^{l (leading)}",
        "\\eta^{l (sub-leading)}",
        "\\pt \\llbar\\",
        "M_{\\llbar}",
        "\\Delta\\ \\Phi\\ (\\llbar)",
        "N_{jets}",
        "\\pt^{b-jet (leading)}",
        "\\pt^{b-jet (sub-leading)}",
        "\\eta^{b-jet (leading)}",
        "\\eta^{b-jet (sub-leading)}",
        "\\pt^{\\bbbar}",
        "M_{\\bbbar}"

    };
    
   // vector<  vector<float>  > chisq;
   // vector<  vector<int>    > ndof;
   // vector<  vector<float>  > pval;
    int nvars = vars.size();
    float   chisq[2][32];
    int     ndof[2][32];
    float   pval[2][32];

    
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


    }
    
    write_latex(mode, modelnames, vars, chisq, ndof, pval);
    
    myfile << "\\end{table}"<<endl;
    myfile.close();
    
    summary_plot("norm_parton.root", "norm_particle.root");
    
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
    
    ndof = (h_model->GetNbinsX() -1); //need to reduce this by 1 for normalised results
    pval = TMath::Prob(chisq_running, ndof);
    
    cout <<"chisq , ndof, pval  =  "<< chisq_running <<"  "<< ndof <<"  "<<  pval << endl;
    
    std::tuple<float, int, float>  results ( chisq_running, ndof, pval);
    
    
    cout <<"returning results" << endl;

    return results;
    
}


void write_latex(string mode, vector<string> model, vector<string> vars, float chisq[2][32], int ndof[2][32], float pval[2][32]){

    string mode_string;
    string mode_rootfile = mode + ".root";

    cout <<"Writing latex "<< mode_rootfile <<endl;
    
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

    myfile << "\\multirow{2}{*}{} &"<<endl;
    myfile << "\\multicolumn{2}{c}{"<<  model[0] <<"} &" <<endl;
    myfile << "\\multicolumn{2}{c|}{"<<  model[1] <<"} \\\\"  <<endl;
    myfile << " &  	$\\chi^{2}$ / ndof  & p-value &   $\\chi^{2}$ / ndof &   p-value \\\\"<<endl;
    myfile << "\\hline"<<endl;

    TFile * f_summary = new TFile(mode_rootfile.c_str(), "RECREATE");
    TH1F * h_summary_A;
    TH1F * h_summary_B ;
    vector<TH1F*> h_summaries;
    int nhist;

    std::string summary_tag;
    
    for (int summary = 0; summary < model.size(); summary++){
        
        if (model[summary] == "\\Powheg+\\Pythia"){
        summary_tag = "pwhg_p8";
        }else if (model[summary] == "\\Powheg+\\Herwigpp"){
            summary_tag = "pwhg_hpp";
        }

        std::string summary_name_A = mode + summary_tag + "_A";
        std::string summary_name_B = mode + summary_tag + "_B";
        
        h_summary_A = new TH1F(summary_name_A.c_str(), summary_name_A.c_str(), 14, 0, 14);
        h_summary_B = new TH1F(summary_name_B.c_str(), summary_name_B.c_str(), vars.size() - 14, 0, vars.size() - 14);
        
        if (mode == "norm_particle"){
            nhist= 4;
        h_summaries.push_back(h_summary_A);
        h_summaries.push_back(h_summary_B);
        } else{
            nhist = 2;
            h_summaries.push_back(h_summary_A);
        }
        
        
        
    }
    
    cout <<"Histograms created... "<<endl;

    
    vector<string> vars_root = {
        "p_{T} top",
        "p_{T} anti-top",
        "p_{T} top (leading)",
        "p_{T} top (sub-leading)",
        "p_{T} top (ttbar rest frame)",
        "p_{T} anti-top (tt rest frame)",
        "Y top",
        "Y anti-top",
        "Y top (leading)",
        "Y top (sub-leading)",
        "p_{T} tt",
        "Y tt",
        "M_{tt}",
        "#Delta Y tt",
        //    "#Delta# #Phi# #ttbar#",
        "p_{T} l ",
        "p_{T} #bar{l}",
        "p_{T} l (leading)",
        "p_{T} l (sub-leading)",
        "#eta l",
        "#eta #bar{l}",
        "#eta l (leading)",
        "#eta l (sub-leading)",
        "p_{T} l #bar{l}",
        "M_{l #bar{l}}",
        "#Delta #Phi l #bar{l}",
        "N_{jets}",
        "p_{T} b-jet (leading)",
        "p_{T} b-jet (sub-leading)",
        "#eta b-jet (leading)",
        "#eta b-jet (sub-leading)",
        "p_{T} b #bar{b}",
        "M_{b #bar{b}}"
    };
    
    vector<string> vars_A, vars_B ;
    
    for (int i  = 0; i< vars.size() ; i++ ){
        if (i < 14) {
            vars_A.push_back(vars[i]);
        }
        else{
            vars_B.push_back(vars[i]);
        }
    }

    
    if (mode == "norm_particle"){
    for (int i  = 0; i< vars_A.size() ; i++ ){
        h_summaries[0]->SetBinContent(i+1, pval[0][i]);
        h_summaries[2]->SetBinContent(i+1, pval[1][i]);
        h_summaries[0]->GetXaxis()->SetBinLabel(i+1,vars_root[i].c_str());
        h_summaries[2]->GetXaxis()->SetBinLabel(i+1,vars_root[i].c_str());
    }
    for (int i  = 0; i< vars_B.size() ; i++ ){
        h_summaries[1]->SetBinContent(i+1, pval[0][i+14]);
        h_summaries[3]->SetBinContent(i+1, pval[1][i+14]);
        h_summaries[1]->GetXaxis()->SetBinLabel(i+1,vars_root[i+14].c_str());;
        h_summaries[3]->GetXaxis()->SetBinLabel(i+1,vars_root[i+14].c_str());;
    }
    } else{
        for (int i  = 0; i< vars_A.size() ; i++ ){
            h_summaries[0]->SetBinContent(i+1, pval[0][i]);
            h_summaries[1]->SetBinContent(i+1, pval[1][i]);
            h_summaries[0]->GetXaxis()->SetBinLabel(i+1,vars_root[i].c_str());
            h_summaries[1]->GetXaxis()->SetBinLabel(i+1,vars_root[i].c_str());
        }
    }
    
    for (int i  = 0; i< vars.size() ; i++ ){
        
        myfile<<vars[i] <<" & "<<chisq[0][i]<<"/"<< ndof[0][i]<<"& "<<pval[0][i]<<" & "<<chisq[1][i]<<"/"<<ndof[1][i]<<" & "<<pval[1][i]<<" \\\\"<<endl;

    }
    myfile << "\\hline"<<endl;
    myfile << "\\end{tabular}"<<endl;
    myfile << "\\caption{The \\chi^{2}/ndof$ and p values quantifying the agreement between theoretical predictions and data for "<< mode_string <<" measurements are shown.}"<< endl;
    myfile << "\\label{tab:"<< mode << "}" << endl;

    for (int hist = 0; hist< nhist ;hist++){
    h_summaries[hist]->Write();
    }
    
    
    
    f_summary->Close();
    
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


void summary_plot(string f1, string f2){

    TCanvas * c_master = new TCanvas("c_master","c_master");
    c_master->Divide(1,2);

    c_master->cd(1);

    gPad->SetCanvasSize(2200, 1800);
   // gPad->SetWindowSize(2200, 1800);
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetFrameLineColor(0);
    gPad->SetLeftMargin(0.085);
    gPad->SetRightMargin(0.045);
    gPad->SetTopMargin(0.02);
    gPad->SetBottomMargin(0.25);
    
    gStyle->SetOptStat(00000);
    
    cout <<"in summary plot routine"<< endl;
    
    TFile * f_norm_parton = new TFile("norm_parton.root");
    TFile * f_norm_particle = new TFile("norm_particle.root");

    //pwhg_p8
    TH1F * h_norm_parton_pwhg_p8_A =  (TH1F*)f_norm_parton->Get("norm_partonpwhg_p8_A;1");
    TH1F * h_norm_particle_pwhg_p8_A =  (TH1F*)f_norm_particle->Get("norm_particlepwhg_p8_A;1");
    TH1F * h_norm_particle_pwhg_p8_B =  (TH1F*)f_norm_particle->Get("norm_particlepwhg_p8_B;1");

    //pwhg_hpp
    TH1F * h_norm_parton_pwhg_hpp_A =  (TH1F*)f_norm_parton->Get("norm_partonpwhg_hpp_A;1");
    TH1F * h_norm_particle_pwhg_hpp_A =  (TH1F*)f_norm_particle->Get("norm_particlepwhg_hpp_A;1");
    TH1F * h_norm_particle_pwhg_hpp_B =  (TH1F*)f_norm_particle->Get("norm_particlepwhg_hpp_B;1");
    
    TFile * summary = new TFile("summary.root", "RECREATE");
    
    //TCanvas *c_A = new TCanvas("c_A","c_A");
    //c_A->SetCanvasSize(2200, 600);
   // c_A->SetWindowSize(2200, 600);
    
    if (h_norm_particle_pwhg_p8_A){
    h_norm_particle_pwhg_p8_A->SetMinimum(0.0000000000000000000001);
    h_norm_particle_pwhg_p8_A->SetMaximum(1.1);
    h_norm_particle_pwhg_p8_A->GetXaxis()->SetLabelSize(0.06);
    h_norm_particle_pwhg_p8_A->GetXaxis()->CenterLabels(kTRUE);
    h_norm_particle_pwhg_p8_A->GetXaxis()->SetNdivisions(32, kFALSE);
    h_norm_particle_pwhg_p8_A->GetXaxis()->SetLabelOffset(0.02);
    //h_norm_particle_pwhg_p8_A->GetXaxis()->LabelsOption("u");

    h_norm_particle_pwhg_p8_A->GetYaxis()->SetTickLength(0.01);
    h_norm_particle_pwhg_p8_A->GetYaxis()->SetTitleOffset(0.6);
    h_norm_particle_pwhg_p8_A->GetYaxis()->SetTitleSize(0.06);

    h_norm_particle_pwhg_p8_A->SetYTitle("p-value");

    h_norm_particle_pwhg_p8_A->SetTitle(" ");
    h_norm_particle_pwhg_p8_A->Draw("p");
    h_norm_parton_pwhg_hpp_A->Draw("psame");
    h_norm_particle_pwhg_p8_A->Draw("psame");
    //h_norm_particle_pwhg_hpp_B->Draw("psame");
    
    h_norm_parton_pwhg_p8_A->SetMarkerStyle(21);
    h_norm_parton_pwhg_hpp_A->SetMarkerStyle(21);
    h_norm_particle_pwhg_p8_A->SetMarkerStyle(22);
    h_norm_particle_pwhg_hpp_A->SetMarkerStyle(22);
        
    h_norm_parton_pwhg_p8_A->SetMarkerSize(3.0);
    h_norm_parton_pwhg_hpp_A->SetMarkerSize(3.0);
    h_norm_particle_pwhg_p8_A->SetMarkerSize(3.0);
    h_norm_particle_pwhg_hpp_A->SetMarkerSize(3.0);
        
    h_norm_parton_pwhg_p8_A->SetMarkerColor(kRed);
    h_norm_parton_pwhg_hpp_A->SetMarkerColor(kViolet);
    h_norm_particle_pwhg_p8_A->SetMarkerColor(kRed);
    h_norm_particle_pwhg_hpp_A->SetMarkerColor(kViolet);
        
    }
    
    auto legend_A = new TLegend(0.6,0.33,0.95,0.69);
    legend_A->AddEntry(h_norm_parton_pwhg_p8_A,"Powheg + Pythia (parton level) ","p");
    legend_A->AddEntry(h_norm_particle_pwhg_p8_A,"Powheg + Pythia (particle level)","p");
    legend_A->AddEntry(h_norm_parton_pwhg_hpp_A,"Powheg + Herwig++ (parton level)","p");
    legend_A->AddEntry(h_norm_particle_pwhg_hpp_A,"Powheg + Herwig++ (particle level)","p");

    legend_A->SetTextSize(0.035);
    legend_A->Draw();
   // c_A->SetLogy();
   // c_A->SetGridx();
   // c_A->SetFrameLineColor(0);
   // c_A->SetLeftMargin(0.035);
    //c_A->SetRightMargin(0.01);
    //c_A->SetTopMargin(0.06);
    //c_A->SetBottomMargin(0.26);

    //TFrame *fr = (TFrame*)gPad->GetFrame();
    //fr->SetLineColor(kWhite);
    //gPad->Modified();
    

    
//    TFrame *frame = (TFrame*)c->GetPrimitive("TFrame");
 //   frame->SetLineWidth(0);
    
    //c_A->SaveAs("GOF_summary_A.pdf");
   // c_A->Write();
    

    
   // TCanvas *c_B = new TCanvas("c_B","c_B");
   // c_B->SetCanvasSize(2200, 600);
   // c_B->SetWindowSize(2200, 600);
    
    c_master->cd(2);
    gPad->SetCanvasSize(2200, 1800);
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetFrameLineColor(0);
    gPad->SetLeftMargin(0.085);
    gPad->SetRightMargin(0.045);
    gPad->SetTopMargin(0.02);
    gPad->SetBottomMargin(0.2);
    
    if (h_norm_particle_pwhg_p8_B){
        h_norm_particle_pwhg_p8_B->SetMinimum(0.0000000000000000000001);
        h_norm_particle_pwhg_p8_B->SetMaximum(1.1);
        h_norm_particle_pwhg_p8_B->GetXaxis()->SetLabelSize(0.06);
        h_norm_particle_pwhg_p8_B->GetXaxis()->CenterLabels(kTRUE);
        h_norm_particle_pwhg_p8_B->GetXaxis()->SetNdivisions(32, kFALSE);
        h_norm_particle_pwhg_p8_B->GetXaxis()->SetLabelOffset(0.01);
        //h_norm_particle_pwhg_p8_B->GetXaxis()->LabelsOption("u");
        
        h_norm_particle_pwhg_p8_B->GetYaxis()->SetTickLength(0.01);
        h_norm_particle_pwhg_p8_B->GetYaxis()->SetTitleOffset(0.6);
        h_norm_particle_pwhg_p8_B->GetYaxis()->SetTitleSize(0.06);

        h_norm_particle_pwhg_p8_B->SetYTitle("p-value");
        
        h_norm_particle_pwhg_p8_B->SetTitle(" ");
        h_norm_particle_pwhg_p8_B->Draw("p");
        h_norm_particle_pwhg_hpp_B->Draw("psame");
        
        h_norm_particle_pwhg_p8_B->SetMarkerStyle(22);
        h_norm_particle_pwhg_hpp_B->SetMarkerStyle(22);
        
        h_norm_particle_pwhg_p8_B->SetMarkerSize(3.0);
        h_norm_particle_pwhg_hpp_B->SetMarkerSize(3.0);
        
        h_norm_particle_pwhg_p8_B->SetMarkerColor(kRed);
        h_norm_particle_pwhg_hpp_B->SetMarkerColor(kViolet);
        
    }
    
   // auto legend_B = new TLegend(0.78,0.25,0.98,0.55);
  //  legend_B->AddEntry(h_norm_particle_pwhg_p8_B,"Powheg + Pythia (particle level)","p");
  //  legend_B->AddEntry(h_norm_particle_pwhg_hpp_B,"Powheg + Herwig++ (particle level)","p");
    
  //  legend_B->SetTextSize(0.028);
  //  legend_B->Draw();
  //  c_B->SetLogy();
  //  c_B->SetGridx();
  //  c_B->SetFrameLineColor(0);
  //  c_B->SetLeftMargin(0.035);
  ///  c_B->SetRightMargin(0.01);
   // c_B->SetTopMargin(0.06);
   // c_B->SetBottomMargin(0.26);
    
    //TFrame *fr = (TFrame*)gPad->GetFrame();
    //fr->SetLineColor(kWhite);
    //gPad->Modified();
    
    
    
    //    TFrame *frame = (TFrame*)c->GetPrimitive("TFrame");
    //   frame->SetLineWidth(0);
    
   // c_B->SaveAs("GOF_summary_B.pdf");
   // c_B->Write();
    
    c_master->SaveAs("GOF_summary_master.pdf");
    
    return;
    
}

