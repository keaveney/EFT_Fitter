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
    //    "\\Delta\\ \\Phi\\ \\ttbar\\",
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
    myfile << " &  	$\\chi^{2}$ / ndof      & p-value &   $\\chi^{2}$ / ndof &   p-value \\\\"<<endl;
    myfile << "\\hline"<<endl;

    TFile * f_summary = new TFile(mode_rootfile.c_str(), "RECREATE");
    TH1F * h_summary;
    vector<TH1F*> h_summaries;
  //  vector<TFile*> f_summaries;

    std::string summary_tag;
    
    for (int summary = 0; summary < model.size(); summary++){
        
        if (model[summary] == "\\Powheg+\\Pythia"){
        summary_tag = "pwhg_p8";
        }else if (model[summary] == "\\Powheg+\\Herwigpp"){
            summary_tag = "pwhg_hpp";
        }

        std::string summary_name = mode + summary_tag;
        
        h_summary = new TH1F(summary_name.c_str(),summary_name.c_str(), vars.size(), 0, 100);
        h_summaries.push_back(h_summary);
    }
    
    
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
        "p_{T} l#bar{l}",
        "M_{l#bar{l}}",
        "#Delta #Phi l#bar{l}",
        "N_{jets}",
        "p_{T} b-jet (leading)",
        "p_{T} b-jet (sub-leading)",
        "#eta b-jet (leading)",
        "#eta b-jet (sub-leading)",
        "p_{T} b#bar{b}",
        "M_{b#bar{b}}"
    };
    
    
    
    for (int i  = 0; i< vars.size() ; i++ ) {
        
        myfile<<vars[i] <<" & "<<chisq[0][i]<<"/"<< ndof[0][i]<<"& "<<pval[0][i]<<" & "<<chisq[1][i]<<"/"<<ndof[1][i]<<" & "<<pval[1][i]<<" \\\\"<<endl;
    
        h_summaries[0]->SetBinContent(i+1, pval[0][i]);
        h_summaries[1]->SetBinContent(i+1, pval[1][i]);
        h_summaries[0]->GetXaxis()->SetBinLabel(i+1,vars_root[i].c_str());;
        h_summaries[1]->GetXaxis()->SetBinLabel(i+1,vars_root[i].c_str());;

        
    }
    myfile << "\\hline"<<endl;
    myfile << "\\end{tabular}"<<endl;
    myfile << "\\caption{The \\chi^{2}/ndof$ and p values quantifying the agreement between theoretical predictions and data for "<< mode_string <<" measurements are shown.}"<< endl;
   myfile << "\\label{tab:"<< mode << "}" << endl;

    h_summaries[0]->Write();
    h_summaries[1]->Write();
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

    gStyle->SetOptStat(00000);

    cout <<"in summary plot routine"<< endl;
    
    TFile * f_norm_parton = new TFile("norm_parton.root");
    TFile * f_norm_particle = new TFile("norm_particle.root");

    TH1F * h_norm_parton_pwhg_p8 =  (TH1F*)f_norm_parton->Get("norm_partonpwhg_p8;1");
    TH1F * h_norm_parton_pwhg_hpp =  (TH1F*)f_norm_parton->Get("norm_partonpwhg_hpp;1");
    
    TH1F * h_norm_particle_pwhg_p8 =  (TH1F*)f_norm_particle->Get("norm_particlepwhg_p8;1");
    TH1F * h_norm_particle_pwhg_hpp =  (TH1F*)f_norm_particle->Get("norm_particlepwhg_hpp;1");
    
    TFile * summary = new TFile("summary.root", "RECREATE");
    
    TCanvas *c = new TCanvas("c","c");
    c->SetCanvasSize(1800, 600);
    c->SetWindowSize(1800, 600);
    
    if (h_norm_parton_pwhg_p8){
    h_norm_parton_pwhg_p8->SetMinimum(0.0000000000000000000001);
    h_norm_parton_pwhg_p8->SetMaximum(1.1);
    h_norm_parton_pwhg_p8->GetXaxis()->SetLabelSize(0.03);
    h_norm_parton_pwhg_p8->GetXaxis()->CenterLabels(kTRUE);
    h_norm_parton_pwhg_p8->GetXaxis()->SetNdivisions(32, kFALSE);
    h_norm_parton_pwhg_p8->GetYaxis()->SetTickLength(0.01);

    h_norm_parton_pwhg_p8->SetYTitle("p-value");

        
    h_norm_parton_pwhg_p8->SetTitle(" ");
    h_norm_parton_pwhg_p8->Draw("p");
    h_norm_parton_pwhg_hpp->Draw("psame");
    h_norm_particle_pwhg_p8->Draw("psame");
    h_norm_particle_pwhg_hpp->Draw("psame");
    
    h_norm_parton_pwhg_p8->SetMarkerStyle(21);
    h_norm_parton_pwhg_hpp->SetMarkerStyle(21);
    h_norm_particle_pwhg_p8->SetMarkerStyle(22);
    h_norm_particle_pwhg_hpp->SetMarkerStyle(22);
    
    h_norm_parton_pwhg_p8->SetMarkerColor(kRed);
    h_norm_parton_pwhg_hpp->SetMarkerColor(kViolet);
    h_norm_particle_pwhg_p8->SetMarkerColor(kRed);
    h_norm_particle_pwhg_hpp->SetMarkerColor(kViolet);
    }
    
    auto legend = new TLegend(0.78,0.15,0.98,0.45);
    legend->AddEntry(h_norm_parton_pwhg_p8,"Powheg + Pythia (parton level) ","p");
    legend->AddEntry(h_norm_particle_pwhg_p8,"Powheg + Pythia (particle level)","p");
    legend->AddEntry(h_norm_parton_pwhg_hpp,"Powheg + Herwig++ (parton level)","p");
    legend->AddEntry(h_norm_particle_pwhg_hpp,"Powheg + Herwig++ (particle level)","p");

    legend->SetTextSize(0.028);
    legend->Draw();
    c->SetLogy();
    c->SetGridx();
    c->SetFrameLineColor(0);
    c->SetLeftMargin(0.07);
    c->SetRightMargin(0.02);
    c->SetTopMargin(0.06);
    c->SetBottomMargin(0.1);

    //TFrame *fr = (TFrame*)gPad->GetFrame();
    //fr->SetLineColor(kWhite);
    //gPad->Modified();
    

    
//    TFrame *frame = (TFrame*)c->GetPrimitive("TFrame");
 //   frame->SetLineWidth(0);
    
    c->SaveAs("GOF_summary.pdf");
    c->Write();
    return;
    
}

