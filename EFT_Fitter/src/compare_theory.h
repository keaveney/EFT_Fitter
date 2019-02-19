#include <iostream>
#include <iomanip>

#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TFrame.h"
#include "TMath.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TObject.h"
#include "TIterator.h"
#include "TArrayD.h"

#include "TPaveText.h"

#include "TF1.h"
#include <utility>
#include <tuple>

#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include "helper_tools.h"

#include "tdrstyle.C"

TGraphAsymmErrors* read_prediction(std::string, bool, std::string);
void make_plot(std::string, std::string, std::string, std::string, vector<std::string>, vector<std::string> , vector<float>);
TGraphAsymmErrors* make_ratio_graph(TGraphAsymmErrors*,TGraphAsymmErrors*);
TH1D* make_ratio_histo(TGraphAsymmErrors*, TH1F*);
TGraphAsymmErrors* refine_graph(TGraphAsymmErrors*, int, double,std::string);
TGraphAsymmErrors* remove_ex(TGraphAsymmErrors*);
vector<TGraphAsymmErrors*> make_unc_graphs(TGraphAsymmErrors*,TGraphAsymmErrors*,TH1F*);

//top pt
vector<std::string> preds_top_pt = {
    "predictions/NNLO_EW/LHC13-CMS-PTt-mTt2-LUXQED17-pheno.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-PTt-mTt2-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-PTt-mTt2-NNPDF31-pheno.dat",
    "predictions/NNLLprime/LHC13-CMS-PTtav-mTt2-NNPDF31.dat",
    "predictions/NNLLprime_172-5/LHC13-CMS-PTtav-mTt2-NNPDF31.dat",
    "predictions/aNNNLO/LHC13-CMS-PTtav-sqrt-mT-NNPDF30.dat",
    "predictions/aNNLO/LHC13-CMS-PTtav-mT-CT14NNLO.dat"
    //"predictions/aNNLO/LHC13-CMS-PTtav-mT-CT14NNLO.dat"
};

vector<std::string> pred_names_top_pt = {
    "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 173.3 GeV",
    "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 172.5 GeV",
    "NNLO+#alpha_{EW}^{3} (NNPDF3.1) m_{t} = 173.3 GeV",
    "NNLO+NNLL' (NNPDF3.1) m_{t} = 173.3 GeV",
    "NNLO+NNLL' (NNPDF3.1) m_{t} = 172.5 GeV",
    "aN^{3}LO (NNPDF3.0) m_{t} = 172.5 GeV",
    "aNNLO (CT14NNLO) m_{t} = 172.5 GeV"
};

vector<bool> divide_options_top_pt = {
    true,
    true,
    true,
    false,
    false,
    true,
    false};

vector<float> plot_options_top_pt = {
    0.0006, 29.0, 0.81,1.46, 0.2,0.019,0.5,0.5, 0.24,0.75,0.41,0.92};

vector<float> plot_options_top_pt_norm = {
    0.0000006, 0.038, 0.81,1.58, 0.2,0.019,0.5,0.5, 0.24,0.72,0.41,0.93};

//anti-top pt
vector<std::string> preds_antitop_pt = {
    "predictions/NNLO_EW/LHC13-CMS-PTtx-mTtx2-LUXQED17-pheno.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-PTtx-mTtx2-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-PTtx-mTtx2-NNPDF31-pheno.dat"
  //  "predictions/NNLLprime/LHC13-CMS-PTtav-mTt2-NNPDF31.dat"
  //  "predictions/aNNNLO/LHC13-CMS-PTtav-sqrt-mT-NNPDF30.dat",
   // "predictions/aNNLO/LHC13-CMS-PTtav-mT-NNPDF30.dat"
};

vector<std::string> pred_names_antitop_pt = {
    "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 173.3 GeV",
    "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 172.5 GeV",
    "NNLO+#alpha_{EW}^{3} (NNPDF3.1) m_{t} = 173.3 GeV"
   // "NNLO+NNLL' (NNPDF3.1)"
  //  "approx. N^{3}LO (NNPDF3.0)",
  //  "approx. NNLO (CT14NNLO)"
};

vector<bool> divide_options_antitop_pt = {
    true,
    true,
    true
 //  false
 //   true,
 //   false
};

vector<float> plot_options_antitop_pt = {
   // 0.01, 20.0, 0.81,1.56, 0.6,0.45,0.8,0.78, 0.24,0.73,0.41,0.925
    0.012, 16.0, 0.81,1.46, 0.2,0.019,0.5,0.4, 0.24,0.75,0.41,0.92
    
};

vector<float> plot_options_antitop_pt_norm = {
    0.000013, 0.016, 0.81,1.58, 0.2,0.019,0.5,0.4, 0.24,0.73,0.41,0.93
//    0.0001, 0.01, 0.81,1.52, 0.6,0.45,0.8,0.78, 0.24,0.73,0.41,0.92
    
};

//yt
vector<std::string> preds_top_y = {
    "predictions/NNLO_EW/LHC13-CMS-Yt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-Yt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-Yt-HT4-NNPDF31-pheno.dat",
    "predictions/aNNNLO/LHC13-CMS-Yt-NNPDF3.0.dat",
    "predictions/aNNLO/LHC13-CMS-Yt-CT14NNLO.dat"};

vector<std::string> pred_names_top_y = {
    "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 173.3 GeV",
    "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 172.5 GeV",
    "NNLO+#alpha_{EW}^{3} (NNPDF3.1) m_{t} = 173.3 GeV",
    "aN^{3}LO (NNPDF3.0) m_{t} = 172.5 GeV",
    "aNNLO (CT14NNLO) m_{t} = 172.5 GeV"};

vector<bool> divide_options_top_y = {
    true,
    true,
    true,
    true,
    false};

vector<float> plot_options_top_y = {
    45.0, 350.0, 0.87,1.26,0.31,0.02,0.58,0.37, 0.4,0.8,0.6,0.97};

vector<float> plot_options_top_y_norm = {
    0.05, 0.5, 0.87,1.19,0.31,0.02,0.58,0.37, 0.4,0.8,0.6,0.97};

//ytbar
vector<std::string> preds_antitop_y = {
    "predictions/NNLO_EW/LHC13-CMS-Ytx-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-Ytx-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-Ytx-HT4-NNPDF31-pheno.dat"
    //"predictions/aNNNLO/LHC13-CMS-Yt-NNPDF3.0.dat",
   // "predictions/aNNLO/LHC13-CMS-Yt-NNPDF3.0.dat"
};

vector<std::string> pred_names_antitop_y = {
    "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 173.3 GeV",
    "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 172.5 GeV",
    "NNLO+#alpha_{EW}^{3} (NNPDF3.1) m_{t} = 173.3 GeV"
  //  "approx. N^{3}LO (NNPDF3.0)",
  //  "approx. NNLO (CT14NNLO)"
};

vector<bool> divide_options_antitop_y = {
    true,
    true,
    true
  //  true,
 //   false
};

vector<float> plot_options_antitop_y = {
    45.0, 350.0, 0.87,1.22,0.3,0.02,0.6,0.34,  0.4,0.77,0.6,0.94
};

vector<float> plot_options_antitop_y_norm = {
    0.05, 0.5, 0.87,1.19,0.3,0.02,0.6,0.34, 0.4,0.79,0.6,0.96
};

//mtt
vector<std::string> preds_mtt = {
    "predictions/NNLO_EW/LHC13-CMS-Mtt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-Mtt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-Mtt-HT4-NNPDF31-pheno.dat",
    "predictions/NNLLprime/LHC13-CMS-Mtt-HT4-NNPDF31.dat",
    "predictions/NNLLprime_172-5/LHC13-CMS-Mtt-HT4-NNPDF31.dat"

};

vector<std::string> pred_names_mtt = {
    "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 173.3 GeV",
    "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 172.5 GeV",
    "NNLO+#alpha_{EW}^{3} (NNPDF3.1) m_{t} = 173.3 GeV",
    "NNLO+NNLL' (NNPDF3.1) m_{t} = 173.3 GeV",
    "NNLO+NNLL' (NNPDF3.1) m_{t} = 172.5 GeV"

};

vector<bool> divide_options_mtt = {
    true,
    true,
    true,
    false,
    false
};

vector<float> plot_options_mtt = {
    0.0003, 20.0, 0.6,1.55, 0.22,0.03,0.42,0.39, 0.25,0.78,0.45,0.985
};

vector<float> plot_options_mtt_norm = {
    0.0000006, 0.02, 0.6, 1.55, 0.22,0.03,0.42,0.39, 0.25,0.78,0.45,0.985
};

//pttt
vector<std::string> preds_pttt = {
    "predictions/NNLO_EW/LHC13-CMS-PTtt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-PTtt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-PTtt-HT4-NNPDF31-pheno.dat"
};

vector<std::string> pred_names_pttt = {
    "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 173.3 GeV",
    "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 172.5 GeV",
    "NNLO+#alpha_{EW}^{3} (NNPDF3.1) m_{t} = 173.3 GeV",
};

vector<bool> divide_options_pttt = {
    true,
    true,
    true
};

vector<float> plot_options_pttt = {
    0.014, 350.0, 0.75,1.45, 0.38,0.4,0.7,0.71, 0.35,0.77,0.55,0.94
};
vector<float> plot_options_pttt_norm = {
    0.000014, 0.35, 0.76, 1.45, 0.38,0.4,0.7,0.71, 0.35,0.74,0.55,0.91
};

//ytt
vector<std::string> preds_ytt = {
    "predictions/NNLO_EW/LHC13-CMS-Ytt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-Ytt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-Ytt-HT4-NNPDF31-pheno.dat"
    //  "predictions/NNLLprime/LHC13-CMS-Ytt-HT4-NNPDF31.dat"
};

vector<std::string> pred_names_ytt = {
    "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 173.3 GeV",
    "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 172.5 GeV",
    "NNLO+#alpha_{EW}^{3} (NNPDF3.1) m_{t} = 173.3 GeV"
    //   "NNLO+NNLL' (NNPDF3.1)"
};


vector<bool> divide_options_ytt = {
    true,
    true,
    true
    // false
};

vector<float> plot_options_ytt = {
    19.0, 500.0, 0.83,1.28, 0.34,0.02,0.59,0.37, 0.4,0.76,0.6,0.97
};

vector<float> plot_options_ytt_norm = {
    0.025, 0.6, 0.83,1.28, 0.34,0.02,0.59,0.37, 0.4,0.76,0.6,0.97
};

//rapidity differnce
vector<std::string> preds_dytt = {
    //  "predictions/NNLO_EW/LHC13-CMS-Ytt-HT4-LUXQED17-pheno.dat"
    //  "predictions/NNLO_EW/LHC13-CMS-dytt-HT4-LUXQED17-pheno.dat"
    "predictions/NNLO_EW/LHC13-CMS-dytt-HT4-LUXQED17-pheno_rewrite.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-dytt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-dytt-HT4-NNPDF31-pheno.dat"
    //  "predictions/NNLLprime/LHC13-CMS-Ytt-HT4-NNPDF31.dat"
};

vector<std::string> pred_names_dytt = {
    "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 173.3 GeV",
    "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 172.5 GeV",
    "NNLO+#alpha_{EW}^{3} (NNPDF3.1) m_{t} = 173.3 GeV"
    //   "NNLO+NNLL' (NNPDF3.1)"
};

vector<bool> divide_options_dytt = {
    true,
    true,
    true
    // false
};

vector<float> plot_options_dytt = {
    35.0, 550.0, 0.85,1.19, 0.41,0.11,0.6,0.42, 0.45,0.81,0.65,0.97
};

vector<float> plot_options_dytt_norm = {
    0.044, 0.7, 0.91,1.13, 0.41,0.11,0.6,0.42,0.45,0.735,0.65,0.88
};


TH1D* make_ratio_histo(TGraphAsymmErrors* g_den, TH1F* h_num){
    
    const TArrayD* xaxis = h_num->GetXaxis()->GetXbins();
    const Double_t* x_axis_ar = xaxis->GetArray();
    
    int bins = h_num->GetNbinsX();
    TH1D* h_ratio = new TH1D("","",bins, x_axis_ar);
    
    double b_num, b_den;
    double x;
    
    for (int bin = 0; bin< h_num->GetNbinsX(); bin++){
        b_num = h_num->GetBinContent(bin+1);
        g_den->GetPoint(bin, x, b_den);
        h_ratio->SetBinContent(bin+1, b_num/b_den);
        //std::cout <<"make_ratio_histo XXX  " <<b_num/b_den <<std::endl;
    }
    h_ratio->SetLineColor(kRed+1);
    h_ratio->SetLineWidth(2.0);

    
    return h_ratio;
};


TGraphAsymmErrors* remove_ex(TGraphAsymmErrors*  g){
    for (int p = 0; p < g->GetN(); p++){
        g->SetPointEXlow(p, 0.0);
        g->SetPointEXhigh(p,0.0);
    }
    return g;
}

TGraphAsymmErrors* refine_graph(TGraphAsymmErrors*  g_orig, int pred_index, double total_predictions,  std::string pred_name){
    
    //std::cout <<"  " << std::endl;
    //std::cout <<"refining   " << std::endl;

    vector<int> marker_styles = {21,22,23,27,28,29};
    vector<int> marker_colours = {3,4,6,7,8,9};
    
    TGraphAsymmErrors* g = new TGraphAsymmErrors();
    double bin_edge, bin_centre, y, new_centre;
    double x_orig, y_orig, exl_orig, exh_orig, eyl_orig , eyh_orig;
    double bin_offset = 0.05;

    for(int bin = 0; bin < g_orig->GetN(); bin++){
        g_orig->GetPoint(bin,x_orig, y_orig);
        
        exl_orig = g_orig->GetErrorXlow(bin);
        exh_orig = g_orig->GetErrorXhigh(bin);
        eyl_orig = g_orig->GetErrorYlow(bin);
        eyh_orig = g_orig->GetErrorYhigh(bin);
        
        bin_edge = x_orig - exl_orig;
        new_centre = bin_edge + (pred_index+1)*((exl_orig+exh_orig)/(total_predictions+1.0)) + ((bin_offset)*(exl_orig+exh_orig));
        g->SetPoint(bin, new_centre, y_orig);
        g->SetPointError(bin, exl_orig, exh_orig, eyl_orig, eyh_orig);
       // std::cout <<"refining   " <<  bin<<"  "<< new_centre<<" " << y_orig  << std::endl;
    }
    
    if (pred_name == "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 173.3 GeV"){
        //g->SetMarkerStyle(21);
        //g->SetMarkerColor(3);
        g->SetMarkerStyle(21);
        g->SetMarkerColor(3);
    }
    else if (pred_name == "NNLO+#alpha_{EW}^{3} (LUXQED17) m_{t} = 172.5 GeV"){
       // g->SetMarkerStyle(21);
       // g->SetMarkerColor(28);
        g->SetMarkerStyle(21);
        g->SetMarkerColor(28);

    }
    else if (pred_name == "NNLO+#alpha_{EW}^{3} (NNPDF3.1) m_{t} = 173.3 GeV"){
       // g->SetMarkerStyle(21);
       // g->SetMarkerColor(4);
        g->SetMarkerStyle(25);
        g->SetMarkerColor(4);

    } else if (pred_name == "NNLO+NNLL' (NNPDF3.1) m_{t} = 173.3 GeV") {
        //g->SetMarkerStyle(22);
        //g->SetMarkerColor(6);
        g->SetMarkerStyle(20);
        g->SetMarkerColor(6);

    }else if (pred_name == "NNLO+NNLL' (NNPDF3.1) m_{t} = 172.5 GeV") {
       // g->SetMarkerStyle(22);
      //  g->SetMarkerColor(32);
         g->SetMarkerStyle(24);
          g->SetMarkerColor(32);
    
    }
    else if (pred_name == "aN^{3}LO (NNPDF3.0) m_{t} = 172.5 GeV") {
        g->SetMarkerStyle(27);
        g->SetMarkerColor(7);
    
    }else if (pred_name == "aNNLO (CT14NNLO) m_{t} = 172.5 GeV") {
        g->SetMarkerStyle(28);
        g->SetMarkerColor(8);

    }
    
    
    g->SetMarkerSize(0.8);

    return g;
}


vector<TGraphAsymmErrors*> make_unc_graphs(TGraphAsymmErrors* g_data, TGraphAsymmErrors* g_data_stat, TH1F* h_ratio ){
    
    vector<TGraphAsymmErrors*> unc_graphs;
    TGraphAsymmErrors * g_ratio_total = new TGraphAsymmErrors();
    TGraphAsymmErrors * g_ratio_stat = new TGraphAsymmErrors();
    double x_total, y_total, x_stat, y_stat, unc_total, unc_stat;

    for(int p = 0; p < g_data->GetN(); p++){
    
        double bin_error = h_ratio->GetBinWidth(p+1)/2.0;
        g_data->GetPoint(p, x_total, y_total);
        
        unc_total = g_data->GetErrorY(p) / y_total;
        unc_stat = g_data_stat->GetErrorY(p) /  y_total;
        
        g_ratio_total->SetPoint(p,x_total, 1.0);
        g_ratio_total->SetPointError(p, bin_error,bin_error, unc_total, unc_total);
        
        g_ratio_stat->SetPoint(p,x_total, 1.0);
        g_ratio_stat->SetPointError(p, bin_error,bin_error, unc_stat, unc_stat);
    }

    //g_ratio_total->SetFillColor(kOrange-4);
    //g_ratio_stat->SetFillColor(kGray+1);
    g_ratio_total->SetFillColor(796);
    g_ratio_stat->SetFillColor(921);
    g_ratio_total->SetLineWidth(0.0);
    g_ratio_stat->SetLineWidth(0.0);

    unc_graphs.push_back(g_ratio_total);
    unc_graphs.push_back(g_ratio_stat);
    
    return unc_graphs;
};


TGraphAsymmErrors* make_ratio_graph(TGraphAsymmErrors* g_data,TGraphAsymmErrors* g_pred){
    
    TGraphAsymmErrors* g_ratio = new TGraphAsymmErrors();
    
    double x_data, y_data, x_pred, y_pred, pred, ratio;
    
    for (int i = 0; i < g_data->GetN(); i++){
        g_data->GetPoint(i, x_data, y_data);
        g_pred->GetPoint(i, x_pred, y_pred);
        ratio = y_pred/y_data;
        g_ratio->SetPoint(i,x_data, ratio);
        g_ratio->SetPointEXlow(i, g_pred->GetErrorXlow(i));
        g_ratio->SetPointEXhigh(i, g_pred->GetErrorXhigh(i));
        g_ratio->SetPointEYlow(i, g_pred->GetErrorYlow(i)/y_pred);
        g_ratio->SetPointEYhigh(i, g_pred->GetErrorYhigh(i)/y_pred);
        //std::cout <<"ratio  "<< ratio << std::endl;
    }
    return  g_ratio;
}


void make_plot(std::string observable, std::string units, std::string filename,std::string type, vector<std::string> predictions, vector<std::string> prediction_names, vector<bool> divide, vector<float> plot_params){

    gStyle->SetOptStat(00000);
    gStyle->SetEndErrorSize(0);
    //gStyle->SetErrorX(0.00001);
    
    vector<int> colors = {3,4,5,6,7,8,9};
    vector<TGraphAsymmErrors*> pred_graphs;
    
    std::cout <<"***READING FILENAME***   "<<  filename.c_str()  << std::endl;
    TFile * f = new TFile(filename.c_str());
    TGraphAsymmErrors * g = (TGraphAsymmErrors*)f->Get("data");
    TGraphAsymmErrors * g_stat = (TGraphAsymmErrors*)f->Get("data_staterror_only");
    TCanvas * c_plot = (TCanvas*)f->Get("canvas");
    TH1F *  h_mc = (TH1F*)c_plot->GetPrimitive("Nominalplot_copy");
    
    TCanvas * c = new TCanvas("compare","",780,850);
    
    TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
    pad1->SetBottomMargin(0);
    pad1->SetTopMargin(0.18);
    pad1->SetLeftMargin(0.19);
    pad1->SetRightMargin(0.05);
    pad1->SetTicks();

    pad1->Draw();
    pad1->cd();
    std::string y_axis_title;
    std::string first_fragment;
    
    if (type =="abs"){
        first_fragment = "#frac{d#sigma}{d";
        if(units != "") {
            y_axis_title = first_fragment + observable + "} " + "[pb/" + units + "]";
        }else{
            y_axis_title = first_fragment + observable + "} " + "[pb" + units + "]";
        }
    
    }else{
        first_fragment = "#frac{1}{#sigma} #frac{d#sigma}{d";
    
        if(units != "") {
            y_axis_title = first_fragment + observable + "} " + "[" + units + "^{-1}]";
        }else{
            y_axis_title = first_fragment + observable + "} ";
        }
    }
    
    h_mc->SetYTitle(y_axis_title.c_str());
    h_mc->GetYaxis()->SetTitleSize(0.07);
    h_mc->GetYaxis()->SetLabelSize(0.06);
    h_mc->GetYaxis()->SetTitleOffset(1.2);
    h_mc->SetTitle("");
    h_mc->SetLineColor(kRed+1);

    
    double integral =0.0;
    double b_mc, b_w;
/*    if (type=="norm"){
        for (int bin = 0 ; bin < h_mc->GetNbinsX(); bin++){
            b_mc = h_mc->GetBinContent(bin+1);
            b_w = h_mc->GetBinWidth(bin+1);
            integral = integral + (b_mc*b_w);
        }
        for (int bin = 0 ; bin < h_mc->GetNbinsX(); bin++){
            b_mc = h_mc->GetBinContent(bin+1);
            h_mc->SetBinContent(bin+1, b_mc/integral);
        }
    }
 */

    h_mc->Draw();
    h_mc->SetMarkerSize(0);
    h_mc->SetMarkerColor(kRed+1);

    h_mc->GetYaxis()->SetRangeUser(plot_params[0], plot_params[1]);
    
    c->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.4);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.28);
    pad2->SetTopMargin(0.0);
    pad2->SetLeftMargin(0.19);
    pad2->SetRightMargin(0.05);
    pad2->SetGridy();

    pad2->SetTicks();
    pad2->Draw();
    pad2->cd();
    
    TH1F *h;
    TH1F * h_ratio;
    TGraphAsymmErrors * g_pred;
    TGraphAsymmErrors * g_ratio;
    vector<TGraphAsymmErrors*> g_ratio_uncs;
    
    TH1F* h_ratio_base = (TH1F*)h_mc->Clone();
    h_ratio_base->Reset();
    
    for (int bin = 0 ; bin < h_ratio_base->GetNbinsX(); bin++) h_ratio_base->SetBinContent(bin+1, 1.0);
    h_ratio_base->SetLineWidth(1.0);
    h_ratio_base->SetLineColor(1.0);
    h_ratio_base->Draw();
    h_ratio_base->GetYaxis()->SetRangeUser(plot_params[2], plot_params[3]);
    h_ratio_base->GetYaxis()->SetNdivisions(505);
    h_ratio_base->GetXaxis()->SetNdivisions(506);
    h_ratio_base->GetYaxis()->SetLabelSize(0.09);
    h_ratio_base->GetXaxis()->SetLabelSize(0.1);
    h_ratio_base->GetXaxis()->SetTitleSize(0.11);
    h_ratio_base->GetXaxis()->SetTitleOffset(1.0);
    h_ratio_base->GetYaxis()->SetTitleOffset(0.63);
    h_ratio_base->GetYaxis()->SetTitleSize(0.1);
    h_ratio_base->SetYTitle("#frac{Theory}{Data}");
    g_ratio_uncs = make_unc_graphs(g, g_stat, h_ratio_base);
    g_ratio_uncs[0]->Draw("E2SAME");
    g_ratio_uncs[1]->Draw("E2SAME");
    h_ratio_base->Draw("SAME");
    TH1D * h_mc_ratio = make_ratio_histo(g,h_mc);
    h_mc_ratio->Draw("SAME");
    
   // float bins[7] = {0.0, 65.00, 125.00, 200.00, 290.00, 400.00, 550.00};
   // TH1F* test = new TH1F("","",6,bins);
   // double data_x, data_y, pred_x, pred_y;
   //    for (int p = 0; p < h_mc_ratio->GetNbinsX(); p++){
   //       pred_y = h_mc->GetBinContent(p+1);
   //       g->GetPoint(p, data_x, data_y);
   //       test->SetBinContent(p+1, pred_y/data_y);
   //   }
   //   test->Draw("SAME");
    

    for (int p = 0; p < predictions.size(); p++){
       // std::cout <<"***LOOPING ON PREDICTIONS**** =   "  << std::endl;

        g_pred = read_prediction(predictions[p], divide[p], type);
        g_pred = refine_graph(g_pred, p, predictions.size(), prediction_names[p]);

        g_ratio = make_ratio_graph(g, g_pred);
        g_ratio = refine_graph(g_ratio, p, predictions.size(), prediction_names[p]);

        g_pred = remove_ex(g_pred);
        g_ratio = remove_ex(g_ratio);
        
        pred_graphs.push_back(g_pred);
        pad1->cd();
        
        g_pred->Draw("E0PSAME");
        pad2->cd();
        g_ratio->Draw("E0PSAME");
    }
    
    std::string x_axis_title;
    if(units == ""){
        x_axis_title = observable;
    }else{
        x_axis_title = observable + " [" + units + "]";
    }
    std::cout <<"***OBSERVABLE****"  << x_axis_title <<"***" <<std::endl;
    h_ratio_base->SetXTitle(x_axis_title.c_str());

    pad2->SetGridy();

    
    pad1->cd();
    pad1->SetLogy();
    g->SetMarkerSize(0.8);
    g->SetLineWidth(2.0);

    
    TH1D* data_histo;
    data_histo = graph2histo(g);
    data_histo->SetMarkerSize(0.8);
    data_histo->SetLineWidth(2.0);
    data_histo->SetLineColor(1);
    data_histo->SetMarkerStyle(20);
    g->Draw("E0PSAME");

    vector<TH1F*> dummy_hists;
    for (int p = 0; p< prediction_names.size(); p++){
        TH1F * dummy_hist = new TH1F("","", 1, 0.0, 1.0);
        dummy_hist->SetMarkerSize(pred_graphs[p]->GetMarkerSize());
        dummy_hist->SetMarkerStyle(pred_graphs[p]->GetMarkerStyle());
        dummy_hist->SetMarkerColor(pred_graphs[p]->GetMarkerColor());
        dummy_hist->SetLineColor(1);
        dummy_hists.push_back(dummy_hist);
    }
    
    TLegend *leg = new TLegend(plot_params[4], plot_params[5], plot_params[6], plot_params[7]);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(data_histo ,"Data","e0p");
    leg->AddEntry(h_mc,"POWHEGV2 + PYTHIA8","l");
    
    for (int p = 0; p< prediction_names.size(); p++){
       // leg->AddEntry(pred_graphs[p], prediction_names[p].c_str(),"p");
        leg->AddEntry(dummy_hists[p], prediction_names[p].c_str(),"e0p");
    }
    
    

    leg->SetTextSize(0.041);
    leg->Draw();
    
    pad2->cd();
    TLegend *leg_2 = new TLegend(plot_params[8], plot_params[9], plot_params[10], plot_params[11]);
    leg_2->SetBorderSize(0);
    leg_2->SetFillStyle(0);
    leg_2->SetTextSize(0.09);
    leg_2->AddEntry(g_ratio_uncs[0], "Stat #oplus Syst","f");
    leg_2->AddEntry(g_ratio_uncs[1], "Stat","f");
    leg_2->Draw();
    pad2->RedrawAxis();
    pad2->SetTicky(1);
    pad2->SetGrid(0,1);
    
    c->cd();
    c->SetTopMargin(0.1);
    c->SetBottomMargin(0.15);
    c->SetLeftMargin(0.0995);
    c->SetRightMargin(0.09);
    
    if (x_axis_title == "m_{t#bar{t}} [GeV]"){
        std::cout <<"found mtt"  << std::endl;
        pad1->SetLogx();
        pad2->SetLogx();
    }
    
    float H = c->GetWh();
    float W = c->GetWw();
    float l = c->GetLeftMargin();
    float t = c->GetTopMargin();
    float r = c->GetRightMargin();
    float b = c->GetBottomMargin();
    float extraOverCmsTextSize  = 0.89;
    
    TString cmsText, extraText, lumiText;
    cmsText += "CMS";
    extraText += "";
    lumiText += "35.9 fb^{-1} (13 TeV)";
    
    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextSize(0.43*t);
    latex.SetTextColor(kBlack);
    latex.SetTextFont(61);
    latex.SetTextAlign(31);
    latex.DrawLatex(0.29,0.905,cmsText);
    
    latex.SetTextFont(52);
    latex.SetTextSize(0.43*t*extraOverCmsTextSize);
    latex.DrawLatex(0.582,0.905,extraText);
    
    latex.SetTextFont(42);
    latex.SetTextSize(0.39*t);
    latex.DrawLatex(0.94,0.905,lumiText);
    
    latex.SetTextSize(0.39*t);
    latex.DrawLatex(0.56,0.84,"Dilepton, parton level");
    
    std::string pdf_title = observable +"_" + type + ".pdf";
    c->SaveAs(pdf_title.c_str());

}

TGraphAsymmErrors* read_prediction(std::string filename, bool divide_binwidth, std::string type){
    //loop over text file
    //first fill vectors from text file
    std::ifstream readFile(filename.c_str());
    std::vector<std::vector<std::string>> matrix;
    std::vector<std::string> tokens;
    std::string line;
    std::string delimiter = " ";
    
    while(getline(readFile,line)){
        tokens.clear();
        std::istringstream iss(line);
        size_t pos = 0;
        std::string token;
        
        while ((pos = line.find(delimiter)) != std::string::npos){
            token = line.substr(0, pos);
            //std::cout <<" mini-loop token "<< token << std::endl;
            line.erase(0, pos + delimiter.length());
            tokens.push_back(token);
        }
        matrix.push_back(tokens);
    }
    readFile.close();

    int nbins = matrix.size();

    //std::cout <<" N bins =  " <<  nbins <<std::endl;
    float bins[nbins+1];
    
    for (int i = 0; i < nbins; i++){
        bins[i] = std::stof(matrix[i][0]) - std::stof(matrix[i][8])/2.0;
    }

    bins[nbins] = std::stof(matrix[nbins-1][0]) + std::stof(matrix[nbins-1][8])/2.0;
    //std::cout <<" nbins " <<  nbins <<" bins  =  "<<  bins << std::endl;

    TGraphAsymmErrors * g = new TGraphAsymmErrors();
    double central_value, bin_centre, bin_width, unc_hi, unc_lo, exl,exh,eyl,eyh;
    float integral = 0.0;
    
    for (int i = 0; i < nbins; i++){
        bin_centre = std::stof(matrix[i][0]);
        central_value = std::stof(matrix[i][1]);
        bin_width = std::stof(matrix[i][8]);
        unc_hi = std::stof(matrix[i][7]);
        unc_lo = std::stof(matrix[i][6]);

        if(divide_binwidth == true){
        integral = integral + central_value;
        g->SetPoint(i,bin_centre, central_value / bin_width);
        g->SetPointError(i,bin_width/2.0,bin_width/2.0,(central_value - unc_lo)/bin_width ,(unc_hi - central_value)/bin_width);
        }else{
        integral = integral + (central_value*bin_width);
        g->SetPoint(i,bin_centre, central_value);
        g->SetPointError(i,bin_width/2.0,bin_width/2.0,(central_value-unc_lo),(unc_hi - central_value));
        }
        //std::cout<<"centre = " << bin_centre   << " width  = " << bin_width <<" unc_lo " <<  central_value-unc_lo  << " unc_hi " << unc_hi - central_value  << std::endl;
        }
    
    if (type == "norm"){
        for (int i = 0; i < nbins; i++){
            g->GetPoint(i,bin_centre,central_value);
            g->SetPoint(i,bin_centre, central_value/integral);
            exl = g->GetErrorXlow(i);
            exh = g->GetErrorXhigh(i);
            eyl = g->GetErrorYlow(i);
            eyh = g->GetErrorYhigh(i);
            g->SetPointError(i,exl,exh,eyl/integral,eyh/integral);
        }
    }
    //std::cout <<"***READING PREDICTION***   "<<  filename.c_str()  <<" sigma_tt (incl.) = "<< integral <<" pb"<<std::endl;

    return g;
}
