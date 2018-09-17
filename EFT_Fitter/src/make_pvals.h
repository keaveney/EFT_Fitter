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

#include "TPaveText.h"
#include "TF1.h"
#include <utility>
#include <tuple>

#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
//#include "helper_tools.h"
#include "compare_theory.h"

//#include "tdrstyle.C"
//#include "CMS_lumi.C"

//to run on beyond-nlo predictions just change n_models and change list of models in .cc
// summary plot does not yet work for b-nlo models.

using namespace std;
using namespace boost;

//const int n_models = 5;
//const int n_vars = 8;

const int n_models = 3;
const int n_vars = 34;

int make_table(string);
void write_latex(string, vector<string>, vector<string>, float[n_models][n_vars] , int[n_models][n_vars], double[n_models][n_vars]);

TH1D* graph2histo(TGraphAsymmErrors*);

void write_results_table(string, vector<string>);
void write_hepdata_tables(string, vector<string>, vector<string>);
void summary_plot(string, string);
std::tuple <float, int, float > process_result(string,string, string, int, int);
std::tuple<float, int, float>  results;

TH1F * h_model, h_data;

TH1F * h_norm_bnlo_nnlo_ew_lux_A;
TH1F * h_norm_bnlo_nnlo_ew_lux_1725_A;
TH1F * h_norm_bnlo_nnlo_ew_nnpdf_A;
TH1F * h_norm_bnlo_nnlo_nnllprime_nnpdf_A;
TH1F * h_norm_bnlo_nnlo_nnllprime_nnpdf_1725_A;
TH1F * h_norm_bnlo_an3lo_nnpdf_A;
TH1F * h_norm_bnlo_annlo_ct14_A;

TH1F * h_norm_bnlo_nnlo_ew_lux_A_clone;
TH1F * h_norm_bnlo_nnlo_ew_lux_1725_A_clone;
TH1F * h_norm_bnlo_nnlo_ew_nnpdf_A_clone;
TH1F * h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone;
TH1F * h_norm_bnlo_nnlo_nnllprime_nnpdf_1725_A_clone;
TH1F * h_norm_bnlo_an3lo_nnpdf_A_clone;
TH1F * h_norm_bnlo_annlo_ct14_A_clone;

TH1F * h_norm_parton_pwhg_p8_A;
TH1F * h_norm_particle_pwhg_p8_A;
TH1F * h_norm_particle_pwhg_p8_B;

TH1F * h_norm_parton_pwhg_hpp_A;
TH1F * h_norm_particle_pwhg_hpp_A;
TH1F * h_norm_particle_pwhg_hpp_B;

TH1F * h_norm_parton_amc_p8_A;
TH1F * h_norm_particle_amc_p8_A;
TH1F * h_norm_particle_amc_p8_B;

TH1F *h_norm_particle_pwhg_p8_A_clone;
TH1F *h_norm_parton_pwhg_p8_A_clone;
TH1F *h_norm_parton_pwhg_hpp_A_clone;
TH1F *h_norm_particle_pwhg_hpp_A_clone;
TH1F *h_norm_parton_amc_p8_A_clone;
TH1F *h_norm_particle_amc_p8_A_clone ;

TH1F *h_norm_particle_pwhg_p8_B_clone;
TH1F *h_norm_parton_pwhg_p8_B_clone;
TH1F *h_norm_parton_pwhg_hpp_B_clone;
TH1F *h_norm_particle_pwhg_hpp_B_clone;
TH1F *h_norm_parton_amc_p8_B_clone;
TH1F *h_norm_particle_amc_p8_B_clone ;

ofstream myfile;
ofstream myfile_2;

double low_lim = 0.001;

TH1D* graph2histo(TGraphAsymmErrors* g){
    double x, y;
    double nbins = g->GetN();
    TH1D * h = new TH1D("", "", nbins, 0.0, nbins);
    for(int i=0; i< nbins; i++){
        g->GetPoint(i, x, y);
        h->SetBinContent(i+1, y);
    }
    return h;
};

std::string bnlo_matrix[8][5] = {
    //pt top
    {   "predictions/NNLO_EW/LHC13-CMS-PTt-mTt2-LUXQED17-pheno.dat",
        "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-PTt-mTt2-LUXQED17-pheno.dat",
        "predictions/NNLO_EW/LHC13-CMS-PTt-mTt2-NNPDF31-pheno.dat",
        "predictions/NNLLprime/LHC13-CMS-PTtav-mTt2-NNPDF31.dat",
        "predictions/NNLLprime_172-5/LHC13-CMS-PTtav-mTt2-NNPDF31.dat"
      //  "predictions/aNNNLO/LHC13-CMS-PTtav-sqrt-mT-NNPDF30.dat",
      //  "predictions/aNNLO/LHC13-CMS-PTtav-mT-NNPDF30.dat"
    },
    //pt antitop

    {
        "predictions/NNLO_EW/LHC13-CMS-PTtx-mTtx2-LUXQED17-pheno.dat",
        "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-PTtx-mTtx2-LUXQED17-pheno.dat",
        "predictions/NNLO_EW/LHC13-CMS-PTtx-mTtx2-NNPDF31-pheno.dat",
        "predictions/NNLLprime/LHC13-CMS-PTtav-mTt2-NNPDF31.dat",
       "predictions/NNLLprime_172-5/LHC13-CMS-PTtav-mTt2-NNPDF31.dat"
   //     "predictions/aNNNLO/LHC13-CMS-PTtav-sqrt-mT-NNPDF30.dat",
    //    "predictions/aNNLO/LHC13-CMS-PTtav-mT-NNPDF30.dat"
    },
    //y top
    {
    "predictions/NNLO_EW/LHC13-CMS-Yt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-Yt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-Yt-HT4-NNPDF31-pheno.dat",
    "NA",
    "NA"
//    "predictions/aNNNLO/LHC13-CMS-Yt-NNPDF3.0.dat",
 //   "predictions/aNNLO/LHC13-CMS-Yt-NNPDF3.0.dat"
    },
    //y antitop
    {
        "predictions/NNLO_EW/LHC13-CMS-Ytx-HT4-LUXQED17-pheno.dat",
        "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-Ytx-HT4-LUXQED17-pheno.dat",
        "predictions/NNLO_EW/LHC13-CMS-Ytx-HT4-NNPDF31-pheno.dat",
        "NA",
        "NA"
      //  "predictions/aNNNLO/LHC13-CMS-Yt-NNPDF3.0.dat",
      //  "predictions/aNNLO/LHC13-CMS-Yt-NNPDF3.0.dat"
    },
    //pt tt
    {
    "predictions/NNLO_EW/LHC13-CMS-PTtt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-PTtt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-PTtt-HT4-NNPDF31-pheno.dat",
    "NA",
    "NA"
    //"NA",
    //    "NA"
    },
    // Ytt

    {
    "predictions/NNLO_EW/LHC13-CMS-Ytt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-Ytt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-Ytt-HT4-NNPDF31-pheno.dat",
    "NA",
    "NA"
    //    "NA",
    //    "NA"
    },
    // Mtt
    {
    "predictions/NNLO_EW/LHC13-CMS-Mtt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-Mtt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-Mtt-HT4-NNPDF31-pheno.dat",
    "predictions/NNLLprime/LHC13-CMS-Mtt-HT4-NNPDF31.dat",
    "predictions/NNLLprime_172-5/LHC13-CMS-Mtt-HT4-NNPDF31.dat"
  //      "NA",
  //      "NA"
    },
    //dytt
    {
    "predictions/NNLO_EW/LHC13-CMS-dytt-HT4-LUXQED17-pheno_rewrite.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-dytt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-dytt-HT4-NNPDF31-pheno.dat",
    "NA",
    "NA"
  //  "NA",
  //  "NA"
    }
};


std::vector<string> units = {
    "\\GeV",
    "\\GeV",
    "\\GeV",
    "\\GeV",
    "\\GeV",
    "\\GeV",
    "",
    "",
    "",
    "",
    "",
    "",
    "\\GeV",
    "",
    "\\GeV",
    "\\GeV",
    "\\GeV",
    "\\GeV",
    "\\GeV",
    "",
    "",
    "",
    "",
    "\\GeV",
    "\\GeV",
    "",
    "",
    "",
    "\\GeV",
    "\\GeV",
    "",
    "",
    "\\GeV",
    "\\GeV"
};



std::vector<string> vars = {
    "$\\pt^{\\text{t}}$",
    "$\\pt^{\\text{\\tbar}}$",
    "$\\pt^{\\text{t}} \\mathrm{(leading)} $",
    "$\\pt^{\\text{t}} \\mathrm{(trailing)} $",
    "$\\pt^{\\text{t}} \\mathrm{(} \\ttbar\\ \\mathrm{RF)}$",
    "$\\pt^{\\tbar} \\mathrm{(} \\ttbar\\ \\mathrm{RF)}$",
    "$y_{\\text{t}}$",
    "$y_{\\tbar}$",
    "$y_{\\text{t}} (\\mathrm{leading})$",
    "$y_{\\text{t}} (\\mathrm{trailing})$",
    "$\\Delta|y| ( \\mathrm{t,\\tbar})$",
    "$\\Delta\\phi (\\mathrm{ t,\\tbar})$",
    "$\\pt^{\\ttbar}$",
    "$y_{\\ttbar}$",
    "$m_{\\ttbar}$",
    "$\\pt^{\\ell}$",
    "$\\pt^{\\lbar}$",
    "$\\pt^{\\ell} \\mathrm{(leading)} $",
    "$\\pt^{\\ell} \\mathrm{(trailing)} $",
    "$\\eta_{\\ell}$",
    "$\\eta_{\\lbar}$",
    "$\\eta_{\\ell} \\mathrm{(leading)}$",
    "$\\eta_{\\ell} \\mathrm{(trailing)}$",
    "$\\pt^{\\llbar}$",
    "$m_{\\llbar}$",
    "$\\Delta\\phi\\ ($\\ell$, \\lbar)$",
    "$\\Delta\\ \\eta ($\\ell$,\\lbar})$",
    "$\Nj$",
    "$\\pt^{\\mathrm{b}} \\mathrm{(leading)}$",
    "$\\pt^{\\mathrm{b}} \\mathrm{(trailing)}$",
    "$\\eta_{\\mathrm{b}} \\mathrm{(leading)}$",
    "$\\eta_{\\mathrm{b}} \\mathrm{(trailing)}$",
    "$\\pt^{\\bbbar}$",
    "$m_{\\bbbar}$"
};

std::vector<string> vars_root = {
    "p_{T}^{t}",
    "p_{T}^{#bar{t}}",
    "p_{T}^{t} (leading)",
    "p_{T}^{t} (trailing)",
    "p_{T}^{t} (t#bar{t} RF)",
    "p_{T}^{#bar{t}}(t#bar{t} RF)",
    "y_{t}",
    "y_{#bar{t}}",
    "y_{t} (leading)",
    "y_{t} (trailing)",
    "#Delta|y|(t,#bar{t})",
    "#Delta#phi(t,#bar{t})",
    "p_{T}^{t#bar{t}}",
    "y_{t#bar{t}}",
    "m_{t#bar{t}}",
    "p_{T}^{l}",
    "p_{T}^{#bar{l}}",
    "p_{T}^{l} (leading)",
    "p_{T}^{l} (trailing)",
    "#eta_{l}",
    "#eta_{#bar{l}}",
    "#eta_{l} (leading)",
    "#eta_{l} (trailing)",
    "p_{T}^{l#bar{l}}",
    "m_{l#bar{l}}",
    "#Delta#phi(l,#bar{l})",
    "#Delta#eta(l,#bar{l})",
    "N_{jets}",
    "p_{T}^{b} (leading)",
    "p_{T}^{b} (trailing)",
    "#eta^{b} (leading)",
    "#eta^{b} (trailing)",
    "p_{T}^{b#bar{b}}",
    "m_{b#bar{b}}"
};




std::vector<string> vars_txt = {
    "pt_t",
    "pt_tbar",
    "pt_t_leading",
    "pt_t_trailing",
    "pt_t_ttrf",
    "pt_tbar_ttrf",
    "y_t",
    "y_tbar",
    "y_t_leading",
    "y_t_trailing",
    "dabsy_tt",
    "dphi_tt",
    "pt_tt",
    "y_tt",
    "m_tt",
    "pt_l",
    "pt_lbar",
    "pt_l_leading",
    "pt_l_trailing",
    "eta_l",
    "eta_lbar",
    "eta_l_leading",
    "eta_l_trailing",
    "pt_ll",
    "m_ll",
    "dphi_ll",
    "deta_ll",
    "Njets",
    "pt_b_leading",
    "pt_b_trailing",
    "eta_b_leading",
    "eta_b_trailing",
    "pt_bb",
    "m_bb"
};

std::vector<string> vars_obskey = {
    "PTT",
    "PTTBAR",
    "PTTLEADING",
    "PTTTRAILING",
    "PTTTTRF",
    "PTTBARTTRF",
    "YT",
    "YTBAR",
    "YTLEADING",
    "YTTRAILING",
    "DYTT",
    "DPHITT",
    "PTTT",
    "YTT",
    "MTT",
    "PTL",
    "PTLBAR",
    "PTLLEADING",
    "PTLTRAILING",
    "ETAL",
    "ETALBAR",
    "ETALLEADING",
    "ETALTRAILING",
    "PTLL",
    "MLL",
    "DPHILL",
    "DETALL",
    "NJETS",
    "PTBLEADING",
    "PTBTRAILING",
    "ETABLEADING",
    "ETABTRAILING",
    "PTBB",
    "MBB"
};



std::vector<string> vars_bnlo = {
    "$\\pt^{\\text{t}}$",
    "$\\pt^{\\tbar}$",
    "$y_{\\text{t}}$",
    "$y_{\\text{\\tbar}}$",
    "$\\pt^{\\ttbar}$",
    "$y^{\\ttbar}$",
    "$m_{\\ttbar}$",
    "$\\Delta|y|$ (t,\\tbar) "
};

std::vector<string> vars_root_bnlo = {
    "p_{T}^{t}",
    "p_{T}^{#bar{t}}",
    "y_{t}",
    "y_{#bar{t}}",
    "p_{T}^{t#bar{t}}",
    "y_{t#bar{t}}",
    "m_{t#bar{t}}",
    "#Delta|y|(t,#bar{t})"
};


void write_hepdata_tables(string mode, vector<string> filenames, vector<string> filenames_cov){
    cout<<"writing hepdata table,  mode =  "<< mode << endl;
    cout<<"writing hepdata table,  filenames[0] =  "<< filenames[0] << endl;

    ofstream HDfile;
    ofstream HDfile_cov;
    string filename;
    string filename_cov;
    string hepdata_filename = "";
    string hepdata_filename_cov = "";
    
    char_separator<char> sep("#");
    char_separator<char> sep_2("_");
    
    //extract data and bins
    TFile * f_var, *f_var_cov;
    TH1F * h_mc;
    TH2D * h_cov;

    TGraphAsymmErrors *g_data, *g_data_stat_unc_only;
    double bin_centre, bin_width, bin_edge, data_x, data_y, stat_unc, sys_unc, tot_unc;
    vector<double> bins;
    string varname;
    string xtitle;
    string xsec_type_text, xsec_type_text_2;
    vector<std::string> tokens_vec;
    vector<std::string> mode_tokens_vec;

    tokenizer< char_separator<char> > mode_tokens(mode, sep_2);
    BOOST_FOREACH (const string& t, mode_tokens) {
        mode_tokens_vec.push_back(t);
    }
    
    cout<<" tokens vec = " << mode_tokens_vec[0] << " " << mode_tokens_vec[1] <<endl;

    if (mode_tokens_vec[0] == "abs"){
        xsec_type_text = "Absolute";
        xsec_type_text_2 = "absolute";
    }else{
        xsec_type_text = "Normalized";
        xsec_type_text = "normalized";
    }
    
    for (int f = 0; f < filenames.size(); f++){
        filename = filenames[f];
        filename_cov = filenames_cov[f];
        hepdata_filename = "";
        hepdata_filename += "HEPData/";
        hepdata_filename += mode;
        hepdata_filename += "/";
        varname = vars_txt[f];
        string unit = units[f];
        hepdata_filename += varname;
        hepdata_filename += "_";
        hepdata_filename += "_hepdata_table";
        hepdata_filename_cov = hepdata_filename;
        hepdata_filename_cov += "_cov";
        hepdata_filename += ".oldhepdata";
        hepdata_filename_cov += ".oldhepdata";

        HDfile.open (hepdata_filename);
        
        bins.clear();
        tokens_vec.clear();
        string f_str = std::to_string(f);
        
        tokenizer< char_separator<char> > tokens_2(filenames[f], sep_2);
        BOOST_FOREACH (const string& t, tokens_2) {
            tokens_vec.push_back(t);
        }
        

        tokens_vec.clear();
        
        //get data files
        f_var = new TFile(filename.c_str());
        h_mc = (TH1F*)f_var->Get("mc");
        g_data = (TGraphAsymmErrors*)f_var->Get("data");
        g_data_stat_unc_only = (TGraphAsymmErrors*)f_var->Get("data_staterror_only");
        
        
        xtitle = h_mc->GetXaxis()->GetTitle();
        tokenizer< char_separator<char> > tokens(xtitle, sep);
        BOOST_FOREACH (const string& t, tokens){
            tokens_vec.push_back(t);
        }


        
        //write data
        for (int bin = 1; bin <= h_mc->GetNbinsX(); bin++){
            bin_centre = h_mc->GetBinCenter(bin);
            bin_width = h_mc->GetBinWidth(bin);
            bin_edge = bin_centre - (bin_width/2.0);
            bins.push_back(bin_edge);
        }
        
        bins.push_back(bin_centre + (bin_width/2.0));
        
        //first write header metadata
        HDfile<<"*dataset: "<<endl;
        HDfile<<"*location: Fig. X, Tab. Y: "<<endl; //does this matter?
        HDfile<<"*dscomment: "<< xsec_type_text <<" cross section at parton level."<<endl; //needs to be changed according to mode variable
        HDfile<<"*reackey: P P --> TOP TOPBAR X"<<endl;
        HDfile<<"*obskey: DSIG/D"<<  vars_obskey[f]  <<endl; //needs to be changed according to var variable
        HDfile<<"*qual: SQRT(S) IN GEV : 13000.0"<<endl;
        //HDfile<<"*yheader: $\\frac{d\\sigma}{d" << vars_root[f] << "}$ [pb GeV$^{-1}$]"<<endl; //this needs to be dynamically constructed
        string hd_unit;
        cout <<"HD UNIT   == = = = = = = "<< units[f] << endl;

        if (units[f] == "\\GeV"){
            cout <<" GeV units " << endl;
            hd_unit = "GeV";
        }else{
            cout <<" no units " << endl;
            hd_unit =  "";
        }
        
        HDfile<<"*yheader: $\\frac{d\\sigma}{d " << vars_root[f] << "}$ [" << hd_unit << "]"<<endl;
        HDfile<<"*xheader: Bin : $" << vars_root[f]  << "$ ["<< hd_unit <<"]"<<endl;
        HDfile<<"*data: x : x : y"<<endl;

        for (int bin = 1; bin < bins.size(); bin++){
            g_data->GetPoint(bin-1, data_x, data_y);
            
            stat_unc = g_data_stat_unc_only->GetErrorY(bin-1);
            tot_unc = g_data->GetErrorY(bin-1);
            sys_unc = pow( (  (pow(tot_unc,2)) -  (pow(stat_unc,2))    ), 0.5     );
 
            HDfile<< bin <<"; "<< std::defaultfloat << std::setprecision(4) << bins[bin-1]<<" TO "<< bins[bin] << "; " << data_y << " +- "<< stat_unc << " (DSYS="<<  sys_unc<< ");"<< endl;
        }
        HDfile<<"*dataend"<<endl;

        HDfile.close();
        
        
        //write covariances
        HDfile_cov.open (hepdata_filename_cov);
        //get covariance matrix
        f_var_cov = new TFile(filename_cov.c_str());
        h_cov = (TH2D*)f_var_cov->Get("cov");
        
        
        //first write header metadata
        HDfile_cov<<"*dataset: "<<endl;
        HDfile_cov<<"*location: Fig. X, Tab. Y: "<<endl;
        HDfile_cov<<"*dscomment:Covariance matrix of "<< xsec_type_text_2 <<" cross section at parton level."<<endl;
        HDfile_cov<<"*reackey: P P --> TOP TOPBAR X"<<endl;
        HDfile_cov<<"*obskey: DSIG/D"<<  vars_obskey[f]  <<endl;
        HDfile_cov<<"*qual: COVARIANCE MATRIX" <<endl;
        HDfile_cov<<"*yheader";
        for (int xbin = 1; xbin <= h_cov->GetNbinsX(); xbin++)HDfile_cov<<": "<< xbin;
        HDfile_cov<<""<<endl;
        HDfile_cov<<"*xheader: Bin"<<endl;
        HDfile_cov<<"*data: x";

        for (int xbin = 1; xbin <= h_cov->GetNbinsX(); xbin++)HDfile_cov<<": y";
        HDfile_cov<<""<<endl;

        for (int xbin = 1; xbin <= h_cov->GetNbinsX(); xbin++){
            HDfile_cov<<xbin<<"; ";

            for (int ybin = 1; ybin <= h_cov->GetNbinsY(); ybin++){
                HDfile_cov<< std::defaultfloat << std::setprecision(15) << h_cov->GetBinContent(xbin,ybin) <<  "; ";
            }
            HDfile_cov<<""<<endl;

        }
        HDfile_cov<<"*dataend"<<endl;
    
        HDfile_cov.close();


    }
}



void write_results_table(string mode, vector<string> filenames){
    
    cout<<"writing results table"<< endl;
    
    string tex_filename;
    string filename_norm;
    string filename_abs;

    char_separator<char> sep("#");
    char_separator<char> sep_2("/");

    tex_filename += mode;
    tex_filename += "_results_table";
    tex_filename += ".tex";
    cout<<"tex filename = " << tex_filename <<endl;
    
    //extract data and bins
    TFile * f_var_norm, * f_var_abs;
    TH1F * h_mc_norm, * h_mc_abs;
    TGraphAsymmErrors * g_data_norm, * g_data_abs,  * g_data_norm_stat_unc_only,  * g_data_abs_stat_unc_only;
    double bin_centre, bin_width, bin_edge, data_x, data_norm_y, data_abs_y, stat_unc_norm, sys_unc_norm, stat_unc_abs, sys_unc_abs, tot_unc_norm, tot_unc_abs;
    vector<double> bins;
    string varname;
    string xtitle;
    vector<std::string> tokens_vec;
    myfile_2.open (tex_filename);
    
    if (mode == "norm_parton"){
        myfile_2 << "All the measured differential cross sections at the parton level are tabulated in Tables \\ref{tab:norm_parton0}--\\ref{tab:norm_parton14}. The statistical and systematic uncertainties are quoted seperately for each bin."<<endl;
    }else{
        myfile_2 << "All the measured differential cross sections at the particle level are tabulated in Tables \\ref{tab:norm_particle0}--\\ref{tab:norm_particle33}. The statistical and systematic uncertainties are quoted seperately for each bin."<<endl;
    }
    
    for (int f = 0; f < filenames.size(); f++){
        bins.clear();
        tokens_vec.clear();
        string f_str = std::to_string(f);
        
        tokenizer< char_separator<char> > tokens_2(filenames[f], sep_2);
        BOOST_FOREACH (const string& t, tokens_2) {
            tokens_vec.push_back(t);
        }
        
        filename_norm = tokens_vec[0] + "/" + tokens_vec[1] + "/"  + tokens_vec[2] + "/normalised/" +tokens_vec[4];
        filename_abs = tokens_vec[0] + "/" + tokens_vec[1] + "/"  + tokens_vec[2] + "/absolute/" +tokens_vec[4];
        tokens_vec.clear();
        
        f_var_norm = new TFile(filename_norm.c_str());
        h_mc_norm = (TH1F*)f_var_norm->Get("mc");
        g_data_norm = (TGraphAsymmErrors*)f_var_norm->Get("data");
        g_data_norm_stat_unc_only = (TGraphAsymmErrors*)f_var_norm->Get("data_staterror_only");
        
        f_var_abs = new TFile(filename_abs.c_str());
        h_mc_abs = (TH1F*)f_var_abs->Get("mc");
        g_data_abs = (TGraphAsymmErrors*)f_var_abs->Get("data");
        g_data_abs_stat_unc_only = (TGraphAsymmErrors*)f_var_abs->Get("data_staterror_only");
        
        xtitle = h_mc_norm->GetXaxis()->GetTitle();
        tokenizer< char_separator<char> > tokens(xtitle, sep);
        BOOST_FOREACH (const string& t, tokens){
            tokens_vec.push_back(t);
        }
        varname   = vars[f];
        string unit = units[f];

        myfile_2 << "\\begin{table}"<<endl;
        //myfile_2 << "\\small"<<endl;
        myfile_2 << "\\centering"<<endl;
        
        if (mode == "norm_parton"){
            myfile_2 << "\\caption{The measured differential cross section and bin boundaries for each bin of the normalised and absolute measurements of the \\ttbar\\ differential cross section at parton level in the full phase space as a function of $" << varname <<"$ are tabulated.}"<<endl;
        }else{
                myfile_2 << "\\caption{The measured differential cross section and bin boundaries for each bin of the normalised and absolute measurements of the \\ttbar\\ differential cross section at particle level in the fiducial phase space as a function of $" << varname <<"$ are tabulated.}"<<endl;
        }
            
        myfile_2 << "\\begin{tabular}{| c | c | c |}"<<endl;
        myfile_2 << "\\hline"<<endl;
        
        cout <<" xtitle = "<< xtitle << endl;
        cout <<" varname  = "<< varname << endl;
        cout <<" "<< endl;
        
        if (unit == ""){
            myfile_2<< "$" << varname << " $ & $\\frac{1}{\\sigma}$ $\\frac{\\mathrm{d}\\sigma}{\\mathrm{d}"<< varname <<"}$ & $\\frac{\\mathrm{d}\\sigma}{\\mathrm{d}"<< varname <<"}$ [pb] \\\\ "<<endl;
        }else{
            myfile_2  <<"$"<< varname << " $ ["<< unit << "]  & $\\frac{1}{\\sigma}$ $\\frac{\\mathrm{d}\\sigma}{\\mathrm{d}"<< varname <<"}$ ["<< unit <<"$^{-1}$] & $\\frac{\\mathrm{d}\\sigma}{\\mathrm{d}"<< varname <<"}$ [pb/"<< unit <<"] \\\\ "<<endl;
        }
        
        myfile_2 << "\\hline"<<endl;
        string exp_str_norm = "";
        string exp_str_abs = "";

        for (int bin = 1; bin <= h_mc_norm->GetNbinsX(); bin++){
            bin_centre = h_mc_norm->GetBinCenter(bin);
            bin_width = h_mc_norm->GetBinWidth(bin);
            bin_edge = bin_centre - (bin_width/2.0);
            bins.push_back(bin_edge);
        }
        
        bins.push_back(bin_centre + (bin_width/2.0));
        
        bool scify_norm = true;
        bool scify_abs = true;

        for (int bin = 1; bin < bins.size(); bin++){
            
            g_data_norm->GetPoint(bin-1, data_x, data_norm_y);
            g_data_abs->GetPoint(bin-1, data_x, data_abs_y);

            stat_unc_norm = g_data_norm_stat_unc_only->GetErrorY(bin-1);
            tot_unc_norm = g_data_norm->GetErrorY(bin-1);
            stat_unc_abs = g_data_abs_stat_unc_only->GetErrorY(bin-1);
            tot_unc_abs = g_data_abs->GetErrorY(bin-1);

            sys_unc_norm = pow( (  (pow(tot_unc_norm,2)) -  (pow(stat_unc_norm,2))    ), 0.5     );
            sys_unc_abs = pow( (  (pow(tot_unc_abs,2)) -   (pow(stat_unc_abs,2))   ), 0.5         );
            
            if (data_norm_y < 0.01 && data_norm_y > 0.001){
                exp_str_norm = " $ \\times 10^{-3}$ ";
                data_norm_y = (data_norm_y*1000.0);
                stat_unc_norm = (stat_unc_norm*1000.0);
                sys_unc_norm = (sys_unc_norm*1000.0);
            } else if (data_norm_y < 0.001 && data_norm_y > 0.0001) {
                exp_str_norm = " $ \\times 10^{-4}$ ";
                data_norm_y = (data_norm_y*10000.0);
                stat_unc_norm = (stat_unc_norm*10000.0);
                sys_unc_norm = (sys_unc_norm*10000.0);
            }else if (data_norm_y < 0.0001 && data_norm_y > 0.00001) {
                exp_str_norm = " $ \\times 10^{-5}$ ";
                data_norm_y = (data_norm_y*100000.0);
                stat_unc_norm = (stat_unc_norm*100000.0);
                sys_unc_norm = (sys_unc_norm*100000.0);
            }else if (data_norm_y < 0.00001 && data_norm_y > 0.000001) {
                exp_str_norm = " $ \\times 10^{-6}$ ";
                data_norm_y = (data_norm_y*1000000.0);
                stat_unc_norm = (stat_unc_norm*1000000.0);
                sys_unc_norm = (sys_unc_norm*1000000.0);
            } else {
                scify_norm = false;
            }
            
            if (data_abs_y < 0.01 && data_abs_y > 0.001){
                exp_str_abs = " $ \\times 10^{-3}$ ";
                data_abs_y = (data_abs_y*1000.0);
                stat_unc_abs = (stat_unc_abs*1000.0);
                sys_unc_abs = (sys_unc_abs*1000.0);
            } else if (data_abs_y < 0.001 && data_abs_y > 0.0001) {
                exp_str_abs = " $ \\times 10^{-4}$ ";
                data_abs_y = (data_abs_y*10000.0);
                stat_unc_abs = (stat_unc_abs*10000.0);
                sys_unc_abs = (sys_unc_abs*10000.0);
            }else if (data_abs_y < 0.0001 && data_abs_y > 0.00001) {
                exp_str_abs = " $ \\times 10^{-5}$ ";
                data_abs_y = (data_abs_y*100000.0);
                stat_unc_abs = (stat_unc_abs*100000.0);
                sys_unc_abs = (sys_unc_abs*100000.0);
            }else if (data_abs_y < 0.00001 && data_abs_y > 0.000001) {
                exp_str_abs = " $ \\times 10^{-6}$ ";
                data_abs_y = (data_abs_y*1000000.0);
                stat_unc_abs = (stat_unc_abs*1000000.0);
                sys_unc_abs = (sys_unc_abs*1000000.0);
            } else {
                scify_abs = false;
            }

            myfile_2<< std::defaultfloat << std::setprecision(4) << bins[bin-1]<<"  --  "<< bins[bin];
            
            if (scify_norm && scify_abs){
                myfile_2 << std::fixed << std::setprecision(4) << " & ("<< data_norm_y  <<" $\\pm$ "<<  stat_unc_norm <<" $\\pm$ "<< sys_unc_norm << ") "<<exp_str_norm << " & ("<< data_abs_y <<" $\\pm$ "<<  stat_unc_abs << " $\\pm$ " << sys_unc_abs <<") "<< exp_str_abs << " \\\\"<<endl;
            }else if (scify_norm && !scify_abs) {
                myfile_2 << std::fixed << std::setprecision(4) << " & ("<< data_norm_y  <<" $\\pm$ "<<  stat_unc_norm <<" $\\pm$ "<< sys_unc_norm << ") "<<exp_str_norm << " & "<< data_abs_y <<" $\\pm$ "<<  stat_unc_abs << " $\\pm$ " << sys_unc_abs  << " \\\\"<<endl;
            }else if (!scify_norm && scify_abs) {
                myfile_2 << std::fixed << std::setprecision(4) << " & "<< data_norm_y  <<" $\\pm$ "<<  stat_unc_norm <<" $\\pm$ "<< sys_unc_norm << " & ("<< data_abs_y <<" $\\pm$ "<<  stat_unc_abs << " $\\pm$ " << sys_unc_abs <<") "<< exp_str_abs << " \\\\"<<endl;
            }
            else if (!scify_norm && !scify_abs) {
                myfile_2 << std::fixed<< std::setprecision(4) << " & "<< data_norm_y  <<" $\\pm$ "<<  stat_unc_norm <<" $\\pm$ "<< sys_unc_norm << " & " << data_abs_y  <<" $\\pm$ "<<  stat_unc_abs << " $\\pm$ " << sys_unc_abs <<" \\\\" <<endl;
            }
        }
        
        string lab_str = mode + f_str;
        myfile_2 << "\\hline"<<endl;
        myfile_2 << "\\end{tabular}"<<endl;
        myfile_2 << "\\label{tab:"<< lab_str <<"}"<<endl;
        myfile_2 << "\\end{table}"<<endl;
        myfile_2 << "  "<<endl;
    //    myfile_2 << "\\clearpage"<<endl;
        myfile_2 << "  "<<endl;
    }
    
    myfile_2.close();

    /* string filename = mode + ".tex";
    
     myfile_2.precision(2);
     myfile_2.open (filename);
     myfile_2 << "\\begin{table}"<<endl;
     myfile_2 << "\\small"<<endl;
     myfile_2 << "\\centering"<<endl;
     myfile_2 << "\\begin{tabular}{| l | c | c | c | c | c | c | c | c | c | c | c |}"<<endl;
     myfile_2 << "\\hline"<<endl;
     */

}


void plot_remaker(string filename, string mode){
    
    string target_path = "/Users/keaveney/Desktop/testplots/";
    
    gStyle->SetOptStat(00000);
    gStyle->SetEndErrorSize(0);

    TFile * f =  new TFile(filename.c_str());
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
    

    h_mc->GetYaxis()->SetTitleSize(0.07);
    h_mc->GetYaxis()->SetTitleOffset(0.95);
    h_mc->SetTitle("");
    h_mc->SetLineColor(kRed+1);
    
    h_mc->Draw();
    //h_mc->GetYaxis()->SetRangeUser(plot_params[0], plot_params[1]);
    
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
  //  h_ratio_base->GetYaxis()->SetRangeUser(plot_params[2], plot_params[3]);
    h_ratio_base->GetYaxis()->SetNdivisions(505);
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
    
    std::string x_axis_title;
   // if(units == ""){
  //      x_axis_title = observable;
  //  }else{
  //      x_axis_title = observable + " [" + units + "]";
  //  }
    std::cout <<"***OBSERVABLE****   "  <<    x_axis_title <<std::endl;
    h_ratio_base->SetXTitle(x_axis_title.c_str());
    
    pad1->cd();
    pad1->SetLogy();
    g->SetMarkerSize(1.3);
    g->Draw("E0PSAME");
    
   // TLegend *leg = new TLegend(plot_params[4], plot_params[5], plot_params[6], plot_params[7]);
    TLegend *leg = new TLegend(0.4, 0.4, 0.6, 0.6);
    leg->SetBorderSize(0);
    leg->AddEntry(g ,"Data","P");
    leg->AddEntry(h_mc,"Powheg v2+Pythia8","l");
    leg->SetTextSize(0.033);
    leg->Draw();
    
    pad2->cd();
    TLegend *leg_2 = new TLegend(0.4, 0.4, 0.6, 0.6);
    leg_2->SetBorderSize(0);
    leg_2->SetTextSize(0.09);
    leg_2->AddEntry(g_ratio_uncs[0], "Stat. #oplus Syst.","f");
    leg_2->AddEntry(g_ratio_uncs[1], "Stat.","f");
    leg_2->Draw();
    pad2->RedrawAxis();
    pad2->SetTicky(1);
    pad2->SetGrid(0,1);
    
    c->cd();
    c->SetTopMargin(0.1);
    c->SetBottomMargin(0.16);
    c->SetLeftMargin(0.0995);
    c->SetRightMargin(0.09);
    
    float H = c->GetWh();
    float W = c->GetWw();
    float l = c->GetLeftMargin();
    float t = c->GetTopMargin();
    float r = c->GetRightMargin();
    float b = c->GetBottomMargin();
    float extraOverCmsTextSize  = 0.89;
    
    TString cmsText, extraText, lumiText;
    cmsText += "CMS";
    extraText += "Preliminary";
    lumiText += "35.9 fb^{-1} (13 TeV)";
    
    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextSize(0.43*t);
    latex.SetTextColor(kBlack);
    latex.SetTextFont(61);
    latex.SetTextAlign(31);
    latex.DrawLatex(0.34,0.905,cmsText);
    
    latex.SetTextFont(52);
    latex.SetTextSize(0.43*t*extraOverCmsTextSize);
    latex.DrawLatex(0.582,0.905,extraText);
    
    latex.SetTextFont(42);
    latex.SetTextSize(0.39*t);
    latex.DrawLatex(0.915,0.905,lumiText);
    
    latex.SetTextSize(0.36*t);
    latex.DrawLatex(0.47,0.83,"Dilepton: parton");
    
    std::string pdf_title = target_path +"_" + mode + ".pdf";
    c->SaveAs(pdf_title.c_str());
    
    
}








