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
void summary_plot(string, string);
std::tuple <float, int, float > process_result(string,string, string, int, int);
std::tuple<float, int, float>  results;

TH1F * h_model, h_data;

TH1F * h_norm_bnlo_nnlo_ew_lux_A;
TH1F * h_norm_bnlo_nnlo_ew_nnpdf_A;
TH1F * h_norm_bnlo_nnlo_nnllprime_nnpdf_A;
TH1F * h_norm_bnlo_an3lo_nnpdf_A;
TH1F * h_norm_bnlo_annlo_ct14_A;

TH1F * h_norm_bnlo_nnlo_ew_lux_A_clone;
TH1F * h_norm_bnlo_nnlo_ew_nnpdf_A_clone;
TH1F * h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone;
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
        "predictions/NNLO_EW/LHC13-CMS-PTt-mTt2-NNPDF31-pheno.dat",
        "predictions/NNLLprime/LHC13-CMS-PTtav-mTt2-NNPDF31.dat",
        "predictions/aNNNLO/LHC13-CMS-PTtav-sqrt-mT-NNPDF30.dat",
        "predictions/aNNLO/LHC13-CMS-PTtav-mT-NNPDF30.dat"
    },
    //pt antitop

    {
        "predictions/NNLO_EW/LHC13-CMS-PTtx-mTtx2-LUXQED17-pheno.dat",
        "predictions/NNLO_EW/LHC13-CMS-PTtx-mTtx2-NNPDF31-pheno.dat",
        "predictions/NNLLprime/LHC13-CMS-PTtav-mTt2-NNPDF31.dat",
        "predictions/aNNNLO/LHC13-CMS-PTtav-sqrt-mT-NNPDF30.dat",
        "predictions/aNNLO/LHC13-CMS-PTtav-mT-NNPDF30.dat"
    },
    //y top
    {
    "predictions/NNLO_EW/LHC13-CMS-Yt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-Yt-HT4-NNPDF31-pheno.dat",
    "NA",
    "predictions/aNNNLO/LHC13-CMS-Yt-NNPDF3.0.dat",
    "predictions/aNNLO/LHC13-CMS-Yt-NNPDF3.0.dat"
    },
    //y antitop
    {
        "predictions/NNLO_EW/LHC13-CMS-Ytx-HT4-LUXQED17-pheno.dat",
        "predictions/NNLO_EW/LHC13-CMS-Ytx-HT4-NNPDF31-pheno.dat",
        "NA",
        "predictions/aNNNLO/LHC13-CMS-Yt-NNPDF3.0.dat",
        "predictions/aNNLO/LHC13-CMS-Yt-NNPDF3.0.dat"
    },
    //pt tt
    {
        "predictions/NNLO_EW/LHC13-CMS-PTtt-HT4-LUXQED17-pheno.dat",
        "predictions/NNLO_EW/LHC13-CMS-PTtt-HT4-NNPDF31-pheno.dat",
        "NA",
        "NA",
        "NA"
    },
    // Ytt

    {
        "predictions/NNLO_EW/LHC13-CMS-Ytt-HT4-LUXQED17-pheno.dat",
        "predictions/NNLO_EW/LHC13-CMS-Ytt-HT4-NNPDF31-pheno.dat",
        "NA",
        "NA",
        "NA"
    },
    // Mtt
    {
    "predictions/NNLO_EW/LHC13-CMS-Mtt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-Mtt-HT4-NNPDF31-pheno.dat",
    "predictions/NNLLprime/LHC13-CMS-Mtt-HT4-NNPDF31.dat",
        "NA",
        "NA"
    },
    //dytt
    {
    "predictions/NNLO_EW/LHC13-CMS-dytt-HT4-LUXQED17-pheno_rewrite.dat",
    "predictions/NNLO_EW/LHC13-CMS-dytt-HT4-NNPDF31-pheno.dat",
    "NA",
    "NA",
    "NA"
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
    "\\GeV",
    "",
    "\\GeV",
    "",
    "",
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
    "\\pt^{\\text{t}}",
    "\\pt^{\\tbar}",
    "\\pt^{\\text{t}} \\mathrm{(leading)} ",
    "\\pt^{\\text{t}} \\mathrm{(trailing)} ",
    "\\pt^{\\text{t}}\\mathrm{(} \\ttbar\\ \\mathrm{r.f.)}",
    "\\pt^{\\tbar} \\mathrm{(} \\ttbar\\ \\mathrm{r.f.)}",
    "y_{\\text{t}}",
    "y_{\\tbar}",
    "y_{\\text{t}} (\\mathrm{leading})",
    "y_{\\text{t}} (\\mathrm{trailing})",
    "\\pt^{\\ttbar}",
    "y_{\\ttbar}",
    "m_{\\ttbar}",
    "\\Delta|y| ( \\mathrm{t,\\tbar})",
    "\\Delta\\phi (\\mathrm{ t,\\tbar})",
    "\\pt^{\\text{l}}",
    "\\pt^{\\lbar}",
    "\\pt^{\\text{l}} \\mathrm{(leading)} ",
    "\\pt^{\\text{l}} \\mathrm{(trailing)} ",
    "\\eta_{\\text{l}}",
    "\\eta_{\\lbar}",
    "\\eta_{\\text{l}} \\mathrm{(leading)}",
    "\\eta_{\\text{l}} \\mathrm{(trailing)}",
    "\\pt^{\\llbar}",
    "m_{\\llbar}",
    "\\Delta\\phi\\ (\\mathrm{ l,\\lbar})",
    "\\Delta\\ \\eta (\\mathrm{ l,\\lbar})",
    "N_{\\mathrm{jets}}",
    "\\pt^{\\mathrm{b}} \\mathrm{(leading)}",
    "\\pt^{\\mathrm{b}} \\mathrm{(trailing)}",
    "\\eta_{\\mathrm{b}} \\mathrm{(leading)}",
    "\\eta_{\\mathrm{b}} \\mathrm{(trailing)}",
    "\\pt^{\\bbbar}",
    "m_{\\bbbar}"
};


std::vector<string> vars_root = {
    "p_{T}^{t}",
    "p_{T}^{#bar{t}}",
    "p_{T}^{t} \\mathrm{(leading)}",
    "p_{T}^{t} \\mathrm{(trailing)}",
    "p_{T}^{t} (t#bar{t} r.f.)",
    "p_{T}^{#bar{t}}(t#bar{t} r.f.)",
    "y_{t}",
    "y_{#bar{t}}",
    "y_{t} \\mathrm{(leading)}",
    "y_{t} \\mathrm{(trailing)}",
    "p_{T}^{t#bar{t}}",
    "y_{t#bar{t}}",
    "m_{t#bar{t}}",
    "#Delta|y|(t,#bar{t})",
    "#Delta#phi(t,#bar{t})",
    "p_{T}^{l}",
    "p_{T}^{#bar{l}}",
    "p_{T}^{l} \\mathrm{(leading)}",
    "p_{T}^{l} \\mathrm{(trailing)}",
    "#eta_{l}",
    "#eta_{#bar{l}}",
    "#eta_{l} \\mathrm{(leading)}",
    "#eta_{l} \\mathrm{(trailing)}",
    "p_{T}^{l#bar{l}}",
    "m_{l#bar{l}}",
    "#Delta#phi(l,#bar{l})",
    "#Delta#eta(l,#bar{l})",
    "N_{jets}",
    "p_{T}^{b} \\mathrm{(leading)}",
    "p_{T}^{b}  \\mathrm{(trailing)}",
    "#eta^{b}  \\mathrm{(leading)}",
    "#eta^{b} \\mathrm{(trailing)}",
    "p_{T}^{b#bar{b}}",
    "m_{b#bar{b}}"
};

std::vector<string> vars_bnlo = {
    "\\pt^{\\text{t}}",
    "\\pt^{\\tbar}",
    "Y t",
    "Y \\tbar",
    "\\pt \\ttbar\\",
    "y \\ttbar\\",
    "m_{\\ttbar}",
    "\\Delta\\ |y| \\ttbar\\"
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
    double bin_centre, bin_width, bin_edge, data_x, data_norm_y, data_abs_y, stat_unc_norm, sys_unc_norm, stat_unc_abs, sys_unc_abs;
    vector<double> bins;
    string varname;
    string xtitle;
    vector<std::string> tokens_vec;
    myfile_2.open (tex_filename);
    
    for (int f = 0; f < filenames.size(); f++){
        bins.clear();
        tokens_vec.clear();
        
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
      //  myfile_2 << "\\small"<<endl;
        myfile_2 << "\\centering"<<endl;
        
        if (mode == "norm_parton"){
    myfile_2 << "\\caption{The measured differential cross section and bin boundaries for each bin of the normalised and absolute measurements of the tt differential cross section at parton level in the full phase space as a function of $" << varname <<"$ are tabulated.}"<<endl;
        }else{
                myfile_2 << "\\caption{The measured differential cross section and bin boundaries for each bin of the normalised and absolute measurements of the tt differential cross section at particle level in the fiducial phase space as a function of $" << varname <<"$ are tabulated.}"<<endl;
        }
            
        myfile_2 << "\\begin{tabular}{| c | c | c |}"<<endl;
        myfile_2 << "\\hline"<<endl;
        
        cout <<" xtitle = "<< xtitle << endl;
        cout <<" varname  = "<< varname << endl;
        cout <<" "<< endl;
        
        if (unit == ""){
            myfile_2<< "$" << varname << " $ & $\\frac{1}{\\sigma}$ $\\frac{\\mathrm{d}\\sigma}{\\mathrm{d}"<< varname <<"}$ & $\\frac{\\mathrm{d}\\sigma}{\\mathrm{d}"<< varname <<"}$ [pb] \\\\ "<<endl;
        }else{
            myfile_2  <<"$"<< varname << " $ ["<< unit << "] $ & $\\frac{1}{\\sigma}$ $\\frac{\\mathrm{d}\\sigma}{\\mathrm{d}"<< varname <<"}$ ["<< unit <<"$^{-1}$] & $\\frac{\\mathrm{d}\\sigma}{\\mathrm{d}"<< varname <<"}$ [pb/"<< unit <<"] \\\\ "<<endl;
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
            sys_unc_norm = g_data_norm->GetErrorY(bin-1);
            stat_unc_abs = g_data_abs_stat_unc_only->GetErrorY(bin-1);
            sys_unc_abs = g_data_abs->GetErrorY(bin-1);

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
            
            
            myfile_2<< std::defaultfloat << std::setprecision(4) << bins[bin-1]<<"--"<< bins[bin];
            
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
        myfile_2 << "\\hline"<<endl;
        myfile_2 << "\\end{tabular}"<<endl;
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











