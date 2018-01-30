#include <iostream>
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
#include "helper_tools.h"

#include "tdrstyle.C"
//#include "CMS_lumi.C"


using namespace std;
using namespace boost;

int make_table(string);
void write_latex(string, vector<string>, vector<string>, float[3][34] , int[3][34], float[3][34]);
void summary_plot(string, string);
std::tuple <float, int, float > process_result(string,string, string);

TH1F * h_model, h_data;

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

double low_lim = 0.001;

std::vector<string> vars = {
    "\\pt t",
    "\\pt \\tbar",
    "\\pt t (leading)",
    "\\pt t (trailing)",
    "\\pt t (\\ttbar\\ r.f.)",
    "\\pt \\bar{t} (\\ttbar\\ r.f.)",
    "Y t",
    "Y \\tbar",
    "Y t (leading)",
    "Y t (trailing)",
    "\\pt \\ttbar\\",
    "Y \\ttbar\\",
    "M_{\\ttbar}",
    "\\Delta\\ Y \\ttbar\\",
    "\\Delta\\ \\Phi\\ \\ttbar\\",
    "\\pt^{l}",
    "\\pt^{\\lbar}",
    "\\pt^{l (leading)}",
    "\\pt^{l (trailing)}",
    "\\eta^{l}",
    "\\eta^{\\lbar}",
    "\\eta^{l (leading)}",
    "\\eta^{l (trailing)}",
    "\\pt \\llbar\\",
    "M_{\\llbar}",
    "\\Delta\\ \\Phi\\ (\\llbar)",
    "\\Delta\\ \\eta\\ (\\llbar)",
    "N_{jets}",
    "\\pt^{b-jet (leading)}",
    "\\pt^{b-jet (trailing)}",
    "\\eta^{b-jet (leading)}",
    "\\eta^{b-jet (trailing)}",
    "\\pt^{\\bbbar}",
    "M_{\\bbbar}"
};

std::vector<string> vars_root = {
    "p_{T}^{t}",
    "p_{T}^{#bar{t}}",
    "p_{T}^{t} (leading)",
    "p_{T}^{t} (trailing)",
    "p_{T}^{t} (tt r.f.)",
    "p_{T}^{#bar{t}}(tt r.f.)",
    "Y_{t}",
    "Y_{#bar{t}}",
    "Y_{t} (leading)",
    "Y_{t} (trailing)",
    "p_{T}^{tt}",
    "Y_{t#bar{t}}",
    "m_{t#bar{t}}",
    "#Delta Y (t,#bar{t})",
    "#Delta#Phi (t,#bar{t})",
    "p_{T}^{l}",
    "p_{T}^{#bar{l}}",
    "p_{T}^{l} (leading)",
    "p_{T}^{l} (trailing)",
    "#eta_{l}",
    "#eta_{#bar{l}}",
    "#eta_{l} (leading)",
    "#eta_{l} (trailing)",
    "p_{T}^{l #bar{l}}",
    "m_{l#bar{l}}",
    "#Delta#Phi(l,#bar{l})",
    "#Delta#eta(l,#bar{l})",
    "N_{jets}",
    "p_{T}^{b-jet} (leading)",
    "p_{T}^{b-jet} (trailing)",
    "#eta^{b-jet} (leading)",
    "#eta^{b-jet} (trailing)",
    "p_{T}^{b#bar{b}}",
    "m_{b#bar{b}}"
};


