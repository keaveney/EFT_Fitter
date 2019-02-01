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
#include "TPaletteAxis.h"
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
#include "compare_theory.h"

//to run on beyond-nlo predictions just change n_models and change list of models in .cc

using namespace std;
using namespace boost;

//const int n_models = 7;
//const int n_vars = 8;

//const int n_models = 3;
//const int n_vars = 33;

int make_table(string);
void write_latex(string, vector<string>, vector<string>, vector<vector<float>> , vector<vector<int>>, vector<vector<double>>);
vector<std::string> make_covariance_matrices(std::string, vector<std::string>);


TH1D* graph2histo(TGraphAsymmErrors*);

void write_results_table(string, vector<string>);
void write_hepdata_tables(string, vector<string>, vector<string>);
string round_off(double, double);
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


std::string bnlo_matrix[8][7] = {
    //pt top
    {   "predictions/NNLO_EW/LHC13-CMS-PTt-mTt2-LUXQED17-pheno.dat",
        "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-PTt-mTt2-LUXQED17-pheno.dat",
        "predictions/NNLO_EW/LHC13-CMS-PTt-mTt2-NNPDF31-pheno.dat",
        "predictions/NNLLprime/LHC13-CMS-PTtav-mTt2-NNPDF31.dat",
        "predictions/NNLLprime_172-5/LHC13-CMS-PTtav-mTt2-NNPDF31.dat",
        "predictions/aNNNLO/LHC13-CMS-PTtav-sqrt-mT-NNPDF30.dat",
        "predictions/aNNLO/LHC13-CMS-PTtav-mT-NNPDF30.dat"
    },
    //pt antitop
    {
        "predictions/NNLO_EW/LHC13-CMS-PTtx-mTtx2-LUXQED17-pheno.dat",
        "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-PTtx-mTtx2-LUXQED17-pheno.dat",
        "predictions/NNLO_EW/LHC13-CMS-PTtx-mTtx2-NNPDF31-pheno.dat",
        "NA",
        "NA",
        "NA",
        "NA"
    },
    //y top
    {
    "predictions/NNLO_EW/LHC13-CMS-Yt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-Yt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-Yt-HT4-NNPDF31-pheno.dat",
    "NA",
    "NA",
    "predictions/aNNNLO/LHC13-CMS-Yt-NNPDF3.0.dat",
    "predictions/aNNLO/LHC13-CMS-Yt-NNPDF3.0.dat"
    },
    //y antitop
    {
        "predictions/NNLO_EW/LHC13-CMS-Ytx-HT4-LUXQED17-pheno.dat",
        "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-Ytx-HT4-LUXQED17-pheno.dat",
        "predictions/NNLO_EW/LHC13-CMS-Ytx-HT4-NNPDF31-pheno.dat",
        "NA",
        "NA",
        "predictions/aNNNLO/LHC13-CMS-Yt-NNPDF3.0.dat",
        "predictions/aNNLO/LHC13-CMS-Yt-NNPDF3.0.dat"
    },
    //pt tt
    {
    "predictions/NNLO_EW/LHC13-CMS-PTtt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-PTtt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-PTtt-HT4-NNPDF31-pheno.dat",
    "NA",
    "NA",
    "NA",
    "NA"
    },
    // Ytt

    {
    "predictions/NNLO_EW/LHC13-CMS-Ytt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-Ytt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-Ytt-HT4-NNPDF31-pheno.dat",
    "NA",
    "NA",
    "NA",
    "NA"
    },
    // Mtt
    {
    "predictions/NNLO_EW/LHC13-CMS-Mtt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-Mtt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-Mtt-HT4-NNPDF31-pheno.dat",
    "predictions/NNLLprime/LHC13-CMS-Mtt-HT4-NNPDF31.dat",
    "predictions/NNLLprime_172-5/LHC13-CMS-Mtt-HT4-NNPDF31.dat",
        "NA",
        "NA"
    },
    //dytt
    {
    "predictions/NNLO_EW/LHC13-CMS-dytt-HT4-LUXQED17-pheno_rewrite.dat",
    "predictions/NNLO_EW_172-5/LHC13-CMS-mt_172.5-dytt-HT4-LUXQED17-pheno.dat",
    "predictions/NNLO_EW/LHC13-CMS-dytt-HT4-NNPDF31-pheno.dat",
    "NA",
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
    "\\GeV",
    "\\GeV",
    "",
    "",
    "\\GeV",
    "\\GeV",
    ""
};

std::vector<string> vars = {
    "$\\pt^{\\text{t}}$",
    "$\\pt^{\\text{\\tbar}}$",
    "$\\pt^{\\text{t}} \\mathrm{(leading)} $",
    "$\\pt^{\\text{t}} \\mathrm{(trailing)} $",
    "$\\pt^{\\text{t}} \\mathrm{(} \\ttbar\\ \\mathrm{RF)}$",
    "$y_{\\text{t}}$",
    "$y_{\\tbar}$",
    "$y_{\\text{t}} (\\mathrm{leading})$",
    "$y_{\\text{t}} (\\mathrm{trailing})$",
    "$\\pt^{\\ttbar}$",
    "$y_{\\ttbar}$",
    "$m_{\\ttbar}$",
    "$\\Delta|y| ( \\mathrm{t,\\tbar})$",
    "$\\Delta\\phi (\\mathrm{ t,\\tbar})$",
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
    "$\\pt^{\\mathrm{b}} \\mathrm{(leading)}$",
    "$\\pt^{\\mathrm{b}} \\mathrm{(trailing)}$",
    "$\\eta_{\\mathrm{b}} \\mathrm{(leading)}$",
    "$\\eta_{\\mathrm{b}} \\mathrm{(trailing)}$",
    "$\\pt^{\\bbbar}$",
    "$m_{\\bbbar}$",
    "$\\Nj$"
};

std::vector<string> vars_root = {
    "p_{T}^{t}",
    "p_{T}^{#bar{t}}",
    "p_{T}^{t} (leading)",
    "p_{T}^{t} (trailing)",
    "p_{T}^{t} (t#bar{t} RF)",
    "y_{t}",
    "y_{#bar{t}}",
    "y_{t} (leading)",
    "y_{t} (trailing)",
    "p_{T}^{t#bar{t}}",
    "y_{t#bar{t}}",
    "m_{t#bar{t}}",
    "#Delta|y|(t,#bar{t})",
    "#Delta#phi(t,#bar{t})",
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
    "p_{T}^{b} (leading)",
    "p_{T}^{b} (trailing)",
    "#eta^{b} (leading)",
    "#eta^{b} (trailing)",
    "p_{T}^{b#bar{b}}",
    "m_{b#bar{b}}",
    "N_{jets}"
};

std::vector<string> vars_txt = {
    "pt_t",
    "pt_tbar",
    "pt_t_leading",
    "pt_t_trailing",
    "pt_t_ttrf",
    "y_t",
    "y_tbar",
    "y_t_leading",
    "y_t_trailing",
    "pt_tt",
    "y_tt",
    "m_tt",
    "dabsy_tt",
    "dphi_tt",
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
    "pt_b_leading",
    "pt_b_trailing",
    "eta_b_leading",
    "eta_b_trailing",
    "pt_bb",
    "m_bb",
    "Njets"
};

std::vector<string> vars_obskey = {
    "PTT",
    "PTTBAR",
    "PTTLEADING",
    "PTTTRAILING",
    "PTTTTRF",
    "YT",
    "YTBAR",
    "YTLEADING",
    "YTTRAILING",
    "PTTT",
    "YTT",
    "MTT",
    "DYTT",
    "DPHITT",
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
    "PTBLEADING",
    "PTBTRAILING",
    "ETABLEADING",
    "ETABTRAILING",
    "PTBB",
    "MBB"
    "NJETS"
};

std::vector<string> vars_bnlo = {
    "$\\pt^{\\text{t}}$",
    "$\\pt^{\\tbar}$",
    "$y_{\\text{t}}$",
    "$y_{\\text{\\tbar}}$",
    "$\\pt^{\\ttbar}$",
    "$y^{\\ttbar}$",
    "$m_{\\ttbar}$",
    "$\\Delta|y|$ (t,\\tbar)"
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

vector<vector<std::string>> models_files(std::string mode){
    
    vector<vector<std::string>> m_f;
    vector<std::string> modelnames;
    vector<std::string> filenames;
    
    if (mode == "norm_parton"){
        modelnames = {
            "\\Powheg+\\Pythia",
            "\\Powheg+\\Herwigpp",
            "\\MGaMCatNLO+\\Pythia"};
        filenames = {
            "files/Jan18/parton/normalised/DiffXS_HypToppT_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypAntiToppT_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypToppTLead_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypToppTNLead_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypToppTTTRestFrame_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypTopRapidity_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypAntiTopRapidity_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypTopRapidityLead_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypTopRapidityNLead_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypTTBarpT_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypTTBarRapidity_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypTTBarMass_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypTTBarDeltaRapidity_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypTTBarDeltaPhi_source.root"
        };
    }else if (mode == "abs_parton"){
        modelnames = {
            "\\Powheg+\\Pythia",
            "\\Powheg+\\Herwigpp",
            "\\MGaMCatNLO+\\Pythia"};
        filenames = {
            "files/Jan18/parton/absolute/DiffXS_HypToppT_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypAntiToppT_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypToppTLead_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypToppTNLead_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypToppTTTRestFrame_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypTopRapidity_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypAntiTopRapidity_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypTopRapidityLead_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypTopRapidityNLead_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypTTBarpT_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypTTBarRapidity_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypTTBarMass_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypTTBarDeltaRapidity_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypTTBarDeltaPhi_source.root"
        };
    }
    else if (mode == "norm_parton_bnlo"){
        modelnames = {
            "NNLO+\\alphaew^{3} \\\\ (LUXQED17) \\\\ \\mt~=~173.3~\\GeV",
            "NNLO+\\alphaew^{3} \\\\ (LUXQED17) \\\\ \\mt~=~172.5~\\GeV",
            "NNLO+\\alphaew^{3} \\\\ (NNPDF3.1) \\\\ \\mt~=~173.3~\\GeV",
            "NNLO+NNLL' \\\\ (NNPDF3.1) \\\\ \\mt~=~173.3~\\GeV",
            "NNLO+NNLL' \\\\ (NNPDF3.1) \\\\ \\mt~=~172.5~\\GeV",
            "aN^{3}LO \\\\ (NNPDF3.0) \\\\ \\mt~=~172.5~\\GeV",
            "aNNLO \\\\ (CT14NNLO) \\\\ \\mt~=~172.5~\\GeV"
        };
        filenames = {
            "files/Jan18/parton/normalised/DiffXS_HypToppT_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypAntiToppT_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypTopRapidity_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypAntiTopRapidity_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypTTBarpT_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypTTBarRapidity_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypTTBarMass_source.root",
            "files/Jan18/parton/normalised/DiffXS_HypTTBarDeltaRapidity_source.root"
        };
    }else if (mode == "abs_parton_bnlo"){
        modelnames = {
            "NNLO+\\alphaew^{3} \\\\ (LUXQED17) \\\\ \\mt~=~173.3~\\GeV",
            "NNLO+\\alphaew^{3} \\\\ (LUXQED17) \\\\ \\mt~=~172.5~\\GeV",
            "NNLO+\\alphaew^{3} \\\\ (NNPDF3.1) \\\\ \\mt~=~173.3~\\GeV",
            "NNLO+NNLL' \\\\ (NNPDF3.1) \\\\ \\mt~=~173.3~\\GeV",
            "NNLO+NNLL' \\\\ (NNPDF3.1) \\\\ \\mt~=~172.5~\\GeV",
            "aN^{3}LO \\\\ (NNPDF3.0) \\\\ \\mt~=~172.5~\\GeV",
            "aNNLO \\\\ (CT14NNLO) \\\\ \\mt~=~172.5~\\GeV"
        };
        filenames = {
            "files/Jan18/parton/absolute/DiffXS_HypToppT_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypAntiToppT_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypTopRapidity_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypAntiTopRapidity_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypTTBarpT_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypTTBarRapidity_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypTTBarMass_source.root",
            "files/Jan18/parton/absolute/DiffXS_HypTTBarDeltaRapidity_source.root"
        };
    }else if (mode == "norm_particle"){
        modelnames = {
            "\\Powheg+\\Pythia",
            "\\Powheg+\\Herwigpp",
            "\\MGaMCatNLO+\\Pythia"};
        filenames = {
            "files/Jan18/particle/normalised/DiffXS_HypToppT_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypAntiToppT_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypToppTLead_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypToppTNLead_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypToppTTTRestFrame_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypTopRapidity_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypAntiTopRapidity_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypTopRapidityLead_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypTopRapidityNLead_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypTTBarpT_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypTTBarRapidity_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypTTBarMass_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypTTBarDeltaRapidity_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypTTBarDeltaPhi_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypLeptonpT_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypAntiLeptonpT_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypLeptonpTLead_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypLeptonpTNLead_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypLeptonEta_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypAntiLeptonEta_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypLeptonEtaLead_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypLeptonEtaNLead_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypLLBarpT_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypLLBarMass_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypLLBarDPhi_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypLLBarDEta_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypBJetpTLead_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypBJetpTNLead_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypBJetEtaLead_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypBJetEtaNLead_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypBBBarpT_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypBBBarMass_source.root",
            "files/Jan18/particle/normalised/DiffXS_HypJetMultpt30_source.root"
        };
    }
    else if (mode == "abs_particle"){
        modelnames = {
            "\\Powheg+\\Pythia",
            "\\Powheg+\\Herwigpp",
            "\\MGaMCatNLO+\\Pythia"};
        filenames = {
            "files/Jan18/particle/absolute/DiffXS_HypToppT_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypAntiToppT_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypToppTLead_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypToppTNLead_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypToppTTTRestFrame_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypTopRapidity_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypAntiTopRapidity_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypTopRapidityLead_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypTopRapidityNLead_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypTTBarpT_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypTTBarRapidity_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypTTBarMass_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypTTBarDeltaRapidity_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypTTBarDeltaPhi_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypLeptonpT_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypAntiLeptonpT_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypLeptonpTLead_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypLeptonpTNLead_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypLeptonEta_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypAntiLeptonEta_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypLeptonEtaLead_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypLeptonEtaNLead_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypLLBarpT_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypLLBarMass_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypLLBarDPhi_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypLLBarDEta_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypBJetpTLead_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypBJetpTNLead_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypBJetEtaLead_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypBJetEtaNLead_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypBBBarpT_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypBBBarMass_source.root",
            "files/Jan18/particle/absolute/DiffXS_HypJetMultpt30_source.root"
        };
    }

    m_f.push_back(modelnames);
    m_f.push_back(filenames);

    return m_f;    
}


vector<std::string> make_covariance_matrices(std::string mode, vector<std::string> filenames){
    
    vector<string> filenames_cov;
    std::string text_filename;
    std::string root_filename;

    //make covariance matrices
    for (int var = 0; var< filenames.size(); var++){
        
        text_filename  = "files/Nov1/";
        root_filename  = "";
        char_separator<char> sep("_");
        tokenizer< char_separator<char> > tokens(filenames[var], sep);
        
        vector<std::string> tokens_vec;
        BOOST_FOREACH (const string& t, tokens) {
            tokens_vec.push_back(t);
        }
        if (mode == "abs_particle") {
            text_filename += "particle/absolute/covariance/";
            text_filename += tokens_vec[1];
            text_filename += "_totCovEnvXSMtrxFile.txt";
            root_filename = text_filename.substr(0, text_filename.length() - 3) + "root";
        }else if (mode == "norm_particle"){
            text_filename += "particle/normalised/covariance/";
            text_filename += tokens_vec[1];
            text_filename += "_totCovEnvXSMtrxFile.txt";
            root_filename = text_filename.substr(0, text_filename.length() - 3) + "root";
        }else if (mode == "abs_parton" || mode == "abs_parton_bnlo" ){
            text_filename += "parton/absolute/covariance/";
            text_filename += tokens_vec[1];
            text_filename += "_totCovEnvXSMtrxFile.txt";
            root_filename = text_filename.substr(0, text_filename.length() - 3) + "root";
        }else if (mode == "norm_parton" || mode == "norm_parton_bnlo"){
            text_filename += "parton/normalised/covariance/";
            text_filename += tokens_vec[1];
            text_filename += "_totCovEnvXSMtrxFile.txt";
            root_filename = text_filename.substr(0, text_filename.length() - 3) + "root";
        }
        make_covariance_matrix(text_filename, filenames[var], vars_root[var]);
        filenames_cov.push_back(root_filename);
    }
    
    return filenames_cov;
}


///change so that this method runs over all measurements, and makes one set of files and submission
//file with the preferred variable ordering
void write_hepdata_tables(string mode, vector<string> filenames, vector<string> filenames_cov){

    ofstream HDfile_sub;
    ofstream HDfile;
    ofstream HDfile_cov;
    string filename;
    string filename_cov;
    string hepdata_filename = "";
    string hepdata_sub_filename = "";
    string hepdata_filename_cov = "";
    TCanvas * c_cov;
    
    char_separator<char> sep("#");
    char_separator<char> sep_2("_");
    char_separator<char> sep_3("/");
    
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
    vector<std::string> file_tokens_vec;


    tokenizer< char_separator<char> > mode_tokens(mode, sep_2);
    BOOST_FOREACH (const string& t, mode_tokens){
        mode_tokens_vec.push_back(t);
    }
    
    if (mode_tokens_vec[0] == "abs"){
        xsec_type_text = "Absolute";
        xsec_type_text_2 = "absolute";
    }else{
        xsec_type_text = "Normalized";
        xsec_type_text = "normalized";
    }
    
    hepdata_sub_filename = "HEPData/test_submission.yaml";
    HDfile_sub.open(hepdata_sub_filename);
    
    HDfile_sub<<"comment: Absolute and normalised differential cross sections of top quark pair production"<<endl;
    HDfile_sub<<"    at parton level in the full phase space and at particle level in a fiducial phase space."<<endl;
    HDfile_sub<<"    Numbers in yaml files provided by James Keaveney (CMS) and are taken from ROOT files prepared "<<endl;
    HDfile_sub<<"    by Mykola Savitskyi (CMS)."<<endl;
    HDfile_sub<<"date updated: 28/11/2018 XX:XX:XX"<<endl;
    HDfile_sub<<"hepdata_doi: XX.XXXXX/hepdata.85705.v1" << endl;
    HDfile_sub<<"record_ids:"<<endl;
    HDfile_sub<<"- {id: 1703993, type: inspire}"<<endl;
    HDfile_sub<<"- {id: XXXX, type: red}"<<endl;
    
    //////////////////////////////////////////////////////
    //write results   ////////////////////////////////////
    //////////////////////////////////////////////////////
    
    std::vector<std::string> parton_particle_files = {
        "DiffXS_HypToppT_source.root",
        "DiffXS_HypAntiToppT_source.root",
        "DiffXS_HypToppTLead_source.root",
        "DiffXS_HypToppTNLead_source.root",
        "DiffXS_HypToppTTTRestFrame_source.root",
        "DiffXS_HypTopRapidity_source.root",
        "DiffXS_HypAntiTopRapidity_source.root",
        "DiffXS_HypTopRapidityLead_source.root",
        "DiffXS_HypTopRapidityNLead_source.root",
        "DiffXS_HypTTBarpT_source.root",
        "DiffXS_HypTTBarRapidity_source.root",
        "DiffXS_HypTTBarMass_source.root",
        "DiffXS_HypTTBarDeltaRapidity_source.root",
        "DiffXS_HypTTBarDeltaPhi_source.root"
    };
    
    std::vector<std::string> particle_only_files = {
        "DiffXS_HypLeptonpT_source.root",
        "DiffXS_HypAntiLeptonpT_source.root",
        "DiffXS_HypLeptonpTLead_source.root",
        "DiffXS_HypLeptonpTNLead_source.root",
        "DiffXS_HypLeptonEta_source.root",
        "DiffXS_HypAntiLeptonEta_source.root",
        "DiffXS_HypLeptonEtaLead_source.root",
        "DiffXS_HypLeptonEtaNLead_source.root",
        "DiffXS_HypLLBarpT_source.root",
        "DiffXS_HypLLBarMass_source.root",
        "DiffXS_HypLLBarDPhi_source.root",
        "DiffXS_HypLLBarDEta_source.root",
        "DiffXS_HypBJetpTLead_source.root",
        "DiffXS_HypBJetpTNLead_source.root",
        "DiffXS_HypBJetEtaLead_source.root",
        "DiffXS_HypBJetEtaNLead_source.root",
        "DiffXS_HypBBBarpT_source.root",
        "DiffXS_HypBBBarMass_source.root",
        "DiffXS_HypJetMultpt30_source.root"
    };
    
    std::vector<std::string> modes = {"parton_particle", "particle_only"};
    std::vector<std::string> levels;
    std::vector<std::string> types = {"absolute", "normalised"};
    
    int fp = 0 ;
    int var_number =0;
    
    for (int mode = 0; mode < modes.size(); mode++){
            if (modes[mode]=="parton_particle"){
                filenames = parton_particle_files;
                levels = {"parton", "particle"};
            }else{
                filenames = particle_only_files;
                levels = {"particle"};
            }
        
            for (int f = 0; f < filenames.size(); f++){
                for (int l = 0; l < levels.size(); l++){
                    for (int ty = 0; ty < types.size(); ty++){
                                    filename =  "files/Jan18/" + levels[l] + "/" + types[ty] + "/" + filenames[f];
                                    std::cout << "smart filename  = " << filename << std::endl;
                        
                                    file_tokens_vec.clear();
                        
                                    tokenizer< char_separator<char> > file_tokens(filenames[f], sep_2);
                                    BOOST_FOREACH (const string& ft, file_tokens){
                                    file_tokens_vec.push_back(ft);
                                    }
                        
                                    //filename = filenames[f];
                                    //filename_cov = filenames_cov[f];
                                    filename_cov = "files/Nov1/" +  levels[l] + "/" + types[ty] + "/covariance/" + file_tokens_vec[1] + "_totCovEnvXSMtrxFile.root";
                        
                                    std::cout  << "filename = " << filename << " filename_cov = " << filename_cov << std::endl;
                        
                                    hepdata_filename = "";
                                    hepdata_filename += "HEPData/";
                                    varname = vars_txt[f];
                        
                                    fp = fp + 1;
                                    std::string yaml_filename = "d0" + std::to_string(fp) + "-x01-y01.yaml";
                                    hepdata_filename += yaml_filename;
                        
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
                        
                                    string hd_unit;
                                    std::string tex_str = h_mc->GetYaxis()->GetTitle();
                                    std::string x_tex_str = h_mc->GetXaxis()->GetTitle();
                        
                                    char_separator<char> sep_2(" ");
                                    vector<std::string> tokens_vec;
                                    tokenizer< char_separator<char> > tokens_3(tex_str, sep_2);
                                    BOOST_FOREACH (const string& t, tokens_3) {
                                        tokens_vec.push_back(t);
                                    }
                                    tex_str = tokens_vec[0] + tokens_vec[1] ;
                                    tokens_vec.clear();
                        
                                    tokenizer< char_separator<char> > tokens_4(x_tex_str, sep_2);
                                    BOOST_FOREACH (const string& t, tokens_4) {
                                        tokens_vec.push_back(t);
                                    }
                                    x_tex_str = tokens_vec[0];
                                    tokens_vec.clear();
                        
                                    replace(tex_str, "#", "\\");
                                    replace(tex_str, "\\left", "");
                                    replace(tex_str, "\\right", "");
                                    replace(x_tex_str, "#", "\\");

                                    std::string super_script = "^{-1}";
                                    std::string xsec_unit = "pb";

                        /*
                                    if (mode == "norm_parton" || mode == "norm_particle"){
                                        xsec_unit = "";
                                        if (units[f] == "\\GeV"){
                                            hd_unit = "1/GEV";
                                        }else{
                                            hd_unit =  "";
                                        }
                                    }else{
                                        xsec_unit = "";
                                        if (units[f] == "\\GeV"){
                                            hd_unit = "PB/GEV";
                                        }else{
                                            hd_unit =  "PB";
                                        }
                                    }
                         */
                        
                                    std::string level;
                                    std::string type;
                                    std::string sig;
/*
                                    if (mode == "norm_particle"){
                                        level = "particle";
                                        type = "Normalised";
                                        sig = "SIG(fiducial)";
                                    }else if (mode == "abs_particle"){
                                        level = "particle";
                                        type = "Absolute";
                                        sig = "SIG(fiducial)";
                                    }
                                    else if (mode == "norm_parton"){
                                        level = "parton";
                                        type = "Absolute";
                                        sig = "SIG";

                                    }else if (mode == "abs_parton"){
                                        level = "parton";
                                        type = "Absolute";
                                        sig = "SIG";
                                    }
 
 */
                                    HDfile_sub<<"---"<<endl;
                                    HDfile_sub<<"data_file: "<< yaml_filename <<endl;
                                    HDfile_sub<<"description: |"<< endl;
                                    HDfile_sub<<"  "<< types[ty] <<" differential cross section at "<< levels[l] << " level as a function of "<< vars[var_number] <<"."<<endl;
                                    HDfile_sub<<"keywords: "<<endl;
                                    HDfile_sub<<"- name: reactions"<<endl;
                                    HDfile_sub<<"  values: [P P --> TOP TOPBAR X]"<<endl;
                                    HDfile_sub<<"- name: observables"<<endl;
                        
                                   // if (mode == "norm_particle" || mode == "abs_particle" ){
                                   //     HDfile_sub<<"  values: [DSIG_FID/" << vars_obskey[f] << "]"  << endl;
                                   // }
                        
                                    HDfile_sub<<"  name: phrases"<<endl;
                                    HDfile_sub<<"  values: [Top, Differential Cross Section, Proton-Proton Scattering, Top Production]"<<endl;
                                    HDfile_sub<<"- name: cmenergies"<<endl;
                                    HDfile_sub<<"  values: [13000.0]"<<endl;
                                    HDfile_sub<<"location: Data from Fig. "<< f+3 <<endl;
                        
                                    std:string unit = units[f];
                                    replace(unit, "\\", "");
                        
                                    if (unit == "GeV"){
                                        unit = "GEV";
                                    }else{
                                        unit = "";
                                    }
                        
                                    //independent variables
                                    HDfile<<"independent_variables: "<<endl;

                                    if (unit != ""){
                                        HDfile<<"- header: {name: '$"<< x_tex_str <<"$ [" <<  unit << "]'}"<<endl;
                                    }else{
                                        HDfile<<"- header: {name: '$"<< x_tex_str <<"$ '}"<<endl;
                                    }
                        
                                    HDfile<<"  values:"<<endl;
                                    for (int bin = 1; bin < bins.size(); bin++){
                                        HDfile<<"  - {high: "<< bins[bin] <<", low: "<< bins[bin-1] <<"}"<<endl;
                                    }
                        
                                    //dependent variables
                                    HDfile<<"dependent_variables: "<<endl;
                                    HDfile<<"- header: {name: "<<  sig  <<", units: " << hd_unit <<"}" <<endl;
                                    HDfile<<"  qualifiers: "<<endl;
                                    HDfile<<"  - {name: RE, value: P P --> TOP TOPBAR X}"<<endl;
                                    HDfile<<"  - {name: SQRT(S), units: GEV, value: '13000.0'}"<<endl;
                                    HDfile<<"  values:"<<endl;
                        
                                    for (int bin = 1; bin < bins.size(); bin++){
                                        g_data->GetPoint(bin-1, data_x, data_y);
                                        
                                        stat_unc = g_data_stat_unc_only->GetErrorY(bin-1);
                                        tot_unc = g_data->GetErrorY(bin-1);
                                        sys_unc = pow( (  (pow(tot_unc,2)) -  (pow(stat_unc,2))    ), 0.5     );
                                        HDfile<<"  - value: "  << data_y  << endl;
                                        HDfile<<"    errors:"<<endl;
                                        HDfile<<"    - {label: stat, symerror: " << stat_unc << "}" << endl;
                                        HDfile<<"    - {label: sys, symerror: "  << sys_unc  << "}" << endl;
                                    }
                        
                                    //get covariance matrix and draw
                                    gStyle->SetPaintTextFormat("4.8f");
                                    //gStyle->SetPalette(51);
                                    const Int_t Number = 3;
                                    Double_t Red[Number]    = { 1.00, 0.00, 0.00};
                                    Double_t Green[Number]  = { 0.00, 1.00, 0.00};
                                    Double_t Blue[Number]   = { 1.00, 0.00, 1.00};
                                    Double_t Length[Number] = { 0.00, 0.50, 1.00 };
                                    Int_t nb = 500;
                                    //TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);

                                    f_var_cov = new TFile(filename_cov.c_str());
                                    h_cov = (TH2D*)f_var_cov->Get("cov");
                                    h_cov->SetContour(nb);
                                    c_cov = new TCanvas();
                                    h_cov->Draw("colz");
                        
                                    //gPad->Update();
                                    TPaletteAxis *palette = (TPaletteAxis*)h_cov->GetListOfFunctions()->FindObject("palette");
                        
                                    // the following lines moe the paletter. Choose the values you need for the position.
                                    palette->SetX1NDC(0.89);
                                    palette->SetX2NDC(0.94);
                                    palette->SetY1NDC(0.1);
                                    palette->SetY2NDC(0.9);
                                    //gPad->Modified();
                                    //gPad->Update();
                        
                                    std::string hepdata_filename_covplot = "";
                                    hepdata_filename_covplot += "CovMatrices/";
                                    varname = vars_txt[f];
                                    hepdata_filename_covplot += varname;
                                    hepdata_filename_covplot += "_";
                                    hepdata_filename_covplot += mode;
                                    hepdata_filename_covplot += "_";
                                    hepdata_filename_covplot += "_cov.pdf";
                        
                                    c_cov->SetTopMargin(0.1);
                                    c_cov->SetBottomMargin(0.15);
                                    c_cov->SetLeftMargin(0.0935);
                                    c_cov->SetRightMargin(0.12);
                        
                                    float H = c_cov->GetWh();
                                    float W = c_cov->GetWw();
                                    float l = c_cov->GetLeftMargin();
                                    float t = c_cov->GetTopMargin();
                                    float r = c_cov->GetRightMargin();
                                    float b = c_cov->GetBottomMargin();
                                    float extraOverCmsTextSize  = 0.8;
                        
                                    TString cmsText, extraText,extraText2, lumiText;
                                    cmsText += "CMS";
                                    extraText += "Supplementary ";
                        
                        /*
                                    if (mode == "norm_parton"){
                                        extraText2 += " parton level, normalised";
                                    }else if (mode == "abs_parton"){
                                        extraText2 += " parton level, absolute";
                                    }else if (mode == "norm_particle"){
                                        extraText2 += " particle level, normalised";
                                    }else if (mode == "abs_particle"){
                                        extraText2 += " particle level, absolute";
                                    }
                         */
                        
                                    lumiText += "35.9 fb^{-1} (13 TeV)";
                        
                                    TLatex latex;
                                    latex.SetNDC();
                                    latex.SetTextAngle(0);
                                    latex.SetTextSize(0.48*t);
                                    latex.SetTextColor(kBlack);
                                    latex.SetTextFont(61);
                                    latex.SetTextAlign(31);
                                    latex.DrawLatex(0.16,0.916,cmsText);
                                    latex.SetTextFont(52);
                                    latex.SetTextSize(0.49*t*extraOverCmsTextSize);
                                    latex.DrawLatex(0.345,0.915,extraText);
                                    latex.SetTextFont(42);
                                    latex.SetTextSize(0.38*t);
                                    latex.DrawLatex(0.628,0.915,extraText2);
                                    latex.SetTextFont(42);
                                    latex.SetTextSize(0.46*t);
                                    latex.DrawLatex(0.88,0.916,lumiText);
                        
                                    c_cov->SaveAs(hepdata_filename_covplot.c_str());
                        
                                    //////////////////////////////////////////////////////
                                    //write covariance matrices //////////////////////////
                                    //////////////////////////////////////////////////////
                                    std::string yaml_filename_cov = "HEPData/d0" + std::to_string(fp) + "-x01-y01_cov.yaml";

                                    HDfile_cov.open (yaml_filename_cov);
                                    HDfile_cov<<"dependent_variables:" <<endl;
                                    HDfile_cov<<"- header: {name: ''}" <<endl;
                                    HDfile_cov<<"  qualifiers:" <<endl;
                                    HDfile_cov<<"  - {name: '', value: COVARIANCE MATRIX}" <<endl;
                                    HDfile_cov<<"  values:" <<endl;
                        
                                    for (int xbin = 1; xbin <= h_cov->GetNbinsX(); xbin++)
                                        for (int xbin = 1; xbin <= h_cov->GetNbinsX(); xbin++){
                                            for (int ybin = 1; ybin <= h_cov->GetNbinsY(); ybin++){
                                                HDfile_cov<<"  - {value: "  << std::defaultfloat << std::setprecision(15) << h_cov->GetBinContent(xbin,ybin) << "}"<< endl;
                                        }
                                    }
                        
                                    HDfile_cov<<"independent_variables:" <<endl;
                                    HDfile_cov<<"- header: {name: Bin}" <<endl;
                                    HDfile_cov<<"  values:" <<endl;
                                    for (int xbin = 1; xbin <= h_cov->GetNbinsX(); xbin++)
                                    for (int xbin = 1; xbin <= h_cov->GetNbinsX(); xbin++){
                                        for (int ybin = 1; ybin <= h_cov->GetNbinsY(); ybin++){
                                            HDfile_cov<<"  - {value: "<< ybin <<".0}" <<endl;
                                        }
                                    }
                                    HDfile_cov<<"- header: {name: Bin}" <<endl;
                                    HDfile_cov<<"  values:" <<endl;
                                    for (int xbin = 1; xbin <= h_cov->GetNbinsX(); xbin++)
                                    for (int xbin = 1; xbin <= h_cov->GetNbinsX(); xbin++){
                                        for (int ybin = 1; ybin <= h_cov->GetNbinsY(); ybin++){
                                            HDfile_cov<<"  - {value: "<< xbin <<".0}" <<endl;
                                        }
                                    }
                        
                                    HDfile.close();
                                    HDfile_cov.close();
                                    f_var->Close();
                }
            }
                var_number++;
            }
    
    HDfile_sub.close();
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
  //  cout<<"tex filename = " << tex_filename <<endl;
    
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
        myfile_2 << "All the measured differential cross sections at the parton level are tabulated in Tables \\ref{tab:norm_parton0}--\\ref{tab:norm_parton13}. The statistical and systematic uncertainties are quoted separately for each bin."<<endl;
    }else{
        myfile_2 << "All the measured differential cross sections at the particle level are tabulated in Tables \\ref{tab:norm_particle0}--\\ref{tab:norm_particle32}. The statistical and systematic uncertainties are quoted separately for each bin."<<endl;
    }
    
    myfile_2 << "\\clearpage"<<endl;

    
    for (int f = 0; f < filenames.size(); f++){
        bins.clear();
        tokens_vec.clear();
        string f_str = std::to_string(f);
        
       // std::cout <<"filename =  "<< f_str  <<  std::endl;
        
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
            myfile_2 << "\\caption{The measured differential cross section and bin boundaries for each bin of the normalized and absolute measurements of the \\ttbar\\ differential cross section at parton level in the full phase space as a function of " << varname <<" are tabulated.}"<<endl;
        }else{
                myfile_2 << "\\caption{The measured differential cross section and bin boundaries for each bin of the normalized and absolute measurements of the \\ttbar\\ differential cross section at particle level in the fiducial phase space as a function of " << varname <<" are tabulated.}"<<endl;
        }
            
        myfile_2 << "\\begin{tabular}{ c | c | c }"<<endl;
        
        varname = varname.substr(1,varname.size() - 2);

        
       // cout <<" xtitle = "<< xtitle << endl;
       // cout <<" varname  = "<< varname << endl;
       // cout <<" "<< endl;
        
        if (unit == ""){
            myfile_2<< " $" << varname << "$  & $\\frac{1}{\\sigma}$ $\\frac{\\mathrm{d}\\sigma}{\\mathrm{d}"<< varname <<"}$ & $\\frac{\\mathrm{d}\\sigma}{\\mathrm{d}"<< varname <<"}$ [pb] \\\\ "<<endl;
        }else{
            myfile_2  <<"$ "<< varname << "$ ["<< unit << "]  & $\\frac{1}{\\sigma}$ $\\frac{\\mathrm{d}\\sigma}{\\mathrm{d}"<< varname <<"}$ ["<< unit <<"$^{-1}$] & $\\frac{\\mathrm{d}\\sigma}{\\mathrm{d}"<< varname <<"}$ [pb/"<< unit <<"] \\\\ "<<endl;
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
            scify_norm = true;
            scify_abs = true;
            
            g_data_norm->GetPoint(bin-1, data_x, data_norm_y);
            g_data_abs->GetPoint(bin-1, data_x, data_abs_y);

            stat_unc_norm = g_data_norm_stat_unc_only->GetErrorY(bin-1);
            tot_unc_norm = g_data_norm->GetErrorY(bin-1);
            stat_unc_abs = g_data_abs_stat_unc_only->GetErrorY(bin-1);
            tot_unc_abs = g_data_abs->GetErrorY(bin-1);

            sys_unc_norm = pow( (  (pow(tot_unc_norm,2)) -  (pow(stat_unc_norm,2))    ), 0.5 );
            sys_unc_abs = pow( (  (pow(tot_unc_abs,2)) -   (pow(stat_unc_abs,2))   ), 0.5    );
            
            if (filename_abs == "files/Jan18/particle/absolute/DiffXS_HypToppTTTRestFrame_source.root"){
                //std::cout << "bin = = =  "<< data_abs_y << "  scify abs = "<< scify_abs  <<std::endl;
            }
            
            if (data_norm_y < 10000.0 && data_norm_y > 1000.0){
                exp_str_norm = " $ \\times 10^{3}$ ";
                data_norm_y = (data_norm_y*(0.001));
                stat_unc_norm = (stat_unc_norm*(0.001));
                sys_unc_norm = (sys_unc_norm*(0.001));
            }
           else if (data_norm_y < 1000.0 && data_norm_y > 100.0){
                exp_str_norm = " $ \\times 10^{2}$ ";
                data_norm_y = (data_norm_y*(0.01));
                stat_unc_norm = (stat_unc_norm*(0.01));
                sys_unc_norm = (sys_unc_norm*(0.01));
            }  else if (data_norm_y < 100.0 && data_norm_y > 10.0){
                exp_str_norm = " $ \\times 10$ ";
                data_norm_y = (data_norm_y*(0.1));
                stat_unc_norm = (stat_unc_norm*(0.1));
                sys_unc_norm = (sys_unc_norm*(0.1));
            }
            
            else if (data_norm_y < 0.1 && data_norm_y > 0.01){
                exp_str_norm = " $ \\times 10^{-2}$ ";
                data_norm_y = (data_norm_y*100.0);
                stat_unc_norm = (stat_unc_norm*100.0);
                sys_unc_norm = (sys_unc_norm*100.0);
            }
            else if (data_norm_y < 0.01 && data_norm_y > 0.001){
                exp_str_norm = " $ \\times 10^{-3}$ ";
                data_norm_y = (data_norm_y*1000.0);
                stat_unc_norm = (stat_unc_norm*1000.0);
                sys_unc_norm = (sys_unc_norm*1000.0);
            } else if (data_norm_y < 0.001 && data_norm_y > 0.0001){
                exp_str_norm = " $ \\times 10^{-4}$ ";
                data_norm_y = (data_norm_y*10000.0);
                stat_unc_norm = (stat_unc_norm*10000.0);
                sys_unc_norm = (sys_unc_norm*10000.0);
            }else if (data_norm_y < 0.0001 && data_norm_y > 0.00001){
                exp_str_norm = " $ \\times 10^{-5}$ ";
                data_norm_y = (data_norm_y*100000.0);
                stat_unc_norm = (stat_unc_norm*100000.0);
                sys_unc_norm = (sys_unc_norm*100000.0);
            }else if (data_norm_y < 0.00001 && data_norm_y > 0.000001){
                exp_str_norm = " $ \\times 10^{-6}$ ";
                data_norm_y = (data_norm_y*1000000.0);
                stat_unc_norm = (stat_unc_norm*1000000.0);
                sys_unc_norm = (sys_unc_norm*1000000.0);
            } else {
                scify_norm = false;
            }
            
            
            if (data_abs_y < 10000.0 && data_abs_y > 1000.0){
                exp_str_abs = " $ \\times 10^{3}$ ";
                data_abs_y = (data_abs_y*(0.001));
                stat_unc_abs = (stat_unc_abs*(0.001));
                sys_unc_abs = (sys_unc_abs*(0.001));
            }
            else if (data_abs_y < 1000.0 && data_abs_y > 100.0){
                exp_str_abs = " $ \\times 10^{2}$ ";
                data_abs_y = (data_abs_y*(0.01));
                stat_unc_abs = (stat_unc_abs*(0.01));
                sys_unc_abs = (sys_unc_abs*(0.01));
            }
            else if (data_abs_y < 100.0 && data_abs_y > 10.0){
                exp_str_abs = " $ \\times 10$ ";
                data_abs_y = (data_abs_y*(0.1));
                stat_unc_abs = (stat_unc_abs*(0.1));
                sys_unc_abs = (sys_unc_abs*(0.1));
            }
        
            else if (data_abs_y < 0.1 && data_abs_y > 0.01){
                exp_str_abs = " $ \\times 10^{-2}$ ";
                data_abs_y = (data_abs_y*100.0);
                stat_unc_abs = (stat_unc_abs*100.0);
                sys_unc_abs = (sys_unc_abs*100.0);
            }
             else if (data_abs_y < 0.1 && data_abs_y > 0.01){
                exp_str_abs = " $ \\times 10^{-2}$ ";
                data_abs_y = (data_abs_y*100.0);
                stat_unc_abs = (stat_unc_abs*100.0);
                sys_unc_abs = (sys_unc_abs*100.0);
            }
            else if (data_abs_y < 0.01 && data_abs_y > 0.001){
                exp_str_abs = " $ \\times 10^{-3}$ ";
                data_abs_y = (data_abs_y*1000.0);
                stat_unc_abs = (stat_unc_abs*1000.0);
                sys_unc_abs = (sys_unc_abs*1000.0);
            } else if (data_abs_y < 0.001 && data_abs_y > 0.0001){
                if (filename_abs == "files/Jan18/particle/absolute/DiffXS_HypToppTTTRestFrame_source.root"){
                   //// std::cout << "scifying this bin  ata_abs_y < 0.001 && data_abs_y > 0.0001 "  <<std::endl;
                }
                exp_str_abs = " $ \\times 10^{-4}$ ";
                data_abs_y = (data_abs_y*10000.0);
                stat_unc_abs = (stat_unc_abs*10000.0);
                sys_unc_abs = (sys_unc_abs*10000.0);
            }else if (data_abs_y < 0.0001 && data_abs_y > 0.00001){
                exp_str_abs = " $ \\times 10^{-5}$ ";
                data_abs_y = (data_abs_y*100000.0);
                stat_unc_abs = (stat_unc_abs*100000.0);
                sys_unc_abs = (sys_unc_abs*100000.0);
            }else if (data_abs_y < 0.00001 && data_abs_y > 0.000001){
                exp_str_abs = " $ \\times 10^{-6}$ ";
                data_abs_y = (data_abs_y*1000000.0);
                stat_unc_abs = (stat_unc_abs*1000000.0);
                sys_unc_abs = (sys_unc_abs*1000000.0);
            } else {
                scify_abs = false;
                if (filename_abs == "files/Jan18/particle/absolute/DiffXS_HypToppTTTRestFrame_source.root"){
                    std::cout << "scifying false "  <<std::endl;
                }
            }
            
            int n_sig = 4;
            
            std::string str_data_norm_y;
            std::string str_stat_unc_norm;
            std::string str_sys_unc_norm;
            std::string str_data_abs_y;
            std::string str_stat_unc_abs;
            std::string str_sys_unc_abs;

            //round to n sigificant figures
            str_data_norm_y = round_off(data_norm_y, n_sig);
            str_stat_unc_norm = round_off(stat_unc_norm, n_sig);
            str_sys_unc_norm = round_off(sys_unc_norm, n_sig);
        
            str_data_abs_y = round_off(data_abs_y, n_sig);
            str_stat_unc_abs = round_off(stat_unc_abs, n_sig);
            str_sys_unc_abs = round_off(sys_unc_abs, n_sig);
            
           // std::cout <<" returned strings =  "<< str_data_abs_y << " " << str_stat_unc_abs << " "<<   str_sys_unc_abs<< std::endl;

            
            //myfile_2<< std::defaultfloat << std::setprecision(4) << bins[bin-1]<<"  --  "<< bins[bin];
            myfile_2<< std::defaultfloat << std::setprecision(4) <<"$\\text{[}$"<< bins[bin-1]<<", "<< bins[bin]<<"$\\text{]}$";

            
            if (scify_norm && scify_abs){
           // myfile_2 << std::fixed << std::setprecision(4) << " & ("<< str_data_norm_y  <<" $\\pm$ "<<  str_stat_unc_norm <<" $\\pm$ "<< str_sys_unc_norm << ")"<<exp_str_norm << " & ("<< str_data_abs_y <<" $\\pm$ "<<  str_stat_unc_abs << " $\\pm$ " << str_sys_unc_abs <<") "<< exp_str_abs << " \\\\"<<endl;
            myfile_2 << " & ("<< str_data_norm_y  <<" $\\pm$ "<<  str_stat_unc_norm <<" $\\pm$ "<< str_sys_unc_norm << ")"<<exp_str_norm << " & ("<< str_data_abs_y <<" $\\pm$ "<<  str_stat_unc_abs << " $\\pm$ " << str_sys_unc_abs <<") "<< exp_str_abs << " \\\\"<<endl;
                
            }else if (scify_norm && !scify_abs){
           // myfile_2 << std::fixed << std::setprecision(4) << " & ("<< str_data_norm_y  <<" $\\pm$ "<<  str_stat_unc_norm <<" $\\pm$ "<< str_sys_unc_norm << ")"<<exp_str_norm << " & "<< str_data_abs_y <<" $\\pm$ "<<  str_stat_unc_abs << " $\\pm$ " << str_sys_unc_abs  << " \\\\"<<endl;
            myfile_2 << " & ("<< str_data_norm_y  <<" $\\pm$ "<<  str_stat_unc_norm <<" $\\pm$ "<< str_sys_unc_norm << ")"<<exp_str_norm << " & "<< str_data_abs_y <<" $\\pm$ "<<  str_stat_unc_abs << " $\\pm$ " << str_sys_unc_abs  << " \\\\"<<endl;
                
            }else if (!scify_norm && scify_abs){
            //myfile_2 << std::fixed << std::setprecision(4) << " & "<< str_data_norm_y  <<" $\\pm$ "<<  str_stat_unc_norm <<" $\\pm$ "<< str_sys_unc_norm << " & ("<< data_abs_y <<" $\\pm$ "<<  str_stat_unc_abs << " $\\pm$ " << str_sys_unc_abs <<") "<< exp_str_abs << " \\\\"<<endl;
            myfile_2 << " & "<< str_data_norm_y  <<" $\\pm$ "<<  str_stat_unc_norm <<" $\\pm$ "<< str_sys_unc_norm << " & ("<< str_data_abs_y <<" $\\pm$ "<<  str_stat_unc_abs << " $\\pm$ " << str_sys_unc_abs <<") "<< exp_str_abs << " \\\\"<<endl;
            }
            else if (!scify_norm && !scify_abs){
            //myfile_2 << std::fixed<< std::setprecision(4) << " & "<< str_data_norm_y  <<" $\\pm$ "<<  str_stat_unc_norm <<" $\\pm$ "<< str_sys_unc_norm << " & " << str_data_abs_y  <<" $\\pm$ "<<  str_stat_unc_abs << " $\\pm$ " << str_sys_unc_abs <<" \\\\" <<endl;
            myfile_2 << " & "<< str_data_norm_y  <<" $\\pm$ "<<  str_stat_unc_norm <<" $\\pm$ "<< str_sys_unc_norm << " & " << str_data_abs_y  <<" $\\pm$ "<<  str_stat_unc_abs << " $\\pm$ " << str_sys_unc_abs <<" \\\\" <<endl;
            }
        }
        
        string lab_str = mode + f_str;
        myfile_2 << "\\end{tabular}"<<endl;
        myfile_2 << "\\label{tab:"<< lab_str <<"}"<<endl;
        myfile_2 << "\\end{table}"<<endl;
        myfile_2 << "  "<<endl;
        
        if ((f % 3 == 0) && (f != 0)){
        myfile_2 << "\\clearpage"<<endl;
    }
    
        myfile_2 << "  "<<endl;

    /* string filename = mode + ".tex";
    
     myfile_2.precision(2);
     myfile_2.open (filename);
     myfile_2 << "\\begin{table}"<<endl;
     myfile_2 << "\\small"<<endl;
     myfile_2 << "\\centering"<<endl;
     myfile_2 << "\\begin{tabular}{| l | c | c | c | c | c | c | c | c | c | c | c |}"<<endl;
     myfile_2 << "\\hline"<<endl;
     */

        f_var_norm->Close();
        f_var_abs->Close();

}
    myfile_2.close();

    
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
    std::string root_title = mode + ".root";

    TFile *f_canvas = new TFile(root_title.c_str(),"RECREATE");
    c->Write();
    c->SaveAs(pdf_title.c_str());
    f_canvas->Close();
    f->Close();
    
    
}


// Function to round - off the number
string round_off(double N, double n)
{
//first do the rounding
    
   // std::cout <<"before rounding = "<< N << std::endl;
    
    //n represents the desired number of significant digits
    if (N < 0.1){
        n = n - 1;
    } else if ((N < 1.0) && (N > 0.1)){
        n = n - 1;
    }else if ((N >= 1.0) && (N < 10.0)){
        n = n - 1;
    }else if((N >= 10.0) && (N < 100.0)){
        n = n - 2;
    }else if((N >= 100.0) && (N < 1000.0)){
        n = n - 3;
    }
    else if((N >= 1000.0) && (N < 10000.0)){
        n = n - 4;
    }else{
        n = n - 4;
    }
    
    N = (N)*(pow(10, n));
    N = roundf(N);
    N = N/(pow(10, n));
    
   // std::cout << " after rounding  = "<< N << " n = " <<  n <<std::endl;

    
    stringstream stream;
    stream << fixed << setprecision(n) << N;
    string str_j = stream.str();
    
    //std::cout << " after string = "<< str_j <<std::endl;

//    std::cout << " " <<  std::endl;

    
    
    //remove trailing zeros
//    std::string str_j;
 //   str_j = std::to_string(j);
    str_j.erase ( str_j.find_last_not_of('0') + 1, std::string::npos );

    if (N < 0.001){
        str_j = "$<10^{-3}$";
    }
    
    
    return str_j;
}








