//make pvals
// script to write tables of p-values, chi2/ndof between// differential
// cross section results and theory/MC predictions.
#include "make_pvals.h"
using namespace std;

int main(int argc, const char * argv[]){
 
   // make_table("norm_parton_bnlo");
    // make_table("abs_parton_bnlo");
    make_table("norm_parton");
    make_table("norm_particle");
    make_table("abs_particle");
    make_table("abs_parton");

}

int make_table(std::string mode){
    
    //string mode = "norm_parton";
    
    vector<string> filenames;
    vector<string> filenames_cov;

    myfile.precision(2);
    string filename = mode + ".tex";
    myfile.open (filename);
    myfile << "\\begin{table}"<<endl;
    myfile << "\\small"<<endl;
    myfile << "\\centering"<<endl;
    myfile << "\\begin{tabular}{| l | c | c | c | c | c | c | c | c | c | c | c |}"<<endl;
    myfile << "\\hline"<<endl;
    
    vector<string> modelnames = {
        "\\Powheg+\\Pythia",
        "\\Powheg+\\Herwigpp",
        "\\MGaMCatNLO+\\Pythia"
    };
    
//  vector<string> modelnames = {
//        "NNLO+\\alphaew^{3} \\\\ (LUXQED17) \\\\ \\mt~=~173.3~\\GeV",
//         "NNLO+\\alphaew^{3} \\\\ (LUXQED17) \\\\ \\mt~=~172.5~\\GeV",
//          "NNLO+\\alphaew^{3} \\\\ (NNPDF3.1) \\\\ \\mt~=~173.3~\\GeV",
//          "NNLO+NNLL' \\\\ (NNPDF3.1) \\\\ \\mt~=~173.3~\\GeV",
//          "NNLO+NNLL' \\\\ (NNPDF3.1) \\\\ \\mt~=~172.5~\\GeV",
//          "aN^{3}LO \\\\ (NNPDF3.0) \\\\ \\mt~=~172.5~\\GeV",
//          "aNNLO \\\\ (CT14NNLO) \\\\ \\mt~=~172.5~\\GeV"
//     };
    

    if (mode == "norm_parton"){
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
    } else if (mode == "abs_parton"){
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
    } else if (mode == "abs_parton_bnlo"){
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
    }
    else if (mode == "norm_particle"){
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

    int nvars = vars.size();
    int nmodels = modelnames.size();

    float   chisq[n_models][n_vars];
    int     ndof[n_models][n_vars];
    double   pval[n_models][n_vars];
    
    std::string text_filename;
    std::string root_filename;
    
    //first write latex tables of differential results
    write_results_table(mode, filenames);
    
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
    
    std::tuple <float, int, float > gof;
    vector<std::string> preds;
    if (mode == "norm_parton" || mode == "abs_parton" ||mode == "norm_particle" || mode == "abs_particle"){
    for (int m = 0; m < modelnames.size(); m++){
        for (int f = 0; f < filenames.size(); f++){
            gof = process_result(modelnames[m], filenames[f], filenames_cov[f], f, m);
            chisq[m][f] = std::get<0>(gof);
            ndof[m][f] = std::get<1>(gof);
            pval[m][f] = std::get<2>(gof);
          //  cout <<"RETURNED chisq , ndof, pval  =  "<< chisq[m][f] <<"  "<< ndof[m][f] <<"  "<<  pval[m][f]  << endl;
        }
    }
    }else if (mode == "norm_parton_bnlo" || mode == "abs_parton_bnlo"){
        for (int m = 0; m < modelnames.size(); m++){
            for (int f = 0; f < filenames.size(); f++){
             //   cout << "bnlo matrix element = " << bnlo_matrix[f][m] << endl;
                if (bnlo_matrix[f][m] != "NA"){
                gof = process_result(modelnames[m], filenames[f], filenames_cov[f], f, m);
                chisq[m][f] = std::get<0>(gof);
                ndof[m][f] = std::get<1>(gof);
                pval[m][f] = std::get<2>(gof);
               // cout <<"RETURNED chisq , ndof, pval  =  "<< chisq[m][f] <<"  "<< ndof[m][f] <<"  "<<  pval[m][f]  << endl;
                }
                else{
                    chisq[m][f] = -1.0;
                    ndof[m][f] = -1.0;
                    pval[m][f] = -1.0;
                }
            }
            }
    }
    
    //remake plots...
   // for (int f = 0; f < filenames.size(); f++){
   //call to plot remaker function...
   //     plot_remaker(filenames[f], mode);
   // }

    write_latex(mode, modelnames, vars, chisq, ndof, pval);
    myfile << "\\end{table}"<<endl;
    myfile.close();
    
   //summary_plot("norm_parton.root", "norm_particle.root");
    summary_plot("norm_parton_bnlo.root", "norm_particle.root");
    
    //second write HEPData tables of differential results
    write_hepdata_tables(mode, filenames, filenames_cov);

    return 1;
}

std::tuple <float, int, float > process_result(string modelname, string filename, string filename_cov, int f, int m){
    
  //  cout <<"processing result 1 "<< "  file  = "<< filename<< " model = " << modelname <<"  filename_cov = "<<filename_cov<< endl;
    std::string mode;
    int ndof;
    TFile * f_results = new TFile(filename.c_str());
    TCanvas * c = (TCanvas*)f_results->Get("canvas;1");
    //default setting for h_model to get ndof (changed later)
    h_model = (TH1F*)c->GetPrimitive("MCATNLOplot");
    TGraphAsymmErrors * g_model;

    //trick to fined out if we are running on abs/norm/div-bw/nodiv-bw
    char_separator<char> sep("/");
    tokenizer< char_separator<char> > tokens(filename, sep);
    vector<std::string> tokens_vec;
    BOOST_FOREACH (const string& t, tokens) {
        tokens_vec.push_back(t);
    }
    if (tokens_vec[3] == "normalised"){
        mode = "norm";
        ndof = (h_model->GetNbinsX() - 1 ); //need to reduce this by 1 for normalised results
    }else if (tokens_vec[3] == "absolute"){
        ndof = (h_model->GetNbinsX() );
        mode = "abs";
    }
    
    string modelhistoname, dataname;
    double chisq_running = 0.0, pval = 0.0;
    double bin_x, bin_data, e_data;
    dataname = "Graph";

    
    if (modelname == "\\Powheg+\\Pythia"){
        modelhistoname = "Nominalplot_copy";
    } else if (modelname == "\\Powheg+\\Herwigpp"){
        modelhistoname = "MCATNLOplot";
    }
    else if (modelname == "\\MGaMCatNLO+\\Pythia"){
        modelhistoname = "POWHEGplot";
    }
    else{
        modelhistoname = modelname;
    }
    
    if((modelhistoname == "Nominalplot_copy") || (modelhistoname == "POWHEGplot") || (modelhistoname == "MCATNLOplot")){
    h_model = (TH1F*)c->GetPrimitive(modelhistoname.c_str());
    }else {
        g_model = read_prediction(bnlo_matrix[f][m], divide_options_top_pt[m], mode);
        for (int bin = 0; bin < g_model->GetN(); bin++){
            g_model->GetPoint(bin, bin_x, bin_data);
            //cout <<" bin model  =  "<<  bin_data  <<endl;
            e_data = g_model->GetErrorY(bin);
            h_model->SetBinContent(bin+1,bin_data);
            h_model->SetBinError(bin+1,e_data);
        }
    }

    TGraphAsymmErrors * g_data = (TGraphAsymmErrors*)f_results->Get("data");
    TH1F*  h_data = (TH1F*)h_model->Clone();
    double tot_data = 0.0;
    for (int bin = 0; bin < h_model->GetNbinsX(); bin++){
        g_data->GetPoint(bin, bin_x, bin_data);
        e_data = g_data->GetErrorY(bin);
        h_data->SetBinContent(bin+1,bin_data);
        h_data->SetBinError(bin+1,e_data);
        tot_data = tot_data + bin_data;
       //cout <<" data  =  "<<  bin_data << "  data error = "<< e_data  <<endl;
    }

    chisq_running = calculate_test_statistic(h_data, h_model, g_data, filename_cov);
    
    pval = TMath::Prob(chisq_running, ndof);
    // cout <<"chisq , ndof, pval  =  "<< chisq_running <<"  "<< ndof <<"  "<<  pval << endl;
    results  = std::make_tuple(chisq_running, ndof, pval);
    f_results->Close();
    return results;
}


void write_latex(string mode, vector<string> model, vector<string> vars, float chisq[n_models][n_vars], int ndof[n_models][n_vars],double pval[n_models][n_vars]){
    
    string mode_string;
    string mode_rootfile = mode + ".root";
    int n_vars;
    
    //cout <<"Writing latex "<< mode_rootfile <<endl;
    
    if (mode == "norm_parton"){
        n_vars = vars.size();
        mode_string = "normalised, parton-level";
    }else if (mode == "abs_parton"){
        n_vars = vars.size();
        mode_string = "absolute, parton-level";
    }
    else if (mode == "norm_particle"){
        n_vars = vars.size();
        mode_string = "normalised, particle-level";
    }
    else if (mode == "abs_particle"){
        mode_string = "absolute, particle-level";
        n_vars = vars.size();
    } else if (mode == "norm_parton_bnlo"){
        vars = vars_bnlo;
        n_vars = 8;
        mode_string = "normalised, parton-level";
    }else if (mode == "abs_parton_bnlo"){
        vars = vars_bnlo;
        n_vars = 8;
        mode_string = "absolute, parton-level";
    }

    //cout <<"nvars  "<< n_vars  <<endl;
    
    myfile << "\\multirow{"<< n_models <<"}{*}{} &"<<endl;
    
    for (int i = 0 ; i< n_models; i++){
        if (i < (n_models -1)){
            // myfile << "\\multicolumn{2}{c|}{"<<  model[i] <<"} &" <<endl;
            myfile << "\\multicolumn{2}{c|}{\\parbox{4.3cm}{\\centering "<<  model[i] <<"}} &" <<endl;

        }else{
            myfile << "\\multicolumn{2}{c|}{\\parbox{4.3cm}{\\centering "<<  model[i] <<"}} \\\\" <<endl;
            //myfile << "\\multicolumn{2}{c|}{"<<  model[n_models-1] <<"} \\\\"  <<endl;
        }
    }

    
    for (int m  = 0; m<(n_models-1); m++ ){
        myfile << " & $\\chi^{2}$/ndof & p-val. ";
    }
    
    myfile << " & $\\chi^{2}$/ndof & p-val. \\\\"<< endl;
    myfile << "\\hline"<<endl;


    TFile * f_summary = new TFile(mode_rootfile.c_str(), "RECREATE");
    TH1F * h_summary_A;
    TH1F * h_summary_B ;
    vector<vector<TH1F*>> h_summaries;
    vector<TH1F*> h_summary;
    int nhist;

    std::string summary_tag;
    vector <std::string> temp_bnlo_strings ={"nnlo_ew_lux",
        "nnlo_ew_lux_1725",
        "nnlo_ew_nnpdf",
        "nnlo_nnllprime_nnpdf",
        "nnlo_nnllprime_nnpdf_1725",
        "an3lo__nnpdf",
        "annlo_ct14"
    };
    
    
    for (int summary = 0; summary < model.size(); summary++){
        h_summary.clear();
        
        cout <<"model "<< summary  <<endl;//
        
        
        //"NNLO+\\alphaew^{3} (LUXQED17) m_{t} = 173.3 GeV",
        //"NNLO+\\alphaew^{3} (LUXQED17) m_{t} = 172.5 GeV",
        //"NNLO+\\alphaew^{3} (NNPDF3.1) m_{t} = 173.3 GeV",
       // "NNLO+NNLL' (NNPDF3.1) m_{t} = 173.3 GeV",
       // "NNLO+NNLL' (NNPDF3.1) m_{t} = 172.5 GeV",
       // "aN^{3}LO (NNPDF3.0) m_{t} = 172.5 GeV",
       // "aNNLO (CT14NNLO) m_{t} = 172.5 GeV"

        
        if (model[summary].compare("\\Powheg+\\Pythia") == 0){
        summary_tag = "pwhg_p8";
        }else if (model[summary] == "\\Powheg+\\Herwigpp"){
            summary_tag = "pwhg_hpp";
        }else if (model[summary] == "\\MGaMCatNLO+\\Pythia"){
            summary_tag = "amc_p8";
        }else if (model[summary].compare("NNLO+\alphaew^{3} \\ (LUXQED17) \\ \mt~=~173.3~\GeV") ==0){
            std::cout <<"GOT NNLO LUX  "<< std::endl;
            summary_tag = "nnlo_ew_lux";
        }else if (model[summary] == "NNLO+\alphaew^{3} \\ (LUXQED17) \\ \mt~=~172.5~\GeV"){
            summary_tag = "nnlo_ew_lux_1725";
        }else if (model[summary] == "NNLO+\alphaew^{3} \\ (NNPDF3.1) \\ \mt~=~173.3~\GeV"){
            summary_tag = "nnlo_ew_nnpdf";
        }else if (model[summary] == "NNLO+NNLL' \\ (NNPDF3.1) \\ \mt~=~173.3~\GeV"){
            summary_tag = "nnlo_nnllprime_nnpdf";
        }else if (model[summary] == "NNLO+NNLL' \\ (NNPDF3.1) \\ \mt~=~172.5~\GeV"){
            summary_tag = "nnlo_nnllprime_nnpdf_1725";
        }else if (model[summary] == "aN^{3}LO \\ (NNPDF3.0) \\ \mt~=~172.5~\GeV"){
            summary_tag = "an3lo__nnpdf";
        }else if (model[summary] == "aNNLO \\ (CT14NNLO) \\ \mt~=~172.5~\GeV"){
            summary_tag = "annlo_ct14";
        }
        std::string summary_name_A;
        std::string summary_name_B;
        

        if (mode == "norm_parton_bnlo" || mode == "abs_parton_bnlo"){
            summary_name_A = mode + temp_bnlo_strings[summary] + "_A";
            summary_name_B = mode + summary_tag + "_B";
        }else{
        summary_name_A = mode + summary_tag + "_A";
         summary_name_B = mode + summary_tag + "_B";
        }
       // std::cout <<"model = **"<<  model[summary]  <<"**summary name   = "<< summary_name_A  << std::endl;

        
        double hist_low_lim_A, hist_low_lim_B;
        double hist_high_lim_A,hist_high_lim_B;
        double offset;

        hist_low_lim_A = (summary*0.32);
       // hist_high_lim_A = 15 + (summary*0.32);
         hist_high_lim_A = 14 + (summary*0.32);

        hist_low_lim_B = (summary*0.32);
       // hist_high_lim_B = vars.size() - 15 + (summary*0.32);
         hist_high_lim_B = vars.size() - 14 + (summary*0.32);

        
         //hist_low_lim_B = 0.0;
         //hist_high_lim_B = vars.size() - 15;
        
        if (mode == "norm_particle"){
            hist_low_lim_A = hist_low_lim_A + 0.08;
            hist_low_lim_B = hist_low_lim_A + 0.08;
            hist_high_lim_A = hist_high_lim_A + 0.08;
            hist_high_lim_B = hist_high_lim_B + 0.08;
        }else if (mode == "norm_parton"){
            hist_low_lim_A = hist_low_lim_A - 0.08;
            hist_high_lim_A = hist_high_lim_A - 0.08;
            hist_low_lim_B = hist_low_lim_B - 0.08;
            hist_high_lim_B = hist_high_lim_B - 0.08;
        }
        
       // h_summary_A = new TH1F(summary_name_A.c_str(), summary_name_A.c_str(), 8, 0, 8);
        h_summary_A = new TH1F(summary_name_A.c_str(), summary_name_A.c_str(), 14, hist_low_lim_A, hist_high_lim_A);
        h_summary_B = new TH1F(summary_name_B.c_str(), summary_name_B.c_str(), vars.size() - 14 , hist_low_lim_B , hist_high_lim_B);

        //h_summary_A = new TH1F(summary_name_A.c_str(), summary_name_A.c_str(), 15, 0, 15);
        //h_summary_B = new TH1F(summary_name_B.c_str(), summary_name_B.c_str(), vars.size() - 15 , 0, vars.size() - 15);
        
        if (mode == "norm_particle" || mode == "abs_particle"  ){
           // nhist= 4;
            h_summary.push_back(h_summary_A);
            h_summary.push_back(h_summary_B);
        }else if (mode == "norm_parton" || mode == "abs_parton"){
           // nhist=4;
           // cout <<"in norm parton mode"<< endl;
            h_summary.push_back(h_summary_A);
            h_summary.push_back(h_summary_B);
        }
        else{
           // nhist = 2;
            h_summary.push_back(h_summary_A);
        }
        h_summaries.push_back(h_summary);
       }


    vector<string> vars_A, vars_B ;
    if (mode == "norm_particle" || mode == "norm_parton"){

    for (int i  = 0; i< vars.size() ; i++ ){
        if (i < 14) {
            vars_A.push_back(vars[i]);
        }
        else{
            vars_B.push_back(vars[i]);
        }
    }} else{
        for (int i  = 0; i< vars.size() ; i++)vars_A.push_back(vars[i]);
    }
    
    /*
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
     
     */
    
    //cout <<"Written hist summaries " <<endl;
    stringstream stream;
    std::string pval_str;
    std::string chisq_str;
    std::string ndof_str;
    
    for (int i = 0; i< n_vars; i++){
        myfile<<vars[i];
        for (int m  = 0; m< n_models; m++ ){
        //    cout <<"latex check, var =  "<< vars[i] <<" model =  "<< model[m] << " pval = " << pval[m][i] <<  "  ndof = " <<  ndof[m][i]<<endl;
            pval_str = std::to_string(pval[m][i]);
        if (pval[m][i] < 0.001 && pval[m][i] >= 0.0){
            pval_str = "$< 10^{-3}$";
                pval[m][i] = low_lim*1.01;
                stream << fixed << setprecision(0) << chisq[m][i];
                chisq_str = stream.str();
                stream.str("");
                stream << fixed << setprecision(0) << ndof[m][i];
                ndof_str = stream.str();
                stream.str("");
                myfile<<" & "<< chisq_str <<"/"<< ndof_str <<" & "<< pval_str;
        } else if (pval[m][i] == -1.0){
            pval_str = "-";
            chisq_str = "-";
            ndof_str = "-";
            myfile<<" & - & -";
        }else {
                stream << fixed << setprecision(3) << pval[m][i];
                pval_str = stream.str();
                stream.str("");
                stream << fixed << setprecision(0) << chisq[m][i];
                chisq_str = stream.str();
                stream.str("");
                stream << fixed << setprecision(0) << ndof[m][i];
                ndof_str = stream.str();
                stream.str("");
                myfile<<" & "<< chisq_str <<"/"<< ndof_str <<" & "<< pval_str;
            }
        }
        myfile<<" \\\\"<<endl;
        }
    
        /*
    for (int i  = 0; i< n_vars; i++ ){
        string pval_str_0 = std::to_string(pval[0][i]);
        if (pval[0][i] < 0.001){
            pval_str_0 = "$<$ 0.01";
            pval[0][i] = low_lim*1.01;
        } else {
            stream << fixed << setprecision(2) << pval[0][i] ;
            pval_str_0 = stream.str();
            stream.str("");
        }
        
        string pval_str_1 = std::to_string(pval[1][i]);
        if (pval[1][i] < 0.001){
            pval_str_1 = "$<$ 0.01$";
            pval[1][i] = low_lim*1.01;
        }else {
            stream << fixed << setprecision(2) << pval[1][i] ;
            pval_str_1 = stream.str();
            stream.str("");
        }
        
        string pval_str_2 = std::to_string(pval[2][i]);
        
        if (pval[2][i] < 0.001){
            pval_str_2 = "$<$ 0.01";
            pval[2][i] = low_lim*1.01;
        }else {
            stream << fixed << setprecision(2) << pval[2][i] ;
            pval_str_2 = stream.str();
            stream.str("");
        }
        
        myfile<<vars[i] <<" & "<<chisq[0][i]<<"/"<< ndof[0][i]<<"& "<<pval_str_0 <<" & "<<chisq[1][i]<<"/"<<ndof[1][i]<<" & "<<pval_str_1<< " & "<<chisq[2][i]<<"/"<<ndof[2][i]<<" & "<< pval_str_2 <<" \\\\"<<endl;
        
    }
        */
    myfile << "\\hline"<<endl;
    myfile << "\\end{tabular}"<<endl;
    myfile << "\\label{tab:"<< mode << "}" << endl;
    
   // cout <<"N vars A = "<< vars_A.size() << endl;
    
    for (int m = 0; m< n_models; m++){
       // cout <<"looping models " << endl;

        if (mode == "norm_particle" || mode == "norm_parton" || mode == "abs_particle" || mode == "abs_parton"){
            for (int i  = 0; i< vars_A.size() ; i++ ){
                h_summaries[m][0]->SetBinContent(i+1, pval[m][i]);
                h_summaries[m][0]->SetBinContent(i+1, pval[m][i]);
                h_summaries[m][0]->GetXaxis()->SetBinLabel(i+1,vars_root[i].c_str());
                h_summaries[m][0]->GetXaxis()->SetBinLabel(i+1,vars_root[i].c_str());
            }
            for (int i  = 0; i< vars_B.size() ; i++ ){
                h_summaries[m][1]->SetBinContent(i+1, pval[m][i+14]);
                h_summaries[m][1]->SetBinContent(i+1, pval[m][i+14]);
                h_summaries[m][1]->GetXaxis()->SetBinLabel(i+1,vars_root[i+14].c_str());;
                h_summaries[m][1]->GetXaxis()->SetBinLabel(i+1,vars_root[i+14].c_str());;
            }
            h_summaries[m][0]->Write();
            h_summaries[m][1]->Write();
            
        }else{
            for (int i  = 0; i< vars_A.size() ; i++ ){
           //     cout <<"looping vars a " << endl;

                h_summaries[m][0]->SetBinContent(i+1, pval[m][i]);
                h_summaries[m][0]->SetBinContent(i+1, pval[m][i]);
                h_summaries[m][0]->GetXaxis()->SetBinLabel(i+1,vars_root_bnlo[i].c_str());
                h_summaries[m][0]->GetXaxis()->SetBinLabel(i+1,vars_root_bnlo[i].c_str());
            }
            h_summaries[m][0]->Write();
        }
    }
    
   // cout <<" filled tables" << endl;

    
    f_summary->Close();
    
    myfile << "\\caption{The $\\chi^{2}$/ndof and p values quantifying the agreement between theoretical predictions and data for "<< mode_string <<" measurements are shown.}"<< endl;
    
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
    
  //  cout <<"summary plot 1  " <<endl;

    TCanvas * c_master = new TCanvas("c_master","c_master");
    gStyle->SetOptStat(00000);

    if (f1 == "norm_parton.root" && f2  == "norm_particle.root"){
    c_master->SetCanvasSize(2200, 1800);
    c_master->Divide(1,2);
    c_master->cd(1);
    }

    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetFrameLineColor(0);
    gPad->SetLeftMargin(0.085);
    gPad->SetRightMargin(0.01);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.33);
    
    TFile * f_norm_parton = new TFile(f1.c_str());
    TFile * f_norm_particle = new TFile(f2.c_str());
    
    if (f1 == "norm_parton_bnlo.root"){
    h_norm_bnlo_nnlo_ew_lux_A =  (TH1F*)f_norm_parton->Get("norm_parton_bnlonnlo_ew_lux_A;1");
    h_norm_bnlo_nnlo_ew_lux_1725_A =  (TH1F*)f_norm_parton->Get("norm_parton_bnlonnlo_ew_lux_1725_A;1");
    h_norm_bnlo_nnlo_ew_nnpdf_A =  (TH1F*)f_norm_parton->Get("norm_parton_bnlonnlo_ew_nnpdf_A;1");
    h_norm_bnlo_nnlo_nnllprime_nnpdf_A =  (TH1F*)f_norm_parton->Get("norm_parton_bnlonnlo_nnllprime_nnpdf_A;1");
    h_norm_bnlo_nnlo_nnllprime_nnpdf_1725_A =  (TH1F*)f_norm_parton->Get("norm_parton_bnlonnlo_nnllprime_nnpdf_1725_A;1");
    h_norm_bnlo_an3lo_nnpdf_A =  (TH1F*)f_norm_parton->Get("norm_parton_bnloan3lo__nnpdf_A;1");
    h_norm_bnlo_annlo_ct14_A =  (TH1F*)f_norm_parton->Get("norm_parton_bnloannlo_ct14_A;1");
    }else{
    //pwhg_p8
     h_norm_parton_pwhg_p8_A =  (TH1F*)f_norm_parton->Get("norm_partonpwhg_p8_A;1");
     h_norm_particle_pwhg_p8_A =  (TH1F*)f_norm_particle->Get("norm_particlepwhg_p8_A;1");
     h_norm_particle_pwhg_p8_B =  (TH1F*)f_norm_particle->Get("norm_particlepwhg_p8_B;1");
    //pwhg_hpp
     h_norm_parton_pwhg_hpp_A =  (TH1F*)f_norm_parton->Get("norm_partonpwhg_hpp_A;1");
     h_norm_particle_pwhg_hpp_A =  (TH1F*)f_norm_particle->Get("norm_particlepwhg_hpp_A;1");
     h_norm_particle_pwhg_hpp_B =  (TH1F*)f_norm_particle->Get("norm_particlepwhg_hpp_B;1");
    //amc_p8
     h_norm_parton_amc_p8_A =  (TH1F*)f_norm_parton->Get("norm_partonamc_p8_A;1");
     h_norm_particle_amc_p8_A =  (TH1F*)f_norm_particle->Get("norm_particleamc_p8_A;1");
     h_norm_particle_amc_p8_B =  (TH1F*)f_norm_particle->Get("norm_particleamc_p8_B;1");
    }
    
    TFile * summary = new TFile("summary.root", "RECREATE");
    float x_lim = n_vars;

    if (f1 == "norm_parton_bnlo.root"){
        double unit_shift = 0.56;
        double offset = -0.5;

        // offset bin markers
        h_norm_bnlo_nnlo_ew_lux_A_clone = (TH1F*)h_norm_bnlo_nnlo_ew_lux_A->Clone("h_norm_bnlo_nnlo_ew_lux_A_clone");
        TF1 *f0 = new TF1("f0","1",0,n_vars);
        //h_norm_bnlo_nnlo_ew_lux_A_clone->Multiply(f0,1);
        Double_t shift = (-3.0*unit_shift)*h_norm_bnlo_nnlo_ew_lux_A->GetBinWidth(1);
        h_norm_bnlo_nnlo_ew_lux_A_clone->GetXaxis()->SetLimits(0+shift,x_lim+shift);
        
        h_norm_bnlo_nnlo_ew_lux_1725_A_clone = (TH1F*)h_norm_bnlo_nnlo_ew_lux_A->Clone("h_norm_bnlo_nnlo_ew_lux_1725_A_clone");
        TF1 *f1 = new TF1("f1","1",0,n_vars);
        // h_norm_bnlo_nnlo_ew_lux_1725_A_clone->Multiply(f1,1);
        shift = (-2.0*unit_shift)*h_norm_bnlo_nnlo_ew_lux_1725_A->GetBinWidth(1);
        h_norm_bnlo_nnlo_ew_lux_1725_A_clone->GetXaxis()->SetLimits(0+shift,x_lim+shift);

        h_norm_bnlo_nnlo_ew_nnpdf_A_clone = (TH1F*)h_norm_bnlo_nnlo_ew_nnpdf_A->Clone("h_norm_bnlo_nnlo_ew_nnpdf_A_clone");
      //  TF1 *f2 = new TF1("f1","1",0,n_vars);
      //  h_norm_bnlo_nnlo_ew_nnpdf_A_clone->Multiply(f1,1);
        shift = (-1.0*unit_shift)*h_norm_bnlo_nnlo_ew_nnpdf_A->GetBinWidth(1);
        h_norm_bnlo_nnlo_ew_nnpdf_A_clone->GetXaxis()->SetLimits(0+shift,x_lim+shift);
        
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone = (TH1F*)h_norm_bnlo_nnlo_nnllprime_nnpdf_A->Clone("h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone");
       // TF1 *f2 = new TF1("f2","1",0,n_vars);
       // h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->Multiply(f2,1);
        shift = (0.0*unit_shift)*h_norm_bnlo_nnlo_nnllprime_nnpdf_A->GetBinWidth(1);
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->GetXaxis()->SetLimits(0+shift,x_lim+shift);
        
        h_norm_bnlo_nnlo_nnllprime_nnpdf_1725_A_clone = (TH1F*)h_norm_bnlo_nnlo_nnllprime_nnpdf_1725_A->Clone("h_norm_bnlo_nnlo_nnllprime_nnpdf_1725_A_clone");
       // TF1 *f2 = new TF1("f2","1",0,n_vars);
       // h_norm_bnlo_nnlo_nnllprime_nnpdf_1725_A_clone->Multiply(f2,1);
        shift = (1.0*unit_shift)*h_norm_bnlo_nnlo_nnllprime_nnpdf_1725_A->GetBinWidth(1);
        h_norm_bnlo_nnlo_nnllprime_nnpdf_1725_A_clone->GetXaxis()->SetLimits(0+shift,x_lim+shift);
        
        h_norm_bnlo_an3lo_nnpdf_A_clone = (TH1F*)h_norm_bnlo_an3lo_nnpdf_A->Clone("h_norm_bnlo_an3lo_nnpdf_A_clone");
       // TF1 *f3 = new TF1("f3","1",0,n_vars);
       // h_norm_bnlo_an3lo_nnpdf_A_clone->Multiply(f3,1);
        shift = (2.0*unit_shift)*h_norm_bnlo_an3lo_nnpdf_A->GetBinWidth(1);
        h_norm_bnlo_an3lo_nnpdf_A_clone->GetXaxis()->SetLimits(0+shift,x_lim+shift);
        
        h_norm_bnlo_annlo_ct14_A_clone = (TH1F*)h_norm_bnlo_annlo_ct14_A->Clone("h_norm_bnlo_annlo_ct14_A_clone");
       // TF1 *f4 = new TF1("f4","1",0,n_vars);
       // h_norm_bnlo_annlo_ct14_A_clone->Multiply(f4,1);
        shift = (3.0*unit_shift)*h_norm_bnlo_annlo_ct14_A_clone->GetBinWidth(1);
        h_norm_bnlo_annlo_ct14_A_clone->GetXaxis()->SetLimits(0+shift,x_lim+shift);
        
    //draw
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->Draw("p");
        h_norm_bnlo_nnlo_ew_lux_A_clone->Draw("psame");
        h_norm_bnlo_nnlo_ew_lux_1725_A_clone->Draw("psame");
        h_norm_bnlo_nnlo_ew_nnpdf_A_clone->Draw("psame");
        h_norm_bnlo_nnlo_nnllprime_nnpdf_1725_A_clone->Draw("psame");
        h_norm_bnlo_an3lo_nnpdf_A_clone->Draw("psame");
        h_norm_bnlo_annlo_ct14_A_clone->Draw("psame");
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->Draw("psame");


        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->SetMinimum(low_lim);
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->SetMaximum(1.1);
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->GetXaxis()->SetLabelSize(0.06);
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->GetXaxis()->LabelsOption("d");
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->GetXaxis()->CenterLabels(kTRUE);
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->GetXaxis()->SetNdivisions(32, kFALSE);
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->GetXaxis()->SetLabelOffset(0.02);
        
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->GetYaxis()->SetTickLength(0.01);
       // h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->GetYaxis()->SetTitleOffset(0.6);
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->GetYaxis()->SetTitleSize(0.06);
        
//        vector<int> marker_styles = {21,22,23,27,28,29};
  //      vector<int> marker_colours = {3,4,6,7,8,9};
        h_norm_bnlo_nnlo_ew_lux_A_clone->SetMarkerStyle(21);
        h_norm_bnlo_nnlo_ew_lux_A_clone->SetMarkerColor(3);

        h_norm_bnlo_nnlo_ew_lux_1725_A_clone->SetMarkerStyle(21);
        h_norm_bnlo_nnlo_ew_lux_1725_A_clone->SetMarkerColor(28);
       // h_norm_bnlo_nnlo_ew_lux_1725_A_clone->SetMarkerColor(6);

        h_norm_bnlo_nnlo_ew_nnpdf_A_clone->SetMarkerStyle(25);
        h_norm_bnlo_nnlo_ew_nnpdf_A_clone->SetMarkerColor(4);
        //h_norm_bnlo_nnlo_ew_nnpdf_A_clone->SetMarkerColor(3);

        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->SetMarkerStyle(20);
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->SetMarkerColor(6);
       // h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->SetMarkerColor(3);

        h_norm_bnlo_nnlo_nnllprime_nnpdf_1725_A_clone->SetMarkerStyle(24);
        h_norm_bnlo_nnlo_nnllprime_nnpdf_1725_A_clone->SetMarkerColor(32);
       // h_norm_bnlo_nnlo_nnllprime_nnpdf_1725_A_clone->SetMarkerColor(6);

        h_norm_bnlo_an3lo_nnpdf_A_clone->SetMarkerStyle(27);
        h_norm_bnlo_an3lo_nnpdf_A_clone->SetMarkerColor(7);
       // h_norm_bnlo_an3lo_nnpdf_A_clone->SetMarkerColor(6);

        h_norm_bnlo_annlo_ct14_A_clone->SetMarkerStyle(28);
        h_norm_bnlo_annlo_ct14_A_clone->SetMarkerColor(8);
        //h_norm_bnlo_annlo_ct14_A_clone->SetMarkerColor(6);


        h_norm_bnlo_nnlo_ew_lux_A_clone->SetMarkerSize(1.2);
        h_norm_bnlo_nnlo_ew_lux_1725_A_clone->SetMarkerSize(1.2);
        h_norm_bnlo_nnlo_ew_nnpdf_A_clone->SetMarkerSize(1.2);
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->SetMarkerSize(1.2);
        h_norm_bnlo_nnlo_nnllprime_nnpdf_1725_A_clone->SetMarkerSize(1.2);
        h_norm_bnlo_an3lo_nnpdf_A_clone->SetMarkerSize(1.2);
        h_norm_bnlo_annlo_ct14_A_clone->SetMarkerSize(1.2);

        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->SetYTitle("p-value");
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->GetYaxis()->SetTitleOffset(0.68);
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->SetTitle(" ");
        h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone->SetTitle(" ");
    }
    else {
       // cout <<"summary plot 4  " <<endl;

    double unit_shift = -0.03;
    double offset = -0.5;
    // offset bin markers
    h_norm_particle_pwhg_p8_A_clone = (TH1F*)h_norm_particle_pwhg_p8_A->Clone("h_norm_particle_pwhg_p8_A_clone");
        
   // cout <<"summary plot 41  " <<endl;

    TF1 *f0 = new TF1("f0","1",0,n_vars);
  //  h_norm_particle_pwhg_p8_A_clone->Multiply(f0,1);
    Double_t shift = offset + (-1.6*unit_shift*h_norm_particle_pwhg_p8_A->GetBinWidth(1));
    //Double_t  shift = -0.06;
    cout <<"shift = " << shift <<endl;

//    h_norm_particle_pwhg_p8_A_clone->GetXaxis()->SetLimits(0+shift,x_lim+shift);

    h_norm_parton_pwhg_p8_A_clone = (TH1F*)h_norm_parton_pwhg_p8_A->Clone("h_norm_parton_pwhg_p8_A_clone");
    TF1 *f1 = new TF1("f1","1",0,n_vars);
 //   h_norm_parton_pwhg_p8_A_clone->Multiply(f1,1);
   shift = offset + (1.6*unit_shift*h_norm_parton_pwhg_p8_A->GetBinWidth(1));
    //shift = -0.04;
   // cout <<"shift = " << shift <<endl;

  //  h_norm_parton_pwhg_p8_A_clone->GetXaxis()->SetLimits(0+shift,x_lim+shift);
        
    h_norm_particle_pwhg_hpp_A_clone = (TH1F*)h_norm_particle_pwhg_hpp_A->Clone("h_norm_particle_pwhg_hpp_A_clone");
    TF1 *f3 = new TF1("f3","1",0,n_vars);
  //  h_norm_particle_pwhg_hpp_A_clone->Multiply(f3,1);
   shift = offset + (-2.6*unit_shift*h_norm_particle_pwhg_hpp_A->GetBinWidth(1));
   // shift = -0.02;
   // cout <<"shift = " << shift <<endl;

    //h_norm_particle_pwhg_hpp_A_clone->GetXaxis()->SetLimits(0+shift,x_lim+shift);
        
    h_norm_parton_pwhg_hpp_A_clone = (TH1F*)h_norm_parton_pwhg_hpp_A->Clone("h_norm_parton_pwhg_hpp_A_clone");
    TF1 *f2 = new TF1("f2","1",0,n_vars);
 //   h_norm_parton_pwhg_hpp_A_clone->Multiply(f2,1);
 shift = offset + (2.6*unit_shift*h_norm_parton_pwhg_hpp_A->GetBinWidth(1));
    //shift = 0.02;
    //cout <<"shift = " << shift <<endl;

   // h_norm_parton_pwhg_hpp_A_clone->GetXaxis()->SetLimits(0+shift,x_lim+shift);
        
    h_norm_particle_amc_p8_A_clone = (TH1F*)h_norm_particle_amc_p8_A->Clone("h_norm_particle_amc_p8_A_clone");
    TF1 *f5 = new TF1("f5","1",0,n_vars);
//    h_norm_particle_amc_p8_A_clone->Multiply(f5,1);
    shift =  offset + (-3.9*unit_shift*h_norm_particle_amc_p8_A->GetBinWidth(1));
    //shift = 0.04;
   // cout <<"shift = " << shift <<endl;

  //  h_norm_particle_amc_p8_A_clone->GetXaxis()->SetLimits(0+shift,x_lim+shift);
        
    h_norm_parton_amc_p8_A_clone = (TH1F*)h_norm_parton_amc_p8_A->Clone("h_norm_parton_amc_p8_A_clone");
    TF1 *f4 = new TF1("f4","1",0,n_vars);
  //  h_norm_parton_amc_p8_A_clone->Multiply(f4,1);
    shift = offset + (3.9*unit_shift*h_norm_parton_amc_p8_A->GetBinWidth(1));
    //shift = 0.06;
       // cout <<"shift = " << shift <<endl;

   // h_norm_parton_amc_p8_A_clone->GetXaxis()->SetLimits(0+shift,x_lim+shift);
        
        //Draw
    TH1F * h_base_A = (TH1F*)h_norm_particle_pwhg_hpp_A_clone->Clone();
    h_base_A->Reset();
    h_base_A->GetXaxis()->SetLimits(0.3,14.3);

    h_base_A->Draw();
    h_norm_particle_pwhg_p8_A_clone->Draw("psame");
    h_norm_parton_pwhg_p8_A_clone->Draw("psame");
    h_norm_particle_pwhg_hpp_A_clone->Draw("psame");
    h_norm_parton_pwhg_hpp_A_clone->Draw("psame");
    h_norm_particle_amc_p8_A_clone->Draw("psame");
    h_norm_parton_amc_p8_A_clone->Draw("psame");
        
    h_base_A->SetMinimum(low_lim);
    h_base_A->SetMaximum(1.1);
    h_base_A->GetXaxis()->SetLabelSize(0.07);
    h_base_A->GetXaxis()->LabelsOption("v");
    h_base_A->GetXaxis()->CenterLabels(kTRUE);
    h_base_A->GetXaxis()->SetNdivisions(32, kFALSE);
    h_base_A->GetXaxis()->SetLabelOffset(0.02);
        
    h_base_A->GetYaxis()->SetTickLength(0.01);
    h_base_A->GetYaxis()->SetTitleOffset(0.6);
    h_base_A->GetYaxis()->SetTitleSize(0.07);
    h_base_A->GetYaxis()->SetLabelSize(0.06);

        
    h_norm_parton_pwhg_p8_A_clone->SetMarkerStyle(24);
    h_norm_parton_pwhg_hpp_A_clone->SetMarkerStyle(25);
    h_norm_parton_amc_p8_A_clone->SetMarkerStyle(26);
        
    h_norm_particle_pwhg_p8_A_clone->SetMarkerStyle(20);
    h_norm_particle_pwhg_hpp_A_clone->SetMarkerStyle(21);
    h_norm_particle_amc_p8_A_clone->SetMarkerStyle(22);
        
    h_norm_parton_pwhg_p8_A_clone->SetMarkerSize(2.2);
    h_norm_parton_pwhg_hpp_A_clone->SetMarkerSize(2.2);
    h_norm_parton_amc_p8_A_clone->SetMarkerSize(2.2);
    
    h_norm_particle_pwhg_p8_A_clone->SetMarkerSize(2.2);
    h_norm_particle_pwhg_hpp_A_clone->SetMarkerSize(2.2);
    h_norm_particle_amc_p8_A_clone->SetMarkerSize(2.2);
        
    h_norm_parton_pwhg_p8_A_clone->SetMarkerColor(kRed);
    h_norm_parton_pwhg_hpp_A_clone->SetMarkerColor(kViolet);
    h_norm_parton_amc_p8_A_clone->SetMarkerColor(kGreen);
        
    h_norm_particle_pwhg_p8_A_clone->SetMarkerColor(kRed);
    h_norm_particle_pwhg_hpp_A_clone->SetMarkerColor(kViolet);
    h_norm_particle_amc_p8_A_clone->SetMarkerColor(kGreen);
        
    h_base_A->SetYTitle("p-value");
    h_base_A->SetTitle(" ");
    }
    
    //cout <<"summary plot 2  " <<endl;

    c_master->cd(2);
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetFrameLineColor(0);
    gPad->SetLeftMargin(0.085);
    gPad->SetRightMargin(0.03);
    //gPad->SetTopMargin(0.185);
    //gPad->SetBottomMargin(0.135);
    //gPad->SetTopMargin(0.16);
//    gPad->SetBottomMargin(0.28);
    gPad->SetBottomMargin(0.37);

    
    
    if (h_norm_particle_pwhg_p8_B){
        //h_norm_particle_pwhg_p8_B->SetMinimum(0.000000000000000000000000000000000001);
        h_norm_particle_pwhg_p8_B->SetMinimum(low_lim);
        h_norm_particle_pwhg_p8_B->SetMaximum(1.1);
        h_norm_particle_pwhg_p8_B->GetXaxis()->LabelsOption("v");
        h_norm_particle_pwhg_p8_B->GetXaxis()->CenterLabels(kTRUE);
        h_norm_particle_pwhg_p8_B->GetXaxis()->SetNdivisions(32, kFALSE);
        h_norm_particle_pwhg_p8_B->GetXaxis()->SetLabelOffset(0.01);
        //h_norm_particle_pwhg_p8_B->GetXaxis()->LabelsOption("u");
        
        h_norm_particle_pwhg_p8_B->GetYaxis()->SetTickLength(0.01);
        h_norm_particle_pwhg_p8_B->GetYaxis()->SetTitleOffset(0.6);
        h_norm_particle_pwhg_p8_B->GetYaxis()->SetTitleSize(0.055);
        h_norm_particle_pwhg_p8_B->GetXaxis()->SetLabelSize(0.07);
        
        h_norm_particle_pwhg_p8_B->SetYTitle("p-value");
        
        h_norm_particle_pwhg_p8_B->SetTitle(" ");
       // h_norm_particle_pwhg_p8_B->Draw("p");
       // h_norm_particle_pwhg_hpp_B->Draw("psame");
       // h_norm_particle_amc_p8_B->Draw("psame");

       // h_norm_particle_pwhg_p8_B->SetMarkerStyle(22);
       // h_norm_particle_pwhg_hpp_B->SetMarkerStyle(22);
       // h_norm_particle_amc_p8_B->SetMarkerStyle(22);

        //h_norm_particle_pwhg_p8_B->SetMarkerSize(2.2);
        //h_norm_particle_pwhg_hpp_B->SetMarkerSize(2.2);
        //h_norm_particle_amc_p8_B->SetMarkerSize(2.2);

        //h_norm_particle_pwhg_p8_B->SetMarkerColor(kRed);
        //h_norm_particle_pwhg_hpp_B->SetMarkerColor(kViolet);
        //h_norm_particle_amc_p8_B->SetMarkerColor(kGreen);
        
        double unit_shift = 0.05;
        
        // offset bin markers
        h_norm_particle_pwhg_p8_B_clone = (TH1F*)h_norm_particle_pwhg_p8_B->Clone("h_norm_particle_pwhg_p8_B_clone");
        TF1 *f0 = new TF1("f0","1",0,14);
        h_norm_particle_pwhg_p8_B_clone->Multiply(f0,1);
        Double_t shift = (-1.3*unit_shift)*h_norm_particle_pwhg_p8_B->GetBinWidth(1);
      //  h_norm_particle_pwhg_p8_B_clone->GetXaxis()->SetLimits(0+shift,x_lim+shift);

        h_norm_particle_pwhg_hpp_B_clone = (TH1F*)h_norm_particle_pwhg_hpp_B->Clone("h_norm_particle_pwhg_hpp_B_clone");
        TF1 *f3 = new TF1("f3","1",0,14);
        h_norm_particle_pwhg_hpp_B_clone->Multiply(f3,1);
        shift = -3.3*(unit_shift)*h_norm_particle_pwhg_hpp_B->GetBinWidth(1);
      //  h_norm_particle_pwhg_hpp_B_clone->GetXaxis()->SetLimits(0+shift,x_lim+shift);
    
        h_norm_particle_amc_p8_B_clone = (TH1F*)h_norm_particle_amc_p8_B->Clone("h_norm_particle_amc_p8_B_clone");
        TF1 *f5 = new TF1("f5","1",0,14);
        h_norm_particle_amc_p8_B_clone->Multiply(f5,1);
        shift = -4.7*(unit_shift)*h_norm_particle_amc_p8_B->GetBinWidth(1);
      //  h_norm_particle_amc_p8_B_clone->GetXaxis()->SetLimits(0+shift,x_lim+shift);
        
        TH1F * h_base_B = (TH1F*)h_norm_particle_pwhg_hpp_B_clone->Clone();
        h_base_B->Reset();
        
        h_base_B->SetMinimum(low_lim);
        h_base_B->SetMaximum(1.1);
        h_base_B->GetXaxis()->SetLabelSize(0.07);
        h_base_B->GetXaxis()->LabelsOption("v");
        h_base_B->GetXaxis()->CenterLabels(kTRUE);
        h_base_B->GetXaxis()->SetNdivisions(32, kFALSE);
        h_base_B->GetXaxis()->SetLabelOffset(0.02);
        
        h_base_B->GetYaxis()->SetTickLength(0.01);
        h_base_B->GetYaxis()->SetTitleOffset(0.6);
        h_base_B->GetYaxis()->SetTitleSize(0.07);
        h_base_B->GetYaxis()->SetLabelSize(0.06);

        h_norm_particle_pwhg_p8_B_clone->SetMarkerSize(2.4);
        h_norm_particle_pwhg_hpp_B_clone->SetMarkerSize(2.4);
        h_norm_particle_amc_p8_B_clone->SetMarkerSize(2.4);
        
        h_norm_particle_pwhg_p8_B_clone->SetMarkerStyle(20);
        h_norm_particle_pwhg_hpp_B_clone->SetMarkerStyle(21);
        h_norm_particle_amc_p8_B_clone->SetMarkerStyle(22);
        
        h_norm_particle_pwhg_p8_B_clone->SetMarkerColor(kRed);
        h_norm_particle_pwhg_hpp_B_clone->SetMarkerColor(kViolet);
        h_norm_particle_amc_p8_B_clone->SetMarkerColor(kGreen);
        
        h_base_B->SetYTitle("p-value");
        h_base_B->SetTitle(" ");
        
        h_base_B->GetXaxis()->SetLimits(0.423,19.423);
        h_base_B->Draw();
        h_norm_particle_pwhg_p8_B_clone->Draw("psame");
        h_norm_particle_pwhg_hpp_B_clone->Draw("psame");
        h_norm_particle_amc_p8_B_clone->Draw("psame");
    
    }

    if (f1 == "norm_parton_bnlo.root"){
        
        float gap = 0.035;
        float width = (1.0 - (gap*2.0))/3.0;
        float bottom = 0.0;
        float top = 0.23;
//        float textsize = 0.02155;
        float textsize = 0.0255;
        float offset = 0.005;


        auto legend_A_1 = new TLegend(offset,bottom,width+offset,top);
        legend_A_1->AddEntry(h_norm_bnlo_nnlo_ew_lux_A_clone,"NNLO+#alpha^{3}_{EW} (LUXQED17) m_{t} = 173.3 GeV","p");
        legend_A_1->AddEntry(h_norm_bnlo_nnlo_ew_lux_1725_A_clone,"NNLO+#alpha^{3}_{EW} (LUXQED17) m_{t} = 172.5 GeV","p");
        legend_A_1->AddEntry(h_norm_bnlo_nnlo_ew_nnpdf_A_clone,"NNLO+#alpha^{3}_{EW} (NNPDF3.1) m_{t} = 173.3 GeV","p");
        legend_A_1->SetTextSize(textsize);
        legend_A_1->SetFillColor(0);
        legend_A_1->SetBorderSize(0);
        legend_A_1->SetEntrySeparation(0.05);
        legend_A_1->SetMargin(0.1);
        legend_A_1->Draw();
        
        auto legend_A_2 = new TLegend((width+gap+offset),bottom,(width*2)+(gap)+offset,top);
        legend_A_2->AddEntry(h_norm_bnlo_nnlo_nnllprime_nnpdf_A_clone,"NNLO+NNLL' (NNPDF3.1) m_{t} = 173.3 GeV","p");
        legend_A_2->AddEntry(h_norm_bnlo_nnlo_nnllprime_nnpdf_1725_A_clone,"NNLO+NNLL' (NNPDF3.1) m_{t} = 172.5 GeV","p");
        legend_A_2->SetTextSize(textsize);
        legend_A_2->SetFillColor(0);
        legend_A_2->SetBorderSize(0);
        legend_A_2->SetMargin(0.1);

        legend_A_2->Draw();
        
        auto legend_A_3 = new TLegend((width*2)+(gap*2)+offset,bottom,1.0,top);
        legend_A_3->AddEntry(h_norm_bnlo_an3lo_nnpdf_A_clone,"aN^{3}LO (NNPDF3.0) m_{t} = 172.5 GeV","p");
        legend_A_3->AddEntry(h_norm_bnlo_annlo_ct14_A_clone,"aNNLO (CT14NNLO) m_{t} = 172.5 GeV","p");
        legend_A_3->SetTextSize(textsize);
        legend_A_3->SetFillColor(0);
        legend_A_3->SetBorderSize(0);
        legend_A_3->SetMargin(0.1);

        legend_A_3->Draw();
    
    }else{
    
    c_master->cd(0);
    auto legend_A_1 = new TLegend(0.0,0.46,0.25,0.54);
    legend_A_1->AddEntry(h_norm_parton_pwhg_p8_A_clone,"POWHEG+PYTHIA8 (parton level)","p");
    legend_A_1->AddEntry(h_norm_particle_pwhg_p8_A_clone,"POWHEG+PYTHIA8 (particle level)","p");
    legend_A_1->SetTextSize(0.0193);
    legend_A_1->SetFillColor(0);
    legend_A_1->SetBorderSize(0);
    legend_A_1->Draw();
        
    auto legend_A_2 = new TLegend(0.297,0.46,0.547,0.54);
    legend_A_2->AddEntry(h_norm_parton_pwhg_hpp_A_clone,"POWHEG+HERWIG++ (parton level)","p");
    legend_A_2->AddEntry(h_norm_particle_pwhg_hpp_A_clone,"POWHEG+HERWIG++ (particle level)","p");
    legend_A_2->SetTextSize(0.0193);
    legend_A_2->SetFillColor(0);
    legend_A_2->SetBorderSize(0);
    legend_A_2->Draw();
        
    auto legend_A_3 = new TLegend(0.61,0.46,0.86,0.54);
    legend_A_3->AddEntry(h_norm_parton_amc_p8_A_clone,"MG5_aMC@NLO+PYTHIA8 [FxFx] (parton level)","p");
    legend_A_3->AddEntry(h_norm_particle_amc_p8_A_clone,"MG5_aMC@NLO+PYTHIA8 [FxFx] (particle level)","p");
    legend_A_3->SetTextSize(0.0193);
    legend_A_3->SetFillColor(0);
    legend_A_3->SetBorderSize(0);
    legend_A_3->Draw();
    
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
    
    c_master->cd(1);

    float H = c_master->GetWh();
    float W = c_master->GetWw();
    float l = c_master->GetLeftMargin();
    float t = c_master->GetTopMargin();
    float r = c_master->GetRightMargin();
    float b = c_master->GetBottomMargin();
    float extraOverCmsTextSize  = 0.9;
    
    if (f1 == "norm_parton_bnlo.root"){
        TString cmsText, extraText, lumiText;
        cmsText += "CMS";
        extraText += "";
        lumiText += "35.9 fb^{-1} (13 TeV)";
        TLatex latex;
        latex.SetNDC();
        latex.SetTextAngle(0);
        latex.SetTextSize(0.4*t);
        latex.SetTextColor(kBlack);
        latex.SetTextFont(61);
        latex.SetTextAlign(31);
        latex.DrawLatex(0.15,0.92,cmsText);
        
        latex.SetTextFont(52);
        latex.SetTextSize(0.28*t*extraOverCmsTextSize);
        latex.DrawLatex(0.37,0.92,extraText);
        
        latex.SetTextFont(42);
        latex.SetTextSize(0.4*t);
        latex.DrawLatex(0.97,0.92,lumiText);
        
        latex.SetTextSize(0.28*t);
        latex.DrawLatex(0.034,0.36,"<");
    }else{
    
    TString cmsText, extraText, lumiText;
    cmsText += "CMS";
    extraText += "";
    lumiText += "35.9 fb^{-1} (13 TeV)";
    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextSize(0.75*t);
    latex.SetTextColor(kBlack);
    latex.SetTextFont(61);
    latex.SetTextAlign(31);
    latex.DrawLatex(0.15,0.92,cmsText);
    
    latex.SetTextFont(52);
    latex.SetTextSize(0.75*t);
    latex.DrawLatex(0.31,0.92,extraText);
    
    latex.SetTextFont(42);
    latex.SetTextSize(0.6*t);
    latex.DrawLatex(0.99,0.92,lumiText);

    latex.SetTextSize(0.45*t);
    latex.DrawLatex(0.032,0.31,"<");
    c_master->cd(2);
    latex.DrawLatex(0.032,0.36,"<");

    }

    
    c_master->Write();
    c_master->SaveAs("GOF_summary_master.pdf");
    
    return;
    
}

