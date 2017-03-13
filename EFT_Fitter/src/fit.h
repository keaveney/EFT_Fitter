#include <iostream>
#include <sstream>
#include <string.h>
#include "TH1F.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include <utility>


using namespace std;

bool debug = true;



class Fitter {
    
    public:
        Fitter(string, string);
        void run_extraction(int, float bins[], std::string, std::string, bool, std::string);
        std::pair <TH1F*, vector<TH1F *>> initialise(std::string, std::string, int, float bins[], bool, std::string);
    void scan_couplings(std::string, std::pair <TH1F*, vector<TH1F *>> );
    void make_summary_plot(vector <TGraphErrors*>);

    private:
        double calculate_test_statistic(TH1F*, TH1F*);

};


Fitter* f_EFT;

vector <TGraphErrors*> scans;
vector <std::string> obs_names;
