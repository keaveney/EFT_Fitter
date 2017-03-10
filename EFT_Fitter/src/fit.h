#include <iostream>
#include <sstream>
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
        void scan_couplings(std::pair <TH1F*, vector<TH1F *>> );


    private:
        double calculate_test_statistic(TH1F*, TH1F*);

};


Fitter* f_EFT;
