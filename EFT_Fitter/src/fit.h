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



using namespace std;



class Fitter {
    
    public:
        Fitter(string, string);
        vector <TH1F*> initialise(double);
        void scan_couplings(vector<TH1F *>,TH1F * );

    private:

        double calculate_test_statistic(TH1F*, TH1F*);

};


