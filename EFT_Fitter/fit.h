#include <iostream>
#include <sstream>
#include "TH1F.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"


using namespace std;



class Fitter {
    

    public:
        vector <TH1F*> initialise();
        void scan_couplings(vector<TH1F *>,TH1F * );

    private:

        double calculate_test_statistic(TH1F*, TH1F*);

};


