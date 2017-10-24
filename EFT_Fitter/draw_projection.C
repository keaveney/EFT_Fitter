{
  //TF1 *fa1 = new TF1("fa1","sin(x)/x",0,10);

  //sigma_ttbb_eff = sigma_ttbb * frac_boosted * acc_boosted * BR_ll
  // = 14000 * 0.10 * 0.04

  //N_ttbb_boosted_ll_300 = 2520
  //N_ttH_boosted_ll_300 = 60

  // eff xsec times acc times eff ttbb boosted ll = 8.4
  // eff xsec times acc times eff ttH boosted ll = 0.2
  

//  TF1 *f_ll_boosted = new TF1("f_ll_boosted"," (0.2*x) / ( pow( 8.4*x  +   pow(0.2*8.4*x, 2.0), 0.5 )) ",0,300);
//  TF1 *f_ljets_boosted = new TF1("f_ljets_boosted"," (1.32*x) / ( pow( 56.0 *x  +   pow(0.3* 56.0 *x, 2.0), 0.5 ))  ",0,300);

    //0.1 tth ll boosted per fb
    //0.7 tth ljets boosted per fb
    //5.0 ttbb ll boosted per fb
    //35.0 ttbb ljets boosted per fb

    //0.8 tth ll resolved per fb
    //5.6 tth ljets resolved per fb
    //40.0 ttbb ll resolved per fb
    //280.0 ttbb ljets resolved per fb
    
    gStyle->SetOptStat(0);

 // Frozen projection (systematics fixed to 50%)
//    TF1 *f_ll_resolved = new TF1("f_ll_resolved"," (0.8*x) / ( pow( 40.0*x  +   pow(0.5 * 40.0 * x, 2.0), 0.5 )) ",0,300);
//    TF1 *f_ll_boosted = new TF1("f_ll_boosted"," (0.1*x) / ( pow( 5.0*x  +   pow(0.5 * 5.0 * x, 2.0), 0.5 )) ",0,300);
//    TF1 *f_ljets_resolved =  new TF1("f_ljets_resolved", "(5.6*x) / ( pow( 280.0*x +pow(0.5 * 280.0*x, 2.0), 0.5 ))   " ,0, 300 );
//    TF1 *f_ljets_boosted =   new TF1("f_ljets_boosted", "(0.7*x) / ( pow( 35.0*x +pow(0.5 * 35.0*x, 2.0), 0.5 ))   " ,0, 300 );

   // TF1 *f_ll_resolved    =  new TF1("f_ll_resolved",    " (1.33*x) / ( pow( 1.33*x +  pow(0.0 * 15.0 * x, 2.0), 0.5 )) ",0,300);
   // TF1 *f_ll_boosted     =  new TF1("f_ll_boosted",     " (0.13*x) / ( pow( 0.13*x +   pow(0.0 * 1.5 * x, 2.0), 0.5 )) ",0,300);
   // TF1 *f_ljets_resolved =  new TF1("f_ljets_resolved", " (6.3*x) / ( pow( 6.3*x +    pow(0.0 * 68.0*x, 2.0), 0.5 ))   " ,0, 300 );
   // TF1 *f_ljets_boosted  =  new TF1("f_ljets_boosted",  " (0.81*x) / ( pow( 0.81*x +   pow(0.0 * 6.3*x, 2.0), 0.5 ))   " ,0, 300 );

    TF1 *f_ljets_resolved =  new TF1("f_ljets_resolved", " (6.3*x) / ( pow( 167.0*x +    pow(0.5 * 167*x, 2.0), 0.5 ))   " ,0, 300 );
    TF1 *f_ll_resolved    =  new TF1("f_ll_resolved",    " (1.33*x) / ( pow( 42.2*x +  pow(0.5 * 42.2 * x, 2.0), 0.5 )) ",0,300);

    TF1 *f_ll_boosted     =  new TF1("f_ll_boosted",     " (0.1*x) / ( pow( 1.75*x +   pow(0.5 * 1.75 * x, 2.0), 0.5 )) ",0,300);
    TF1 *f_ljets_boosted  =  new TF1("f_ljets_boosted",  " (0.81*x) / ( pow( 14.4*x +   pow(0.5 * 14.4*x, 2.0), 0.5 ))   " ,0, 300 );
    
    TH1F * h_dummy = new TH1F("","",  1, 1.0, 300);
    
    h_dummy->Draw();
    f_ljets_boosted->Draw("same");
    f_ll_resolved->Draw("same");
    f_ljets_resolved->Draw("same");
    f_ll_boosted->Draw("same");

    
    f_ll_boosted->SetMinimum(0.01);
    f_ll_resolved->SetMinimum(0.01);
    f_ljets_boosted->SetMinimum(0.01);
    f_ljets_resolved->SetMinimum(0.01);


    h_dummy->SetXTitle("Integrated Luminosity (fb^{-1})");
    h_dummy->SetYTitle("N_{ttH} /  #surd(  (#sigma^{stat.}_{ttbb})^{2}  + (#sigma^{sys.}_{ttbb})^{2} )");
    f_ll_boosted->SetTitle(" ");


    f_ljets_boosted->SetLineColor(kBlue);
    f_ll_boosted->SetLineColor(kBlue);

    
    f_ljets_boosted->SetLineStyle(7);
    f_ljets_resolved->SetLineStyle(9);
    f_ll_boosted->SetLineStyle(10);
    f_ll_resolved->SetLineStyle(3);

    leg = new TLegend(0.1,0.7,0.48,0.9);
    leg->AddEntry(f_ljets_resolved,"l+jets (resolved)","l");
    leg->AddEntry(f_ll_resolved,"ll (resolved)","l");
    leg->AddEntry(f_ljets_boosted,"l+jets (boosted) ","l");
    leg->AddEntry(f_ll_boosted,"ll  (boosted)","l");
    leg->Draw();

    c1->SetLogx();

}
