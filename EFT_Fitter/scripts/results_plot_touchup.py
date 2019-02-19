#touch up control plots
from ROOT import *
import os
import subprocess

gStyle.SetOptStat(00000)

directory_base = "files/Jan18"
target_directory_base = "/Users/keaveney/Desktop/testplots/"
#target_directory_base = "/Users/keaveney/docs/TOP-17-014/myDir/papers/TOP-17-014/trunk/fig_logo/diff_results/primary/"

xsec_types = ["absolute", "normalised"]
xsec_levels = ["parton", "particle"]

legend_moves = ["DiffXS_HypTTBarpT_source.root", "DiffXS_HypLLBarMass_source.root", "DiffXS_HypLeptonpT_source.root", "DiffXS_HypAntiLeptonpT_source.root","DiffXS_HypLeptonpTLead_source.root","DiffXS_HypLeptonpTNLead_source.root", "DiffXS_HypLLBarpT_source.root","DiffXS_HypBJetpTLead_source.root", "DiffXS_HypBJetpTNLead_source.root","DiffXS_HypBBBarpT_source.root","DiffXS_HypBBBarMass_source.root", "DiffXS_HypTTBarMass_source.root","DiffXS_HypTopRapidityLead_source.root","DiffXS_HypTopRapidityNLead_source.root" ]

for xsec_type in xsec_types:
    for xsec_level in xsec_levels:
        directory = directory_base + "/" + xsec_level +  "/" + xsec_type
        for filename in os.listdir(directory):
            if filename.startswith("DiffXS") & filename.endswith(".root"):
                filename_path = directory + "/" + filename
                file = TFile(filename_path)
                canvasname = "canvas"
                canvas = file.Get(canvasname)

                h_data = canvas.GetPrimitive("Nominalplot")
                h_data.SetTitle("")
                h_data_1 = canvas.GetPrimitive("Nominalplot_copy")
                
                print "filename  = " + str(filename)
                if(  (filename == "DiffXS_HypTTBarDeltaPhi_source.root") | (filename == "DiffXS_HypTTBarDeltaRapidity_source.root")  | (filename == "DiffXS_HypTTBarRapidity_source.root") | (filename == "DiffXS_HypTopRapidity_source.root") | (filename == "DiffXS_HypAntiTopRapidity_source.root") | (filename == "DiffXS_HypTopRapidityLead_source.root")  | (filename == "DiffXS_HypTopRapidityNLead_source.root") | (filename == "DiffXS_HypBBBarpT_source.root")   | (filename == "DiffXS_HypBBBarMass_source.root") ):
                    canvas.SetLogy()
                if( ((filename == "DiffXS_HypTopRapidity_source.root") | (filename == "DiffXS_HypAntiTopRapidity_source.root"))):
                    if (xsec_level == "parton"):
                        if (xsec_type == "normalised"):
                            h_data.GetYaxis().SetRangeUser(0.051, 0.52)
                        if (xsec_type == "absolute"):
                            h_data.GetYaxis().SetRangeUser(41.0, 390)
                    elif (xsec_level == "particle"):
                        if (xsec_type == "absolute"):
                            h_data.GetYaxis().SetRangeUser(0.21, 7.8)
                        if (xsec_type == "normalised"):
                            h_data.GetYaxis().SetRangeUser(0.018, 0.77)
                if(filename == "DiffXS_HypTTBarDeltaRapidity_source.root"):
                    if (xsec_level == "parton"):
                        if (xsec_type == "normalised"):
                            h_data.GetYaxis().SetRangeUser(0.001, 0.52)
                        if (xsec_type == "absolute"):
                            h_data.GetYaxis().SetRangeUser(1.0, 460)
                    elif (xsec_level == "particle"):
                        if (xsec_type == "absolute"):
                            h_data.GetYaxis().SetRangeUser(0.1, 7.2)
                        if (xsec_type == "normalised"):
                            h_data.GetYaxis().SetRangeUser(0.001, 0.58)
                if( (xsec_type == "absolute") & ((filename == "DiffXS_HypLeptonEtaLead_source.root"))):
                    h_data.GetYaxis().SetRangeUser(0.1, 4.7)
                if (filename == "DiffXS_HypTTBarDeltaRapidity_source.root"):
                    if( (xsec_level == "parton") &(xsec_type == "absolute")):
                        h_data.GetYaxis().SetRangeUser(20.0, 550.0)
                    if((xsec_level == "parton") & (xsec_type == "normalised")):
                        h_data.GetYaxis().SetRangeUser(0.02, 0.9)
                    if((xsec_level == "particle") &(xsec_type == "absolute")):
                        h_data.GetYaxis().SetRangeUser(0.12, 10.5)
                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                        h_data.GetYaxis().SetRangeUser(0.015, 0.9)
                if (filename == "DiffXS_HypTTBarDeltaPhi_source.root"):
                    if( (xsec_level == "parton") &(xsec_type == "absolute")):
                        h_data.GetYaxis().SetRangeUser(30.0, 3500.0)
                    if((xsec_level == "parton") & (xsec_type == "normalised")):
                        h_data.GetYaxis().SetRangeUser(0.04, 5.0)
                    if((xsec_level == "particle") &(xsec_type == "absolute")):
                        h_data.GetYaxis().SetRangeUser(0.4, 50.0)
                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                        h_data.GetYaxis().SetRangeUser(0.04, 5.0)
                if (filename == "DiffXS_HypTTBarRapidity_source.root"):
                    if( (xsec_level == "parton") &(xsec_type == "absolute")):
                        h_data.GetYaxis().SetRangeUser(30.0, 460.0)
                    if((xsec_level == "parton") & (xsec_type == "normalised")):
                        h_data.GetYaxis().SetRangeUser(0.04, 0.6)
                    if((xsec_level == "particle") &(xsec_type == "absolute")):
                        h_data.GetYaxis().SetRangeUser(0.08, 12.0)
                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                        h_data.GetYaxis().SetRangeUser(0.0065, 0.95)
                if (filename == "DiffXS_HypTTBarMass_source.root"):
                    if( (xsec_level == "parton") & (xsec_type == "absolute")):
                        h_data.GetYaxis().SetRangeUser(0.0003, 19.0)
                    if((xsec_level == "parton") & (xsec_type == "normalised")):
                        h_data.GetYaxis().SetRangeUser(0.0000004, 0.03)
                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                        h_data.GetYaxis().SetRangeUser(0.000005, 0.18)
                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                        h_data.GetYaxis().SetRangeUser(0.0000004, 0.02)
                if (filename == "DiffXS_HypBBBarpT_source.root"):
                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                        print "BBAR BBAR BBAR1"
                        h_data.GetYaxis().SetRangeUser(0.0006, 0.385)
                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                        print "BBAR BBAR BBAR2"
                        h_data.GetYaxis().SetRangeUser(0.00006, 0.0385)
                if (filename == "DiffXS_HypLeptonEta_source.root"):
                    if( (xsec_type == "absolute")):
                        h_data.GetYaxis().SetRangeUser(0.31, 4.4)
                if (filename == "DiffXS_HypLLBarDEta_source.root"):
                    if( (xsec_type == "absolute")):
                        h_data.GetYaxis().SetRangeUser(0.1, 5.9)
                    if (xsec_type == "normalised"):
                        h_data.GetYaxis().SetRangeUser(0.007, 0.56)
                if (filename == "DiffXS_HypLLBarpT_source.root"):
                    if( (xsec_type == "absolute")):
                        h_data.GetYaxis().SetRangeUser(0.0008, 0.65)
                    if (xsec_type == "normalised"):
                        h_data.GetYaxis().SetRangeUser(0.00008, 0.055)
                if (filename == "DiffXS_HypLLBarMass_source.root"):
                    if((xsec_type == "absolute")):
                        h_data.GetYaxis().SetRangeUser(0.0008, 0.43)
                    if(xsec_type == "normalised"):
                        h_data.GetYaxis().SetRangeUser(0.00008, 0.025)
                if (filename == "DiffXS_HypBBBarMass_source.root"):
                    if((xsec_type == "absolute")):
                        h_data.GetYaxis().SetRangeUser(0.0026, 0.12)
                    if(xsec_type == "normalised"):
                        h_data.GetYaxis().SetRangeUser(0.00025, 0.015)
                        h_data.GetYaxis().SetTitleOffset(1.6)
                if ( (filename == "DiffXS_HypTopRapidityLead_source.root") | (filename == "DiffXS_HypTopRapidityNLead_source.root")):
                    if( (xsec_type == "absolute") & (xsec_level == "parton")) :
                       h_data.GetYaxis().SetRangeUser(45.0, 350.0)
                    if( (xsec_type == "normalised") & (xsec_level == "parton")) :
                        h_data.GetYaxis().SetRangeUser(0.06, 0.42)
                    if( (xsec_type == "absolute") & (xsec_level == "particle")) :
                        h_data.GetYaxis().SetRangeUser(0.32, 7.0)
                    if( (xsec_type == "normalised") & (xsec_level == "particle")) :
                        h_data.GetYaxis().SetRangeUser(0.032, 0.62)
                if (filename == "DiffXS_HypToppTTTRestFrame_source.root"):
                    if (xsec_type == "normalised"):
                        h_data.GetYaxis().SetTitle("#frac{1}{#sigma} #frac{d#sigma}{p^{t}_{T} (t#bar{t} RF)} [GeV^{-1}]")
                    elif(xsec_type == "absolute"):
                        h_data.GetYaxis().SetTitle("#frac{d#sigma}{p^{t}_{T} (t#bar{t} RF)} [pb/GeV]")

                if (  (filename == "DiffXS_HypLeptonEta_source.root") |  (filename == "DiffXS_HypAntiLeptonEta_source.root") |   (filename == "DiffXS_HypLeptonEtaLead_source.root") | (filename == "DiffXS_HypLeptonEtaNLead_source.root") ):
                    if( (xsec_type == "normalised") & (xsec_level == "particle")) :
                        h_data.GetYaxis().SetRangeUser(0.04, 0.38)

                for prim in canvas.GetListOfPrimitives():
                    if prim.InheritsFrom("TH1"):
                        xmin = prim.GetXaxis().GetXmin()
                        xmax = prim.GetXaxis().GetXmax()
                        prim.SetMarkerColorAlpha(kWhite, 0.0)
                        prim.SetLineColorAlpha(kWhite, 0.0)
                    if prim.GetName() == "POWHEGplot":
                        prim.SetLineColor(kGreen+1)
                    if prim.GetName() == "MCATNLOplot":
                        prim.SetLineColor(kViolet+1)
                    if prim.GetName() == "Nominalplot":
                        prim.SetLineColor(kRed+1)
                    if prim.ClassName() == "TPad":
                        print "padName  = " + str(prim.GetName())
                        #      if ((filename == "DiffXS_HypTopRapidity_source.root") | (filename == "DiffXS_HypAntiTopRapidity_source.root")  | (filename == "DiffXS_HypTopRapidityLead_source.root")  | (filename == "DiffXS_HypTopRapidityNLead_source.root")  | (filename == "DiffXS_HypTTBarDeltaRapidity_source.root")):
                        canvas.cd(2)
                        l = TLine(xmin,1.0,xmax,1.0)
                        l.Draw()
                        for prim_pad in prim.GetListOfPrimitives():
                            if ((filename == "DiffXS_HypTopRapidity_source.root") |  (filename == "DiffXS_HypAntiTopRapidity_source.root") |  (filename == "DiffXS_HypTopRapidityLead_source.root") |  (filename == "DiffXS_HypTopRapidityNLead_source.root") | (filename  == "DiffXS_HypTTBarDeltaRapidity_source.root") | (filename == "DiffXS_HypTTBarDeltaPhi_source.root") ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if(xsec_type == "absolute"):
                                        prim_pad.GetYaxis().SetRangeUser(0.84, 1.35)
                                    if(xsec_type == "normalised"):
                                        prim_pad.GetYaxis().SetRangeUser(0.89, 1.29)
                            if ((filename == "DiffXS_HypTTBarpT_source.root")  ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if(xsec_type == "absolute"):
                                        prim_pad.GetYaxis().SetRangeUser(0.74, 1.64)
                                    if(xsec_type == "normalised"):
                                        prim_pad.GetYaxis().SetRangeUser(0.74, 1.64)
                            if ((filename == "DiffXS_HypTTBarRapidity_source.root")  ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if( (xsec_level == "parton") &(xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.81, 1.33)
                                    if( (xsec_level == "parton") &(xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.86, 1.29)
                                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.81, 1.33)
                                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.86, 1.29)
                            if ((filename == "DiffXS_HypTTBarMass_source.root")  ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if( (xsec_level == "parton") &(xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.71, 1.43)
                                    if( (xsec_level == "parton") &(xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.71, 1.39)
                                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.71, 1.43)
                                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.71, 1.39)
                            if ( (filename == "DiffXS_HypLeptonpT_source.root")  |  (filename == "DiffXS_HypAntiLeptonpT_source.root") | (filename == "DiffXS_HypLeptonpTLead_source.root")  | (filename == "DiffXS_HypLeptonpTNLead_source.root") ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.88, 1.33)
                                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.88, 1.33)
                            if ( (filename == "DiffXS_HypLeptonpTNLead_source.root") ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.89, 1.56)
                                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.89, 1.56)
                            if ( (filename == "DiffXS_HypLeptonEta_source.root") | (filename == "DiffXS_HypAntiLeptonEta_source.root") | (filename == "DiffXS_HypLeptonEtaLead_source.root") | (filename == "DiffXS_HypLeptonEtaNLead_source.root") ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.89, 1.32)
                                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.95, 1.16)
                            if ( (filename == "DiffXS_HypLLBarpT_source.root") ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.87, 1.33)
                                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.91, 1.17)
                            if ( (filename == "DiffXS_HypLLBarMass_source.root") ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.88, 1.33)
                                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.945, 1.16)
                            if ( (filename == "DiffXS_HypLLBarDPhi_source.root") ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.92, 1.22)
                                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.977, 1.08)
                            if ( (filename == "DiffXS_HypLLBarDEta_source.root") ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.91, 1.31)
                                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.961, 1.09)
                            if ( (filename == "DiffXS_HypBJetpTLead_source.root") ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.87, 1.29)
                                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.92, 1.27)
                            if ( (filename == "DiffXS_HypBJetpTNLead_source.root") ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.83, 1.37)
                                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.84, 1.37)
                            if ( (filename == "DiffXS_HypBJetEtaLead_source.root") ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.87, 1.27)
                                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.92, 1.14)
                            if ( (filename == "DiffXS_HypBJetEtaNLead_source.root") ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.86, 1.28)
                                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.91, 1.16)
                            if ( (filename == "DiffXS_HypBBBarpT_source.root") ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.83, 1.31)
                                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.86, 1.28)
                            if ( (filename == "DiffXS_HypBBBarpT_source.root") ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.87, 1.24)
                                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.86, 1.24)
                            if ( (filename == "DiffXS_HypBBBarMass_source.root") ):
                                if (prim_pad.InheritsFrom("TH1")):
                                    if((xsec_level == "particle") & (xsec_type == "absolute")):
                                        prim_pad.GetYaxis().SetRangeUser(0.86, 1.215)
                                    if((xsec_level == "particle") & (xsec_type == "normalised")):
                                        prim_pad.GetYaxis().SetRangeUser(0.87, 1.22)
                            if (filename == "DiffXS_HypToppTTTRestFrame_source.root"):
                                    if (prim_pad.InheritsFrom("TH1")):
                                        prim_pad.GetXaxis().SetTitle("p^{t}_{T} (t#bar{t} RF) [GeV]")

                            if (prim_pad.InheritsFrom("TF1")):
                                print " Setting Alpha White"
                                print "********"
                                prim_pad.SetLineColorAlpha(kWhite, 0.0)

                            if prim_pad.ClassName() == "TLegend":
                                prim_pad.SetFillStyle(0)
                                prim_pad.GetListOfPrimitives()[0].SetLabel("Stat #oplus Syst")
                                prim_pad.GetListOfPrimitives()[1].SetLabel("Stat")
                                prim_pad.SetX1NDC(0.23)
                                prim_pad.SetY1NDC(0.66)
                                if ((filename == "DiffXS_HypTopRapidity_source.root") | (filename == "DiffXS_HypAntiTopRapidity_source.root") ):
                                    if(xsec_type == "normalised"):
                                        prim_pad.SetX1NDC(0.23)
                                        prim_pad.SetY1NDC(0.71)
                                if ((filename == "DiffXS_HypTTBarRapidity_source.root") ):
                                    if( (xsec_level == "parton") &(xsec_type == "absolute")):
                                        prim_pad.SetX1NDC(0.45)
                                        prim_pad.SetX2NDC(0.7)
                                        prim_pad.SetY1NDC(0.72)
                                        prim_pad.SetY2NDC(0.95)
                                    if( (xsec_level == "parton") &(xsec_type == "normalised")):
                                        prim_pad.SetX1NDC(0.45)
                                        prim_pad.SetX2NDC(0.7)
                                        prim_pad.SetY1NDC(0.73)
                                        prim_pad.SetY2NDC(0.96)
                                    if( (xsec_level == "particle") &(xsec_type == "absolute")):
                                        prim_pad.SetX1NDC(0.45)
                                        prim_pad.SetX2NDC(0.7)
                                        prim_pad.SetY1NDC(0.72)
                                        prim_pad.SetY2NDC(0.95)
                                    if( (xsec_level == "particle") &(xsec_type == "normalised")):
                                        prim_pad.SetX1NDC(0.45)
                                        prim_pad.SetX2NDC(0.7)
                                        prim_pad.SetY1NDC(0.73)
                                        prim_pad.SetY2NDC(0.96)
                                if ( filename == "DiffXS_HypLeptonpTNLead_source.root" ):
                                        prim_pad.SetX1NDC(0.21)
                                        prim_pad.SetX2NDC(0.64)
                                if ( (filename == "DiffXS_HypBJetpTLead_source.root")):
                                    if( (xsec_level == "particle") &(xsec_type == "absolute")):
                                        prim_pad.SetX1NDC(0.21)
                                        prim_pad.SetX2NDC(0.39)
                                        prim_pad.SetY1NDC(0.73)
                                        prim_pad.SetY2NDC(0.96)
                                    if( (xsec_level == "particle") &(xsec_type == "normalised")):
                                        prim_pad.SetX1NDC(0.21)
                                        prim_pad.SetX2NDC(0.39)
                                        prim_pad.SetY1NDC(0.74)
                                        prim_pad.SetY2NDC(0.97)
                                if ( (filename == "DiffXS_HypBJetpTNLead_source.root")):
                                    if( (xsec_level == "particle") &(xsec_type == "absolute")):
                                        prim_pad.SetX1NDC(0.21)
                                        prim_pad.SetX2NDC(0.39)
                                        prim_pad.SetY1NDC(0.68)
                                        prim_pad.SetY2NDC(0.91)
                                    if( (xsec_level == "particle") &(xsec_type == "normalised")):
                                        prim_pad.SetX1NDC(0.21)
                                        prim_pad.SetX2NDC(0.39)
                                        prim_pad.SetY1NDC(0.68)
                                        prim_pad.SetY2NDC(0.91)
                                if ( (filename == "DiffXS_HypBJetEtaLead_source.root")):
                                    if( (xsec_level == "particle") &(xsec_type == "absolute")):
                                        prim_pad.SetX1NDC(0.21)
                                        prim_pad.SetX2NDC(0.39)
                                        prim_pad.SetY1NDC(0.75)
                                        prim_pad.SetY2NDC(0.98)
                                    if( (xsec_level == "particle") &(xsec_type == "normalised")):
                                        prim_pad.SetX1NDC(0.21)
                                        prim_pad.SetX2NDC(0.39)
                                        prim_pad.SetY1NDC(0.75)
                                        prim_pad.SetY2NDC(0.98)
                                if ( (filename == "DiffXS_HypBJetEtaNLead_source.root")):
                                    if( (xsec_level == "particle") &(xsec_type == "absolute")):
                                        prim_pad.SetX1NDC(0.21)
                                        prim_pad.SetX2NDC(0.39)
                                        prim_pad.SetY1NDC(0.74)
                                        prim_pad.SetY2NDC(0.97)
                                    if( (xsec_level == "particle") &(xsec_type == "normalised")):
                                        prim_pad.SetX1NDC(0.21)
                                        prim_pad.SetX2NDC(0.39)
                                        prim_pad.SetY1NDC(0.71)
                                        prim_pad.SetY2NDC(0.94)
                                if ( (filename == "DiffXS_HypBBBarpT_source.root")):
                                    if( (xsec_level == "particle") &(xsec_type == "absolute")):
                                        prim_pad.SetX1NDC(0.24)
                                        prim_pad.SetX2NDC(0.42)
                                        prim_pad.SetY1NDC(0.73)
                                        prim_pad.SetY2NDC(0.96)
                                    if( (xsec_level == "particle") &(xsec_type == "normalised")):
                                        prim_pad.SetX1NDC(0.24)
                                        prim_pad.SetX2NDC(0.42)
                                        prim_pad.SetY1NDC(0.75)
                                        prim_pad.SetY2NDC(0.98)
                                if ( (filename == "DiffXS_HypBBBarMass_source.root")):
                                    if( (xsec_level == "particle") &(xsec_type == "absolute")):
                                        prim_pad.SetX1NDC(0.24)
                                        prim_pad.SetX2NDC(0.42)
                                        prim_pad.SetY1NDC(0.755)
                                        prim_pad.SetY2NDC(0.975)
                                    if( (xsec_level == "particle") &(xsec_type == "normalised")):
                                        prim_pad.SetX1NDC(0.24)
                                        prim_pad.SetX2NDC(0.42)
                                        prim_pad.SetY1NDC(0.75)
                                        prim_pad.SetY2NDC(0.97)
                                if ( (filename == "DiffXS_HypJetMultpt30_source.root")):
                                    if( (xsec_level == "particle") &(xsec_type == "absolute")):
                                        prim_pad.SetX1NDC(0.24)
                                        prim_pad.SetX2NDC(0.42)
                                        prim_pad.SetY1NDC(0.65)
                                        prim_pad.SetY2NDC(0.87)
                                    if( (xsec_level == "particle") &(xsec_type == "normalised")):
                                        prim_pad.SetX1NDC(0.24)
                                        prim_pad.SetX2NDC(0.42)
                                        prim_pad.SetY1NDC(0.65)
                                        prim_pad.SetY2NDC(0.87)
                    if prim.ClassName() == "TPaveText":
                        if prim.GetLineWith("").GetTitle() == "CMS":
                            prim.SetX1NDC(0.17)
                        elif prim.GetLineWith("").GetTitle() == "35.9 fb^{-1} (13 TeV)":
                            prim.SetX1NDC(0.87)
                        elif prim.GetLineWith("").GetTitle() == "Dilepton: parton":
                            prim.GetLineWith("").SetTitle("Dilepton, parton level")
                            prim.SetY1NDC(0.84)
                            prim.SetX1NDC(0.232)
                        elif prim.GetLineWith("").GetTitle() == "Dilepton: particle":
                            prim.GetLineWith("").SetTitle("Dilepton, particle level")
                            prim.SetY1NDC(0.84)
                            prim.SetX1NDC(0.232)
                    if (prim.ClassName() == "TLegend"):
                        prim.SetTextSize(0.027)
                        print "*************************** tlegend = = = = =    " + str(prim.GetListOfPrimitives()[0].GetLabel())
                        if (prim.GetListOfPrimitives()[0].GetLabel() == "p_{T}^{jet}> 30 GeV, |#eta^{jet}| < 2.4"):
                            prim.SetTextSize(0.03)
                        if ((prim.GetListOfPrimitives()[0].GetLabel() == "Data") & (filename in legend_moves)):
                            prim.SetX1NDC(0.41)
                            prim.SetX2NDC(0.81)
                        if ((prim.GetListOfPrimitives()[0].GetLabel() == "Data") & ((filename == "DiffXS_HypTopRapidityLead_source.root") | (filename == "DiffXS_HypTopRapidityNLead_source.root")) ):
                            if (xsec_type == "absolute"):
                                if (xsec_level == "parton"):
                                    prim.SetX1NDC(0.34)
                                    prim.SetX2NDC(0.77)
                                    prim.SetY1NDC(0.38)
                                    prim.SetY2NDC(0.56)
                                elif (xsec_level == "particle"):
                                    prim.SetX1NDC(0.42)
                                    prim.SetX2NDC(0.82)
                                    prim.SetY1NDC(0.43)
                                    prim.SetY2NDC(0.6)
                            elif (xsec_type == "normalised"):
                                if (xsec_level == "parton"):
                                    prim.SetX1NDC(0.37)
                                    prim.SetX2NDC(0.8)
                                    prim.SetY1NDC(0.38)
                                    prim.SetY2NDC(0.54)
                                elif (xsec_level == "particle"):
                                    prim.SetX1NDC(0.39)
                                    prim.SetX2NDC(0.82)
                                    prim.SetY1NDC(0.44)
                                    prim.SetY2NDC(0.61)
                        if ((prim.GetListOfPrimitives()[0].GetLabel() == "Data") & (filename == "DiffXS_HypTTBarDeltaRapidity_source.root")):
                            if (xsec_type == "absolute"):
                                if (xsec_level == "parton"):
                                    prim.SetX1NDC(0.37)
                                    prim.SetX2NDC(0.81)
                                    prim.SetY1NDC(0.4)
                                    prim.SetY2NDC(0.55)
                                elif (xsec_level == "particle"):
                                    prim.SetX1NDC(0.39)
                                    prim.SetX2NDC(0.69)
                                    prim.SetY1NDC(0.4)
                                    prim.SetY2NDC(0.57)
                            elif (xsec_type == "normalised"):
                                if (xsec_level == "parton"):
                                    prim.SetX1NDC(0.37)
                                    prim.SetX2NDC(0.8)
                                    prim.SetY1NDC(0.39)
                                    prim.SetY2NDC(0.55)
                                elif (xsec_level == "particle"):
                                    prim.SetX1NDC(0.39)
                                    prim.SetX2NDC(0.69)
                                    prim.SetY1NDC(0.38)
                                    prim.SetY2NDC(0.55)
                        if ((prim.GetListOfPrimitives()[0].GetLabel() == "Data") & (filename == "DiffXS_HypTTBarDeltaPhi_source.root")):
                            if (xsec_type == "absolute"):
                                if (xsec_level == "parton"):
                                    prim.SetX1NDC(0.22)
                                    prim.SetX2NDC(0.65)
                                    prim.SetY1NDC(0.6)
                                    prim.SetY2NDC(0.82)
                                elif (xsec_level == "particle"):
                                    prim.SetX1NDC(0.22)
                                    prim.SetX2NDC(0.65)
                                    prim.SetY1NDC(0.6)
                                    prim.SetY2NDC(0.82)
                            elif (xsec_type == "normalised"):
                                if (xsec_level == "parton"):
                                    prim.SetX1NDC(0.22)
                                    prim.SetX2NDC(0.65)
                                    prim.SetY1NDC(0.6)
                                    prim.SetY2NDC(0.82)
                                elif (xsec_level == "particle"):
                                    prim.SetX1NDC(0.22)
                                    prim.SetX2NDC(0.65)
                                    prim.SetY1NDC(0.6)
                                    prim.SetY2NDC(0.82)
                        if ((prim.GetListOfPrimitives()[0].GetLabel() == "Data") & (filename == "DiffXS_HypTTBarRapidity_source.root")):
                            if (xsec_type == "absolute"):
                                if (xsec_level == "parton"):
                                    prim.SetX1NDC(0.36)
                                    prim.SetX2NDC(0.78)
                                    prim.SetY1NDC(0.39)
                                    prim.SetY2NDC(0.57)
                                elif (xsec_level == "particle"):
                                    prim.SetX1NDC(0.38)
                                    prim.SetX2NDC(0.78)
                                    prim.SetY1NDC(0.39)
                                    prim.SetY2NDC(0.58)
                            elif (xsec_type == "normalised"):
                                if (xsec_level == "parton"):
                                    prim.SetX1NDC(0.35)
                                    prim.SetX2NDC(0.79)
                                    prim.SetY1NDC(0.39)
                                    prim.SetY2NDC(0.59)
                                elif (xsec_level == "particle"):
                                    prim.SetX1NDC(0.36)
                                    prim.SetX2NDC(0.79)
                                    prim.SetY1NDC(0.39)
                                    prim.SetY2NDC(0.58)
                        if ((prim.GetListOfPrimitives()[0].GetLabel() == "Data") & (filename == "DiffXS_HypBBBarpT_source.root")):
                            if (xsec_type == "absolute"):
                                if (xsec_level == "particle"):
                                    prim.SetX1NDC(0.4)
                                    prim.SetX2NDC(0.8)
                                    prim.SetY1NDC(0.47)
                                    prim.SetY2NDC(0.66)
                            elif (xsec_type == "normalised"):
                                if (xsec_level == "particle"):
                                    prim.SetX1NDC(0.4)
                                    prim.SetX2NDC(0.8)
                                    prim.SetY1NDC(0.47)
                                    prim.SetY2NDC(0.66)
                        if ((prim.GetListOfPrimitives()[0].GetLabel() == "Data") & (filename == "DiffXS_HypBBBarMass_source.root")):
                            if (xsec_type == "absolute"):
                                if (xsec_level == "particle"):
                                    prim.SetX1NDC(0.4)
                                    prim.SetX2NDC(0.8)
                                    prim.SetY1NDC(0.49)
                                    prim.SetY2NDC(0.68)
                            elif (xsec_type == "normalised"):
                                if (xsec_level == "particle"):
                                    prim.SetX1NDC(0.4)
                                    prim.SetX2NDC(0.8)
                                    prim.SetY1NDC(0.49)
                                    prim.SetY2NDC(0.68)

                        if ((prim.GetListOfPrimitives()[0].GetLabel() == "Data") & ((filename == "DiffXS_HypTopRapidityLead_source.root") | (filename == "DiffXS_HypTopRapidityNLead_source.root")  | (filename == "DiffXS_HypTopRapidity_source.root")  | (filename == "DiffXS_HypAntiTopRapidity_source.root") )):
                            if (xsec_type == "absolute"):
                                if (xsec_level == "parton"):
                                    prim.SetX1NDC(0.38)
                                    prim.SetX2NDC(0.8)
                                    prim.SetY1NDC(0.38)
                                    prim.SetY2NDC(0.57)
                                elif (xsec_level == "particle"):
                                    prim.SetX1NDC(0.35)
                                    prim.SetX2NDC(0.77)
                                    prim.SetY1NDC(0.38)
                                    prim.SetY2NDC(0.57)
                            elif (xsec_type == "normalised"):
                                if (xsec_level == "parton"):
                                    prim.SetX1NDC(0.38)
                                    prim.SetX2NDC(0.8)
                                    prim.SetY1NDC(0.38)
                                    prim.SetY2NDC(0.57)
                                elif (xsec_level == "particle"):
                                    prim.SetX1NDC(0.36)
                                    prim.SetX2NDC(0.77)
                                    prim.SetY1NDC(0.38)
                                    prim.SetY2NDC(0.57)


                        for prim_leg in prim.GetListOfPrimitives():
                            if prim_leg.GetLabel() == "Data":
                                prim_leg.SetOption("E1p")
                            elif prim_leg.GetLabel() == "Powheg v2+Pythia8":
                                prim_leg.SetLabel("POWHEGV2 + PYTHIA8")
                            elif prim_leg.GetLabel() == "Powheg v2+Herwig++":
                                prim_leg.SetLabel("POWHEGV2 + HERWIG++")
                            elif prim_leg.GetLabel() == "MG5_aMC@NLO+Pythia8 [FxFx]":
                                prim_leg.SetLabel("MG5_aMC@NLO + PYTHIA8 [FxFx]")

                plotname = filename.split("_")
                if len(plotname) < 3:
                    plotname_base = plotname[0]
                else:
                    plotname_base = plotname[0] + "_" + plotname[1]
                if xsec_type == "normalised":
                    target_directory = target_directory_base + "normalized" + "/" + xsec_level + "/"
                else:
                    target_directory = target_directory_base + xsec_type + "/" + xsec_level + "/"
                plotfilename = target_directory + plotname_base + "_rev3.pdf"
                canvas.SaveAs(plotfilename)



#dir_list = []
#dir_list.append(directory)
#insert_eps_text(dir_list)

#r_pad = canvas.GetPrimitive("rPad")
#f_ratio = r_pad.GetPrimitive("f")
#f2_ratio = r_pad.GetPrimitive("f2")
#f_ratio.SetLineColor(0)
#f2_ratio.SetLineColor(0)
#h_ratio = r_pad.GetPrimitive("ratio")
#h_ratio.SetYTitle("#frac{Data}{Pred.}")


def insert_eps_text(path_list):
    
    searchstring = "gsave  2268 2176 0 0 C 691.54 1986 t 0 r /Helvetica-Bold findfont 91.3174 sf 0 0 m (CMS) show NC gr"
    insertstring = "gsave  2268 2176 0 0 C 909.37 1993.76 t 0 r /Helvetica-Oblique findfont 79.9028 sf 0 0 m (Preliminary) show NC gr"
    
    for path in path_list:
        for filename in os.listdir(path):
            if filename.endswith(".eps"):
                targetfile = os.path.join(path, filename)
                print(targetfile)
                inputfile = open(targetfile, 'r').readlines()
                write_file = open(targetfile,'w')
                for line in inputfile:
                    write_file.write(line)
                    if searchstring in line:
                        write_file.write(insertstring + "\n")
                write_file.close()
                subprocess.call(['epstopdf', targetfile])
                continue
            else:
                continue
