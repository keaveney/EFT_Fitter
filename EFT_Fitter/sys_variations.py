#sys variation plots
import ROOT
from ROOT import *
import os
import subprocess
from array import array

gStyle.SetOptStat(00000)

mergeExp = {"absolute":{
"ToppT":["JER", "TOT_TRIG", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"] ,
"AntiToppT":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"ToppTLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"ToppTNLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"ToppTTTRestFrame":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"TopRapidity":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"AntiTopRapidity":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"TopRapidityLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"TopRapidityNLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"TTBarpT":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY"],
"TTBarRapidity":["JER", "TOT_TRIG", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"TTBarMass":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"TTBarDeltaRapidity":["JER", "TOT_TRIG", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"TTBarDeltaPhi":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"LeptonpT":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"AntiLeptonpT":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"LeptonpTLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"LeptonpTNLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"LeptonEta":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"AntiLeptonEta":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"LeptonEtaLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"LeptonEtaNLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"LLBarMass":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"LLBarpT":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"LLBarDPhi":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"LLBarDEta":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"JetMultpt30":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"BJetpTLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"BJetpTNLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"BJetEtaLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"BJetEtaNLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"BBBarpT":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"BBBarMass":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"]
},

"normalised":{
"ToppT":["JER", "TOT_TRIG", "LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"] ,
"AntiToppT":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"ToppTLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"ToppTNLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"ToppTTTRestFrame":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"TopRapidity":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"AntiTopRapidity":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"TopRapidityLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"TopRapidityNLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"TTBarpT":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY"],
"TTBarRapidity":["JER", "TOT_TRIG", "LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"TTBarMass":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"TTBarDeltaRapidity":["JER", "TOT_TRIG", "LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"TTBarDeltaPhi":["JER", "TOT_TRIG", "LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU"],
"LeptonpT":["JER", "TOT_TRIG", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "UNCLUSTERED"],
"AntiLeptonpT":["JER", "TOT_TRIG", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU",  "UNCLUSTERED"],
"LeptonpTLead":["JER", "TOT_TRIG",  "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU",  "UNCLUSTERED"],
"LeptonpTNLead":["JER", "TOT_TRIG", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU",  "UNCLUSTERED"],
"LeptonEta":["JER", "TOT_TRIG", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"AntiLeptonEta":["JER", "TOT_TRIG", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"LeptonEtaLead":["JER", "TOT_TRIG", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"LeptonEtaNLead":["JER", "TOT_TRIG", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"LLBarMass":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "UNCLUSTERED"],
"LLBarpT":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "UNCLUSTERED"],
"LLBarDPhi":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "UNCLUSTERED"],
"LLBarDEta":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "UNCLUSTERED"],
"JetMultpt30":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"BJetpTLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "UNCLUSTERED"],
"BJetpTNLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "UNCLUSTERED"],
"BJetEtaLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "UNCLUSTERED"],
"BJetEtaNLead":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"BBBarpT":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU", "DY", "UNCLUSTERED"],
"BBBarMass":["JER", "TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "PU",  "UNCLUSTERED"]
}
}


directory_base = "/Users/keaveney/Downloads/sys_variation_files/nfs/dust/cms/user/savitsky/13TeV_Mor17/ForResults_8026p1_min/"
directory_tail = "/CMSSW_8_0_26_patch1/src/TopAnalysis/Configuration/analysis/diLeptonic/DIFFXS_PrimaryResults_TOP17014_backup/UnfoldingResults/"

directory_base_results = "/Users/keaveney/EFT_Fitter/EFT_Fitter/files/Nov1/"
directory_tail_results = "/results/"

#xsec_types = ["normalized"]
xsec_levels = ["parton", "pseudo"]
xsec_types = ["absolute", "normalized"]
#xsec_levels = ["pseudo"]

systematic_names =  ["JER", "Trigger","Lepton SF", "Kin Reco", "b tagging (b)", "b tagging (l)", "Drell-Yan","Pile up","Unclustered p^{miss}_{T}","ME scales","m_{t}", "PDF", "PDF #alpha_{S}", "h_{damp}", "Backgrounds", "PS ISR","PS FSR", "UE tune",  "Colour rec", "b fragmentation", "Br(b semileptonic)", "JESAbsoluteStat","JESAbsoluteScale","JESAbsoluteMPFBias","JESFragmentation","JESSinglePionECAL", "JESSinglePionHCAL", "JESFlavorQCD", "JESTimePtEta", "JESRelativeBal", "JESRelativeJEREC1", "JESRelativePtBB", "JESRelativePtEC1", "JESRelativeFSR", "JESRelativeStatFSR", "JESRelativeStatEC", "JESPileUpDataMC", "JESPileUpPtRef", "JESPileUpPtEC1", "JESPileUpPtBB"]

systematics =  ["JER","TOT_TRIG","LEPT", "KIN", "TOT_BTAG", "TOT_BTAG_LJET", "DY", "PU","UNCLUSTERED", "TOT_SCALE", "MASS", "PDF","PDF_ALPHAS","MATCH", "BG", "PSISRSCALE", "PSFSRSCALE", "UETUNE", "TOT_COLORREC", "TOT_BFRAG","BSEMILEP", "JESAbsoluteStat","JESAbsoluteScale","JESAbsoluteMPFBias","JESFragmentation","JESSinglePionECAL", "JESSinglePionHCAL", "JESFlavorQCD", "JESTimePtEta", "JESRelativeBal", "JESRelativeJEREC1", "JESRelativePtBB", "JESRelativePtEC1", "JESRelativeFSR", "JESRelativeStatFSR", "JESRelativeStatEC", "JESPileUpDataMC", "JESPileUpPtRef", "JESPileUpPtEC1", "JESPileUpPtBB" ]

variations = ["UP", "DOWN"]

observables = {
"parton":["ToppT","AntiToppT","ToppTLead", "ToppTNLead", "ToppTTTRestFrame", "TopRapidity","AntiTopRapidity", "TopRapidityLead", "TopRapidityNLead", "TTBarpT", "TTBarRapidity",  "TTBarMass",  "TTBarDeltaRapidity", "TTBarDeltaPhi"],

    "particle":["ToppT","AntiToppT","ToppTLead", "ToppTNLead", "ToppTTTRestFrame", "TopRapidity","AntiTopRapidity", "TopRapidityLead", "TopRapidityNLead", "TTBarpT",  "TTBarRapidity",  "TTBarMass",  "TTBarDeltaRapidity", "TTBarDeltaPhi", "LeptonpT", "AntiLeptonpT", "LeptonpTLead", "LeptonpTNLead", "LeptonEta" , "AntiLeptonEta", "LeptonEtaLead", "LeptonEtaNLead", "LLBarMass", "LLBarpT", "LLBarDPhi", "LLBarDEta", "JetMultpt30", "BJetpTLead", "BJetpTNLead", "BJetEtaLead", "BJetEtaNLead",    "BBBarpT", "BBBarMass"]
}

vars_root = {
"parton":["p_{T}^{t}","p_{T}^{#bar{t}}", "p_{T}^{t} (leading)",  "p_{T}^{t} (trailing)", "p_{T}^{t} (t#bar{t} RF)", "y_{t}","y_{#bar{t}}", "y_{t} (leading)", "y_{t} (trailing)", "p_{T}^{t#bar{t}}","y_{t#bar{t}}","m_{t#bar{t}}","#Delta|y|(t,#bar{t})", "#Delta#phi(t,#bar{t})"],

"particle":["p_{T}^{t}","p_{T}^{#bar{t}}", "p_{T}^{t} (leading)",  "p_{T}^{t} (trailing)", "p_{T}^{t} (t#bar{t} RF)", "y_{t}","y_{#bar{t}}", "y_{t} (leading)", "y_{t} (trailing)", "p_{T}^{t#bar{t}}","y_{t#bar{t}}","m_{t#bar{t}}","#Delta|y|(t,#bar{t})", "#Delta#phi(t,#bar{t})", "p_{T}^{l}", "p_{T}^{#bar{l}}", "p_{T}^{l} (leading)", "p_{T}^{l} (trailing)", "#eta_{l}", "#eta_{#bar{l}}", "#eta_{l} (leading)", "#eta_{l} (trailing)","m_{l#bar{l}}", "p_{T}^{l#bar{l}}", "#Delta#phi(l,#bar{l})", "#Delta#eta(l,#bar{l})",  "N_{jets}",   "p_{T}^{b} (leading)", "p_{T}^{b} (trailing)",  "#eta^{b} (leading)", "#eta^{b} (trailing)", "p_{T}^{b#bar{b}}", "m_{b#bar{b}}" ]}


units = {
"parton":["GeV","GeV", "GeV","GeV","GeV", "","","","","GeV","","GeV","",""],
"particle":["GeV","GeV", "GeV","GeV","GeV", "","","","","GeV","","GeV","","","GeV","GeV", "GeV","GeV","","","","","GeV","GeV", "","","","GeV","GeV","","","GeV","GeV"]
}


systematic_colors = {"JER":1,"Trigger":1,"Lepton SF":46,"Kin Reco":1,"b tagging (b)":1,"b tagging (l)":1,"Drell-Yan":29,"Pile up":1,"Unclustered p^{miss}_{T}":1,"ME scales":3,"m_{t}":4,"PDF":6,"PDF #alpha_{S}":6,"h_{damp}":2,"Backgrounds":5,"PS ISR":28,"PS FSR":28,"UE tune":7,"Colour rec":7,"b fragmentation":38,"Br(b semileptonic)":38,"JESAbsoluteStat":1,"JESAbsoluteScale":1,"JESAbsoluteMPFBias":1,"JESFragmentation":1,"JESSinglePionECAL":1, "JESSinglePionHCAL":1, "JESFlavorQCD":1, "JESTimePtEta":1, "JESRelativeBal":1, "JESRelativeJEREC1":1, "JESRelativePtBB":1, "JESRelativePtEC1":1, "JESRelativeFSR":1, "JESRelativeStatFSR":1, "JESRelativeStatEC":1, "JESPileUpDataMC":1, "JESPileUpPtRef":1, "JESPileUpPtEC1":1, "JESPileUpPtBB":1}

systematic_styles = {"JER":1,"Trigger":1,"Lepton SF":1,"Kin Reco":1,"b tagging (b)":1,"b tagging (l)":1,"Drell-Yan":1,"Pile up":1,"Unclustered p^{miss}_{T}":1,"ME scales":1,"m_{t}":1,"PDF":1,"PDF #alpha_{S}":7,"h_{damp}":1,"Backgrounds":1,"PS ISR":7,"PS FSR":1,"UE tune":1,"Colour rec":7,"b fragmentation":1,"Br(b semileptonic)":7,"JESAbsoluteStat":1,"JESAbsoluteScale":1,"JESAbsoluteMPFBias":1,"JESFragmentation":1,"JESSinglePionECAL":1, "JESSinglePionHCAL":1, "JESFlavorQCD":1, "JESTimePtEta":1, "JESRelativeBal":1, "JESRelativeJEREC1":1, "JESRelativePtBB":1, "JESRelativePtEC1":1, "JESRelativeFSR":1, "JESRelativeStatFSR":1, "JESRelativeStatEC":1, "JESPileUpDataMC":1, "JESPileUpPtRef":1, "JESPileUpPtEC1":1, "JESPileUpPtBB":1}

iobs = 0

def makeOutputFile(fileName):
    file = TFile(fileName, "RECREATE")
    return file

def processNominal(path, histName, minY, maxY):
    inputfile = open(path, 'r').readlines()
    bins = []
    for line in inputfile:
        bins.append(float(line.split( )[3]))
    bins.append(float(line.split( )[5]))
    bins = sorted(bins)
    binsArray = array('f',bins)
    hist = TH1F(histName, histName, len(bins)-1, binsArray)
    ibin = 0

    for line in inputfile:
        hist.SetBinContent(ibin+1, float(line.split( )[7]))
        ibin = ibin + 1

    hist.GetXaxis().SetNdivisions(10)
    hist.GetYaxis().SetRangeUser(minY,maxY)
    hist.SetLineWidth(0)

    return hist

def totalGraph(histTotSysUp, histTotSysDown):
    gr = TGraphAsymmErrors(histTotSysUp.GetNbinsX())
    for ibin in range(0,histTotSysUp.GetNbinsX()):
        gr.SetPoint(ibin,histTotSysUp.GetBinCenter(ibin+1),0.0)
        gr.SetPointEXlow(ibin,histTotSysUp.GetBinWidth(ibin+1)/2.0)
        gr.SetPointEXhigh(ibin,histTotSysUp.GetBinWidth(ibin+1)/2.0)
        gr.SetPointEYhigh(ibin,histTotSysDown.GetBinContent(ibin+1))
        gr.SetPointEYlow(ibin, histTotSysUp.GetBinContent(ibin+1))
        #print "in total graph  = " + str(histTotSysUp.GetBinContent(ibin+1))
            #if ((histTotSysUp.GetBinContent(ibin+1) > 1.0) & (histTotSysDown.GetBinContent(ibin+1) < 1.0)):
            # gr.SetPointEYhigh(ibin,histTotSysUp.GetBinContent(ibin+1)-1.0)
            # gr.SetPointEYlow(ibin,1.0 - histTotSysDown.GetBinContent(ibin+1))
            #elif ((histTotSysUp.GetBinContent(ibin+1) < 1.0) & (histTotSysDown.GetBinContent(ibin+1) > 1.0)):
            #gr.SetPointEYhigh(ibin,histTotSysDown.GetBinContent(ibin+1)-1.0)
            #gr.SetPointEYlow(ibin,1.0 - histTotSysUp.GetBinContent(ibin+1))
    return gr

def centerHistAtZero(hist):
    for iBin in range(0, hist.GetNbinsX()):
        origBin = hist.GetBinContent(iBin+1)
        newBin  = 100.0*(origBin - 1.0)
        hist.SetBinContent(iBin+1, newBin)
    return hist

def centerGraphAtZero(graph):
    n = graph.GetN()
    for iPoint in range(0, n):
        tmpX, tmpY = ROOT.Double(0), ROOT.Double(0)
        graph.GetPoint(iPoint,tmpX, tmpY)
        graph.SetPoint(iPoint, tmpX,0.0)
    return graph

def reNormalise(histNominal, histUp, histDown):
    histVisUp = histNominal.Clone()
    histVisDown = histNominal.Clone()
    histVisUp.Reset()
    histVisDown.Reset()

    for iBin in range(0, histNominal.GetNbinsX()):
        nomDiffXS = histNominal.GetBinContent(iBin+1)
        upDiffXS = (nomDiffXS)*(histUp.GetBinContent(iBin+1))
        downDiffXS = (nomDiffXS)*(histDown.GetBinContent(iBin+1))
        histVisUp.SetBinContent(iBin+1,upDiffXS)
        histVisDown.SetBinContent(iBin+1,downDiffXS)

    targetIntegral = histNominal.Integral("width")
    scaleUp = targetIntegral/histVisUp.Integral("width")
    scaleDown = targetIntegral/histVisDown.Integral("width")
    histVisUp.Scale(scaleUp)
    histVisDown.Scale(scaleDown)

    for iBin in range(0, histNominal.GetNbinsX()):
        nomDiffXS  =  histNominal.GetBinContent(iBin+1)
        upDiffXS   =  histVisUp.GetBinContent(iBin+1)
        downDiffXS =  histVisDown.GetBinContent(iBin+1)
        avSys = abs((upDiffXS - nomDiffXS)/(nomDiffXS*2.0))
        histUp.SetBinContent(iBin+1,1.0 + avSys)
        histDown.SetBinContent(iBin+1,1.0 -avSys)

    return (histUp, histDown)

def processSystematic(observable, xsecType, xsecLevel, systematic, histNominal):
    """
    * basic pre-processing on the raw variatons
    * returns a single histogram for each systematic with bin content represting the symmterised uncertainty
    Pre-processing - top mass variations are scaled by 1/3
               TODO: scale variaitons are 'enveloped' and FSR variations
                     are scaled by 1/ sqrt(2)
                     """
    varHists = []

    linkStr = "_"
    variations = [""]

    for variation in variations:
        if xsecType == "normalised":
            xsecType = "normalized"
        if xsecLevel == "particle":
            xsecLevel = "pseudo"
        path = directory_base + xsecType + "_" + xsecLevel + directory_tail + systematic + linkStr + variation + "/combinedUnfolded/Hyp" + observable + "Results.txt"
        #print "directory  = " + str(path)
        inputfile = open(path, 'r').readlines()
        bins = []
        for line in inputfile:
            bins.append(float(line.split( )[3]))
        bins.append(float(line.split( )[5]))
        bins = sorted(bins)
        binsArray = array('f',bins)
        histNameUp = systematic + "_UP"  
        histNameDown = systematic + "_DOWN" 
        histUp = TH1F(histNameUp, histNameUp, len(bins)-1, binsArray)
        histDown = TH1F(histNameDown, histNameDown, len(bins)-1, binsArray)
        histUpFinal = TH1F("", "", len(bins)-1, binsArray)
        histDownFinal = TH1F("", "", len(bins)-1, binsArray)
        
        ibin = 0

        for line in inputfile:
            nomBin = histNominal.GetBinContent(ibin+1)
            nomBinCenter = histNominal.GetBinCenter(ibin+1)
            unc = float(line.split( )[7])
#            if systematic == "MASS":
#                unc = unc/(3.0)
#            if systematic == "PSFSRSCALE":
#                unc = unc/(sqrt(2.0))

            histUp.SetBinContent(ibin+1, 1.0 + unc)
            histDown.SetBinContent(ibin+1,1.0 - unc)
            ibin = ibin + 1 

        histUpVis = histUp.Clone()
        histDownVis = histDown.Clone()
        histUpFinal = histUp.Clone()
        histDownFinal = histDown.Clone()

    if systematic == "PDF":
        histUpFinal, histDownFinal = reNormalise(histNominal, histUpVis, histDownVis)

    return (histUpFinal, histDownFinal)

def make_systematics_canvas():
    canvas = TCanvas()
    canvas.SetLeftMargin(0.1);
    canvas.SetRightMargin(0.23);
    return canvas

def make_sys_sum(fileName, systematics, variations):
    """
    * takes the pre-processed variation histograms and makes the final
    summed histograms
    - symmeterises effects in bins 
    """
    file = TFile(fileName)
    histTotSysUp = file.Get("Nominal").Clone()
    histTotSysDown = file.Get("Nominal").Clone()
    histTotSysUp.Reset()
    histTotSysDown.Reset()
    histTotSysUp.SetDirectory(0)
    histTotSysDown.SetDirectory(0)

    #print "in make_sys_sum, filename = " + str(fileName)
    for ibin in range(0,histTotSysUp.GetNbinsX()):
        totSys = 0.0
        totSysUp = 0.0
        for systematic in systematics:
                histSysUpName = systematic + "_UP"
                histSysUp = file.Get(histSysUpName)
                sysUp = 100.0*(histSysUp.GetBinContent(ibin+1) - 1.0)
                totSysUp = totSysUp + (sysUp**2.0)
        if fileName.split("_")[1] == "absolute":
            sysUp = 100.0*(0.025)
            totSysUp = totSysUp + (sysUp**2.0)
            if fileName.split("_")[2].split(".")[0] == "parton":
                totSysUp = totSysUp + ((1.5)**2.0)



        totSysUp = totSysUp**0.5

        histTotSysUp.SetBinContent(ibin+1, totSysUp)
        histTotSysDown.SetBinContent(ibin+1, -totSysUp)

    return (histTotSysUp, histTotSysDown)

def make_error_graph(fileName):
    file = TFile(fileName)
    gr = file.Get("data")
    gr_stat = file.Get("data_staterror_only")
    histMC = file.Get("mc")
    gr_n = TGraphAsymmErrors(gr.GetN())
    x, y = Double(0), Double(0)
    x_stat, y_stat = Double(0), Double(0)
    x_sys, y_sys = Double(0), Double(0)
    for ibin in range(0, gr.GetN()):
        gr.GetPoint(ibin, x, y)
        gr_stat.GetPoint(ibin, x_stat, y_stat)
        eyh = gr.GetErrorYhigh(ibin)
        eyl = gr.GetErrorYlow(ibin)
        binError = histMC.GetBinWidth(ibin+1)/2.0
        eyh_stat = gr_stat.GetErrorYhigh(ibin)
        eyl_stat = gr_stat.GetErrorYlow(ibin)
        eyh_sys = ((eyh**2.0) - (eyh_stat**2.0))**0.5
        eyl_sys = ((eyl**2.0) - (eyl_stat**2.0))**0.5

        gr_n.SetPoint(ibin, x, 0.0)
        gr_stat.SetPoint(ibin, x, 0.0)
        gr_n.SetPointEYhigh(ibin, 100.0*(eyh/y))
        gr_n.SetPointEYlow(ibin, 100.0*(eyl/y))
        gr_n.SetPointEXhigh(ibin,binError)
        gr_n.SetPointEXlow(ibin,binError)

        gr_stat.SetPointEYhigh(ibin, 100.0*(eyh_stat/y))
        gr_stat.SetPointEYlow(ibin, 100.0*(eyl_stat/y))
        gr_stat.SetPointEXhigh(ibin,binError)
        gr_stat.SetPointEXlow(ibin,binError)
    return (gr_n, gr_stat)

def histo2graph(hist):
    graph = TGraphErrors()
    for ibin in range(0, hist.GetNbinsX()):
        binX = hist.GetBinCenter(ibin+1)
        binY = hist.GetBinContent(ibin+1)
        binW = hist.GetBinWidth(ibin+1)
        graph.SetPoint(ibin, binX, binY)
        graph.SetPointError(ibin, binW/2.0, 0.0)
    lineColor = hist.GetLineColor()
    lineStyle = hist.GetLineStyle()
    graph.SetLineColor(lineColor)
    graph.SetLineStyle(lineStyle)
    graph.SetLineWidth(2)
    return graph

def make_summary_plot(fileName, systematics, sysHistUp, sysHistDown, gr_er, gr_er_stat):
    observable = fileName.split("_")[0]
    canvas = make_systematics_canvas()
    file = TFile(fileName)
    gr = totalGraph(sysHistUp, sysHistDown)
    gr.SetFillColor(kOrange)
    gr.SetLineColor(kOrange)

    histNominal = file.Get("Nominal")
    print "making summary plot for = " + str(fileName)

    mode = fileName.split("_")[1]
    mode2 = fileName.split("_")[2].split(".")[0]


    # if limits are set at defaut values, try find a "smart range"
    unc_hi_max = 0.0
    x = ROOT.Double(0.0)
    y = ROOT.Double(0.0)
        #if ((histNominal.GetMinimum() == 0.0) & (histNominal.GetMaximum() == 12.5)):
    for n in range(0, gr_er.GetN() ):
        gr_er.GetPoint(n, x, y)
        er = (gr.GetErrorY(n))
        erX = (gr.GetErrorX(n))
        er2 = (gr_er.GetErrorY(n))
        erX2 = (gr_er.GetErrorX(n))
        
        er_stat = (gr_er_stat.GetErrorY(n))
        er = (er**2 + er_stat**2)**0.5
        
        #print "missing uncertainty = " + str(((er2**2)-(er**2))**0.5)

        #print "n x y.  = " + str(n) + " x = " + str(x) + " y = " + str(y) + " error from sum  = " + str(er) + " Xerror from sum  = " + str(erX) + " orig error  = " + str(er2)  + " Xerror from orig  = " + str(erX2)
        gr.SetPoint(n, x, 0.0)
        gr.SetPointEYhigh(n, er)
        #print "rolling uncert.  = " + str(er)
        if (er > unc_hi_max):
            unc_hi_max = er
    histNominal.GetYaxis().SetRangeUser(0.0, (unc_hi_max*1.1))

    #print "max uncert.  = " + str(unc_hi_max)

    histCombJES = histNominal.Clone()
    histCombJES.Reset()
    histCombJES.SetDirectory(0)

    histMergeExp = histNominal.Clone()
    histMergeExp.Reset()
    histMergeExp.SetDirectory(0)

    grModList = {}

    nUnmerged = len(systematics) - len(mergeExp[mode][observable])
        #for i in range(0, nUnmerged):
        #grMod = TGraphErrors()
        #grModList.append(grMod)
        #grModList[]

    histNominal.Reset()
    histNominal.Draw()

    #gr_er.SetFillColor(kOrange)
    #gr_er.SetLineColor(kOrange)
    gr.Draw("E2SAME")
    #gr_er.Draw("SAME")

    #lumi, BR
    gr_lumi = TGraphAsymmErrors(gr_er.GetN())
    gr_br = TGraphAsymmErrors(gr_er.GetN())

    if (mode == "absolute"):
        for p in range(0,gr_er.GetN()):
            gr_er.GetPoint(p, x, y)
            gr_lumi.SetPoint(p, x, 2.5)
            ex = gr.GetErrorX(p)
            gr_lumi.SetPointEXlow(p, ex)
            gr_lumi.SetPointEXhigh(p, ex)
        if (mode2 == "parton"):
            for p in range(0,gr_er.GetN()):
                gr_er.GetPoint(p, x, y)
                gr_br.SetPoint(p, x, 1.5)
                ex = gr.GetErrorX(p)
                gr_br.SetPointEXlow(p, ex)
                gr_br.SetPointEXhigh(p, ex)


    gr_er_stat.SetFillColor(kGray+1)
    gr_er_stat.SetLineColor(kGray+1)
    gr_er_stat.Draw("E2SAME")

    gr_lumi.SetLineColor(kGreen+1)
    gr_lumi.SetLineStyle(3)
    gr_lumi.SetLineWidth(4)
    gr_lumi.Draw("ZSAME")

    gr_br.SetLineColor(kBlue+1)
    gr_br.SetLineStyle(3)
    gr_br.SetLineWidth(4)
    gr_br.Draw("ZSAME")

    isis = 0
    iModSys = 0
    for systematic in systematics:
        histSysUpName = systematic + "_UP" 
        histSysDownName = systematic + "_DOWN" 
        histSysUp = file.Get(histSysUpName)
        histSysDown = file.Get(histSysDownName)
        color = histSysUp.GetLineColor()
        if "JES" not in systematic:
            if (systematic not in mergeExp[mode][observable]):
                histSysUp =  centerHistAtZero(histSysUp)
                grMod = TGraphErrors()
                grModList[systematic] = histo2graph(histSysUp)
                grModList[systematic].Draw("ZSAME")
                iModSys = iModSys + 1
            elif (systematic in mergeExp[mode][observable]):
                for iBin in range(0,histMergeExp.GetNbinsX()):
                    runningMergeExpBin = histMergeExp.GetBinContent(iBin+1) +  ((histSysUp.GetBinContent(iBin+1)-1.0)**2)
                    histMergeExp.SetBinContent(iBin+1,runningMergeExpBin)
        else:
           for iBin in range(0,histCombJES.GetNbinsX()):
              runningJESBin = histCombJES.GetBinContent(iBin+1) +  ((histSysUp.GetBinContent(iBin+1)-1.0)**2)
              histCombJES.SetBinContent(iBin+1,runningJESBin)
        isis = isis + 1

    for iBin in range(0,histCombJES.GetNbinsX()):
        sqrtJESBin = (histCombJES.GetBinContent(iBin+1))**0.5
        histCombJES.SetBinContent(iBin+1,1.0 + sqrtJESBin)
        sqrtMergeExpBin = (histMergeExp.GetBinContent(iBin+1))**0.5
        histMergeExp.SetBinContent(iBin+1,1.0 + sqrtMergeExpBin)

    # add merged JES contribution
    histCombJES = centerHistAtZero(histCombJES)
    grCombJES = histo2graph(histCombJES) 
    histCombJES.SetLineColor(kRed)
    histCombJES.SetLineStyle(7)
    histCombJES.SetLineWidth(4)
    histCombJES.Draw("HISTSAME")

    # add merged exp. contribution
    histMergeExp = centerHistAtZero(histMergeExp)
    grMergeExp = histo2graph(histMergeExp) 
    histMergeExp.SetLineColor(kBlack)
    histMergeExp.SetLineStyle(7)
    histMergeExp.SetLineWidth(4)
    histMergeExp.Draw("HISTSAME")

#####legends
    legend_split = 5
    n_legend_entries = len(systematics) - 19 - len(mergeExp)
    n_leg = (n_legend_entries)/(legend_split)
    legendWidth = 0.18
#    legendsStart = (1.0/n_leg) - (legendWidth/1.5)
    legendsStart = 0.78
    legendY1 = 0.78
    legendY2 = 0.91

    leg = TLegend(legendsStart, legendY1, legendsStart+legendWidth, legendY2)
    leg.AddEntry(gr,"Total", "f")
    leg.AddEntry(gr_er_stat, "Stat", "f")
    
    if (mode == "normalised"):
        leg_1 = TLegend(legendsStart+0.1, legendY1+0.065, legendsStart+0.27, legendY2)
    if ((mode == "absolute") & (mode2 == "particle")):
        leg_1 = TLegend(legendsStart+0.1, legendY1+0.065, legendsStart+0.27, legendY2)
        leg_1.AddEntry(gr_lumi,"Lumi", "l")
    if ((mode == "absolute") & (mode2 == "parton")):
        leg_1 = TLegend(legendsStart+(legendWidth*0.5), legendY1, legendsStart+(legendWidth*1.5), legendY2)
        leg_1.AddEntry(gr_lumi,"Lumi", "l")
        leg_1.AddEntry(gr_br,"Br(W#rightarrowl#nu)", "l")

    leg_2 = TLegend(legendsStart, 0.1, legendsStart+legendWidth, (legendY2-0.13))
    leg_2.AddEntry(histCombJES,"JES", "l")
    leg_2.AddEntry(histMergeExp, "Other exp syst", "l")

    #print "there are "  +  str(len(grModList)) + " mod sys."
    #print "mod list "  +  str(grModList)

    for iG in range (0,len(systematics)):
        if ((systematics[iG] not in mergeExp[mode][observable]) & ("JES" not in systematic_names[iG])) :
            leg_2.AddEntry(grModList[systematics[iG]],systematic_names[iG], "l")
        """    if iG < 3:
                leg_2.AddEntry(grModList[iG],systematic_names[iG+8], "l")
            elif ((iG >= 3) & (iG <= 7)):
                leg_3.AddEntry(grModList[iG],systematic_names[iG+8], "l")
            elif (iG > 7):
                leg_4.AddEntry(grModList[iG],systematic_names[iG+8], "l")
                """

    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.029)
    leg.Draw()

    leg_1.SetBorderSize(0)
    leg_1.SetFillStyle(0)
    leg_1.SetTextSize(0.029)
    leg_1.Draw()

    leg_2.SetBorderSize(0)
    leg_2.SetFillStyle(0)
    leg_2.SetTextSize(0.029)
    leg_2.Draw()

    canvas.RedrawAxis()
    canvas.SetTicky(1)
    canvas.SetTickx(0)
#    canvas.SetLogy()
    plotName = fileName.split(".")[0] + ".pdf"

    H = canvas.GetWh()
    W = canvas.GetWw()
    l = canvas.GetLeftMargin()
    t = canvas.GetTopMargin()
    r = canvas.GetRightMargin()
    b = canvas.GetBottomMargin()
    extraOverCmsTextSize  = 0.89

    latex = TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextSize(0.43*t)
    latex.SetTextColor(kBlack)
    latex.SetTextFont(61)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.165,0.913,"CMS")

    latex.SetTextFont(42)
    latex.SetTextSize(0.39*t)
    latex.DrawLatex(0.77,0.913,"35.9 fb^{-1} (13 TeV)");

    latex.SetTextSize(0.34*t)
    resultLevel = fileName.split("_")[2].split(".")[0] 
    resultType  =  fileName.split("_")[1]
    if resultLevel == "pseudo":
        resultLevel = "particle"
    if resultType == "normalized":
        resultType = "normalised"


    plotString = "Dilepton, " + resultLevel +" level, \n" + resultType
    latex.DrawLatex(0.45,0.85,plotString);
    plotPath = "/Users/keaveney/Desktop/TOP-17-014-SupplementaryMaterial/UncBreakdown_" + plotName
    canvas.SaveAs(plotPath)
####MAIN PROGRAM#############
for xsec_type in xsec_types:
    for xsec_level in xsec_levels:
        iobs = 0
        if xsec_level == "pseudo":
            xsec_level = "particle"
        if xsec_type == "normalized":
            xsec_type = "normalised"
        for observable in observables[xsec_level]:
            resultsFileName = directory_base_results + xsec_level + "/" + xsec_type + directory_tail_results + "DiffXS_Hyp"+ observable+ "_source.root"
            gr_er, gr_er_stat = make_error_graph(resultsFileName)
            fileName = observable + "_" + xsec_type + "_" + xsec_level + ".root"
            canvasName = observable + "_" + xsec_type + "_" + xsec_level + ".pdf"
            file = makeOutputFile(fileName)
            if xsec_level == "particle":
                xsec_level = "pseudo"
            if xsec_type == "normalised":
                xsec_type = "normalized"
            nominal = directory_base + xsec_type + "_" + xsec_level + directory_tail + "Nominal" + "/combinedUnfolded/Hyp" + observable + "Results.txt"
            if observable == "TTBarMass":
                histNominal = processNominal(nominal, "Nominal", 0.0, 24.0)
            elif observable == "ToppT":
                histNominal = processNominal(nominal, "Nominal", 0.0, 8.5)
            elif observable == "AntiToppT":
                histNominal = processNominal(nominal, "Nominal", 0.0, 14.0)
            elif observable == "TopRapidity":
                histNominal = processNominal(nominal, "Nominal", 0.0, 6.5)
            elif observable == "AntiTopRapidity":
                histNominal = processNominal(nominal, "Nominal", 0.0, 6.5)
            elif observable == "TTBarpT":
                histNominal = processNominal(nominal, "Nominal", 0.0, 14.0)
            elif observable == "LLBarDPhi":
                histNominal = processNominal(nominal, "Nominal", 0.0, 8.0)
            elif observable == "TTBarDeltaRapidity":
                histNominal = processNominal(nominal, "Nominal", 0.0, 8.5)
            else:
                histNominal = processNominal(nominal, "Nominal", 0.0, 12.5)
            histNominal.SetTitle("")
            if xsec_level == "pseudo":
                xsec_level = "particle"
            if xsec_type == "normalized":
                xsec_type = "normalised"
            if  units[xsec_level][iobs] != "":
                xTitle = vars_root[xsec_level][iobs] + " [" +  units[xsec_level][iobs] + "]"
            else:
                xTitle = vars_root[xsec_level][iobs]
            iobs = iobs + 1

            histNominal.GetXaxis().SetTitleOffset(1.2)
            histNominal.SetXTitle(xTitle)
            yTitle = "Uncertainty (%)"      
            histNominal.SetYTitle(yTitle)
            histNominal.Write()
            icolor = 0
            for systematic in systematic_names:
                color = systematic_colors[systematic]
                style = systematic_styles[systematic]
                histSysUp, histSysDown = processSystematic(observable, xsec_type, xsec_level, systematics[icolor], histNominal)
                histSysUp.SetLineColor(color) 
                histSysDown.SetLineColor(color) 
                histSysUp.SetLineStyle(style) 
                histSysDown.SetLineStyle(style)
                histSysUp.Write()         
                histSysDown.Write()
                icolor = icolor + 1

            file.Close()
            totSysUp = TH1F()
            totSysDown = TH1F()
            totSysUp, totSysDown = make_sys_sum(fileName, systematics, variations)
            make_summary_plot(fileName, systematics, totSysUp, totSysDown, gr_er,gr_er_stat)
