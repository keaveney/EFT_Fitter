#sys variation plots
import ROOT
from ROOT import *
import os
import subprocess
from array import array

gStyle.SetOptStat(00000)

vars_root = [
    "p_{T}^{t}",
    "p_{T}^{#bar{t}}",
    "y_{t}",
    "y_{#bar{t}}",
    "p_{T}^{t#bar{t}}",
    "y_{t#bar{t}}",
    "m_{t#bar{t}}",
    "#Delta|y|(t,#bar{t})"]

units = [
    "GeV",
    "GeV",
    "",
    "",
    "GeV",
    "",
    "GeV",
    ""]

directory_base = "/nfs/dust/cms/user/savitsky/13TeV_Mor17/ForResults_8026p1_min/"
directory_tail = "/CMSSW_8_0_26_patch1/src/TopAnalysis/Configuration/analysis/diLeptonic/DIFFXS_PrimaryResults_TOP17014_backup/UnfoldingResults/"

directory_base_results = "/nfs/dust/cms/user/savitsky/RESULTS_MOR17/PlotsFor17014Paper/"
directory_tail_results = "_FinalResults_PlotsForPhD_CMSLogo/combinedUnfolded/" 

xsec_types = ["normalized"]
xsec_levels = ["parton"]

#systematic_names = ["ERD On", "Gluon move", "b-tagging", "b-tagging (pt)","b-tagging (eta)", "b-tagging (ljet)","b-tagging ljet (pt)","b-tagging ljet (eta)" , "m_{t}", "PS FSR Scale", "PS ISR Scale", "M.E. Ren. scale", "M.E. Fac. scale", "PDF", "JER", "Lepton SF","Backgrounds", "Pile Up","h_damp", "b frag.",  "b frag. (P)", "Br(b semiltonic) ", "Kin. Reco.", "Generator tunes", "Unclustered MET", "Trigger", "Trigger (eta)", "JESAbsoluteMPFBias", "JESAbsoluteScale", "JESAbsoluteStat", "JESFlavorQCD","JESFragmentation", "JESPileUpDataMC", "JESPileUpPtBB", "JESPileUpPtEC1", "JESPileUpPtEC2", "JESPileUpPtHF","JESPileUpPtRef","JESRelativeBal", "JESRelativeFSR","JESRelativeJEREC1","JESRelativeJEREC2", "JESRelativeJERHF","JESRelativePtBB","JESRelativePtEC1", "JESRelativePtEC2", "JESRelativePtHF","JESRelativeStatEC","JESRelativeStatFSR", "JESRelativeStatHF","JESSinglePionECAL","JESSinglePionHCAL","JESTimePtEta" ]

#systematics = ["ERDON","GLUONMOVETUNE", "BTAG", "BTAG_PT", "BTAG_ETA","BTAG_LJET", "BTAG_LJET_PT", "BTAG_LJET_ETA",  "MASS", "PSFSRSCALE", "PSISRSCALE",  "MERENSCALE", "MEFACSCALE", "PDF_ALPHAS","JER","LEPT","BG","PU","MATCH","BFRAG", "BFRAG_PETERSON", "BSEMILEP","KIN","UETUNE","UNCLUSTERED","TRIG","TRIG_ETA", "JESAbsoluteMPFBias", "JESAbsoluteScale", "JESAbsoluteStat", "JESFlavorQCD", "JESFragmentation","JESPileUpDataMC", "JESPileUpPtBB", "JESPileUpPtEC1", "JESPileUpPtEC2","JESPileUpPtHF","JESPileUpPtRef","JESRelativeBal","JESRelativeFSR","JESRelativeJEREC1","JESRelativeJEREC2","JESRelativeJERHF","JESRelativePtBB","JESRelativePtEC1", "JESRelativePtEC2","JESRelativePtHF","JESRelativeStatEC","JESRelativeStatFSR","JESRelativeStatHF","JESSinglePionECAL","JESSinglePionHCAL","JESTimePtEta"]

systematic_names =  ["JER", "Trigger","Lepton SF","Unclustered MET", "Kin. Reco", "b tagging (b)", "b tagging (l)", "Drell-Yan","Pile up", "h_{damp}", "Underlying event", "Br(b semileptonic)", "Backgrounds",  "PS FSR","PS ISR", "PDF", "PDF #alpha_{S}",  "U.E.Tune", "m_{t}", "Color Rec.", "ME Scale", "b fragmentation", "JESAbsoluteStat","JESAbsoluteScale","JESAbsoluteMPFBias","JESFragmentation","JESSinglePionECAL", "JESSinglePionHCAL", "JESFlavorQCD", "JESTimePtEta", "JESRelativeBal", "JESRelativeJEREC1", "JESRelativePtBB", "JESRelativePtEC1", "JESRelativeFSR", "JESRelativeStatFSR", "JESRelativeStatEC", "JESPileUpDataMC", "JESPileUpPtRef", "JESPileUpPtEC1", "JESPileUpPtBB"]


systematics =  ["JER","TOT_TRIG","LEPT", "UNCLUSTERED","KIN", "TOT_BTAG", "TOT_BTAG_LJET", "DY", "PU", "MATCH", "UETUNE", "BSEMILEP", "BG", "PSFSRSCALE", "PSISRSCALE", "PDF","PDF_ALPHAS","UETUNE", "MASS", "TOT_COLORREC", "TOT_SCALE", "TOT_BFRAG", "JESAbsoluteStat","JESAbsoluteScale","JESAbsoluteMPFBias","JESFragmentation","JESSinglePionECAL", "JESSinglePionHCAL", "JESFlavorQCD", "JESTimePtEta", "JESRelativeBal", "JESRelativeJEREC1", "JESRelativePtBB", "JESRelativePtEC1", "JESRelativeFSR", "JESRelativeStatFSR", "JESRelativeStatEC", "JESPileUpDataMC", "JESPileUpPtRef", "JESPileUpPtEC1", "JESPileUpPtBB" ]


#systematic_names = ["PS FSR Scale", "PS ISR Scale"]
#systematics = [ "PSFSRSCALE", "PSISRSCALE"]

#systematic_names = [ "JESAbsoluteMPFBias", "JESAbsoluteScale", "JESAbsoluteStat", "JESFlavorQCD","JESFragmentation", "JESPileUpDataMC", "JESPileUpPtBB", "JESPileUpPtEC1", "JESPileUpPtEC2", "JESPileUpPtHF","JESPileUpPtRef","JESRelativeBal", "JESRelativeFSR","JESRelativeJEREC1","JESRelativeJEREC2", "JESRelativeJERHF","JESRelativePtBB","JESRelativePtEC1", "JESRelativePtEC2", "JESRelativePtHF","JESRelativeStatEC","JESRelativeStatFSR", "JESRelativeStatHF","JESSinglePionECAL","JESSinglePionHCAL","JESTimePtEta"]

#systematics = ["JESAbsoluteMPFBias", "JESAbsoluteScale", "JESAbsoluteStat", "JESFlavorQCD", "JESFragmentation","JESPileUpDataMC", "JESPileUpPtBB", "JESPileUpPtEC1", "JESPileUpPtEC2","JESPileUpPtHF","JESPileUpPtRef","JESRelativeBal","JESRelativeFSR","JESRelativeJEREC1","JESRelativeJEREC2","JESRelativeJERHF","JESRelativePtBB","JESRelativePtEC1", "JESRelativePtEC2","JESRelativePtHF","JESRelativeStatEC","JESRelativeStatFSR","JESRelativeStatHF","JESSinglePionECAL","JESSinglePionHCAL","JESTimePtEta"]

variations = ["UP", "DOWN"]
observables = ["ToppT","AntiToppT","TopRapidity","AntiTopRapidity", "TTBarpT",  "TTBarRapidity",  "TTBarMass",  "TTBarDeltaRapidity"]
systematic_colors = [i for i in range (3, len(systematics)+3)]
#print "sys colors = " + str(systematic_colors)
iobs = 0

def makeOutputFile(fileName):
    file = TFile(fileName, "RECREATE")
    return file

def processNominal(path, histName):
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
        #print str(line.split( )[7])
        hist.SetBinContent(ibin+1, float(line.split( )[7]))
        ibin = ibin + 1
    hist.GetYaxis().SetRangeUser(0.6,1.6)
    hist.SetLineWidth(0)

    return hist

def totalGraph(histTotSysUp, histTotSysDown):
    gr = TGraphAsymmErrors(histTotSysUp.GetNbinsX())
    for ibin in range(0,histTotSysUp.GetNbinsX()):
        gr.SetPoint(ibin,histTotSysUp.GetBinCenter(ibin+1),1.0)
        gr.SetPointEXlow(ibin,histTotSysUp.GetBinWidth(ibin+1)/2.0)
        gr.SetPointEXhigh(ibin,histTotSysUp.GetBinWidth(ibin+1)/2.0)

        if ((histTotSysUp.GetBinContent(ibin+1) > 1.0) & (histTotSysDown.GetBinContent(ibin+1) < 1.0)):
            gr.SetPointEYhigh(ibin,histTotSysUp.GetBinContent(ibin+1)-1.0)
            gr.SetPointEYlow(ibin,1.0 - histTotSysDown.GetBinContent(ibin+1))
        elif ((histTotSysUp.GetBinContent(ibin+1) < 1.0) & (histTotSysDown.GetBinContent(ibin+1) > 1.0)):    
            gr.SetPointEYhigh(ibin,histTotSysDown.GetBinContent(ibin+1)-1.0)
            gr.SetPointEYlow(ibin,1.0 - histTotSysUp.GetBinContent(ibin+1))
    return gr

def reNormalise(histNominal, histUp, histDown):
    histVisUp = histNominal.Clone()
    histVisDown = histNominal.Clone()
    histVisUp.Reset()
    histVisDown.Reset()

    print "in reNormalise..."
    print "int. histNom  = "  + str(histNominal.Integral("width"))
    print "int. histUp  = "  + str(histVisUp.Integral("width"))
    print "int. histDown  = "  + str(histVisDown.Integral("width"))

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
    linkStr = ""
    singlePointSystematics = ["ERDON", "ERDONRETUNE", "GLUONMOVETUNE", "BFRAG_PETERSON"]

    sPS = 0

    if any(singlePointSystematic in systematic for singlePointSystematic in singlePointSystematics):
        sPS = 1

    linkStr = "_"
    variations = [""]
    for variation in variations:
        path = directory_base + xsec_type + "_" + xsec_level + directory_tail + systematic + linkStr + variation + "/combinedUnfolded/Hyp" + observable + "Results.txt"
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
            if systematic == "DY":
                print "DY UP = " + str(1.0 + unc)
                print "DY DOWN = " + str(1.0 - unc)


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

    for ibin in range(0,histTotSysUp.GetNbinsX()):
        totSys = 0.0
        totSysUp = 0.0
        for systematic in systematics:
                histSysUpName = systematic + "_UP"
                histSysUp = file.Get(histSysUpName)
                sysUp = 1.0 - histSysUp.GetBinContent(ibin+1)
                totSysUp = totSysUp + (sysUp**2.0)
   
        totSysUp = totSysUp**0.5
        
        histTotSysUp.SetBinContent(ibin+1, 1.0 + totSysUp)
        histTotSysDown.SetBinContent(ibin+1, 1.0 - totSysUp)

    return (histTotSysUp, histTotSysDown)

def make_error_graph(fileName):
    file = TFile(fileName)
    gr = file.Get("data")
    gr_stat = file.Get("data_staterror_only")
    gr_n = TGraphAsymmErrors(gr.GetN())
    x, y = Double(0), Double(0)
    x_stat, y_stat = Double(0), Double(0)
    x_sys, y_sys = Double(0), Double(0)
    for ibin in range(0, gr.GetN()):
        gr.GetPoint(ibin, x, y)
        gr_stat.GetPoint(ibin, x_stat, y_stat)
        eyh = gr.GetErrorYhigh(ibin)
        eyl = gr.GetErrorYlow(ibin)
        eyh_stat = gr_stat.GetErrorYhigh(ibin)
        eyl_stat = gr_stat.GetErrorYlow(ibin)
        eyh_sys = ((eyh**2.0) - (eyh_stat**2.0))**0.5
        eyl_sys = ((eyl**2.0) - (eyl_stat**2.0))**0.5
        gr_n.SetPoint(ibin, x ,1.0)
        gr_n.SetPointEYhigh(ibin, (eyh_sys/y))
        gr_n.SetPointEYlow(ibin, (eyl_sys/y))
        #print "data errors up (total, stat, sys) = " + str(eyh) + ", " + str(eyh_stat) + " , " + str(eyh_sys)
    return gr_n

def histo2graph(hist):
    graph = TGraphErrors()
    for ibin in range(0, hist.GetNbinsX()):
        binX = hist.GetBinCenter(ibin+1)
        binY = hist.GetBinContent(ibin+1)
        binW = hist.GetBinWidth(ibin+1)
        graph.SetPoint(ibin, binX, binY)
        graph.SetPointError(ibin, binW/2.0, 0.0)
    return graph

def make_summary_plot(fileName, systematics, sysHistUp, sysHistDown, gr_er):
    canvas = make_systematics_canvas()
    file = TFile(fileName)
    sysHistUp.SetLineWidth(3)
    sysHistDown.SetLineWidth(3)
    gr = totalGraph(sysHistUp, sysHistDown)
    gr.SetFillColor(kOrange)
    gr.SetLineColor(kOrange)

    histNominal = file.Get("Nominal")
 
    histCombJESUp = histNominal.Clone()
    histCombJESUp.Reset()
    histCombJESUp.SetDirectory(0)

    histCombJESDown = histNominal.Clone()
    histCombJESDown.Reset()
    histCombJESDown.SetDirectory(0)

    runningJESBin = 0.0
    grListUp = []
    grListDown = []

    for i in range(0, len(systematics)-19):
        grNJES = TGraphErrors()
        grListUp.append(grNJES)
        grListDown.append(grNJES)

    histNominal.Draw()
    gr.Draw("E2SAME")

    isis = 0
    for systematic in systematics:
        histSysUpName = systematic + "_UP" 
        histSysDownName = systematic + "_DOWN" 
        histSysUp = file.Get(histSysUpName)
        histSysDown = file.Get(histSysDownName)
        color = histSysUp.GetLineColor() 
        if "JES" not in systematic:
            grListUp[isis] = histo2graph(histSysUp)
            grListDown[isis] = histo2graph(histSysDown)
            grListUp[isis].SetLineColor(color) 
            grListDown[isis].SetLineColor(color) 
            grListUp[isis].Draw("ZSAME")
            grListDown[isis].Draw("ZSAME")
        else:
           for iBin in range(0,histCombJESDown.GetNbinsX()):
              runningJESBin = histCombJESDown.GetBinContent(iBin+1) +  ((histSysUp.GetBinContent(iBin+1)-1.0)**2)
              histCombJESDown.SetBinContent(iBin+1,runningJESBin)
        isis = isis + 1

    for iBin in range(0,histCombJESDown.GetNbinsX()):
        sqrtJESBin = (histCombJESDown.GetBinContent(iBin+1))**0.5
        histCombJESUp.SetBinContent(iBin+1,1.0 + sqrtJESBin)
        histCombJESDown.SetBinContent(iBin+1,1.0 - sqrtJESBin)

    grCombJESUp = histo2graph(histCombJESUp) 
    grCombJESDown = histo2graph(histCombJESDown)

    grCombJESUp.Draw("ZSAME")
    grCombJESDown.Draw("ZSAME")
    gr_er.Draw("E2SAME")

    grCombJESUp.SetLineColor(kRed)
    grCombJESDown.SetLineColor(kRed)


#####legends

    legend_split = 5
    n_leg = (len(systematics)-18)/(legend_split)
    legendWidth = 0.14
    legendsStart = (1.0/n_leg) - (legendWidth/1.5)
    legendY1 = 0.68
    legendY2 = 0.89

    leg = TLegend(legendsStart,  legendY1,legendsStart+legendWidth, legendY2)
    leg_2 = TLegend(legendsStart+legendWidth, legendY1 ,legendsStart+(legendWidth*2.0), legendY2)
    leg_3 = TLegend(legendsStart+(legendWidth*2.0), legendY1 ,legendsStart+(legendWidth*3.0), legendY2)
    leg_4 = TLegend(legendsStart+(legendWidth*3.0), legendY1 ,legendsStart+(legendWidth*4.0), legendY2)
    leg_5 = TLegend(legendsStart+(legendWidth*4.0), legendY1 ,legendsStart+(legendWidth*5.0), legendY2)

    leg.AddEntry(gr,"Total Syst", "f")
    leg.AddEntry(grCombJESUp,"JES", "l")

    iG = 0
    for g in grListUp:
        if iG < (legend_split - 2):
            leg.AddEntry(g,systematic_names[iG], "l")
        elif(iG >= (legend_split - 2)) & (iG < (legend_split*2  -2) ):
            leg_2.AddEntry(g,systematic_names[iG], "l")
        elif(iG >= ((legend_split*2) -2)) & (iG < (legend_split*3 -2 ) ):
            leg_3.AddEntry(g,systematic_names[iG], "l")
        elif(iG >= ((legend_split*3) -2)) & (iG < (legend_split*4 -2 ) ):
            leg_4.AddEntry(g,systematic_names[iG], "l")
        else:
            leg_5.AddEntry(g,systematic_names[iG], "l")
        iG = iG + 1

    plotName = fileName.split(".")[0] + ".pdf"
    leg.SetBorderSize(0)
    leg_2.SetBorderSize(0)
    leg_3.SetBorderSize(0)
    leg_4.SetBorderSize(0)
    leg_5.SetBorderSize(0)
    leg.SetTextSize(0.02)
    leg_2.SetTextSize(0.02)
    leg_3.SetTextSize(0.02)
    leg_4.SetTextSize(0.02) 
    leg_5.SetTextSize(0.02)
    leg.Draw()

    if len(systematics) > (legend_split -1):
        leg_2.Draw()
    if len(systematics) > (legend_split*2 -1):
        leg_2.Draw()
        leg_3.Draw()
    if len(systematics) > (legend_split*3 -1) :
        leg_2.Draw() 
        leg_3.Draw()
        leg_4.Draw()
    if len(systematics) > (legend_split*4 -1) :
        leg_2.Draw() 
        leg_3.Draw()
        leg_4.Draw()
        leg_5.Draw()

    canvas.RedrawAxis()
    canvas.SetTicky(1)

    H = canvas.GetWh();
    W = canvas.GetWw();
    l = canvas.GetLeftMargin();
    t = canvas.GetTopMargin();
    r = canvas.GetRightMargin();
    b = canvas.GetBottomMargin();
    extraOverCmsTextSize  = 0.89;

    latex = TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextSize(0.43*t)
    latex.SetTextColor(kBlack)
    latex.SetTextFont(61)
    latex.SetTextAlign(31)
    latex.DrawLatex(0.18,0.912,"CMS")

    latex.SetTextFont(42)
    latex.SetTextSize(0.39*t)
    latex.DrawLatex(0.9,0.912,"35.9 fb^{-1} (13 TeV)");

    latex.SetTextSize(0.32*t)
    resultLevel = fileName.split("_")[2].split(".")[0] 
    resultType  =  fileName.split("_")[1]
    if resultType == "normalized":
       resultType = "normalised"
    plotString = "Dilepton, " + resultLevel +", " + resultType
    latex.DrawLatex(0.4,0.6,plotString);
    canvas.SaveAs(plotName)

####MAIN PROGRAM#############

for observable in observables:                    
    for xsec_type in xsec_types:
        for xsec_level in xsec_levels: 
            resultsFileName = directory_base_results + xsec_type + "_" + xsec_level + directory_tail_results + "DiffXS_Hyp"+ observable+ "_source.root"
            gr_er = make_error_graph(resultsFileName)
            fileName = observable + "_" + xsec_type + "_" + xsec_level + ".root"
            canvasName = observable + "_" + xsec_type + "_" + xsec_level + ".pdf"
            file = makeOutputFile(fileName) 
            nominal = directory_base + xsec_type + "_" + xsec_level + directory_tail + "Nominal" + "/combinedUnfolded/Hyp" + observable + "Results.txt"
            #print "nominal = " + str(nominal)
            histNominal = processNominal(nominal, "Nominal")
            histNominal.SetTitle("")
            if  units[iobs] != "":
                xTitle = vars_root[iobs] + " [" +  units[iobs] + "]"
            else:
                xTitle = vars_root[iobs]

            histNominal.SetXTitle(xTitle)
#            if  units[iobs] != "":
#                yTitle = "#frac{d#sigma}{d" + vars_root[iobs] + "} [" + units[iobs] + "]"
#            else:
#                yTitle = "#frac{d#sigma}{d" + vars_root[iobs] + "}"
         #   yTitle = "#frac{d#sigma^{varied}}{d" + vars_root[iobs] + "} / " + "#frac{d#sigma^{nominal}}{d" + vars_root[iobs] + "}"
            yTitle = "Relative uncertainty (%)"      
            histNominal.SetYTitle(yTitle)
            histNominal.Write()
            icolor = 0
            for systematic in systematics:
                color = systematic_colors[icolor]
                icolor = icolor + 1
                histSysUp, histSysDown = processSystematic(observable, xsec_type, xsec_level, systematic, histNominal)
                histSysUp.SetLineColor(color) 
                histSysDown.SetLineColor(color) 
                histSysUp.Write()         
                histSysDown.Write()         
            file.Close()
            totSysUp = TH1F()
            totSysDown = TH1F()
            totSysUp, totSysDown = make_sys_sum(fileName, systematics, variations)
            make_summary_plot(fileName, systematics, totSysUp, totSysDown, gr_er)
    iobs = iobs + 1
