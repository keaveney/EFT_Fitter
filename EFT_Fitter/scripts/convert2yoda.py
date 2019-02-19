from ROOT import *
import os


filesPath = "/Users/keaveney/EFT_Fitter/EFT_Fitter/files/Jan18/particle/normalised/"
targetFile = "particle_level.yoda"

files = [
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
"DiffXS_HypTTBarDeltaPhi_source.root",
"DiffXS_HypTTBarDeltaRapidity_source.root",
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
"DiffXS_HypBBBarMass_source.root",
"DiffXS_HypBBBarpT_source.root",
"DiffXS_HypJetMultpt30_source.root"
]


allFiles = []

for filename in os.listdir(filesPath):
    if ((filename.endswith(".root")) & (filename.split("_")[0] == "DiffXS" )):
        allFiles.append(filename)
        if (filename not in files):
            print filename

print "There is a total of " + str(len(allFiles)) + " suitable files. We will analsyse " + str(len(files))

y = open(targetFile,"w+")

i = 1

for filename in files:
    filePath = filesPath + filename
    f =  TFile(filePath)
    h = f.Get("mc")
    if (i < 10):
        histoName = "d0" + str(i) + "-x01-y01"
    else:
        histoName = "d" + str(i) + "-x01-y01"

    y.write("BEGIN YODA_SCATTER2D /REF/CMS_1D_DIFF/")
    y.write(histoName)
    y.write("\n")
        
    y.write("Path=/REF/CMS_1D_DIFF/")
    y.write(histoName)
    y.write("\n")
        
    y.write("Title=\n")
    y.write("Type=Scatter2D\n")

    y.write("# xval     xerr-     xerr+     yval     yerr-     yerr+\n")

    for bin in range(1,h.GetSize()-1):
        #print "bin = " + str(bin) + " bin centre = " + str(h.GetBinCenter(bin))
        binString =  str(h.GetBinCenter(bin))  + " "  + str(h.GetBinWidth(bin)/2.0) + " " + str(h.GetBinWidth(bin)/2.0) + " 0.0 0.0 0.0 \n"
        y.write(binString)

    y.write("END YODA_SCATTER2D\n")
    y.write("    \n")
    i = i+1
