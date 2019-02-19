#touch up control plots
from ROOT import *
import os
import subprocess

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

directory = "/Users/keaveney/EFT_Fitter/EFT_Fitter/control_plot_files"
gStyle.SetOptStat(00000)
gStyle.SetErrorX(0)

for filename in os.listdir(directory):
    if filename.endswith(".root"):
        #print(os.path.join(directory, filename))
        file = TFile(filename)
        plotname = filename.split("_")
        if len(plotname) < 3:
            plotname_base = plotname[0]
        else:
            plotname_base = plotname[0] + "_" + plotname[1]
        print("Accessing " + plotname_base)
        canvasname = plotname_base + "_canvas"
        canvas = file.Get(canvasname)

        h_data = canvas.GetPrimitive(plotname_base)
        h_data.SetTitle("")
        
        print list(canvas.GetListOfPrimitives())

        for prim in canvas.GetListOfPrimitives():
            if prim.ClassName() == "TPaveText":
                if prim.GetLineWith("").GetTitle() == "CMS":
                    prim.SetX1NDC(0.17)
                if "35.9" in prim.GetLineWith("").GetTitle():
                    prim.SetX1NDC(0.79)
        r_pad = canvas.GetPrimitive("rPad")
        #print list(r_pad.GetListOfPrimitives())
        f_ratio = r_pad.GetPrimitive("f")
        f2_ratio = r_pad.GetPrimitive("f2")
        #f_ratio.SetLineColor(0)
        f2_ratio.SetLineColor(0)
        h_ratio = r_pad.GetPrimitive("ratio")
        h_ratio.SetYTitle("#frac{Data}{Pred.}")
        h_data.SetDrawOption("E0p")
        
        if (plotname_base == "HypBjetMulti_noBTag"):
            h_data.GetYaxis().SetRangeUser(2.0, 1000000000)
            h_ratio.GetYaxis().SetRangeUser(0.1, 3.5)
            #TAxis::SetBinLabel(bin, label)
            h_ratio.GetXaxis().SetBinLabel(7, "#geq 6")
            canvas.SetLogy()
        elif (plotname_base == "HypLeptonpT"):
            h_data.GetYaxis().SetRangeUser(10.0, 10000000)
            h_ratio.GetYaxis().SetRangeUser(0.8, 1.14)
            bin_width = int(h_data.GetBinWidth(1))
            bin_width_str = "Leptons / " + str(bin_width) + " GeV"
            h_data.SetYTitle(bin_width_str)
            canvas.SetLogy()
        elif (plotname_base == "HypJetMultpt30"):
            h_data.GetYaxis().SetRangeUser(10.0, 10000000)
            h_ratio.GetYaxis().SetRangeUser(0.5, 1.5)
            h_ratio.SetXTitle("N_{jets}")
            h_ratio.GetXaxis().SetBinLabel(8, "#geq 9")
            canvas.SetLogy()
        elif (plotname_base == "HypTTBarRapidity"):
            h_data.GetYaxis().SetRangeUser(10.0, 100000000)
            h_ratio.GetYaxis().SetRangeUser(0.05, 1.98)
            bin_width = int(h_data.GetBinWidth(1))
            bin_width_str = "Top quark pairs"
            h_data.SetYTitle(bin_width_str)
            canvas.SetLogy()
        elif (plotname_base == "HypTopRapidity"):
            h_data.GetYaxis().SetRangeUser(10.0, 100000000)
            h_ratio.GetYaxis().SetRangeUser(0.45, 1.65)
            bin_width = int(h_data.GetBinWidth(1))
            bin_width_str = "Top quarks"
            h_data.SetYTitle(bin_width_str)
            canvas.SetLogy()
        elif (plotname_base == "HypToppT"):
            h_data.GetYaxis().SetRangeUser(10.0, 10000000)
            h_ratio.GetYaxis().SetRangeUser(0.65, 1.19)
            bin_width = int(h_data.GetBinWidth(1))
            bin_width_str = "Top quarks / " + str(bin_width) + " GeV"
            h_data.SetYTitle(bin_width_str)
            canvas.SetLogy()
        elif (plotname_base == "HypTTBarpT"):
            h_data.GetYaxis().SetRangeUser(10.0, 10000000)
            h_ratio.GetYaxis().SetRangeUser(0.62, 1.28)
            bin_width = int(h_data.GetBinWidth(1))
            bin_width_str = "Top quark pairs / " + str(bin_width) + " GeV"
            h_data.SetYTitle(bin_width_str)
            canvas.SetLogy()
        elif (plotname_base == "HypTTBarMass"):
            h_data.GetYaxis().SetRangeUser(10.0, 10000000)
            h_ratio.GetYaxis().SetRangeUser(0.25, 1.75)
            bin_width = int(h_data.GetBinWidth(1))
            bin_width_str = "Top quark pairs / " + str(bin_width) + " GeV"
            h_data.SetYTitle(bin_width_str)
            canvas.SetLogy()
        elif (plotname_base == "HypBJetpT"):
            h_data.GetYaxis().SetRangeUser(10.0, 1000000)
            h_ratio.GetYaxis().SetRangeUser(0.35, 1.75)
            h_data.SetYTitle("b jets")
            h_ratio.SetXTitle("p^{b-jet}_{T} [GeV]")
            bin_width = int(h_data.GetBinWidth(1))
            bin_width_str = "b jets / " + str(bin_width) + " GeV"
            h_data.SetYTitle(bin_width_str)
            canvas.SetLogy()
        plotfilename = plotname_base + "_rev3.pdf"
        canvas.SaveAs(plotfilename)
        continue
    else:
        continue

#dir_list = []
#dir_list.append(directory)
#insert_eps_text(dir_list)
