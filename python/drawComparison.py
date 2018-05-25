#  Compare different plots from the same
#  OR
#  compare the same plot from different files

#!/bin/python
import glob, os, ROOT, math, sys
from ROOT import TCanvas, TGraph, TLegend, TF1, TH1, TH1F, TLatex
from ROOT import gROOT, gStyle
from array import array
from math import tan
import yaml # WARNING use this environment ~/env/setupPyTools27.env to enable yaml package!!!
ROOT.gROOT.LoadMacro("CLIC_style/CLICdpStyle/rootstyle/CLICdpStyle.C+")
ROOT.CLICdpStyle()

#  yamlFile = "/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/python/config/drawComparison/efficiency_vs_energy.yml"
#  yamlFile = "/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/python/config/drawComparison/recoTheta_thetaBins.yml"
#  yamlFile = "/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/python/config/drawComparison/efficiency_vs_theta_superimpose.yml"
yamlFile = "/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/python/config/drawComparison/energy.yml"

histColor = [ROOT.kBlack, ROOT.kRed-7, ROOT.kBlue, ROOT.kGreen+2, ROOT.kCyan+1, ROOT.kRed+2, ROOT.kOrange, ROOT.kViolet+2, ROOT.kGray]
markerStyle = [ROOT.kOpenCircle, ROOT.kOpenSquare, ROOT.kOpenTriangleUp, ROOT.kOpenDiamond, ROOT.kOpenCross, ROOT.kOpenCircle, ROOT.kOpenCircle]

if __name__ == "__main__":

        print ("")
	if (len(sys.argv)>=2):
		yamlFile = sys.argv[1]

	gStyle.SetOptStat(0)

	with open(yamlFile, 'r') as ymlfile:
		globalCfg = yaml.load(ymlfile)
                print ("[info]\t Read yaml file: \n\t\t%s" % (yamlFile))

        defaultCfg = None
        if (globalCfg.get("default") is not None):
            defaultCfg = globalCfg.get("default")
            for cfgIterator in globalCfg:
                cfg = globalCfg[cfgIterator]
                if (cfg == defaultCfg):
                    continue
                for cfgIt in defaultCfg:
                    if (cfg.get(cfgIt) is None):
                        cfg[cfgIt] = defaultCfg.get(cfgIt)


        for cfgIterator in globalCfg:
                cfg = globalCfg[cfgIterator]
                if (cfg == defaultCfg):
                    continue
                #  print (cfgIterator)
                #  print (cfg)

	for cfgIterator in globalCfg:
		cfg = globalCfg[cfgIterator]
                if (cfg == defaultCfg):
                    continue

		if ("dirPrefix" in cfg):
			dirPrefix = cfg['dirPrefix']
		if ("legTitle" in cfg):
			legTitle = cfg['legTitle']
		if ("fileName" in cfg):
			fileName = cfg['fileName']
		if ("histName" in cfg):
			histName = cfg['histName']
		if ("histColor" in cfg):
			histColor = cfg['histColor']
		if ("markerStyle" in cfg):
			markerStyle = cfg['markerStyle']
		if ("legPos" in cfg):
			leg = TLegend(cfg["legPos"][0],cfg["legPos"][1],cfg["legPos"][2],cfg["legPos"][3])

		hists = []
		nEntries = []
                histCounter = -1
		for i in range(0,len(fileName)):
			myFile = ROOT.TFile.Open(dirPrefix+fileName[i],"read")
                        print ("[info]\t Open file: \n\t\t%s" % (dirPrefix+fileName[i]))
			ROOT.SetOwnership(myFile,False)
			for k in range(0,len(histName)):
                                print ("[info]\t Getting hist: \n\t\t%s" % (histName[k]))
                                if not myFile.Get(histName[k]):
                                    print ('[ERROR]\t Hist "%s" is not found in file "%s"!' % (histName[k],fileName[i]))
                                    if ("skipMissingHists" in cfg):
                                        continue
                                    else:
                                        print ('[FATAL]\t Terminating...')
                                        sys.exit()
                                hist = myFile.Get(histName[k])
                                histCounter = histCounter + 1
                                #  print (hist)
                                #  print ("i: %i, k: %i; k*len(fileName)+i: %i" % (i,k,k*len(fileName)+i))
				hist.SetLineColor(histColor[histCounter])
				hist.SetMarkerColor(histColor[histCounter])
				hist.SetMarkerStyle(markerStyle[histCounter])
				hist.SetLineWidth(2)
                                #  hist.Sumw2()
				if ("lineWidth" in cfg):
					hist.SetLineWidth(cfg["lineWidth"])
				if ("markerSize" in cfg):
					hist.SetMarkerSize(cfg["markerSize"])
				if ("legCaptionMode" in cfg):
					nEntries.append(hist.GetEntries())
				if ("histTitle" in cfg):
					hist.SetTitle(cfg['histTitle'])
				if ("xAxisTitleOffset" in cfg):
					hist.GetXaxis().SetTitleOffset(cfg['xAxisTitleOffset'])
				if ("yAxisTitleOffset" in cfg):
					hist.GetYaxis().SetTitleOffset(cfg['yAxisTitleOffset'])
				if ("rebinFactor" in cfg):
					hist.Rebin(cfg['rebinFactor'])
				if ("scaleFactor" in cfg):
					hist.Scale(cfg['scaleFactor'])
                                #  if ("drawOption" in cfg):
                                #      if ("HIST" in cfg['drawOption']):
                                #          for j in range(0,hist.GetNbinsX()):
                                #              hist.SetBinError(j+1,0)
				if ("sortXaxisInGraph" in cfg):
					if (cfg["sortXaxisInGraph"]==True):
						histPoints = []
						for j in range(0,hist.GetN()):
							histPoints.append([hist.GetX()[j],hist.GetY()[j],hist.GetEX()[j],hist.GetEY()[j]])
						histPoints = sorted(histPoints, key=lambda x: x[0])
						for j in range(0,hist.GetN()):
							hist.SetPoint(j+1,histPoints[j][0],histPoints[j][1])
							hist.SetPointError(j+1,histPoints[j][2],histPoints[j][3])
				hists.append(hist)
		totalEntries = sum(nEntries)

		if ("canPos" in cfg):
			c1 = TCanvas( 'c1', 'A Simple Graph Example',  cfg["canPos"][0],cfg["canPos"][1],cfg["canPos"][2],cfg["canPos"][3])
		else:
			c1 = TCanvas( 'c1', 'A Simple Graph Example', 0, 0, 800, 700 )
		for i in range(0,len(hists)):
			if ("legPos" in cfg):
				if ("legCaptionMode" in cfg):
					if (cfg["legCaptionMode"]=="fraction"):
						leg.AddEntry(hists[i],"%s (%i%%)" % (legTitle[i], round(100.0*nEntries[i]/totalEntries)),"lp")
					elif (cfg["legCaptionMode"]=="entries"):
						leg.AddEntry(hists[i],"%s (%i)" % (legTitle[i], nEntries[i]),"lp")
					elif (cfg["legCaptionMode"]=="mean"):
						leg.AddEntry(hists[i],"%s (%.2f GeV)" % (legTitle[i], hists[i].GetMean()),"lp")
				else:
					leg.AddEntry(hists[i],"%s" % (legTitle[i]),"lp")

			if (i==0):
				if ("xAxisRange" in cfg):
					hists[i].GetXaxis().SetRangeUser(cfg['xAxisRange'][0],cfg['xAxisRange'][1])
				if ("yAxisRange" in cfg):
					hists[i].GetYaxis().SetRangeUser(cfg['yAxisRange'][0],cfg['yAxisRange'][1])
				if ("xAxisLabel" in cfg):
					hists[i].GetXaxis().SetTitle(cfg['xAxisLabel'])
				if ("yAxisLabel" in cfg):
					hists[i].GetYaxis().SetTitle(cfg['yAxisLabel'])
                                drawOption = ""
				if ("drawOption" in cfg):
					drawOption = cfg['drawOption']
                                        sameDrawOption = drawOption
				if ("gridX" in cfg):
					c1.SetGridx(cfg['gridX'])
				if ("gridY" in cfg):
					c1.SetGridy(cfg['gridY'])
				if ("logX" in cfg):
					c1.cd(1).SetLogx(cfg['logX'])
				if ("logY" in cfg):
					c1.cd(1).SetLogy(cfg['logY'])
                                
                                print ("[info]\t drawOption: \n\t\t%s" % (drawOption))
				hists[i].Draw(drawOption)
			else:
                                sameDrawOption = ""
				if ("drawOption" in cfg):
                                        sameDrawOption = cfg['drawOption']
				if ("sameDrawOption" in cfg):
					sameDrawOption = cfg['sameDrawOption']
                                print ("[info]\t sameDrawOption: \n\t\t%s" % (sameDrawOption))
				hists[i].Draw(sameDrawOption+"same")
		if ("legPos" in cfg):
			leg.Draw("same")
		if ("additionalLabel" in cfg):
			l = ROOT.TLatex()
			l.SetNDC()
			l.SetTextFont(42)
			l.SetTextColor(ROOT.kBlack)
			l.SetTextSize(0.045)
			label = cfg['additionalLabel']
			x = 0.17
			y = 0.932
			l.DrawLatex(x,y,label);

                print ("")
                ROOT.gPad.SetPhi(150);
		c1.SaveAs("pictures/"+cfgIterator+".png")
		c1.SaveAs("pictures/"+cfgIterator+".pdf")
		c1.SaveAs("pictures/"+cfgIterator+".C")
