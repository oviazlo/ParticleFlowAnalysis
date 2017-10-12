#  Compare different plots from the same
#  OR
#  compare the same plot from different files

#!/bin/python
import glob, os, ROOT, math, sys
from ROOT import TCanvas, TGraph, TLegend, TF1, TH1, TH1F
from ROOT import gROOT, gStyle
from array import array
from math import tan
import yaml # WARNING use this environment ~/env/setupPyTools27.env to enable yaml package!!!

#  yamlFile = "/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/python/config/drawComparison/efficiency_vs_energy.yml"
#  yamlFile = "/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/python/config/drawComparison/recoTheta_thetaBins.yml"
yamlFile = "/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/python/config/drawComparison/efficiency_vs_theta_superimpose.yml"

if __name__ == "__main__":

	if (len(sys.argv)>=2):
		yamlFile = sys.argv[1]

	gStyle.SetOptStat(0)

	with open(yamlFile, 'r') as ymlfile:
		globalCfg = yaml.load(ymlfile)

	for cfgIterator in globalCfg:
		cfg = globalCfg[cfgIterator]

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
		if ("legPos" in cfg):
			leg = TLegend(cfg["legPos"][0],cfg["legPos"][1],cfg["legPos"][2],cfg["legPos"][3])

		hists = []
		nEntries = []
		for i in range(0,len(fileName)):
			myFile = ROOT.TFile.Open(dirPrefix+fileName[i],"read")
			ROOT.SetOwnership(myFile,False)
			for k in range(0,len(histName)):
				hist = myFile.Get(histName[k])
				hist.SetLineColor(histColor[k*len(fileName)+i])
				hist.SetMarkerColor(histColor[k*len(fileName)+i])
				hist.SetLineWidth(2)
				if ("legCaptionMode" in cfg):
					if (cfg["legCaptionMode"]=="entries"):
						nEntries.append(hist.GetEntries())
				if ("histTitle" in cfg):
					hist.SetTitle(cfg['histTitle'])
				if ("rebinFactor" in cfg):
					hist.Rebin(cfg['rebinFactor'])
				if ("sortXaxisInGraph" in cfg):
					if (cfg["sortXaxisInGraph"]==True):
						histPoints = []
						for j in range(0,hist.GetN()):
							histPoints.append([hist.GetX()[j],hist.GetY()[j]])
						histPoints = sorted(histPoints, key=lambda x: x[0])
						for j in range(0,hist.GetN()):
							hist.SetPoint(j+1,histPoints[j][0],histPoints[j][1])
				hists.append(hist)
		totalEntries = sum(nEntries)

		c1 = TCanvas( 'c1', 'A Simple Graph Example', 0, 0, 800, 600 )
		for i in range(0,len(hists)):
			if ("legPos" in cfg):
				if ("legCaptionMode" in cfg):
					if (cfg["legCaptionMode"]=="entries"):
						leg.AddEntry(hists[i],"%s (%i%%)" % (legTitle[i], round(100.0*nEntries[i]/totalEntries)),"lp")
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

				hists[i].Draw(drawOption)
			else:
				sameDrawOption = ""
				if ("sameDrawOption" in cfg):
					sameDrawOption = cfg['sameDrawOption']
				hists[i].Draw(sameDrawOption+"same")
		if ("legPos" in cfg):
			leg.Draw("same")
		c1.SaveAs("pictures/"+cfgIterator+".png")
		c1.SaveAs("pictures/"+cfgIterator+".pdf")
