#!/bin/python
import glob, os, ROOT, math, sys
from ROOT import TCanvas, TGraph, TLegend, TF1
from ROOT import gROOT, gStyle
from array import array
from math import tan

def getThetaFromFileName(fileName):
	# example of input file name:
	#  ECAL_photonGun_E55_theta40.root
	return int(round(float(('.'.join(fileName.split(".")[0:-1])).split("_")[-1].split("theta")[-1])))

def getFiles(dirName, regexString):
	os.chdir(dirName)
	return glob.glob(regexString)

def styleGraph(inGr, iColor):
	inGr.SetMarkerStyle(34)
	inGr.SetMarkerSize(1.2)
	inGr.SetMarkerColor(iColor)

if __name__ == "__main__":

	if (len(sys.argv)<=1):
		print ("Specify input file!")
		sys.exit()

	gStyle.SetOptStat(0)

	myFile = ROOT.TFile.Open(sys.argv[1],"read")
	numerator = [myFile.Get("photonEfficiency/matchedMC_vs_E"),myFile.Get("photonEfficiency/matchedMC_vs_theta")]
	denominator = [myFile.Get("photonEfficiency/findableMC_vs_E"),myFile.Get("photonEfficiency/findableMC_vs_theta")]
	ratioPlots = []
	for i in range(0,len(numerator)):
		numerator[i].Sumw2()
		denominator[i].Sumw2()
		tmpHist = numerator[i].Clone("ratio")
		tmpHist.GetYaxis().SetTitle("Efficiency")
		if (i==1):
			tmpHist.GetXaxis().SetTitle("Theta [degree]")
		tmpHist.SetTitle("")
		tmpHist.Divide(denominator[i])
		ratioPlots.append(tmpHist)

	c1 = TCanvas( 'c1', 'A Simple Graph Example', 0, 0, 800, 600 )
	c1.Divide(1,len(numerator))

	for i in range(0,len(numerator)):
		c1.cd(i+1)
		ratioPlots[i].Draw()

	c1.SaveAs("temp.png")

