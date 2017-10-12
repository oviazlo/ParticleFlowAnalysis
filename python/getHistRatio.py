#!/bin/python
import glob, os, ROOT, math, sys
from ROOT import TCanvas, TGraph, TLegend, TF1, TH1, TH1F
from ROOT import gROOT, gStyle

#  numeratorDir = "PandoraPFOs_22_2112"
numeratorDir = "PandoraPFOs_2112"
numeratorHist =   [numeratorDir+"/truthParticle_onlyOneRecoPFO_Theta", numeratorDir+"/truthParticle_twoOrMoreRecoPFO_Theta"]
denominatorHist = "MCParticlesSkimmed/truthParticle_Theta"

if __name__ == "__main__":

	gStyle.SetOptStat(0)
	os.chdir("./")

	outRootFile = ROOT.TFile.Open("ratioFile.root","RECREATE")

	for iFile in glob.glob("particleGun*.root"):
		myFile = ROOT.TFile.Open(iFile,"read")
		num = myFile.Get(numeratorHist[0])
		for i in range(1,len(numeratorHist)):
			num.Add(myFile.Get(numeratorHist[i]))
		den = myFile.Get(denominatorHist)
		num.Divide(den)
		num.SetName(iFile.split(".")[0])
		outRootFile.cd()
		num.Write()
	outRootFile.Close()
				
