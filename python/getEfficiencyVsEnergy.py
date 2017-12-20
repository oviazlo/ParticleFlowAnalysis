#!/bin/python
import glob, os, ROOT, math, sys
from ROOT import TCanvas, TGraph, TLegend, TF1, TH1, TH1F
from ROOT import gROOT, gStyle

#  numeratorDir = "PandoraPFOs_22_2112"
#  absPath = "/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/outData/FCCee_singleParticle_performace/gamma/FCCee_testConvertionDetection/"
#  absPath = "/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/outData/FCCee_singleParticle_performace/gamma/FCCee_testConvertionDetection/e-/theta10_170/"
absPath = "./"
#  rootDir = "eventHists_photonAndNeutralRecl"
rootDir = ["eventHists","eventHists_photonRecl","eventHists_noConv","eventHists_photonAndNeutralRecl"]
#  rootDir = "eventHists"
#  rootDir = "eventHists_noConv"
#  rootDir = "eventHists_conv"
histName = "efficiencyVsEnergy"
fileNamePrefix = "particleGun_E"
fileNameIndex = ["1","2","5","10","20","50","100"]
fileNamePostfix = ".root"

def styleGraph(inGr, iColor):
    inGr.SetMarkerStyle(34)
    inGr.SetMarkerSize(1.2)
    inGr.SetMarkerColor(iColor)

if __name__ == "__main__":

    for iDir in rootDir:
        gStyle.SetOptStat(0)
        os.chdir("./")

        outRootFile = ROOT.TFile.Open(iDir+".root","RECREATE")

        #  tmpGr = TGraph(len(fileNameIndex))
        tmpGr = TGraph(len(fileNameIndex)-1)
        tmpGr.SetTitle("")
        tmpGr.GetXaxis().SetTitle("Energy [GeV]")
        tmpGr.GetYaxis().SetTitle("Efficiency")
        styleGraph(tmpGr,1)
        for i in range(0,len(fileNameIndex)):
            fileName = absPath + fileNamePrefix + fileNameIndex[i] + fileNamePostfix
            myFile = ROOT.TFile.Open(fileName,"read")
            iHist = myFile.Get(iDir+"/"+histName)
            maxBin = iHist.GetMaximumBin()
            energy = iHist.GetBinCenter(maxBin)
            eff = iHist.GetBinContent(maxBin)
            tmpGr.SetPoint(i,energy,eff)
            print("E:%f, eff:%f" % (energy, eff))

        
        outRootFile.cd()
        tmpGr.Write()
        outRootFile.Close()
        c1 = TCanvas( 'c1', 'A Simple Graph Example', 0, 0, 800, 600 )
        c1.cd()
        tmpGr.Draw("ALP")
        c1.SaveAs(iDir+".png")
                               

