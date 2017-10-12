#!/bin/python
import glob, os, ROOT, math
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
	fileList = getFiles("./","*.root")
	print (fileList)

	resultList = []

	for iFile in fileList:
		myFile = ROOT.TFile.Open(iFile,"read")
		barrelHist = myFile.Get("ECALBarrel/ECAL_Energy")
		endcapHist = myFile.Get("ECALEndcap/ECAL_Energy")
		barrelMeanE = barrelHist.GetMean()
		endcapMeanE = endcapHist.GetMean()
		#  print (iFile)
		largerMean = 0.0
		if (barrelMeanE>endcapMeanE):
			#  print ("Barrel E: %f" % (barrelMeanE) )
			largerMean = barrelMeanE	
		else:
			#  print ("Endcap E: %f" % (endcapMeanE) )
			largerMean = endcapMeanE / 1.022345
	
		theta = getThetaFromFileName(iFile)
		if (theta>45):
			# barrel
			largerMean = largerMean
		else:
			# endcap
			largerMean = largerMean
		resultList.append( [theta, largerMean] )
		#  print ("Theta: %d; E: %f" % (theta, largerMean))

	dict_p10 = dict()
	dict_p55 = dict()
	dict_p100 = dict()

	for iEl in resultList:
		if iEl[1]<15:
			dict_p10[iEl[0]] = iEl[1]
		elif iEl[1]<60:
			dict_p55[iEl[0]] = iEl[1]
		else:
			dict_p100[iEl[0]] = iEl[1]

	ratioDict = dict()
	iAngle = 20
	ratioDict[iAngle] = [ float(dict_p10[90-iAngle]/dict_p10[iAngle]), float(dict_p55[90-iAngle]/dict_p55[iAngle]), float(dict_p100[90-iAngle]/dict_p100[iAngle]) ]
	iAngle = 30
	ratioDict[iAngle] = [ float(dict_p10[90-iAngle]/dict_p10[iAngle]), float(dict_p55[90-iAngle]/dict_p55[iAngle]), float(dict_p100[90-iAngle]/dict_p100[iAngle]) ]
	iAngle = 40
	ratioDict[iAngle] = [ float(dict_p10[90-iAngle]/dict_p10[iAngle]), float(dict_p55[90-iAngle]/dict_p55[iAngle]), float(dict_p100[90-iAngle]/dict_p100[iAngle]) ]
	
	
	#  dummyArr = [20,30,40]
	for kk in [20,30,40]:
		print ("Theta: %d - E: %f; Theta: %d - E: %f; ratio: %f" % (kk, dict_p10[kk], 90-kk, dict_p10[90-kk], dict_p10[90-kk]/dict_p10[kk]))
	for kk in [20,30,40]:
		print ("Theta: %d - E: %f; Theta: %d - E: %f; ratio: %f" % (kk, dict_p55[kk], 90-kk, dict_p55[90-kk], dict_p55[90-kk]/dict_p55[kk]))
	for kk in [20,30,40]:
		print ("Theta: %d - E: %f; Theta: %d - E: %f; ratio: %f" % (kk, dict_p100[kk], 90-kk, dict_p100[90-kk], dict_p100[90-kk]/dict_p100[kk]))

	energyArr = array ('d')
	energyArr.append(10.0)
	energyArr.append(55.0)
	energyArr.append(100.0)

	dummyDoubleArr = array ('d')
	iAngle = 20
	dummyDoubleArr.append(ratioDict[iAngle][0])
	dummyDoubleArr.append(ratioDict[iAngle][1])
	dummyDoubleArr.append(ratioDict[iAngle][2])

	gr_d20 = TGraph(3, energyArr, dummyDoubleArr)

	iAngle = 30
	dummyDoubleArr = array ('d')
	dummyDoubleArr.append(ratioDict[iAngle][0])
	dummyDoubleArr.append(ratioDict[iAngle][1])
	dummyDoubleArr.append(ratioDict[iAngle][2])

	gr_d30 = TGraph(3, energyArr, dummyDoubleArr)

	iAngle = 40
	dummyDoubleArr = array ('d')
	dummyDoubleArr.append(ratioDict[iAngle][0])
	dummyDoubleArr.append(ratioDict[iAngle][1])
	dummyDoubleArr.append(ratioDict[iAngle][2])

	gr_d40 = TGraph(3, energyArr, dummyDoubleArr)

	styleGraph(gr_d20,1)
	styleGraph(gr_d30,2)
	styleGraph(gr_d40,3)

        c1 = TCanvas( 'c1', 'A Simple Graph Example', 0, 0, 800, 600 )
        myFrame = c1.DrawFrame(0,1,110,1.03)
        myFrame.GetXaxis().SetLabelSize(0.03)
        myFrame.GetXaxis().SetTitle("Energy [GeV]")
        myFrame.GetYaxis().SetTitle("Correction")


	gr_d20.Draw("P2 same")
	gr_d30.Draw("P2 same")
	gr_d40.Draw("P2 same")

	leg = TLegend(0.35, 0.14, 0.64, 0.36)
	leg.SetBorderSize(0)
	leg.SetTextSize(0.03)
	leg.AddEntry(gr_d20,"20 degree","p")
	leg.AddEntry(gr_d30,"30 degree","p")
	leg.AddEntry(gr_d40,"40 degree","p")

	leg.Draw("same")

	c1.SaveAs("test.png")

	tmpSum = 0.0
	for it in [20,30,40]:
		for k in [0,1,2]:
			tmpSum = tmpSum + ratioDict[it][k]
	print ("Average correction: %f" % (tmpSum/9.0))
	print ("tan(TMath::Pi()*30./180.): %f" % (tan(3.14159*30./180.)))

