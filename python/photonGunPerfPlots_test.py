#!/bin/python
import glob, os, ROOT, math, sys
from ROOT import TCanvas, TGraph, TLegend, TF1, TH1, TH1F
from ROOT import gROOT, gStyle
from array import array
from math import tan, sqrt

def getRMS90Resolution(pTH1F, resolution, resolutionError):
	

#	FLOAT_MAX = sys.float_info.max
	FLOAT_MAX = 9999.0

	sum = 0.
	total = 0.
	sx = 0.
	sxx = 0.
	nbins = pTH1F.GetNbinsX()

	for i in range(0,nbins):

		binx = pTH1F.GetBinLowEdge(i) + 0.5 * pTH1F.GetBinWidth(i)
		yi = pTH1F.GetBinContent(i)
		sx = sx + yi * binx
		sxx = sxx + yi * binx * binx
		total = total+ yi
    

	rawMean = sx / total
	rawMeanSquared = sxx / total

	rawRms = sqrt(rawMeanSquared - rawMean * rawMean)
	print ("rawRms: %f" % (rawRms))

	sum = 0.
	is0 = 0

    
	for i in range(0,nbins+1):
		if (sum < total / 10.):
			sum += pTH1F.GetBinContent(i)
			is0 = i
		else:
			break
	print ("sum: %f" % (sum))
	print ("total: %f" % (total))
	print ("is0: %d" % (is0))

	rmsmin = FLOAT_MAX
	sigma = FLOAT_MAX
	sigmasigma = FLOAT_MAX
	frac = FLOAT_MAX
	efrac = FLOAT_MAX
	mean = FLOAT_MAX
	low = FLOAT_MAX
	rms = FLOAT_MAX
	high = 0.0

	for istart in range(0,is0+1):
		sumn = 0.
		csum = 0.
		sumx = 0.
		sumxx = 0.
		iend = 0

		for istart in range(istart,nbins+1):
			if (csum < 0.9 * total):
				binx = pTH1F.GetBinLowEdge(i) + (0.5 * pTH1F.GetBinWidth(i))
				yi = pTH1F.GetBinContent(i)
				csum = csum + yi

				if (sumn < 0.9 * total):
					sumn = sumn + yi
					sumx = sumx + yi * binx
					sumxx = sumxx + yi * binx * binx
					iend = i
			else:
				break

		print ("iend: %d" % (iend))

		localMean = sumx / sumn
		localMeanSquared = sumxx / sumn
		localRms = sqrt(localMeanSquared - localMean * localMean)

        	if (localRms < rmsmin):
			mean = localMean
			rms = localRms
			low = pTH1F.GetBinLowEdge(istart)
			high = pTH1F.GetBinLowEdge(iend)
			rmsmin = localRms

		sigma = rms
		sigmasigma = sigma / sqrt(total)
    
	print ("mean: %f" % (mean))
	print ("rms: %f" % (rms))
	print ("rmsmin: %f" % (rmsmin))

	resolution = frac
	resolutionError = efrac

	print ("resolution: %f" % (resolution))
	print ("resolutionError: %f" % (resolutionError))

def getFiles(dirName, regexString):
	os.chdir(dirName)
	return glob.glob(regexString)

def addPathToFileName(path, fileList):
	for i in range(0, len(fileList)):
		fileList[i] = path + "/" + fileList[i]

def styleGraph(inGr, iColor):
	inGr.SetMarkerStyle(34)
	inGr.SetMarkerSize(1.2)
	inGr.SetMarkerColor(iColor)

def getTruthInfo(inFile):
	scale = 0.1
	energyHistName = "MCParticlesSkimmed/truthParticle/Energy"
	thetaHistName = "MCParticlesSkimmed/truthParticle/Theta"
	energyHist = inFile.Get(energyHistName)
	thetaHist = inFile.Get(thetaHistName)
	if not energyHist:
		print ("[ERROR]\tno hist found --> inFile: %s, hist: %s" % (inFile,energyHistName))
		sys.exit()
	if not energyHist:
		print ("[ERROR]\tno hist found --> inFile: %s, hist: %s" % (inFile,thetaHistName))
		sys.exit()
	energy = energyHist.GetMean()
	theta = thetaHist.GetMean()

	energy = math.floor( energy / scale + 0.5)*scale
	theta = math.floor( theta / scale + 0.5)*scale
	#  print ("%f, %f" % (energy, theta))
	return [energy, theta]

#  def getRAWEcalEnergy(inFile):
#          return [inFile.Get("ECalBarrelCollection/ECAL_Energy").GetMean(), inFile.Get("ECalEndcapCollection/ECAL_Energy").GetMean()]
#
#  def getRAWHcalEnergy(inFile):
#          return [inFile.Get("HCalBarrelCollection/ECAL_Energy").GetMean(), inFile.Get("HCalEndcapCollection/ECAL_Energy").GetMean()]
#
#  def getEcalEnergy(inFile):
#          return [inFile.Get("ECALBarrel/ECAL_Energy").GetMean(), inFile.Get("ECALEndcap/ECAL_Energy").GetMean()]
#
#  def getHcalEnergy(inFile):
#          return [inFile.Get("HCALBarrel/ECAL_Energy").GetMean(), inFile.Get("HCALEndcap/ECAL_Energy").GetMean()]
#
def getRAWEcalEnergy(inFile):
	return [0,0];

def getRAWHcalEnergy(inFile):
	return [0,0];

def getEcalEnergy(inFile):
	return [0,0];

def getHcalEnergy(inFile):
	return [0,0];


def getEfficiency(inFile):
	numerator = inFile.Get("photonEfficiency/matchedMC_vs_theta").GetEntries()
	denominator = inFile.Get("photonEfficiency/findableMC_vs_theta").GetEntries()
	return numerator/denominator

def getThetaEffPlot(inFile):
	numerator = inFile.Get("photonEfficiency/matchedMC_vs_theta")
	denominator = inFile.Get("photonEfficiency/findableMC_vs_theta")
	ratioPlot = numerator.Clone("thetaEff")
	ratioPlot.Divide(denominator)
	return ratioPlot


def getFitResult(inFile):
	dEHist = inFile.Get("photonEfficiency/dE_matched")
	f1 = TF1("f1", "gaus", -50, 50);
	dEHist.Fit("f1", "Rq");
	print ("FIT: inFile: %s, par1: %f #pm %f, par2: %f #pm %f" % (inFile.GetName(), f1.GetParameter(1), f1.GetParError(1), f1.GetParameter(2), f1.GetParError(2)))
	resolution = 0
	resolutionError = 0
	getRMS90Resolution(inFile.Get("photonEfficiency/PFO_E"), resolution, resolutionError)
	print ("RMS90: res: %f #pm %f" % (resolution, resolutionError))
	return [f1.GetParameter(1),f1.GetParameter(2)]

def createGraphs(inArr, index):
	uniqueEnergy = []
	uniqueTheta = []
	for arrElement in inArr:
		energy = arrElement[0]
		theta = arrElement[1]
		if not energy in uniqueEnergy:
			uniqueEnergy.append(energy)
		if not theta in uniqueTheta:
			uniqueTheta.append(theta)
	#  print (uniqueEnergy)
	#  print (uniqueTheta)

	energyGraph = []
	thetaGraph = []
	for i in range(0,len(uniqueEnergy)):
		tmpGr = TGraph(len(uniqueTheta))
		tmpGr.SetTitle("Energy: %i" % (int(uniqueEnergy[i])))
		styleGraph(tmpGr,1)
		thetaGraph.append(tmpGr)
	for i in range(0,len(uniqueTheta)):
		tmpGr = TGraph(len(uniqueEnergy))
		tmpGr.SetTitle("Theta: %i" % (int(uniqueTheta[i])))
		styleGraph(tmpGr,1)
		energyGraph.append(tmpGr)

	for kk in range(0,len(uniqueEnergy)):
		targetEnergy = uniqueEnergy[kk]
		targetGraph = thetaGraph[kk]
		pointCounter = 0
		
		for arrElement in inArr:
			energy = arrElement[0]
			theta = arrElement[1]
			if (energy==targetEnergy):
				targetGraph.SetPoint(pointCounter,theta,arrElement[index])
				#  print ("targetEnergy: %f, theta: %f, mean: %f" % (targetEnergy,theta,arrElement[index]))
				pointCounter = pointCounter + 1

	for kk in range(0,len(uniqueTheta)):
		targetTheta = uniqueTheta[kk]
		targetGraph = energyGraph[kk]
		pointCounter = 0
		
		for arrElement in inArr:
			energy = arrElement[0]
			theta = arrElement[1]
			if (theta==targetTheta):
				targetGraph.SetPoint(pointCounter,energy,arrElement[index])
				#  print ("targetTheta: %f, energy: %f, mean: %f" % (targetTheta,energy,arrElement[index]))
				pointCounter = pointCounter + 1
	return [energyGraph, thetaGraph]


if __name__ == "__main__":

	if (len(sys.argv)!=2):
		print ("[ERROR]\tSpecify input directory!")
		sys.exit()

	fileList = getFiles(sys.argv[1],"*.root")
	if (len(fileList)==0):
		print ("[ERROR]\tNo input files found... terminate!")
		sys.exit()
	
	addPathToFileName(sys.argv[1], fileList)
	print (("[INFO]\t%i file found. Processing...")%(len(fileList)))
	gStyle.SetOptStat(0)


	globalArr = []
	for iFile in fileList:
		print ("Read file: %s" % (iFile))
		myFile = ROOT.TFile.Open(iFile,"read")
		truthInfo = getTruthInfo(myFile)
		eff = getEfficiency(myFile)
		fitResult = getFitResult(myFile)
		EcalRAWEnergy = getRAWEcalEnergy(myFile)[0] + getRAWEcalEnergy(myFile)[1]
		HcalRAWEnergy = getRAWHcalEnergy(myFile)[0] + getRAWHcalEnergy(myFile)[1]

		EcalEnergy = getEcalEnergy(myFile)[0] + getEcalEnergy(myFile)[1]
		HcalEnergy = getHcalEnergy(myFile)[0] + getHcalEnergy(myFile)[1]

		#  print ("[INFO]\tFile %s" % (iFile))
		#  print ("[INFO]\t Truth Energy: %f, theta: %f" % (truthInfo[0],truthInfo[1]))
		#  print ("[INFO]\t Efficiency: %f" % (eff))
		#  print ("[INFO]\t Fit info, mean: %f, sigma: %f" % (fitResult[0],fitResult[1]))

		tmpArr = []
		tmpArr.append(round(truthInfo[0]))
		tmpArr.append(round(truthInfo[1]))
		#  tmpArr.append(int(truthInfo[0]))
		#  tmpArr.append(int(truthInfo[1]))
		tmpArr.append(100.0*fitResult[0]/truthInfo[0])
		#  tmpArr.append(fitResult[0])
		tmpArr.append(100.0*fitResult[1]/truthInfo[0])
		tmpArr.append(eff)
		tmpArr.append(EcalRAWEnergy)
		tmpArr.append(HcalRAWEnergy)
		tmpArr.append(getRAWEcalEnergy(myFile)[0])
		tmpArr.append(getRAWHcalEnergy(myFile)[0])
		tmpArr.append(EcalEnergy)
		tmpArr.append(EcalEnergy+HcalEnergy)
		tmpArr.append(EcalRAWEnergy+HcalRAWEnergy)
		tmpArr.append((fitResult[0]+truthInfo[0]))
		#  tmpArr.append(getThetaEffPlot(myFile))
		#  print (tmpArr)
		globalArr.append(tmpArr)

		#  histToDraw.append([myFile.Get("ECALBarrel/ECAL_EnergyPerLayers"),myFile.Get("ECALEndcap/ECAL_EnergyPerLayers")])
		


	outFileNamePrefix = ["energyScaleBias","energyResolution","efficiency","EcalRawEnergy","HcalRawEnergy","BarrelEcalRawEnergy","BarrelHcalRawEnergy","EcalEnergy","EcalHcalEnergy","TotalRawEnergy","PFOEnergy","thetaEfficiency"]
	yAxisTitle = ["(E_{PFO}-E_{truth})/E_{truth} [%]", "#sigma(E_{PFO})/E_{truth} [%]", "Efficiency","ECAL Raw Energy [GeV]","HCAL Raw Energy [GeV]","Barrel ECAL Raw Energy [GeV]","Barrel HCAL Raw Energy [GeV]","ECAL Energy [GeV]","ECAL+HCAL Energy [GeV]", "ECAL+HCAL Raw Energy [GeV]","PFO Energy [GeV]", "Efficiency"]
	xAxisTitle = ["Energy [GeV]","Theta"]
	postFix = ["energy","theta"]
	yAxisScale = [ [-1.0,3.0],[0.5,6.0],[0.7,1.1],[0.28,0.3],[0,0.04], [0,0.4],[0,0.04], [9.5,10.2],[9.5,10.2], [0.27,0.29], [0.5,51.0], [0.5,1.1] ]
	legendPosition = [ [0.6,0.75,0.85,0.9],[0.6,0.65,0.85,0.9],[0.6,0.15,0.85,0.4],[0.6,0.65,0.85,0.9],[0.6,0.15,0.85,0.4],[0.6,0.15,0.85,0.4],[0.6,0.15,0.85,0.4],[0.6,0.15,0.85,0.4],[0.6,0.15,0.85,0.4],[0.6,0.15,0.85,0.4], [0.6,0.15,0.85,0.4], [0.6,0.15,0.85,0.4] ]

	outRootDir = "out/"
	if not os.path.exists(outRootDir):
		os.makedirs(outRootDir)
	outRootFile = ROOT.TFile.Open(outRootDir+"outFile.root","RECREATE")

	#  for jj in range(1,2):
	#          for ii in {0,3,4,5,6,7,8,9,10}:
	for jj in range(0,1):
		for ii in {0,1,2,10}:
		#  for ii in range(0,9):
			graphs = createGraphs(globalArr,ii+2) # +2 is needed since two first elements in the array correspond to the truth information
			c1 = TCanvas( 'c1', 'A Simple Graph Example', 0, 0, 800, 600 )
			c1.cd()
			leg = TLegend(legendPosition[ii][0],legendPosition[ii][1],legendPosition[ii][2],legendPosition[ii][3])
			leg.SetBorderSize(0)
			leg.SetTextSize(0.03)

			for kk in range(0,len(graphs[jj])):
				iGr = graphs[jj][kk]
				styleGraph(iGr,kk+1)
				leg.AddEntry(iGr,iGr.GetTitle(),"p")
				if kk==0:
					iGr.SetTitle("")
					iGr.GetYaxis().SetTitle(yAxisTitle[ii])
					iGr.GetXaxis().SetTitle(xAxisTitle[jj])
					iGr.GetYaxis().SetTitleOffset(1.3)
					iGr.Draw("AP")
					iGr.GetYaxis().SetRangeUser(yAxisScale[ii][0],yAxisScale[ii][1])
					iGr.Draw("AP")
				else:
					iGr.Draw("P same")
				
				outRootFile.cd()
				iGr.SetName(postFix[jj]+"_"+outFileNamePrefix[ii])
				iGr.Write()
				
			#  leg.Draw("same")

			c1.SaveAs("%s_%s.png" % (outFileNamePrefix[ii], postFix[jj]))
			#  c1.SaveAs("%s_%s.root" % (outFileNamePrefix[ii], postFix[jj]))

	for iFile in fileList:
		myFile = ROOT.TFile.Open(iFile,"read")
		tmpPlot = getThetaEffPlot(myFile)
		truthInfo = getTruthInfo(myFile)
		tmpPlot.SetTitle(("Energy: %i GeV; Theta; Efficiency" % (truthInfo[0])))
		tmpPlot.SetName(("eff_E_%i" % (truthInfo[0])))
		pfoEPlot = myFile.Get("photonEfficiency/PFO_E")
		pfoEPlot.SetTitle(("Energy: %i GeV; Theta; Efficiency" % (truthInfo[0])))
		pfoEPlot.SetName(("PFO_E_%i" % (truthInfo[0])))
		outRootFile.cd()
		tmpPlot.Write()
		pfoEPlot.Write()

	outRootFile.Close()
#################################################################################

	c2 = TCanvas( 'c2', 'A Simple Graph Example', 0, 0, 800, 600 )
	c2.cd()
	xShift = 0
	yShift = 0.55
	leg2 = TLegend(0.65+xShift, 0.14+yShift, 0.94+xShift, 0.36+yShift)
	leg2.SetBorderSize(0)
	leg2.SetTextSize(0.03)
	#  grList = []
	#  thetaToDraw = [20,70]
	#  momentumToDraw = [10,100]
	thetaToDraw = []
	momentumToDraw = [10]

	i=0
	for iFile in fileList:
		myFile = ROOT.TFile.Open(iFile,"read")
		ROOT.SetOwnership(myFile,False)
		truthInfo = getTruthInfo(myFile)

		thetaInt = int(truthInfo[1])
		momentumInt = int(truthInfo[0])

		if ((not (thetaInt in thetaToDraw)) and (not len(thetaToDraw)==0)):
			continue
		if ((not (momentumInt in momentumToDraw)) and (not len(momentumToDraw)==0) ):
			continue

		iGr_barrel = myFile.Get("ECALBarrel/ECAL_AverageEnergyPerLayers")
		iGr_endcap = myFile.Get("ECALEndcap/ECAL_AverageEnergyPerLayers")
		if truthInfo[1]<45:
			iGr_barrel = iGr_endcap 
			
		iGr_barrel.SetLineColor(i+1)
		iGr_barrel.SetLineWidth(2)
		iGr_endcap.SetLineColor(i+1)
		iGr_endcap.SetLineWidth(2)
		if i==0:
			iGr_barrel.GetYaxis().SetRangeUser(0.0,0.07)
			iGr_barrel.Draw("H")
		else:
			iGr_barrel.Draw("Hsame")
		#  iGr_endcap.Draw("Hsame")
		leg2.AddEntry(iGr_barrel, ( "E: %d, theta: %d" % (truthInfo[0],truthInfo[1]) ) ,"l")
		i = i+1

	#  print ("OUT_OF_LOOP: %f" % (grList[0].GetMean()))

	leg2.Draw()
	c2.Update()
	#  raw_input()
	c2.SaveAs("E_vs_layerNumber.png")
	
