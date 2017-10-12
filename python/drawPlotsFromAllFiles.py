#  Draw plots specific plots for all files in the directory.
#  One TPad can contain many plots from ONE file

#!/bin/python
import glob, os, ROOT, math, sys
from ROOT import TCanvas, TGraph, TLegend, TF1, TH1, TH1F
from ROOT import gROOT, gStyle
import yaml # WARNING use this environment ~/env/setupPyTools27.env to enable yaml package!!!

#  rebinFactor = 2
#  yAxisRange = [0, 600]
#  yAxisRange = []
#  xAxisRange = [0, 110]
#  histNamePrefix = "truthParticle_"
#  histNamePostfix = "_CosTheta"
#  histNamePrefix = "firstEnergeticPartOfType_"
#  histNamePostfix = "Energy"

#  histInfo = [
#  ["onlyOneRecoPFO",1,"one reconstructed PFO"],
#  ["twoOrMoreRecoPFO",2,"2 or more reconstructed PFOs"]
#  ["noAdditionalPFOs",1,"one reconstructed PFO"],
#  ["thereAreAdditionalPFOs",2,"2 or more reconstructed PFOs"]
#  ]

# string format: particleGun_E<>_Theta<>_Phi<>.root
def getPhaseSpacePoint(inStr):
	tmpStrArr = ('.'.join(inStr.split(".")[:-1])).split("_")
	E = float( tmpStrArr[1].split("E")[1] )
	Theta = float( tmpStrArr[2].split("Theta")[1] )
	Phi = float( tmpStrArr[3].split("Phi")[1] )
	return [E, Theta, Phi]

if __name__ == "__main__":

	with open("/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/python/config/drawPlotsFromAllFiles/globalConfig.yml", 'r') as ymlfile:
		tmpCfg = yaml.load(ymlfile)
		pfoTypeMap = tmpCfg["pfo_type_map"]

	with open("/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/python/config/drawPlotsFromAllFiles/plotConfigExample.yml", 'r') as ymlfile:
		globalCfg = yaml.load(ymlfile)
	

	if (len(sys.argv)==1):
		pfoType = "Photon"
	else:
		pfoType = sys.argv[1]
	
	fileCounter = 0

	gStyle.SetOptStat(0)
	os.chdir("./")
	for iFile in glob.glob("*.root"):
		myFile = ROOT.TFile.Open(iFile,"read")
		phaceSpacePoint = getPhaseSpacePoint(iFile)
		print (iFile)
		print (phaceSpacePoint)
		iE = phaceSpacePoint[0]
		iTheta = phaceSpacePoint[1]
		iPhi = phaceSpacePoint[2]
		counter = 0
		for cfgIterator in globalCfg:
			cfg = globalCfg[cfgIterator]
			histNamePrefix = ""
			histNamePostfix = ""
			if ("histNamePrefix" in cfg):
				histNamePrefix = cfg['histNamePrefix']
			if ("histNamePostfix" in cfg):
				histNamePostfix = cfg['histNamePostfix']
			histNameBase = cfg['histNameBase']
			inFileDir = "PandoraPFOs_" + str(pfoTypeMap[pfoType]) + "/"
			if ("customInFileDirOrHistName" in cfg):
				inFileDir = ""

			hists = []
			nEntries = []
			for histIt in range(0,len(histNameBase)):
				histName = inFileDir+histNamePrefix+histNameBase[histIt]+histNamePostfix
				print (histName)

				hist = myFile.Get(histName)
				hist.SetTitle(pfoType+" Gun, E=%i GeV, Theta=%i, Phi=%i" % (iE,iTheta,iPhi))
				nEntries.append(hist.GetEntries())

				if ("histColor" in cfg):
					hist.SetLineColor(cfg['histColor'][histIt])
				else:
					hist.SetLineColor(histIt+1)
				if ("rebinFactor" in cfg):
					hist.Rebin(cfg['rebinFactor'])
				if ("xAxisRange" in cfg):
					hist.GetXaxis().SetRangeUser(cfg["xAxisRange"][0],cfg["xAxisRange"][1])
				if ("yAxisRange" in cfg):
					hist.GetYaxis().SetRangeUser(cfg["yAxisRange"][0],cfg["yAxisRange"][1])

				hists.append(hist)

			totalEntries = sum(nEntries)
			can = TCanvas( 'c1_'+cfgIterator, 'A Simple Graph Example', 0, 0, 800, 600 )
			if ("legPos" in cfg):
				leg = TLegend(cfg["legPos"][0],cfg["legPos"][1],cfg["legPos"][2],cfg["legPos"][3])
			else:
				leg = TLegend(0.25,0.3,0.75,0.55)

			for i in range(0,len(hists)):
				hist = hists[i]
				if ("histLegend" in cfg):
					leg.AddEntry(hist,"%s (%i%%)" % (cfg["histLegend"][i],round(100.0*nEntries[i]/totalEntries)),"l")
				if (i==0):
					hist.Draw()
				else:
					hist.Draw("same")

			if ("histLegend" in cfg):
				leg.Draw("same")
			#  can.SetLogy()
			#  can.SaveAs("E%i_Theta%i_Phi%i" % (iE, iTheta, iPhi) +cfgIterator+".png"  )
			if (counter == 0):
				can.Print("E%i_Theta%i_Phi%i" % (iE, iTheta, iPhi) + ".pdf(","pdf")
			else:
				can.Print("E%i_Theta%i_Phi%i" % (iE, iTheta, iPhi) + ".pdf","pdf")
			counter = counter + 1

			if (fileCounter == 0):
				can.Print(cfgIterator + ".pdf(","pdf")
			elif (fileCounter==(len(glob.glob("*.root"))-1)):
				can.Print(cfgIterator + ".pdf)","pdf")
			else:
				can.Print(cfgIterator + ".pdf","pdf")

		can = TCanvas( 'c2_'+cfgIterator, 'A Simple Graph Example', 0, 0, 800, 600 )
		can.Print("E%i_Theta%i_Phi%i" % (iE, iTheta, iPhi) + ".pdf)","pdf")

		
		fileCounter = fileCounter + 1

	print ("****************************************************************")
	print ("[INFO]\tRUN OVER PFO TYPE: %s" % (pfoType))
