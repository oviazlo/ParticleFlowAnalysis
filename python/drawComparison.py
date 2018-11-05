#  Compare different plots from the same
#  OR
#  compare the same plot from different files

#!/bin/python
import glob, os, ROOT, math, sys
from ROOT import TCanvas, TGraph, TLegend, TF1, TH1, TH1F, TLatex, TText
from ROOT import gROOT, gStyle
from array import array
from math import tan
import yaml # WARNING use this environment ~/env/setupPyTools27.env to enable yaml package!!!
ROOT.gROOT.LoadMacro("CLIC_style/CLICdpStyle/rootstyle/CLICdpStyle.C+")
ROOT.CLICdpStyle()

from classesForDrawComparison import histList


#  def processSedMe(globalCfg):
#      sedMePairs = []
#      for cfgIterator in globalCfg:
#          print ("")
#          cfg = globalCfg[cfgIterator]
#          strToPrintInBeginning = "********** " + cfgIterator + " **********"
#          print (strToPrintInBeginning )
#          #  print ("")
#          #  print (cfg)
#          #  print ("")
#          yamlKey = 'fileName'
#          if (yamlKey in cfg):
#                  dummy = cfg[yamlKey]
#                  print (yamlKey + ":")
#                  print (dummy)
#                  print ("")
#
#          yamlKey = 'legTitle'
#          if (yamlKey in cfg):
#                  dummy = cfg[yamlKey]
#                  print (yamlKey + ":")
#                  print (dummy)
#                  print ("")
#
#          for innerIter in cfg:
#              if "sedMe" in innerIter:
#                  sedMePairs.append([innerIter, cfg[innerIter]])
#
#
#          print ("*"*len(strToPrintInBeginning))
#
#      print ("")
#      print ("sedMePairs: ")
#      print (sedMePairs)
#      sys.exit()

def main(yamlFile):

        print ("")
	gStyle.SetOptStat(0)

	with open(yamlFile, 'r') as ymlfile:
		globalCfg = yaml.load(ymlfile)
                print ("[info]\t Read yaml file: \n\t\t%s" % (yamlFile))

        #  processSedMe(globalCfg)

        # copy default settings to every plot entry:
        if (globalCfg.get("default") is not None):
            defaultCfg = globalCfg.get("default")
            for cfgIterator in globalCfg:
                cfg = globalCfg[cfgIterator]
                if (cfg == defaultCfg):
                    continue
                for cfgIt in defaultCfg:
                    if (cfg.get(cfgIt) is None):
                        cfg[cfgIt] = defaultCfg.get(cfgIt)
            globalCfg.pop("default")

        #  DEBUGGING: print dict
        #  for cfgIterator in globalCfg:
        #          cfg = globalCfg[cfgIterator]
        #          if (cfg == defaultCfg):
        #              continue
        #          print (cfgIterator)
        #          print (cfg)

	for cfgIterator in globalCfg:
		cfg = globalCfg[cfgIterator]
                print('[INFO]\tProcess plot: %s' % (cfgIterator))
               
                classHelper = histList()
                hists = classHelper.getProcessedHists(cfg)
                leg = classHelper.getLegend(hists, cfg)
                canvas = classHelper.getCanvas(cfg)

		for i in range(0,len(hists)):
                    if (i==0):
                        drawOption = ""
                        if ("drawOption" in cfg):
                            drawOption = cfg['drawOption']
                            sameDrawOption = drawOption
                        
                        print ("[info]\t drawOption: \n\t\t%s" % (drawOption))
                        hists[i].Draw(drawOption)
                        ttext = TText(0.05,hists[i].GetMaximum()*1.03,"CLD")
                    else:
                        sameDrawOption = ""
                        if ("drawOption" in cfg):
                            sameDrawOption = cfg['drawOption']
                        if ("sameDrawOption" in cfg):
                            sameDrawOption = cfg['sameDrawOption']
                        print ("[info]\t sameDrawOption: \n\t\t%s" % (sameDrawOption))
                        hists[i].Draw(sameDrawOption+"same")


                ttext.Draw("same")
		if (leg):
			leg.Draw("same")
		if ("additionalLabel" in cfg):
			l = ROOT.TLatex()
			l.SetNDC()
			l.SetTextFont(42)
			l.SetTextColor(ROOT.kBlack)
			l.SetTextSize(0.045)
			label = cfg['additionalLabel']
			x = 0.41
			y = 0.42
                        if ("legPos" in cfg):
                            xShift = 0.01
                            yShift = 0.0
                            if ("additionalLabel_xShift" in cfg):
                                xShift = cfg['additionalLabel_xShift']
                            x = cfg["legPos"][0]+xShift
                            y = cfg["legPos"][3]+yShift
			l.DrawLatex(x,y,label);

                print ("")
                ROOT.gPad.SetPhi(150);
                outHistName = cfgIterator
                if ("histNamePrefix" in cfg):
                    outHistName = cfg['histNamePrefix'] + outHistName
                if ("savePictDir" in cfg):
                    pythonFileDir = os.path.dirname(os.path.realpath(__file__))
                    fullPath_savePictDir = pythonFileDir + '/pictures/' + cfg['savePictDir'] + '/' + outHistName
                    if not os.path.exists(fullPath_savePictDir):
                        os.makedirs(fullPath_savePictDir)
                    outHistName = cfg['savePictDir'] + '/' + outHistName
		canvas.SaveAs("pictures/"+outHistName+".png")
		canvas.SaveAs("pictures/"+outHistName+".pdf")
		canvas.SaveAs("pictures/"+outHistName+".C")



if __name__ == "__main__":
    if len(sys.argv)==2:
        print ("[INFO]\tRead config from: {0}".format(sys.argv[1]))
        main(sys.argv[1])
    else:
        print ("[ERROR]\tPass config file to read!!! Terminating...")
        sys.exit()
