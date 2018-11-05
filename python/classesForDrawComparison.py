import ROOT

class nLegendCaptions(Exception):
    pass

class objectPosition(Exception):
    pass

################################################################################
################################################################################
class histList:

################################################################################
    def __init__(self):
        self.histColor = [ROOT.kBlack, ROOT.kRed-7, ROOT.kBlue, ROOT.kGreen+2, ROOT.kCyan+1, ROOT.kRed+2, ROOT.kOrange, ROOT.kViolet+2, ROOT.kGray]
        self.markerStyle = [ROOT.kOpenCircle, ROOT.kOpenSquare, ROOT.kOpenTriangleUp, ROOT.kOpenDiamond, ROOT.kOpenCross, ROOT.kOpenCircle, ROOT.kOpenCircle]
        self.leg = None
        self.canvas = None

################################################################################
    def getCanvas(self, cfg):
        if (self.canvas == None):
            c1 = None
            if ("canPos" in cfg):
                if (len(cfg["canPos"])!=4):
                    raise objectPosition("Canvas position should have 4 coordinates!")
                else:
                    c1 = ROOT.TCanvas( 'c1', 'A Simple Graph Example',  cfg["canPos"][0],cfg["canPos"][1],cfg["canPos"][2],cfg["canPos"][3])
                cfg.pop('canPos')
            else:
                c1 = ROOT.TCanvas( 'c1', 'A Simple Graph Example', 0, 0, 800, 700 )

            if ("gridX" in cfg):
                c1.SetGridx(cfg['gridX'])
                cfg.pop('gridX')
            if ("gridY" in cfg):
                c1.SetGridy(cfg['gridY'])
                cfg.pop('gridY')
            if ("logX" in cfg):
                c1.cd(1).SetLogx(cfg['logX'])
                cfg.pop('logX')
            if ("logY" in cfg):
                c1.cd(1).SetLogy(cfg['logY'])
                cfg.pop('logY')
            self.canvas = c1

        return self.canvas


################################################################################
    def getLegend(self, hists, cfg):
        if self.leg is None:
            leg = None
            legTitle = []
            if ("legTitle" in cfg):
                legTitle = cfg['legTitle']
                cfg.pop('legTitle')
            if ("legPos" in cfg):
                leg = ROOT.TLegend(cfg["legPos"][0],cfg["legPos"][1],cfg["legPos"][2],cfg["legPos"][3])
                cfg.pop('legPos')
            if ("legTextSize" in cfg):
                leg.SetTextSize(cfg["legTextSize"])
                cfg.pop('legTextSize')
            if ("legDrawOption" in cfg):
                legDrawOption = cfg["legDrawOption"]
                cfg.pop('legDrawOption')
            else:
                legDrawOption = "lp"

            try:
                for i in range(len(hists)):
                    if ("legCaptionMode" in cfg):
                        if (cfg["legCaptionMode"]=="fraction"):
                            totalEntries = sum( list(x.GetEntries() for x in hists) )
                            leg.AddEntry(hists[i],"%s (%i%%)" % (legTitle[i], round(100.0*hists[i].GetEntries()/totalEntries)),legDrawOption)
                        elif (cfg["legCaptionMode"]=="entries"):
                            leg.AddEntry(hists[i],"%s (%i)" % (legTitle[i], hists[i].GetEntries()),legDrawOption)
                        elif (cfg["legCaptionMode"]=="mean"):
                            leg.AddEntry(hists[i],"%s (mean: %.2f)" % (legTitle[i], hists[i].GetMean()),legDrawOption)
                        elif (cfg["legCaptionMode"]=="RMS"):
                            leg.AddEntry(hists[i],"%s (RMS: %.2f)" % (legTitle[i], hists[i].GetRMS()),legDrawOption)
                        cfg.pop('legCaptionMode')
                    else:
                        leg.AddEntry(hists[i],"%s" % (legTitle[i]),legDrawOption)
            except IndexError:
                raise nLegendCaptions("number of captions in the legend IS NOT EQUAL to number of histograms to draw")

            self.leg = leg
        return self.leg

################################################################################
    def getProcessedHists(self, cfg):
        if ("histColor" in cfg):
            self.histColor = cfg['histColor'] + self.histColor
            cfg.pop('histColor')
        if ("markerStyle" in cfg):
            self.markerStyle = cfg['markerStyle'] + self.markerStyle
            cfg.pop('markerStyle')

        if ("fileName" in cfg):
            self.fileName = cfg['fileName']
            cfg.pop('fileName')
        if ("dirPrefix" in cfg):
            dirPrefix = cfg['dirPrefix']
            cfg.pop('dirPrefix')
        if ("histName" in cfg):
            histName = cfg['histName']
            cfg.pop('histName')

        self.hists = []
        for i in range(0,len(self.fileName)):
            myFile = ROOT.TFile.Open(dirPrefix+self.fileName[i],"read")
            #  print ("[info]\t Open file: \n\t\t%s" % (dirPrefix+self.fileName[i]))
            ROOT.SetOwnership(myFile,False)
            for k in range(0,len(histName)):
                #  print ("[info]\t Getting hist: \n\t\t%s" % (histName[k]))
                if not myFile.Get(histName[k]):
                    print ('[ERROR]\t Hist "%s" is not found in file "%s"!' % (histName[k],self.fileName[i]))
                    if ("skipMissingHists" in cfg):
                        cfg.pop('skipMissingHists')
                        continue
                    else:
                        print ('[FATAL]\t Terminating...')
                        sys.exit()
                hist = myFile.Get(histName[k])

                hist.SetLineColor(self.histColor[0])
                hist.SetMarkerColor(self.histColor[0])
                self.histColor.pop(0)

                hist.SetMarkerStyle(self.markerStyle[0])
                self.markerStyle.pop(0)

                hist.SetLineWidth(2)

                if ("lineWidth" in cfg):
                    hist.SetLineWidth(cfg["lineWidth"])
                    cfg.pop('lineWidth')
                if ("markerSize" in cfg):
                    hist.SetMarkerSize(cfg["markerSize"])
                    cfg.pop('markerSize')
                if ("histTitle" in cfg):
                    hist.SetTitle(cfg['histTitle'])
                    cfg.pop('histTitle')
                if ("rebinFactor" in cfg):
                    hist.Rebin(cfg['rebinFactor'])
                    cfg.pop('rebinFactor')
                if ("scaleFactor" in cfg):
                    hist.Scale(cfg['scaleFactor'])
                    cfg.pop('scaleFactor')
                if ("normalize" in cfg):
                    hist.Scale(1.0/hist.Integral())
                    cfg.pop('normalize')

                if ("sortXaxisInGraph" in cfg):
                    if (cfg["sortXaxisInGraph"]==True):
                        histPoints = []
                        for j in range(0,hist.GetN()):
                            histPoints.append([hist.GetX()[j],hist.GetY()[j],hist.GetEX()[j],hist.GetEY()[j]])
                        histPoints = sorted(histPoints, key=lambda x: x[0])
                        for j in range(0,hist.GetN()):
                            hist.SetPoint(j+1,histPoints[j][0],histPoints[j][1])
                            hist.SetPointError(j+1,histPoints[j][2],histPoints[j][3])
                    cfg.pop('sortXaxisInGraph')

                if ("xAxisRange" in cfg):
                    hist.GetXaxis().SetRangeUser(cfg['xAxisRange'][0],cfg['xAxisRange'][1])
                    cfg.pop('xAxisRange')
                if ("yAxisRange" in cfg):
                    hist.GetYaxis().SetRangeUser(cfg['yAxisRange'][0],cfg['yAxisRange'][1])
                    cfg.pop('yAxisRange')
                if ("xAxisLabel" in cfg):
                    hist.GetXaxis().SetTitle(cfg['xAxisLabel'])
                    cfg.pop('xAxisLabel')
                if ("yAxisLabel" in cfg):
                    hist.GetYaxis().SetTitle(cfg['yAxisLabel'])
                    cfg.pop('yAxisLabel')
                if ("xAxisTitleOffset" in cfg):
                    hist.GetXaxis().SetTitleOffset(float(cfg['xAxisTitleOffset']))
                    cfg.pop('xAxisTitleOffset')
                if ("yAxisTitleOffset" in cfg):
                    hist.GetYaxis().SetTitleOffset(float(cfg['yAxisTitleOffset']))
                    cfg.pop('yAxisTitleOffset')
                if ("yAxisTitleSize" in cfg):
                    hist.GetYaxis().SetTitleSize(float(cfg['yAxisTitleSize']))
                    cfg.pop('yAxisTitleSize')
                if ("xAxisTitleSize" in cfg):
                    hist.GetXaxis().SetTitleSize(float(cfg['xAxisTitleSize']))
                    cfg.pop('xAxisTitleSize')

                self.hists.append(hist)         
        return self.hists

