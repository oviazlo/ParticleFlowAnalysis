#ifndef objectFill_H
#define objectFill_H
 
//ROOT
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TList.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TGaxis.h>
#include <TMath.h>
#include <TVector3.h>
#include <TVectorD.h>

//LCIO
#include <IOIMPL/LCFactory.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Cluster.h>
#include <Exceptions.h>

//STD
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>

#include <serviceFunctions.h>
// #include <globalConfig.h>

struct histStruct{
	string title;
	unsigned int nBins;
	double xLow;
	double xHigh;
	double yLow;
	double yHigh;
	unsigned int ynBins;
	string histType;
	histStruct(){}
	histStruct( string _title, unsigned int _nBins, double _xLow, double _xHigh, string _histType = "TH1D", unsigned int _ynBins = 0, double _yLow = 0.0, double _yHigh = 0.0 ) : title(_title), nBins(_nBins), xLow(_xLow), xHigh(_xHigh), histType(_histType), ynBins(_ynBins), yLow(_yLow), yHigh(_yHigh) {}
};


class objectFill{
	public:
		objectFill(const string _outDirName);
		~objectFill();
		virtual int writeToFile(TFile* outFile);
		virtual int init(){return 0;}
		virtual int fillEvent(const EVENT::LCEvent*){return 0;}
		void setDebugFlag(const bool inFlag){debugFlag = inFlag;}
		// template <class T> vector<T> getObjVecFromCollection(EVENT::LCCollection* inCollection);
		vector<EVENT::ReconstructedParticle*> getObjVecFromCollection(const EVENT::LCCollection* inCollection);
		void createHistsFromMap(const map<string,histStruct> inHistStructMap, const string prefix);

	protected:
		std::map<std::string, TH1*> histMap;
		string outDirName;
		bool debugFlag;
};

#endif
