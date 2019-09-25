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
#include <TEfficiency.h>

//LCIO
#include <IOIMPL/LCFactory.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
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
		// vector<EVENT::ReconstructedParticle*> getObjVecFromCollection(const EVENT::LCCollection* inCollection);
		void createHistsFromMap(const map<string,histStruct> inHistStructMap, const string prefix);
		TH1* getHistFromMap(string histID);
		IMPL::ReconstructedParticleImpl* CopyReconstructedParticle (const EVENT::ReconstructedParticle* const pfo_orig );

		// template <class T> vector<T> getObjVecFromCollection(EVENT::LCCollection* inCollection);
		template <class T> vector<T> getObjVecFromCollection(EVENT::LCCollection* inCollection){
			int nElements = inCollection->getNumberOfElements();
			vector<T> outVec;
			for(int j=0; j < nElements; j++) {
				auto part = dynamic_cast<T>(inCollection->getElementAt(j));
				outVec.push_back(part);
			}
			return outVec;
		}
		void DeleteHists();
		void createTH1I(string histName, string histTitle, unsigned int nBins, double leftRange, double rightRange); 
		void createTH1D(string histName, string histTitle, unsigned int nBins, double leftRange, double rightRange);
		void createTEff(string numeratorHistName, string denominatorHistName){
			TEfficiency* tmpTEff = new TEfficiency(*getHistFromMap(numeratorHistName),*getHistFromMap(denominatorHistName));	
			tmpTEff->SetName(numeratorHistName.c_str());
			tmpTEff->SetDirectory(0);
			tEffMap[numeratorHistName] = tmpTEff;
		}

	protected:
		std::map<std::string, TH1*> histMap;
		std::map<std::string, TEfficiency*> tEffMap;
		string outDirName;
		bool debugFlag;
		double get_dPhi(double phi_reco, double phi_truth);
		string className = "objectFill";
};

#endif
