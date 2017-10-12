#ifndef truthParticleSelector_H
#define truthParticleSelector_H

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

#include <energyFill.h>
#include <particleFill.h>
#include <photonEffCalculator.h>
#include <serviceFunctions.h>

class truthParticleSelector{
	public:
		truthParticleSelector(string _mcTruthCollection);
		~truthParticleSelector();
		
		void setDebugFlag(bool inFlag){debugFlag = inFlag;}
		void setEfficiencyCollection(string _effCollection){effCollection=_effCollection;}
		void setEfficiencyPFOType(int pfoType){efficiencyPFOType=pfoType;}
		void setParticleFillCollections(vector<string> _particleFillCollections){particleFillCollections = _particleFillCollections;}
		void setPFOTypes(vector<vector<int> > inVec){PFOTypes = inVec;}
		void setEnergyFillCollections(vector<string> _energyFillCollections){energyFillCollections = _energyFillCollections;}
		void init(); 
	
		void setEnergyRange(double min, double max){energyRange = make_pair(min,max);}
		void setThetaRange(double min, double max){thetaRange = make_pair(min,max);}
		void setPhiRange(double min, double max){phiRange = make_pair(min,max);}
		string getPostFixString(){return "E"+DoubToStr((energyRange.first+energyRange.second)/2.0)+"_Theta"+DoubToStr((thetaRange.first+thetaRange.second)/2.0)+"_Phi"+DoubToStr((phiRange.first+phiRange.second)/2.0);  }

		bool selectEvent(EVENT::LCEvent*);

		void writeToFile(TFile *outFile); 
		void setDiscardFSREvents(bool inBool){discardFSREvents = inBool;}
		void setDPhiMergeValue(double inVal){dPhiMergeValue = inVal;}

	private:
		string mcTruthCollection;
		string effCollection;
		bool debugFlag;
		map<string, objectFill*> objFillMap;
		pair<double,double> energyRange;
		pair<double,double> thetaRange;
		pair<double,double> phiRange;
		vector<string> particleFillCollections;
		vector<string> energyFillCollections;
		vector<vector<int> > PFOTypes;
		int efficiencyPFOType;
		bool discardFSREvents;
		double dPhiMergeValue;
};

#endif
