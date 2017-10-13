#ifndef particleFill_H
#define particleFill_H

#include <objectFill.h>

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

class particleFill : public objectFill{

	public:
		particleFill(string _outDirName) : objectFill(_outDirName) {}
		~particleFill(){}

		int init();
		int fillEvent(EVENT::LCEvent*);
		double getMeanEnergy();
		double getMeanTheta();
		void setCollectionName(string _collectionName){collectionName = _collectionName;}
		void setReconstructedParticleType(int _partTypeToSelect){partTypeToSelect.push_back(_partTypeToSelect);}
		void setReconstructedParticleType(vector<int> _partTypeToSelect){partTypeToSelect = _partTypeToSelect;}
		void updateRootDirName(string inStr){outDirName = inStr;}
		void setDPhiMergeValue(double inVal){dPhiMergeValue = inVal;}
		int writeToFile(TFile* outFile);

	private:
		void initHistStructs();
		void createSingleRecoParticleClustersHists(string prefix);
		void createSingleParticleHists(string prefix);
		void createTwoParticleCorrelationHists(string prefix);
		void createHistsFromMap(map<string,histStruct> inHistStructMap, string prefix);
		int fillHist(double inVal, string baseString, string prefix);
		int fillHist(double inVal1, double inVal2, string baseString, string prefix);
		void dumpTruthPart(EVENT::MCParticle* part, int counter = 0);
		void dumpReconstructedPart(EVENT::ReconstructedParticle* part, int counter = 0);
		int fillPart (EVENT::MCParticle* inPart, string prefix="truthParticle_");
		int fillPart (EVENT::ReconstructedParticle* inPart, string prefix);
		int fillParticleCorrelations (EVENT::ReconstructedParticle* inPart1, EVENT::ReconstructedParticle* inPart2, string prefix);
		int fillClusterInfo (EVENT::ReconstructedParticle* inPart, string prefix);
		// FIXME hardcoded default name of the truth collection
		vector<EVENT::MCParticle*> getTruthMCParticlesFromCollection(EVENT::LCEvent* event, string truthCollectionName = "MCParticlesSkimmed");

		EVENT::LCCollection *collection;
		string collectionName;
		vector<int> partTypeToSelect;
		map<string,histStruct> singleParticleHistStructMap;
		map<string,histStruct> singleRecoParticleClustersHistStructMap;
		map<string,histStruct> twoParticleCorrelationHistStructMap;
		double dPhiMergeValue;
		vector<bool> boolVecDefaultFalse;

};

#endif
