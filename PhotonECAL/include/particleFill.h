#ifndef particleFill_H
#define particleFill_H

#include <objectFill.h>
#include <globalConfig.h>

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
		int fillEvent(const EVENT::LCEvent*);
		double getMeanEnergy();
		double getMeanTheta();
		void setCollectionName(const string _collectionName){collectionName = _collectionName;}
		void setReconstructedParticleType(const int _partTypeToSelect){partTypeToSelect.push_back(_partTypeToSelect);}
		void setReconstructedParticleType(const vector<int> _partTypeToSelect){partTypeToSelect = _partTypeToSelect;}
		void updateRootDirName(const string inStr){outDirName = inStr;}
		void setDPhiMergeValue(const double inVal){dPhiMergeValue = inVal;}
		int writeToFile(TFile* outFile);

	private:
		void initHistStructs();
		void createSingleRecoParticleClustersHists(const string prefix);
		void createSingleParticleHists(const string prefix);
		void createTwoParticleCorrelationHists(const string prefix);
		void createHistsFromMap(const map<string,histStruct> inHistStructMap, const string prefix);
		int fillHist(const double inVal, const string baseString, const string prefix);
		int fillHist(const double inVal1, const double inVal2, const string baseString, const string prefix);
		void dumpTruthPart(const EVENT::MCParticle* part, const int counter = 0);
		void dumpReconstructedPart(const EVENT::ReconstructedParticle* part, const int counter = 0);
		int fillPart (const EVENT::MCParticle* inPart, const string prefix="truthParticle_");
		int fillPart (const EVENT::ReconstructedParticle* inPart, const string prefix);
		void fillRecoPhoton(const EVENT::ReconstructedParticle* inPart, const EVENT::MCParticle* mcPart, const string prefix);
		int fillParticleCorrelations (const EVENT::ReconstructedParticle* inPart1, const EVENT::ReconstructedParticle* inPart2, const string prefix);
		int fillParticleCorrelations (const EVENT::ReconstructedParticle* inPart1, const EVENT::MCParticle* inPart2, const string prefix);
		int fillClusterInfo (const EVENT::ReconstructedParticle* inPart, const string prefix);
		// FIXME hardcoded default name of the truth collection
		vector<EVENT::MCParticle*> getTruthMCParticlesFromCollection(const EVENT::LCEvent* event, const string truthCollectionName = "MCParticlesSkimmed");

		EVENT::LCCollection *collection;
		string collectionName;
		vector<int> partTypeToSelect;
		map<string,histStruct> singleParticleHistStructMap;
		map<string,histStruct> singleRecoParticleClustersHistStructMap;
		map<string,histStruct> twoParticleCorrelationHistStructMap;
		double dPhiMergeValue;
		vector<bool> boolVecDefaultFalse;
		map<string,int> intMap;
		vector< pair<int, double> > PFOtypeAndEnergyVec;
		void clasiffyPFO(EVENT::ReconstructedParticle* inPFO);



};

#endif
