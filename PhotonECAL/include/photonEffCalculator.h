#ifndef photonEffCalculator_H
#define photonEffCalculator_H

#include <objectFill.h>

class photonEffCalculator : public objectFill{

	public:
		photonEffCalculator(string _outDirName) : objectFill(_outDirName) {PFOPartType = 0;}
		~photonEffCalculator(){}

		int init();
		int fillEvent(const EVENT::LCEvent*);
		void setPFOCollection(const string _collectionName){PFOCollectionName = _collectionName;}
		void setPFOType(const int pfoType){PFOPartType = pfoType;}
		void setMCTruthCollection(const string _collectionName){MCTruthCollectionName = _collectionName;}
		void setDPhiMergeValue(const double inVal){dPhiMergeValue = inVal;}
		int writeToFile(TFile* outFile);
		void setEfficiencyOneClusterRequirement(bool inVal){onlyOneRecoClusterPerEvent = inVal;}


	private:
		string PFOCollectionName;
		string MCTruthCollectionName;
		EVENT::LCCollection *PFOCollection;
		EVENT::LCCollection *MCTruthCollection;
		EVENT::ReconstructedParticle* getMatchedPFO(const EVENT::MCParticle* inMCPart, const vector<EVENT::ReconstructedParticle*> findablePFO);
		int PFOPartType;
		double dPhiMergeValue;
		bool onlyOneRecoClusterPerEvent;
};
#endif
