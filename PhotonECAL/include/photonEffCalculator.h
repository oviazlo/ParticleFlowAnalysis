#ifndef photonEffCalculator_H
#define photonEffCalculator_H

#include <objectFill.h>

class photonEffCalculator : public objectFill{

	public:
		photonEffCalculator(string _outDirName) : objectFill(_outDirName) {PFOPartType = 0;}
		~photonEffCalculator(){}

		int init();
		int fillEvent(EVENT::LCEvent*);
		void setPFOCollection(string _collectionName){PFOCollectionName = _collectionName;}
		void setPFOType(int pfoType){PFOPartType = pfoType;}
		void setMCTruthCollection(string _collectionName){MCTruthCollectionName = _collectionName;}
		void setDPhiMergeValue(double inVal){dPhiMergeValue = inVal;}
		int writeToFile(TFile* outFile);

	private:
		string PFOCollectionName;
		string MCTruthCollectionName;
		EVENT::LCCollection *PFOCollection;
		EVENT::LCCollection *MCTruthCollection;
		EVENT::ReconstructedParticle* getMatchedPFO(EVENT::MCParticle* inMCPart, vector<EVENT::ReconstructedParticle*> findablePFO);
		int PFOPartType;
		double dPhiMergeValue;
};
#endif
