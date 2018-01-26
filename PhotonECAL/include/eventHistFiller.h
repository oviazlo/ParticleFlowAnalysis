/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/include/eventHistFiller.h
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	05th Dec 2017
 * 	Last Update:	05th Dec 2017
 */
#ifndef EVENTHISTFILLER_H
#define EVENTHISTFILLER_H

#include <objectFill.h>
#include <globalConfig.h>
#include "truthCondition.h"

struct PFOMergeSettings{
	unsigned int pfoTypeToMerge;
	double thetaCone;
	double phiCone;
	double phiConeMomentumDep = 0.0;
};

// struct mergedEnergyContainer{
//         double mergedEnergy;
//         mergedEnergyContainer(): mergedEnergy(0.0){}
//         mergedEnergyContainer(double partE): mergedEnergy(partE){}
//         vector<PFOMergeSettings> > PFOmergeMap;
//         bool tryToMergeParticle(EVENT::ReconstructedParticle* inPart){
//
//         }
// };

class eventHistFiller : public objectFill{
	public:
		eventHistFiller(string _outDirName, string _PFOCollectionName = "PandoraPFO") : objectFill(_outDirName) {PFOCollectionName = _PFOCollectionName; mergeTag = "nominal";}
		~eventHistFiller(){}

		int init();
		int fillEvent(const EVENT::LCEvent*);
		int writeToFile(TFile* outFile);

		void setClusterMerging(string _mergeTag);

	private:
		unsigned int nSelectecTruthParticles;
		string PFOCollectionName;
		EVENT::LCCollection *PFOCollection;
		int checkPfoType(vector<unsigned int> inVec);
		// map<string,histStruct> singleParticleHistStructMap;
		// vector<string> effType = {"nominal","photonMerge","photonAndNeutralMerge","photonAndNeutralLooseMerge","bremRecovery"};
		vector<string> effType = {"nominal","photonMerge","photonAndNeutralMerge","photonAndNeutralLooseMerge","photonMergeMomentumDep"};
		// vector<string> effType = {"nominal","photonMerge"};
		map<string,vector<PFOMergeSettings> > PFOmergeMap;
		string mergeTag;
		double mergedEnergy;
};
#endif // EVENTHISTFILLER_H
