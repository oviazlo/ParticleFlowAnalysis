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

class eventHistFiller : public objectFill{
	public:
		eventHistFiller(string _outDirName, string _PFOCollectionName = "PandoraPFO") : objectFill(_outDirName) {PFOCollectionName = _PFOCollectionName;}
		~eventHistFiller(){}

		int init();
		int fillEvent(const EVENT::LCEvent*);
		int writeToFile(TFile* outFile);

	private:
		unsigned int nSelectecTruthParticles;
		string PFOCollectionName;
		EVENT::LCCollection *PFOCollection;
		int checkPfoType(vector<unsigned int> inVec);
		// map<string,histStruct> singleParticleHistStructMap;
};
#endif // EVENTHISTFILLER_H
