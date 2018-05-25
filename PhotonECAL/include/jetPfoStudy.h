/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/include/jetPfoStudy.h
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
// #include "truthCondition.h"

class jetPfoStudy : public objectFill{
	public:
		jetPfoStudy(string _outDirName, string _mcCollectionName = "MCParticle", string _PFOCollectionName = "PandoraPFO") : objectFill(_outDirName) {MCCollectionName = _mcCollectionName; PFOCollectionName = _PFOCollectionName;}
		~jetPfoStudy(){}

		int init();
		int fillEvent(const EVENT::LCEvent*);
		int writeToFile(TFile* outFile);
		string getPFOCollectionName(){return PFOCollectionName;}
		void setDebugFlag(bool _debugFlag){debugFlag=_debugFlag;}

	private:
		unsigned int nSelectecTruthParticles;
		string PFOCollectionName;
		string MCCollectionName;
		EVENT::LCCollection *PFOCollection;
		int checkPfoType(vector<unsigned int> inVec);
		vector<const EVENT::MCParticle*> mcQuarkVector;
		double jetTheta, jetCosTheta;
		bool debugFlag;
};
#endif // EVENTHISTFILLER_H
