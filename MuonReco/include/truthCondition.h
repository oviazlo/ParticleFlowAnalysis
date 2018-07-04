/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/include/truthCondition.h
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	15th Dec 2017
 * 	Last Update:	15th Dec 2017
 */
#ifndef TRTUTHCONDITION_H
#define TRTUTHCONDITION_H

// ROOT
#include <TVector3.h>
#include <TMath.h>

//LCIO
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <Exceptions.h>
#include <globalConfig.h>

class truthCondition
{

public:
	static truthCondition *instance(){
		if (!s_instance){
			s_instance = new truthCondition;
			s_instance->initDefault();
		}
		return s_instance;
	}
	static truthCondition& get(){
		static truthCondition instance;
		return instance;
	}

	// Set functions
	void setEvent(const EVENT::LCEvent* _event){event = _event;}
	void setMCTruthCollectionName(std::string inStr){MCTruthCollectionName = inStr;}
	void setDebugFlag(bool _debugFlag){debugFlag=_debugFlag;}

	// Dump function
	void dumpTruthPart(const EVENT::MCParticle* part, const int counter = 0);
	void dumpTruthCondition();
	
	unsigned int getnTruthParticles(){return nTruthParticles;}
	int getpartGun_stablePartType(){return partGun_stablePartType;}
	int get_partGun_isStablePartDecayedInTracker(){return partGun_isStablePartDecayedInTracker;}
	bool get_simFSRPresent(){return simFSRPresent;}

	// Main functions
	void processEvent();
	EVENT::MCParticle* getGunParticle(){return static_cast<EVENT::MCParticle*>(MCTruthCollection->getElementAt(partGun_idStableGenPart));}

protected:

private:
	truthCondition(){};
	truthCondition(const truthCondition&){};
	truthCondition& operator=(const truthCondition&){};
	static truthCondition* s_instance;

	void initDefault(){debugFlag = config::vm.count("debug"); simFSRPresent=false;}

	const EVENT::LCEvent* event;
	unsigned int nTruthParticles;
	unsigned int nStableGenParticles;
	unsigned int partGun_idStableGenPart;
	bool partGun_isStablePartDecayedInTracker;
	bool simFSRPresent;
	int partGun_stablePartType;
	EVENT::LCCollection *MCTruthCollection;
	std::string MCTruthCollectionName;
	bool debugFlag;
};


#endif // TRTUTHCONDITION_H
