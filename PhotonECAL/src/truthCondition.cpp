/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/src/truthCondition.cpp
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	15th Dec 2017
 * 	Last Update:	15th Dec 2017
 */

#include "truthCondition.h"

truthCondition* truthCondition::s_instance = NULL;  

/*===========================================================================*/
/*===============================[ implementation ]===============================*/
/*===========================================================================*/

void truthCondition::processEvent(){

	nTruthParticles = 0;
	nStableGenParticles = 0;
	simFSRPresent = false;

	try {
		MCTruthCollection = event->getCollection(MCTruthCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		std::cout << "[ERROR|truthCondition::setEvent]\tCan't find collection: " << MCTruthCollectionName << std::endl;
	}

	nTruthParticles = MCTruthCollection->getNumberOfElements();
	if (config::vm.count("debug"))
		std::cout << "Truth particles:" << std::endl;
	for(int j=0; j < nTruthParticles; j++) {
		auto part = static_cast<EVENT::MCParticle*>(MCTruthCollection->getElementAt(j));
		if (config::vm.count("debug"))
			dumpTruthPart(part,j);
		if (part->getGeneratorStatus()==1){
			nStableGenParticles++;
			partGun_idStableGenPart = j;
			partGun_stablePartType = part->getPDG();
			if (part->isDecayedInTracker())
				partGun_isStablePartDecayedInTracker = true;
			else
				partGun_isStablePartDecayedInTracker = false;
		}
		else{
			if (part->vertexIsNotEndpointOfParent()==true)
				simFSRPresent = true;
		}
	}
	// if (config::vm.count("debug"))
	//         dumpTruthCondition();
	// if (simFSRPresent==false && nTruthParticles>1)
	//         dumpTruthCondition();

}

void truthCondition::dumpTruthCondition(){

	std::cout << "Event\t" << event->getEventNumber() << "; nTruthParticles: " << nTruthParticles << "; nStableGenParticles: " << nStableGenParticles << "; partGun_stablePartType: " << partGun_stablePartType << "; partGun_isStablePartDecayedInTracker: " << partGun_isStablePartDecayedInTracker << "; partGun_idStableGenPart: " << partGun_idStableGenPart << "; simFSRPresent: " << simFSRPresent << std::endl << std::endl;	

}

void truthCondition::dumpTruthPart(const EVENT::MCParticle* part, const int counter){
	const double *partMom = part->getMomentum();
	TVector3 v1;
	v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
	double partPhi = 180.*v1.Phi()/TMath::Pi();
	double partTheta = 180.*v1.Theta()/TMath::Pi();
	bool inTracker = part->isDecayedInTracker();
	bool inCal = part->isDecayedInCalorimeter();
	int genStatus = part->getGeneratorStatus();
	int pdgId = part->getPDG();
	bool vertexIsNotEndpointOfParent = part->vertexIsNotEndpointOfParent();
	const double *vertexPos = part->getVertex();
	double vertexR = sqrt(vertexPos[0]*vertexPos[0]+vertexPos[1]*vertexPos[1]);
	double vertexZ = vertexPos[2];
	std::cout << "t" << counter << ": pdg: " << std::setw(5) << pdgId << ": E: " << std::setw(6) << (round(100*part->getEnergy())/100.0) << ": pT: " << std::setw(6) << (round(100*v1.Pt())/100.0) << "; theta: " << std::setw(6) << round(100*partTheta)/100.0 << "; phi: " << std::setw(6) << round(100*partPhi)/100.0 << "; inTracker: " << inTracker << "; inCal: " << inCal << "; genStatus: " << genStatus << "; isRadiation: " << vertexIsNotEndpointOfParent << "; vertexR: " << vertexR << "; vertexZ: " << vertexZ << std::endl;
}
/*===========================================================================*/
/*===============================[ implementation ]===============================*/
/*===========================================================================*/



