/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/src/truthZWCondition.cpp
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	15th Dec 2017
 * 	Last Update:	15th Dec 2017
 */

#include "truthZWCondition.h"

truthZWCondition* truthZWCondition::s_instance = NULL;  

/*===========================================================================*/
/*===============================[ implementation ]===============================*/
/*===========================================================================*/

void truthZWCondition::processEvent(const EVENT::LCEvent* _event){
	event = _event;

	nTruthParticles = 0;
	nStableGenParticles = 0;
	particlesOfInterest.clear();

	try {
		MCTruthCollection = event->getCollection(MCTruthCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		std::cout << "[ERROR|truthZWCondition::setEvent]\tCan't find collection: " << MCTruthCollectionName << std::endl;
	}

	nTruthParticles = MCTruthCollection->getNumberOfElements();
	if (config::vm.count("debug"))
		std::cout << "Truth particles:" << std::endl;
	std::cout << "Truth particles:" << std::endl;
	for(int j=0; j < nTruthParticles; j++) {
		auto part = static_cast<EVENT::MCParticle*>(MCTruthCollection->getElementAt(j));
		TVector3 v1(part->getMomentum());
		double partTheta = v1.Theta()*TMath::RadToDeg();
		auto parents = part->getParents();
		int parentPDG = 0;
		if (parents.size()==1)
			parentPDG = parents[0]->getPDG();

		if (config::vm.count("debug") )
			dumpTruthPart(part,j);
		if (
				part->getGeneratorStatus()>=1
				&& abs(part->getPDG())==partTypeToSelect 
				&& (part->getEnergy()>particleEnergyCut)
				&& (partTheta>particleThetaCut) && (partTheta<(180.0-particleThetaCut))
				// && (abs(parentPDG)==abs(motherPDG))
		   )
		{
			dumpTruthPart(part,j);
			auto daughters = part->getDaughters();
			for (auto iPartInnerLoop: daughters){
				dumpTruthPart(iPartInnerLoop,j);
			}
			// EVENT::MCParticle* loopPart = nullptr;
			// bool noStableLeptonsOfInterest = false;
			// while (noStableLeptonsOfInterest == false){
			//         for (auto iPartInnerLoop: daughters){
			//                 dumpTruthPart(iPartInnerLoop,j);
			//                 // if abs(iPartInnerLoop->getPDG())==partTypeToSelect &&
			//         }
			//         noStableLeptonsOfInterest = true;
			// }
			nStableGenParticles++;
			particlesOfInterest.push_back(part);
		}
	}
	// if (config::vm.count("debug"))
	//         dumpTruthCondition();
	// if (simFSRPresent==false && nTruthParticles>1)
	//         dumpTruthCondition();

}

void truthZWCondition::dumpTruthCondition(){

	std::cout << "Event\t" << event->getEventNumber() << "; nTruthParticles: " << nTruthParticles << "; nStableGenParticles: " << nStableGenParticles << std::endl << std::endl;	

}

// bool hasParentOfType(const EVENT::MCParticle* part, unsigned int parentPDG){
//
// }

void truthZWCondition::dumpTruthPart(const EVENT::MCParticle* part, const int counter){
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
	unsigned int nParents = part->getParents().size();
	unsigned int nDaughters = part->getDaughters().size();
	std::cout << "t" << counter << ": pdg: " << std::setw(5) << pdgId << ": E: " << std::setw(6) << (round(100*part->getEnergy())/100.0) << ": pT: " << std::setw(6) << (round(100*v1.Pt())/100.0) << "; theta: " << std::setw(6) << round(100*partTheta)/100.0 << "; phi: " << std::setw(6) << round(100*partPhi)/100.0 << "; inTracker: " << inTracker << "; inCal: " << inCal << "; genStatus: " << genStatus << "; isRadiation: " << vertexIsNotEndpointOfParent << "; vertexR: " << vertexR << "; vertexZ: " << vertexZ << "; nParents: " << nParents << "; nDaughters: " << nDaughters << std::endl;
}
/*===========================================================================*/
/*===============================[ implementation ]===============================*/
/*===========================================================================*/



