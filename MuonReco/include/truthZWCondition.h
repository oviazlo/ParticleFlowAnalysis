/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/include/truthZWCondition.h
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	15th Dec 2017
 * 	Last Update:	15th Dec 2017
 */
#ifndef TRUTHZWCONDITION_H
#define TRUTHZWCONDITION_H

// ROOT
#include <TVector3.h>
#include <TMath.h>

//LCIO
#include <EVENT/LCEvent.h>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <Exceptions.h>
#include <globalConfig.h>
#include <vector>

class truthZWCondition
{

public:
	static truthZWCondition *instance(){
		if (!s_instance){
			s_instance = new truthZWCondition;
			s_instance->initDefault();
		}
		return s_instance;
	}
	static truthZWCondition& get(){
		static truthZWCondition instance;
		return instance;
	}

	// Set functions
	void setMCTruthCollectionName(std::string inStr){MCTruthCollectionName = inStr;}
	void setDebugFlag(bool _debugFlag){debugFlag=_debugFlag;}
	void setPartTypeToSelect(unsigned int _partType){partTypeToSelect = _partType;}
	void setEnergySelectionCut(double _energyCut){particleEnergyCut = _energyCut;}
	void setThetaSelectionCut(double _thetaCut){particleThetaCut = _thetaCut;}
	void setMotherPDG(double _motherPDG){motherPDG = _motherPDG;}

	unsigned int getPartTypeToSelect(){return partTypeToSelect;}
	double getEnergySelectionCut(){return particleEnergyCut;}
	double getThetaSelectionCut(){return particleThetaCut;}

	// Dump function
	void dumpTruthPart(const EVENT::MCParticle* part, const int counter = 0);
	void dumpTruthCondition();
	
	unsigned int getnTruthParticles(){return nTruthParticles;}
	unsigned int getnStableGenParticles(){return nStableGenParticles;}

	// Main functions
	void processEvent(const EVENT::LCEvent* _event);
	std::vector<EVENT::MCParticle*> getGenParticles(){return particlesOfInterest;}

protected:

private:
	truthZWCondition(){};
	truthZWCondition(const truthZWCondition&){};
	truthZWCondition& operator=(const truthZWCondition&){};
	static truthZWCondition* s_instance;

	void initDefault(){debugFlag = config::vm.count("debug"); partTypeToSelect = 13; particleEnergyCut = 0.0; particleThetaCut = 0.0; }

	unsigned int partTypeToSelect = 0;
	const EVENT::LCEvent* event;
	unsigned int nTruthParticles;
	unsigned int nStableGenParticles;
	double particleEnergyCut = 0.0;
	double particleThetaCut = 0.0;
	double motherPDG = 0;
	std::vector<EVENT::MCParticle*> particlesOfInterest;
	EVENT::LCCollection *MCTruthCollection;
	std::string MCTruthCollectionName;
	bool debugFlag;
};


#endif // TRUTHZWCONDITION_H
