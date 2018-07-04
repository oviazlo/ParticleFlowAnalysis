/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/src/eventHistFiller.cpp
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	05th Dec 2017
 * 	Last Update:	05th Dec 2017
 */

#include "eventHistFiller.h"

/*===========================================================================*/
/*===============================[ function implementations ]================*/
/*===========================================================================*/

// double round(double value, unsigned int prec = 1)
// {
//             return floor( value*pow(10,prec) + 0.5 )/pow(10,prec);
// }


int eventHistFiller::fillEvent(const EVENT::LCEvent* event){
	try {
		PFOCollection = event->getCollection(PFOCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|eventHistFiller]\tCan't find collection: " << PFOCollectionName << endl;
		return -1;
	
	}
	nPFOs = PFOCollection->getNumberOfElements();
	// cout << "nPFOs: " << nPFOs << endl;

	// if (config::vm.count("debug"))
	//         cout << "[INFO]	eventHistFiller::fillEvent: " << event->getEventNumber() << endl;

	initTruthInfoAndFillIt();
	// return 0;

	vector<EVENT::ReconstructedParticle*> recoPFOs = getObjVecFromCollection<EVENT::ReconstructedParticle*>(PFOCollection);
	getHistFromMap("nPFOs")->Fill(recoPFOs.size());
	if (recoPFOs.size()==0) 
		return 0; // no reco PFOs

	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++)
		pfoCounter[it->second] = 0;

	totalRecoEnergy = 0.0;
	bool efficiencyHistWasAlreadyFilled = false;

	bool passTypeCut = false;
	bool passAngularCut = false;
	bool passEnergyCut = false;

	double tmpEnergy = std::numeric_limits<double>::min();

	if (config::vm.count("debug"))
		cout << "Reconstructed particles:" << endl;
	for (int i=0; i<recoPFOs.size(); i++){

		fillPfoCounterAndFillGeneralPfoInfo(i);

		if (recoPFOs[i]->getType()==truthType){
			passTypeCut = true;
			if (recoPFOs[i]->getEnergy()>tmpEnergy){
				tmpEnergy = recoPFOs[i]->getEnergy();
				idOfPartCandidate = i;
			}
		}
	}


	if (passTypeCut==false){
		fillEventsFailedSelection();
		if (config::vm.count("debug"))
			cout << "***Do Not Pass Type Cut***" << endl;
	}
	else{
		partCandidate = static_cast<IMPL::ReconstructedParticleImpl*>(PFOCollection->getElementAt(idOfPartCandidate));
		getHistFromMap("efficiencyVsEnergy_onlyType")->Fill(truthEnergy);

		//##### Cluster merging
		vector<unsigned int> mergedParticleIds = mergeClusters();


		bool isMergedCandidate = (mergedParticleIds.size()>0);

		if (config::vm.count("debug"))
			std::cout << endl << "Main Truth and Reconstructed Candidates:" << endl << "t:  pdg: " << std::setw(5) << truthType << ": E: " << std::setw(6) << (round(100*truthEnergy)/100.0) << ": pT: " << std::setw(6) << (round(100*truthPt)/100.0) << "; theta: " << std::setw(6) << round(100*truthTheta)/100.0 << "; phi: " << std::setw(6) << round(100*truthPhi)/100.0 << std::endl;

		//##### Angular and energy matching
		angularAndEnergyMatching();

		//##### Histogram filling
		fillOtherHists();


		if (config::vm.count("debug")){
			cout << endl << "Merged particles: "; 
			for(auto const& value: mergedParticleIds) 
			     std::cout << value << " ";
			cout << endl << endl;
		}

	}

}



