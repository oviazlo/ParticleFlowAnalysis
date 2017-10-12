#include <truthParticleSelector.h>

bool truthParticleSelector::selectEvent(EVENT::LCEvent* event){
	EVENT::LCCollection *MCTruthCollection = NULL;
	try {
		MCTruthCollection = event->getCollection(mcTruthCollection);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|particleFill]\tCan't find collection: " << mcTruthCollection << endl;
		return -1;
	}
	
	// find primary generated MC particle which is the only findable in the event by definition [particle gun]
	EVENT::MCParticle* genPart = NULL;
	int nElements = MCTruthCollection->getNumberOfElements();
	if (discardFSREvents && nElements!=1) return 0;
	for(int j=0; j < nElements; j++) {
		auto part = static_cast<EVENT::MCParticle*>(MCTruthCollection->getElementAt(j));
		if (part->getGeneratorStatus()==1){
			const double *partMom = part->getMomentum();
			TVector3 v1;
			v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
			double partTheta = 180.*v1.Theta()/TMath::Pi();
			if (partTheta<8 || partTheta>172) return false;
			genPart = part;
			double partMomMag = v1.Mag();
			double partPhi = 180.*v1.Phi()/TMath::Pi();
			if ((partMomMag<energyRange.first) || (partMomMag>energyRange.second))
				return false;
			if ((partTheta<thetaRange.first) || (partTheta>thetaRange.second))
				return false;
			if ((partPhi<phiRange.first) || (partPhi>phiRange.second))
				return false;
			
			for(auto const &mapElement : objFillMap){
				mapElement.second->fillEvent(event);
			}
			return true;
		}
	}
}

truthParticleSelector::truthParticleSelector(string _mcTruthCollection){
	mcTruthCollection=_mcTruthCollection;
	energyRange = make_pair(0.0,std::numeric_limits<double>::max());
	thetaRange = make_pair(-180.0,180.);
	phiRange = make_pair(-360.0,360.);

	efficiencyPFOType = 0;	

}

truthParticleSelector::~truthParticleSelector(){
	for(auto const &mapElement : objFillMap){
		delete mapElement.second;
	} 
}

void truthParticleSelector::init(){
	for (int i; i<particleFillCollections.size(); i++){
		particleFill* tmpPartFill = new particleFill(particleFillCollections[i]);
		tmpPartFill->setCollectionName(particleFillCollections[i]);
		tmpPartFill->setReconstructedParticleType(PFOTypes[i]);
		string postfixRootDirName = "";
		for (auto j=0; j<PFOTypes[i].size(); j++)
			postfixRootDirName += "_" + DoubToStr(PFOTypes[i][j]);
		tmpPartFill->updateRootDirName(particleFillCollections[i]+postfixRootDirName);
		if (dPhiMergeValue > 0.0)
			tmpPartFill->setDPhiMergeValue(dPhiMergeValue);
		objFillMap[particleFillCollections[i]+postfixRootDirName] = tmpPartFill;
	}
	for (int i; i<energyFillCollections.size(); i++){
		energyFill* tmpEnergyFill = new energyFill(energyFillCollections[i]);
		tmpEnergyFill->setCollectionName(energyFillCollections[i]); 
		objFillMap[energyFillCollections[i]] = tmpEnergyFill;
	}

	photonEffCalculator* effCalculator = new photonEffCalculator("photonEfficiency");
	effCalculator->setPFOCollection(effCollection);
	effCalculator->setMCTruthCollection(mcTruthCollection);
	effCalculator->setPFOType(efficiencyPFOType);
	if (dPhiMergeValue > 0.0)
		effCalculator->setDPhiMergeValue(dPhiMergeValue);
	objFillMap["photonEfficiency"] = effCalculator;

	for(auto const &mapElement : objFillMap){
		mapElement.second->init();
	}

	discardFSREvents = false;
	dPhiMergeValue = 0;
}

void truthParticleSelector::writeToFile(TFile *outFile){
	for(auto const &mapElement : objFillMap){
		mapElement.second->writeToFile(outFile);
	}
}

