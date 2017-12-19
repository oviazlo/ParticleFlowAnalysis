#include <truthParticleSelector.h>

bool truthParticleSelector::selectEvent(const EVENT::LCEvent* event){
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
	if (discardFSREvents && nElements!=1) 
		return 0;
	for(int j=0; j < nElements; j++) {
		auto part = dynamic_cast<EVENT::MCParticle*>(MCTruthCollection->getElementAt(j));
		if (part->getGeneratorStatus()==1){
			const double *partMom = part->getMomentum();
			TVector3 v1;
			v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
			double partTheta = 180.*v1.Theta()/TMath::Pi();
			if (partTheta<8 || partTheta>172) 
				return false;
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

truthParticleSelector::truthParticleSelector(const string _mcTruthCollection){
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
	discardFSREvents = false;
	dPhiMergeValue = 0;
	onlyOneRecoClusterPerEvent = false;
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

	int pfoTypeToUse = 22;
	eventHistFiller* eventFill = NULL;
        //
	eventFill = new eventHistFiller("eventHists",pfoTypeToUse);
	eventFill->setPFOCollection(effCollection);
	eventFill->setMCTruthCollection(mcTruthCollection);
	objFillMap["eventHists"] = eventFill;
        //
	eventFill = new eventHistFiller("eventHists_noConv",pfoTypeToUse);
	eventFill->setPFOCollection(effCollection);
	eventFill->setMCTruthCollection(mcTruthCollection);
	eventFill->setDiscardConvertions(true);
	objFillMap["eventHists_noConv"] = eventFill;
        //
	eventFill = new eventHistFiller("eventHists_photonRecl",pfoTypeToUse);
	eventFill->setPFOCollection(effCollection);
	eventFill->setMCTruthCollection(mcTruthCollection);
	eventFill->setMergePfoType(22);
	objFillMap["eventHists_photonRecl"] = eventFill;
        //
	eventFill = new eventHistFiller("eventHists_photonAndNeutralRecl",pfoTypeToUse);
	eventFill->setPFOCollection(effCollection);
	eventFill->setMCTruthCollection(mcTruthCollection);
	eventFill->setMergePfoType({22,2112});
	objFillMap["eventHists_photonAndNeutralRecl"] = eventFill;
	//
	eventFill = new eventHistFiller("eventHists_photonAndNeutralRecl_looseThetaCut",pfoTypeToUse);
	eventFill->setPFOCollection(effCollection);
	eventFill->setMCTruthCollection(mcTruthCollection);
	eventFill->setMergePfoType({22,2112});
	eventFill->setThetaMergingCut(2.0);
	objFillMap["eventHists_photonAndNeutralRecl_looseThetaCut"] = eventFill;
        //
	eventFill = new eventHistFiller("eventHists_conv",pfoTypeToUse);
	eventFill->setPFOCollection(effCollection);
	eventFill->setMCTruthCollection(mcTruthCollection);
	eventFill->setSelectConvertions(true);
	objFillMap["eventHists_conv"] = eventFill;
        //
	eventFill = new eventHistFiller("eventHistsPhotonsAndNeutrals_photonAndNeutralRecl",{22,2112});
	eventFill->setPFOCollection(effCollection);
	eventFill->setMCTruthCollection(mcTruthCollection);
	eventFill->setMergePfoType({22,2112});
	objFillMap["eventHistsNeutrals_photonAndNeutralRecl"] = eventFill;

	// photonEffCalculator* effCalculator = new photonEffCalculator("photonEfficiency");
	// effCalculator->setPFOCollection(effCollection);
	// effCalculator->setMCTruthCollection(mcTruthCollection);
	// effCalculator->setPFOType(efficiencyPFOType);
	// effCalculator->setEfficiencyOneClusterRequirement(onlyOneRecoClusterPerEvent);
	// if (dPhiMergeValue > 0.0)
	//         effCalculator->setDPhiMergeValue(dPhiMergeValue);
	// objFillMap["photonEfficiency"] = effCalculator;

	for(auto const &mapElement : objFillMap){
		mapElement.second->init();
	}

}

void truthParticleSelector::writeToFile(TFile *outFile){
	for(auto const &mapElement : objFillMap){
		mapElement.second->writeToFile(outFile);
	}
}

string truthParticleSelector::getPostFixString(){
	string postFix = "E"+DoubToStr((energyRange.first+energyRange.second)/2.0);
	if (thetaRange != make_pair(-180.0,180.))
		postFix += "_Theta" + DoubToStr(thetaRange.first)+"_"+DoubToStr(thetaRange.second);
	if (phiRange != make_pair(-360.0,360.))
		postFix += "_Phi"+ DoubToStr(phiRange.first)+"_"+DoubToStr(phiRange.second);	
	return postFix;
}
