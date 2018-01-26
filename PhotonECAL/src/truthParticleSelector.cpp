#include <truthParticleSelector.h>

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

	eventHistFiller* eventFill = NULL;
	string mergeTag = "";

	eventFill = new eventHistFiller("eventHists",effCollection);
	objFillMap["eventHists"] = eventFill;

	eventFill = new eventHistFiller("eventHists_noConv",effCollection);
	objFillMap["eventHists_noConv"] = eventFill;

	eventFill = new eventHistFiller("eventHists_conv",effCollection);
	objFillMap["eventHists_conv"] = eventFill;

	eventFill = new eventHistFiller("eventHists_noFSR",effCollection);
	objFillMap["eventHists_noFSR"] = eventFill;

	eventFill = new eventHistFiller("eventHists_FSR",effCollection);
	objFillMap["eventHists_FSR"] = eventFill;

	mergeTag = "photonRecl";
	eventFill = new eventHistFiller("eventHists_"+mergeTag,effCollection);
	eventFill->setClusterMerging("photonMerge");
	objFillMap["eventHists_"+mergeTag] = eventFill;

	mergeTag = "photonReclMomDep";
	eventFill = new eventHistFiller("eventHists_"+mergeTag,effCollection);
	eventFill->setClusterMerging("photonMergeMomentumDep");
	objFillMap["eventHists_"+mergeTag] = eventFill;

	mergeTag = "photonAndNeutralRecl";
	eventFill = new eventHistFiller("eventHists_"+mergeTag,effCollection);
	eventFill->setClusterMerging("photonAndNeutralMerge");
	objFillMap["eventHists_"+mergeTag] = eventFill;

	mergeTag = "photonAndNeutralRecl_looseThetaCut";
	eventFill = new eventHistFiller("eventHists_"+mergeTag,effCollection);
	eventFill->setClusterMerging("photonAndNeutralLooseMerge");
	objFillMap["eventHists_"+mergeTag] = eventFill;

	// vector<int> pfoTypeVec = {11,22};
	// vector<int> pfoTypeVec = {11};
	// for (int ii=0; ii<pfoTypeVec.size(); ii++){
	//         int pfoTypeToUse = pfoTypeVec[ii];
	//         cout << "pfoTypeToUse: " << pfoTypeToUse << endl;
	//         string histDirPrefix = config::pfoTypeIntStringMap[pfoTypeToUse]+"_eventHist";
	//         cout << "histDirPrefix: " << histDirPrefix << endl;
	//         //
	//         eventFill = new eventHistFiller(histDirPrefix+"",pfoTypeToUse);
	//         eventFill->setPFOCollection(effCollection);
	//         objFillMap[histDirPrefix+""] = eventFill;
	//         //
	//         eventFill = new eventHistFiller(histDirPrefix+"_noConv",pfoTypeToUse);
	//         eventFill->setPFOCollection(effCollection);
	//         eventFill->setDiscardConvertions(true);
	//         objFillMap[histDirPrefix+"_noConv"] = eventFill;
	//         //
	//         eventFill = new eventHistFiller(histDirPrefix+"_photonRecl",pfoTypeToUse);
	//         eventFill->setPFOCollection(effCollection);
	//         eventFill->setMergePfoType(22);
	//         objFillMap[histDirPrefix+"_photonRecl"] = eventFill;
	//         //
	//         eventFill = new eventHistFiller(histDirPrefix+"_photonAndNeutralRecl",pfoTypeToUse);
	//         eventFill->setPFOCollection(effCollection);
	//         eventFill->setMergePfoType({22,2112});
	//         objFillMap[histDirPrefix+"_photonAndNeutralRecl"] = eventFill;
	//         //
	//         eventFill = new eventHistFiller(histDirPrefix+"_photonAndNeutralRecl_looseThetaCut",pfoTypeToUse);
	//         eventFill->setPFOCollection(effCollection);
	//         eventFill->setMergePfoType({22,2112});
	//         eventFill->setThetaMergingCut(2.0);
	//         objFillMap[histDirPrefix+"_photonAndNeutralRecl_looseThetaCut"] = eventFill;
	//         //
	//         eventFill = new eventHistFiller(histDirPrefix+"_conv",pfoTypeToUse);
	//         eventFill->setPFOCollection(effCollection);
	//         eventFill->setSelectConvertions(true);
	//         objFillMap[histDirPrefix+"_conv"] = eventFill;
	// }
        // //
	// eventFill = new eventHistFiller("Photons_Neutral_eventHists_photonAndNeutralRecl",{22,2112});
	// eventFill->setPFOCollection(effCollection);
	// eventFill->setMergePfoType({22,2112});
	// objFillMap["Photons_Neutral_eventHists_photonAndNeutralRecl"] = eventFill;

	// photonEffCalculator* effCalculator = new photonEffCalculator("photonEfficiency");
	// effCalculator->setPFOCollection(effCollection);
	// effCalculator->setPFOType(efficiencyPFOType);
	// effCalculator->setEfficiencyOneClusterRequirement(onlyOneRecoClusterPerEvent);
	// if (dPhiMergeValue > 0.0)
	//         effCalculator->setDPhiMergeValue(dPhiMergeValue);
	// objFillMap["photonEfficiency"] = effCalculator;

	for(auto const &mapElement : objFillMap){
		cout << "Init mapElement: " << mapElement.first << endl;
		mapElement.second->init();
	}

}

bool truthParticleSelector::selectEvent(const EVENT::LCEvent* event){
	
	EVENT::MCParticle* part = truthCondition::instance()->getGunParticle();
	if (discardFSREvents && truthCondition::instance()->getnTruthParticles()!=1) 
		return 0;
	const double *partMom = part->getMomentum();
	TVector3 v1;
	v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
	double partTheta = 180.*v1.Theta()/TMath::Pi();
	if (partTheta<8 || partTheta>172) 
		return false;
	double partMomMag = v1.Mag();
	double partPhi = 180.*v1.Phi()/TMath::Pi();
	if ((partMomMag<energyRange.first) || (partMomMag>energyRange.second))
		return false;
	if ((partTheta<thetaRange.first) || (partTheta>thetaRange.second))
		return false;
	if ((partPhi<phiRange.first) || (partPhi>phiRange.second))
		return false;

	if (config::vm.count("debug"))
		cout << "[INFO]\t *****EVENT: " << event->getEventNumber() << " *****" <<endl;
	
	for(auto const &mapElement : objFillMap){
		if (IsInWord(mapElement.first,"_noFSR")){
			if (truthCondition::instance()->get_simFSRPresent()==false)
				mapElement.second->fillEvent(event);
		}
		else if (IsInWord(mapElement.first,"_FSR")){
			if (truthCondition::instance()->get_simFSRPresent()==true)
				mapElement.second->fillEvent(event);
		}
		else if (IsInWord(mapElement.first,"_conv")){
			if (truthCondition::instance()->get_partGun_isStablePartDecayedInTracker()==true)
				mapElement.second->fillEvent(event);
		}
		else if (IsInWord(mapElement.first,"_noConv")){
			if (truthCondition::instance()->get_partGun_isStablePartDecayedInTracker()==false)
				mapElement.second->fillEvent(event);
		}
		else{
			mapElement.second->fillEvent(event);
		}
	}
	return true;
}

truthParticleSelector::truthParticleSelector(){
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
