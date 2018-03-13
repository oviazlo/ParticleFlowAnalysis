#include <photonEffCalculator.h>
#include <IMPL/ReconstructedParticleImpl.h>

int photonEffCalculator::init(){

	TH1D* tmpHist;
	tmpHist = new TH1D("dAngle","(PFO,MCTruth) angle; angle [degree]; Counts",1800,0.0,90.0);
	histMap["dAngle"] = tmpHist;
	tmpHist = new TH1D("dE","E_{PFO}-E_{MC}; Energy [GeV]; Counts",100,-25.0,25.0);
	histMap["dE"] = tmpHist;
	tmpHist = new TH1D("dE_matched","E_{PFO}-E_{MC}^{matched}; Energy [GeV]; Counts",1000,-25.0,25.0);
	histMap["dE_matched"] = tmpHist;
	tmpHist = new TH1D("PFO_E","E_{PFO}; Energy [GeV]; Counts",750,0.0,150.0);
	histMap["PFO_E"] = tmpHist;
	tmpHist = new TH1D("findableMC_vs_theta","Number of findable MC particle vs theta; Theta [degree]; Counts",180,0.0,180.0);
	histMap["findableMC_vs_theta"] = tmpHist;
	tmpHist = new TH1D("matchedMC_vs_theta","Number of matched MC particle vs theta; Theta [degree]; Counts",180,0.0,180.0);
	histMap["matchedMC_vs_theta"] = tmpHist;
	tmpHist = new TH1D("findableMC_vs_cosTheta","Number of findable MC particle vs cos(theta); Cos(Theta); Counts",2*180,-1.0,1.0);
	histMap["findableMC_vs_cosTheta"] = tmpHist;
	tmpHist = new TH1D("matchedMC_vs_cosTheta","Number of matched MC particle vs cos(theta); Cos(Theta); Counts",2*180,-1.0,1.0);
	histMap["matchedMC_vs_cosTheta"] = tmpHist;
	tmpHist = new TH1D("findableMC_vs_E","Number of findable MC particle vs energy; Energy [GeV]; Counts",19,7.5,102.5);
	histMap["findableMC_vs_E"] = tmpHist;
	tmpHist = new TH1D("matchedMC_vs_E","Number of matched MC particle vs energy; Energy [GeV]; Counts",19,7.5,102.5);
	histMap["matchedMC_vs_E"] = tmpHist;
	debugFlag = false;
	dPhiMergeValue = 0;
	onlyOneRecoClusterPerEvent = false;
	return 0;

}

EVENT::ReconstructedParticle* photonEffCalculator::getMatchedPFO(const EVENT::MCParticle* inMCPart, const vector<EVENT::ReconstructedParticle*> findablePFO){

	EVENT::ReconstructedParticle* matchedPFO = NULL;

	// get MC part vector
	TVector3 v0;
	const double *partMom = inMCPart->getMomentum();
	v0.SetXYZ(partMom[0],partMom[1],partMom[2]);
	double inPartTheta = 180.*v0.Theta()/TMath::Pi();
	double inPartEnergy = inMCPart->getEnergy();
	histMap["findableMC_vs_theta"]->Fill(inPartTheta);
	histMap["findableMC_vs_cosTheta"]->Fill(cos(inPartTheta/180.0*TMath::Pi()));
	histMap["findableMC_vs_E"]->Fill(inPartEnergy);
	
	// check matching with each PFO
	for (unsigned int i=0; i<findablePFO.size(); i++){
		auto part = findablePFO[i];
		const double* partMom = part->getMomentum();
		TVector3 v1;
		v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
		double partTheta = 180.*v1.Theta()/TMath::Pi();
		// Direction requirement:
		// within a cone of 1 degree
		histMap["dAngle"]->Fill(180.0*v0.Angle(v1)/TMath::Pi());
		// if (180.0*v0.Angle(v1)/TMath::Pi()>1.0) continue;
		// Energy requirement:
		// energy difference less than resolution
		double eRes = 2*TMath::Sqrt(inPartEnergy)+0.5;
		double dE = part->getEnergy()- inPartEnergy;
		histMap["PFO_E"]->Fill(part->getEnergy());
		histMap["dE"]->Fill(dE);
		// if ( abs(dE) > eRes ) continue;
		histMap["dE_matched"]->Fill(dE);
		histMap["matchedMC_vs_theta"]->Fill(inPartTheta);
		histMap["matchedMC_vs_cosTheta"]->Fill(cos(inPartTheta/180.0*TMath::Pi()));
		histMap["matchedMC_vs_E"]->Fill(inPartEnergy);
		matchedPFO = part;
		break;
	}

	return matchedPFO;
}

int photonEffCalculator::fillEvent(const EVENT::LCEvent* event){
	// read collections
	try {
		PFOCollection = event->getCollection(PFOCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|particleFill]\tCan't find collection: " << PFOCollectionName << endl;
		return -1;
	}
	try {
		MCTruthCollection = event->getCollection(MCTruthCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|particleFill]\tCan't find collection: " << MCTruthCollectionName << endl;
		return -1;
	}
	
	// find primary generated MC particle which is the only findable in the event by definition [particle gun]
	EVENT::MCParticle* genPart = NULL;
	int nElements = MCTruthCollection->getNumberOfElements();
	for(int j=0; j < nElements; j++) {
		auto part = dynamic_cast<EVENT::MCParticle*>(MCTruthCollection->getElementAt(j));
		if (part->getGeneratorStatus()==1){
			const double *partMom = part->getMomentum();
			TVector3 v1;
			v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
			double partTheta = 180.*v1.Theta()/TMath::Pi();
			if (partTheta<8 || partTheta>172) 
				return 0;
			if (part->isDecayedInTracker()) 
				return 0;
			genPart = part;
			break;
		}
	}

	if (PFOPartType==0){
		PFOPartType = genPart->getPDG();
		// WARNING FIXME according to Mathias Padora can return only these types:
		vector<int> PandoraAllowedTypes = {11,13,22,211,2112};
		if (find(PandoraAllowedTypes.begin(), PandoraAllowedTypes.end(), abs(PFOPartType))==PandoraAllowedTypes.end() )
			PFOPartType = 2112; // neutron
		cout << "[INFO]\t[photonEffCalculator] Efficiency PFO type is assigned to: " << PFOPartType << endl;
	}

	vector<EVENT::ReconstructedParticle*> recoPFOs = getObjVecFromCollection<EVENT::ReconstructedParticle*>(PFOCollection);
	// cout << "[DEBUG]\trecoPFOs.size():" << recoPFOs.size() << endl;
	if (recoPFOs.size()==0) 
		return 0; // no reco PFOs
	// FIXME hardcoded type!!!
	// if (PFOPartType==22 && onlyOneRecoClusterPerEvent && recoPFOs.size()!=1) return 0;
	int itMostEnergeticOfType = -1;
	for (int i=0; i<recoPFOs.size(); i++){
		// cout << "[DEBUG]\tPFOPartType: " << PFOPartType << "\trecoPFOs[i]->getType(): " << recoPFOs[i]->getType() << endl;
		if (recoPFOs[i]->getType()!=PFOPartType) 
			continue;
		if (itMostEnergeticOfType==-1){
			itMostEnergeticOfType = i;
			break;
		}
		if (recoPFOs[i]->getEnergy()>recoPFOs[itMostEnergeticOfType]->getEnergy())
			itMostEnergeticOfType = i;
	}

	// if (itMostEnergeticOfType == -1) return 0; // no PFO of needed type
	// cout << "[DEBUG]\titMostEnergeticOfType: " << itMostEnergeticOfType << endl;

	// FIXME WARNING **********************************
	// create modifiable PFO
	if (itMostEnergeticOfType != -1 && dPhiMergeValue>0.0){
		IMPL::ReconstructedParticleImpl* partMod = new IMPL::ReconstructedParticleImpl();
		partMod->setEnergy(recoPFOs[itMostEnergeticOfType]->getEnergy());
		partMod->setMomentum(recoPFOs[itMostEnergeticOfType]->getMomentum());
		partMod->setType(recoPFOs[itMostEnergeticOfType]->getType());

		// delete non-modifiable one
		recoPFOs.erase(recoPFOs.begin() + itMostEnergeticOfType);

		// merge allPFOs within criteria
		TVector3 vecToMerge;
		const double *partMomMod = partMod->getMomentum();
		vecToMerge.SetXYZ(partMomMod[0],partMomMod[1],partMomMod[2]);
		for (int i=0; i<recoPFOs.size(); i++){
			const double *partMom2 = recoPFOs[i]->getMomentum();
			TVector3 v2;
			v2.SetXYZ(partMom2[0],partMom2[1],partMom2[2]);
			double dPhi = vecToMerge.DeltaPhi(v2)*180./TMath::Pi();
			if (abs(dPhi)<dPhiMergeValue)
				partMod->setEnergy(partMod->getEnergy()+recoPFOs[i]->getEnergy());
		}
		recoPFOs.push_back(partMod);
	}
	// FIXME WARNING **********************************

	auto maxEnergeticRecoPFO = recoPFOs.back();
	for (int j=0; j<(recoPFOs.size()-1); j++){
		auto part = recoPFOs[j];
		if (part->getEnergy()>maxEnergeticRecoPFO->getEnergy())
			maxEnergeticRecoPFO = part;
	}

	vector<EVENT::ReconstructedParticle*> findablePFO;
	if (maxEnergeticRecoPFO->getType() == PFOPartType)  
		findablePFO.push_back(maxEnergeticRecoPFO);

	EVENT::ReconstructedParticle* matchedPFO = getMatchedPFO(genPart,findablePFO);
	if (matchedPFO==NULL)
		return 0;
	// TODO continue logic...
	// fill some histograms for number of findable and number of matched
	return 0;

}

int photonEffCalculator::writeToFile(TFile* outFile){
	vector<string> postFix = {"E", "theta", "cosTheta"};
	for (unsigned int i=0; i<postFix.size(); i++) {
		if (histMap.find("findableMC_vs_" + postFix[i])!=histMap.end() && histMap.find("matchedMC_vs_" + postFix[i])!=histMap.end()){
			TH1D* tmpHist = dynamic_cast<TH1D*>(histMap["matchedMC_vs_"+postFix[i]]->Clone(("efficiency_vs_"+postFix[i]).c_str()));
			tmpHist->Divide(histMap["findableMC_vs_" + postFix[i]]);
			tmpHist->GetYaxis()->SetTitle("Efficiency");
			histMap["efficiency_vs_"+postFix[i]] = tmpHist;
		}
	}
	return objectFill::writeToFile(outFile);
}


