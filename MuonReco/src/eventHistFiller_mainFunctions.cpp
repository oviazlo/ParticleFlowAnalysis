#include "eventHistFiller.h"

void eventHistFiller::initTruthInfoAndFillIt(){
	genPart = truthCondition::instance()->getGunParticle();

	truthType = genPart->getPDG();
	const double *truthMom = genPart->getMomentum(); 
	TVector3 vTruthMom(truthMom[0],truthMom[1],truthMom[2]);
	truthTheta = vTruthMom.Theta()*TMath::RadToDeg();
	cosTruthTheta = TMath::Cos(truthTheta*TMath::DegToRad());
	truthPhi = vTruthMom.Phi()*TMath::RadToDeg();
	truthPt = vTruthMom.Pt();
	truthEnergy = genPart->getEnergy();

	getHistFromMap("truthParticle_isDecayedInTracker")->Fill(genPart->isDecayedInTracker());
	getHistFromMap("nTruthPartsVsCosTheta")->Fill(cosTruthTheta);
	getHistFromMap("nTruthPartsVsTheta")->Fill(truthTheta);
	getHistFromMap("nTruthPartsVsEnergy")->Fill(truthEnergy);
}

void eventHistFiller::fillPfoCounterAndFillGeneralPfoInfo(unsigned int partId){
	EVENT::ReconstructedParticle* recoPart = static_cast<IMPL::ReconstructedParticleImpl*>(PFOCollection->getElementAt(partId));

	const int pfoType = abs(recoPart->getType());
	const double *partMom = recoPart->getMomentum();
	TVector3 vPartMom(partMom[0],partMom[1],partMom[2]);
	const double partTheta = vPartMom.Theta()*TMath::RadToDeg();
	const double partPhi = vPartMom.Phi()*TMath::RadToDeg();
	const double cosPartTheta = TMath::Cos(partTheta*TMath::DegToRad());
	const double partEnergy = recoPart->getEnergy();
	const double dPhiPartTruth = (partPhi-truthPhi)*TMath::DegToRad();
	const double dThetaPartTruth = (partTheta-truthTheta)*TMath::DegToRad();

	getHistFromMap("PFOType")->Fill(pfoType);
	getHistFromMap("nPFOsVsCosTheta_all")->Fill(cosTruthTheta);
	getHistFromMap("nPFOsVsTheta_all")->Fill(truthTheta);
	getHistFromMap("totalEnergyVsTheta")->Fill(truthTheta,recoPart->getEnergy());
	totalRecoEnergy += recoPart->getEnergy();

	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
		if (abs(pfoType)==it->first) {
			getHistFromMap("nPFOsVsCosTheta_"+it->second)->Fill(cosTruthTheta);
			getHistFromMap("nPFOsVsTheta_"+it->second)->Fill(truthTheta);
			getHistFromMap("thetaResolution_"+it->second)->Fill(dThetaPartTruth);
			getHistFromMap("phiResolution_"+it->second)->Fill(dPhiPartTruth);
			getHistFromMap("energyResolution_"+it->second)->Fill(partEnergy);
			pfoCounter[it->second]++;
		}
	}

	if (config::vm.count("debug"))
		cout << "r" << partId << ": pdg: " << setw(5) << pfoType << "; E: " << setw(6) << round(100*partEnergy)/100.0 << ": pT: " << std::setw(6) << (round(100*vPartMom.Pt())/100.0) << "; theta: " << setw(6) << round(100*partTheta)/100.0 << "; phi: " << setw(6) << round(100*partPhi)/100.0 << "; dPhi: " <<  setw(6) << round(100*dPhiPartTruth*TMath::RadToDeg())/100.0 << "; dTheta: " << setw(6) << round(100*dThetaPartTruth*TMath::RadToDeg())/100.0 << endl;
}

void eventHistFiller::fillEventsFailedSelection(){

	for (int kkk=0; kkk<nPFOs; kkk++){
		for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
			if (abs(static_cast<EVENT::ReconstructedParticle*>(PFOCollection->getElementAt(kkk))->getType())==it->first) {
				getHistFromMap("nPFOsVsCosThetaFailType_"+it->second)->Fill(cosTruthTheta);
			}
		}
	}
	getHistFromMap("efficiencyVsCosThetaFailType_all")->Fill(cosTruthTheta);
	if (pfoCounter["Electron"]==0 && pfoCounter["Muon"]==0 && pfoCounter["Pion"]==0)
		getHistFromMap("efficiencyVsCosThetaFailType_noChargedParts")->Fill(cosTruthTheta);
	else if (pfoCounter["Electron"]==0 && pfoCounter["Muon"]>0 && pfoCounter["Pion"]==0)
		getHistFromMap("efficiencyVsCosThetaFailType_onlyMuon")->Fill(cosTruthTheta);
	else if (pfoCounter["Electron"]==0 && pfoCounter["Muon"]==0 && pfoCounter["Pion"]>0) 
		getHistFromMap("efficiencyVsCosThetaFailType_onlyPion")->Fill(cosTruthTheta);
	else if (pfoCounter["Electron"]>0 && pfoCounter["Muon"]==0 && pfoCounter["Pion"]==0) 
		getHistFromMap("efficiencyVsCosThetaFailType_onlyElectron")->Fill(cosTruthTheta);
	else
		getHistFromMap("efficiencyVsCosThetaFailType_chargePartsOfTwoOrMoreTypes")->Fill(cosTruthTheta);

}

vector <unsigned int> eventHistFiller::mergeClusters(){
	vector<unsigned int> idOfMergedParticles;
	if (config::vm.count("debug")){
		cout << "[DEBUG]\tpartCandidate pointer: " << partCandidate << endl;
		cout << "[DEBUG]\tpartCandidate->getClusters().size(): " << partCandidate->getClusters().size() << endl;
	}

	if (partCandidate->getClusters().size()==0){
		cout << "[ERROR]\tin eventHistFiller::mergeClusters(): partCandidate->getClusters().size()==0" << endl;
		return idOfMergedParticles;
	}

	double tmpPartCandidateClusterEnergy = partCandidate->getClusters()[0]->getEnergy();

	if (config::vm.count("debug"))
		cout << "[DEBUG]\tpartCandEnergy: " << partCandidate->getEnergy() << "; partCandClusterEnergy: " << tmpPartCandidateClusterEnergy << endl;
	partCandidate = CopyReconstructedParticle(partCandidate);

	if (useCaloInfoDuringEnergyMerging){
		partCandidate->setEnergy(tmpPartCandidateClusterEnergy);
	}

	const double* candMom = partCandidate->getMomentum();
	TVector3 vCandMom(candMom[0],candMom[1],candMom[2]);
	const double candTheta = vCandMom.Theta()*TMath::RadToDeg();
	const double candPhi = vCandMom.Phi()*TMath::RadToDeg();

	double tmpMomentum[3];
	tmpMomentum[0] = partCandidate->getMomentum()[0]; 
	tmpMomentum[1] = partCandidate->getMomentum()[1]; 
	tmpMomentum[2] = partCandidate->getMomentum()[2]; 

	if (config::vm.count("debug"))
		cout << "[DEBUG]\tnPFOs: " << nPFOs << endl;
	if (PFOmergeMap[mergeTag].size()>0 ) {
		for (int j=0; j<nPFOs; j++){
			if (config::vm.count("debug"))
				cout << "[DEBUG]\tj: " << j << "; idOfPartCandidate: " << idOfPartCandidate << endl;
			if (j==idOfPartCandidate) 
				continue;
			EVENT::ReconstructedParticle* recoPart = static_cast<EVENT::ReconstructedParticle*>(PFOCollection->getElementAt(j));
			const unsigned int partTypeLoop = recoPart->getType();
			const double* partMomLoop = recoPart->getMomentum();
			TVector3 vPartMomLoop(partMomLoop[0],partMomLoop[1],partMomLoop[2]);
			const double partThetaLoop = vPartMomLoop.Theta()*TMath::RadToDeg();
			const double partPhiLoop = vPartMomLoop.Phi()*TMath::RadToDeg();
			const double partEnergyLoop = recoPart->getEnergy();

			for (int iMerge=0; iMerge<PFOmergeMap[mergeTag].size(); iMerge++){
				PFOMergeSettings mergeSettings = PFOmergeMap[mergeTag][iMerge];
				if (abs(partTypeLoop)!=mergeSettings.pfoTypeToMerge)
					continue;

				bool passTheta = abs(candTheta-partThetaLoop)<mergeSettings.thetaCone;
				bool passPhi = abs(candPhi-partPhiLoop)<mergeSettings.phiCone;
			
				if (config::vm.count("debug"))
					cout << "[INFO] passTheta: " << passTheta << "; passPhi: " << passPhi << "; E: " << partEnergyLoop << " GeV; thetaCone: " << mergeSettings.thetaCone << "; dTheta: " << abs(candTheta-partThetaLoop) << "; phiCone: " << mergeSettings.phiCone << "; dPhi: " << abs(candPhi-partPhiLoop) << endl;

				if (passTheta && passPhi) {
				
					tmpMomentum[0] += recoPart->getMomentum()[0]; 
					tmpMomentum[1] += recoPart->getMomentum()[1]; 
					tmpMomentum[2] += recoPart->getMomentum()[2]; 
					// TODO implement recalculation of momentum based on Cluster from charged particle not from ReconstructedParticle
					// if (useCaloInfoDuringEnergyMerging){
					//         tmpMomentum[0] += partCandidate->getMomentum()[0];
					//         tmpMomentum[1] += partCandidate->getMomentum()[1];
					//         tmpMomentum[2] += partCandidate->getMomentum()[2];
					// }
					// else{
					//         tmpMomentum[0] += partCandidate->getClusters()[0]->getMomentum()[0];
					// }

					double tmpEnergy = recoPart->getEnergy();
					if (useCaloInfoDuringEnergyMerging){
						if (config::vm.count("debug"))
							cout << "[SAHSA]\tnClusters: " << recoPart->getClusters().size() << endl;
						tmpEnergy = recoPart->getClusters()[0]->getEnergy();
						// tmpEnergy = partCandidate->getEnergy();
					}

					tmpEnergy += partCandidate->getEnergy();

					partCandidate->setMomentum(tmpMomentum);
					partCandidate->setEnergy(tmpEnergy);

					idOfMergedParticles.push_back(j);
				}
			}
		}
	}

	return idOfMergedParticles;
}

void eventHistFiller::angularAndEnergyMatching(){

	// Here we already have most energetic candidate of the type
	bool passAngularMatching = false;
	bool passTransMomMatching = false;
	bool passEnergyMatching = false;
	double dPhi_angularMatchingCut = 0.02;
	double dTheta_angularMatchingCut = 0.01;
	double caloCut_energyMatching = 0.75;
	double transMomCut_energyMatching = 0.05;

	const double *partMom = partCandidate->getMomentum();
	TVector3 vPartMom(partMom[0],partMom[1],partMom[2]);
	const double partTheta = vPartMom.Theta()*TMath::RadToDeg();
	const double partPhi = vPartMom.Phi()*TMath::RadToDeg();
	const double cosPartTheta = TMath::Cos(partTheta*TMath::DegToRad());
	const double partEnergy = partCandidate->getEnergy();
	const double partPt = vPartMom.Pt();
	const double dPhiPartTruth = (partPhi-truthPhi)*TMath::DegToRad();
	const double dThetaPartTruth = (partTheta-truthTheta)*TMath::DegToRad();

	if (abs(dPhiPartTruth)<dPhi_angularMatchingCut && abs(dThetaPartTruth)<dTheta_angularMatchingCut)
		passAngularMatching = true;

	if (abs(partEnergy-truthEnergy)<caloCut_energyMatching*sqrt(truthEnergy))
		passEnergyMatching = true;

	if ((abs(truthPt-partPt))<transMomCut_energyMatching*truthPt)
		passTransMomMatching = true;

	passAngularMatching = passAngularMatching || (!applyAngularMatching);
	passEnergyMatching = passEnergyMatching || (!applyEnergyMatching);
	passTransMomMatching = passTransMomMatching || (!applyEnergyMatching);

	bool passFinalEnergyMomMatching = ((passEnergyMatching && (truthType==22 || useCaloCutInsteadMomentum==true)) || (passTransMomMatching==true && truthType!=22));

	if (config::vm.count("debug"))
		std::cout << "r:  pdg: " << std::setw(5) << partCandidate->getType()  << ": E: " << std::setw(6) << (round(100*partEnergy)/100.0) << ": pT: " << std::setw(6) << (round(100*partPt)/100.0) << "; theta: " << std::setw(6) << round(100*partTheta)/100.0 << "; phi: " << std::setw(6) << round(100*partPhi)/100.0 << "; AngularMatching: " << passAngularMatching << "; EnergyMatching: " << passEnergyMatching << "; TransMomMatching: " << passTransMomMatching << "; passFinalEnMom: " << passFinalEnergyMomMatching << std::endl;

	if (passFinalEnergyMomMatching == true && passAngularMatching == true){
		getHistFromMap("efficiencyVsTheta")->Fill(truthTheta);
		getHistFromMap("efficiencyVsCosTheta")->Fill(cosTruthTheta);
		getHistFromMap("efficiencyVsEnergy")->Fill(truthEnergy);

		getHistFromMap("PFO_passed_eff_E")->Fill(partEnergy);
		getHistFromMap("PFO_passed_eff_Pt")->Fill(partPt);
		getHistFromMap("PFO_passed_eff_dE")->Fill((partEnergy-truthEnergy)/truthEnergy);
		getHistFromMap("PFO_passed_eff_dPt")->Fill((truthPt-partPt)/truthPt);
		getHistFromMap("PFO_passed_eff_dTheta")->Fill(dThetaPartTruth);
		getHistFromMap("PFO_passed_eff_dPhi")->Fill(dPhiPartTruth);

	}
	else if (passAngularMatching==false && passFinalEnergyMomMatching == true){
		getHistFromMap("efficiencyVsCosThetaFailAngularMatching")->Fill(cosTruthTheta);
	}
	else if (passFinalEnergyMomMatching==false){
		getHistFromMap("efficiencyVsCosThetaFailEnergyMatching")->Fill(cosTruthTheta);
	}

}

void eventHistFiller::fillOtherHists(){

	getHistFromMap("totalRecoEnergy")->Fill(totalRecoEnergy);

	if (pfoCounter["Electron"]==0 && pfoCounter["Pion"]==0)
		getHistFromMap("efficiencyVsCosThetaCat1")->Fill(cosTruthTheta);
	else if (pfoCounter["Electron"]==1 && pfoCounter["Pion"]==0)
		getHistFromMap("efficiencyVsCosThetaCat2")->Fill(cosTruthTheta);
	else if (pfoCounter["Pion"]==1 && pfoCounter["Electron"]==0)
		getHistFromMap("efficiencyVsCosThetaCat3")->Fill(cosTruthTheta);
	else
		getHistFromMap("efficiencyVsCosThetaCat4")->Fill(cosTruthTheta);
	if (pfoCounter["Electron"]==1 && pfoCounter["Pion"]==0)
		getHistFromMap("efficiencyVsCosThetaCat5")->Fill(cosTruthTheta);


}

void eventHistFiller::setClusterMerging(string _mergeTag){
	if (!IsInVector<string>(_mergeTag,effType))
		cout << "[eventHistFiller::setClusterMerging]\tERROR mergeTag <" << _mergeTag << "> is not supported!!!" << endl; 
	else
		mergeTag = _mergeTag;
}


