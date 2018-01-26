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

int eventHistFiller::init(){

	if (config::vm.count("debug"))
		cout << "[INFO]\teventHistFiller::init()" << endl;

	TH1* tmpHist;
	tmpHist = new TH1I("nPFOs","Number of PFOs in event; Number of PFOs; Counts",5,0,5);
	histMap["nPFOs"] = tmpHist;
	tmpHist= new TH1I("PFOType","PFO particle type; Type; Counts",2200,0,2200); // max part.type = 2112 (neutron)
	histMap["PFOType"] = tmpHist;
	
	tmpHist = new TH1D("nPFOsVsTheta_all","nPFOs vs Theta; Theta; Counts per Event",180*2,0,180);
	histMap["nPFOsVsTheta_all"] = tmpHist;
	tmpHist = new TH1D("nPFOsVsCosTheta_all","nPFOs vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["nPFOsVsCosTheta_all"] = tmpHist;

	tmpHist = new TH1D("totalEnergyVsTheta","Sum of PFOs Energy vs Theta; Theta; Energy [GeV]",180*2,0,180);
	histMap["totalEnergyVsTheta"] = tmpHist;
	tmpHist = new TH1D("matchedEnergyVsTheta","Sum of matched PFOs Energy vs Theta; Theta; Energy [GeV]",180*2,0,180);
	histMap["matchedEnergyVsTheta"] = tmpHist;


	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
		tmpHist = new TH1D(("nPFOsVsTheta_"+it->second).c_str(),("n"+it->second+"s vs Theta; Theta; Counts per Event").c_str(),180*2,0,180);
		histMap["nPFOsVsTheta_"+it->second] = tmpHist;
		tmpHist = new TH1D(("nPFOsVsCosTheta_"+it->second).c_str(),("n"+it->second+"s vs Cos(#Theta); Cos(#Theta); Counts per Event").c_str(),180*2,-1,1);
		histMap["nPFOsVsCosTheta_"+it->second] = tmpHist;
		tmpHist = new TH1D(("nPFOsVsCosThetaFailType_"+it->second).c_str(),("n"+it->second+"s vs Cos(#Theta); Cos(#Theta); Counts per Event").c_str(),180*2,-1,1);
		histMap["nPFOsVsCosThetaFailType_"+it->second] = tmpHist;

	}
	
	tmpHist = new TH1D("nTruthPartsVsTheta","nTruthParts vs Theta; Theta; Counts per Event",180*2,0,180);
	histMap["nTruthPartsVsTheta"] = tmpHist;
	tmpHist = new TH1D("nTruthPartsVsCosTheta","nTruthParts vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["nTruthPartsVsCosTheta"] = tmpHist;
	tmpHist = new TH1D("nTruthPartsVsEnergy","nTruthParts vs Energy ;Energy [GeV]; Counts per Event",100,0.5,100.5); 
	histMap["nTruthPartsVsEnergy"] = tmpHist;

	tmpHist = new TH1D("efficiencyVsTheta","efficiency vs Theta; Theta; Counts per Event",180*2,0,180);
	histMap["efficiencyVsTheta"] = tmpHist;
	tmpHist = new TH1D("efficiencyVsCosTheta","efficiency vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["efficiencyVsCosTheta"] = tmpHist;
	tmpHist = new TH1D("efficiencyVsEnergy","efficiency vs Energy; Energy [GeV]; Counts per Event",100,0.5,100.5);
	histMap["efficiencyVsEnergy"] = tmpHist;


	tmpHist = new TH1D("efficiencyVsCosThetaFailType","efficiency vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["efficiencyVsCosThetaFailType"] = tmpHist;
	tmpHist = new TH1D("efficiencyVsCosThetaFailAngularMatching","efficiency vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["efficiencyVsCosThetaFailAngularMatching"] = tmpHist;
	tmpHist = new TH1D("efficiencyVsCosThetaFailEnergyMatching","efficiency vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["efficiencyVsCosThetaFailEnergyMatching"] = tmpHist;

	tmpHist = new TH1D("efficiencyVsCosThetaSum","efficiency vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["efficiencyVsCosThetaSum"] = tmpHist;

	for (int iCount=0; iCount<=9; iCount++){
		tmpHist = new TH1D(("efficiencyVsCosThetaCat"+DoubToStr(iCount)).c_str(),"efficiency vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
		histMap["efficiencyVsCosThetaCat"+DoubToStr(iCount)] = tmpHist;
	}

	tmpHist = new TH1I("truthParticle_isDecayedInTracker","isDecayedInTracker; isDecayedInTracker flag; Counts",2,-0.5,1.5);
	histMap["truthParticle_isDecayedInTracker"] = tmpHist;

	// cout << "debug1" << endl;

	tmpHist = new TH1D("mergedEnergy","Merged energy; E [GeV]; Counts",1250,0,125);
	histMap["mergedEnergy"] = tmpHist;
	tmpHist = new TH1D("totalRecoEnergy","Total Reconstructed energy; E [GeV]; Counts",1250,0,125);
	histMap["totalRecoEnergy"] = tmpHist;

	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++){
		tmpHist = new TH1D(("phiResolution_"+it->second).c_str(),(it->second+" Phi resolution; dPhi [rad]; Counts").c_str(),20000,-0.2,0.2);
		histMap["phiResolution_"+it->second] = tmpHist;
		tmpHist = new TH1D(("thetaResolution_"+it->second).c_str(),(it->second+" Theta resolution; dTheta [rad]; Counts").c_str(),400,-0.01,0.01);
		histMap["thetaResolution_"+it->second] = tmpHist;
		tmpHist = new TH1D(("energyResolution_"+it->second).c_str(),(it->second+" Energy resolution; E [GeV]; Counts").c_str(),1250,0,125);
		histMap["energyResolution_"+it->second] = tmpHist;
	}

	// cout << "debug2" << endl;

	vector<PFOMergeSettings> tmpVec;
	PFOMergeSettings tmpPFOMergeSettings;

	PFOmergeMap["nominal"] = tmpVec;
	tmpVec.clear();

	tmpPFOMergeSettings.pfoTypeToMerge = 22;
	tmpPFOMergeSettings.thetaCone = 0.01*180.0/TMath::Pi();
	tmpPFOMergeSettings.phiCone = 0.2*180.0/TMath::Pi();
	tmpVec.push_back(tmpPFOMergeSettings);
	PFOmergeMap["photonMerge"] = tmpVec;

	tmpPFOMergeSettings.pfoTypeToMerge = 2112;
	tmpPFOMergeSettings.thetaCone = 0.01*180.0/TMath::Pi();
	tmpPFOMergeSettings.phiCone = 0.2*180.0/TMath::Pi();
	tmpVec.push_back(tmpPFOMergeSettings);
	PFOmergeMap["photonAndNeutralMerge"] = tmpVec;
	tmpVec.clear();

	tmpPFOMergeSettings.pfoTypeToMerge = 22;
	tmpPFOMergeSettings.thetaCone = 0.01*180.0/TMath::Pi();
	tmpPFOMergeSettings.phiCone = 0.2*180.0/TMath::Pi();
	tmpVec.push_back(tmpPFOMergeSettings);
	tmpPFOMergeSettings.pfoTypeToMerge = 2112;
	tmpPFOMergeSettings.thetaCone = 0.035*180.0/TMath::Pi();
	tmpPFOMergeSettings.phiCone = 0.2*180.0/TMath::Pi();
	tmpVec.push_back(tmpPFOMergeSettings);
	PFOmergeMap["photonAndNeutralLooseMerge"] = tmpVec;
	tmpVec.clear();

	tmpPFOMergeSettings.pfoTypeToMerge = 22;
	tmpPFOMergeSettings.thetaCone = 0.01*180.0/TMath::Pi();
	tmpPFOMergeSettings.phiCone = 0.2*180.0/TMath::Pi();
	tmpPFOMergeSettings.phiConeMomentumDep = 60.0;
	tmpVec.push_back(tmpPFOMergeSettings);
	PFOmergeMap["photonMergeMomentumDep"] = tmpVec;
	tmpVec.clear();

	nSelectecTruthParticles = 0;
	return 0;

}

int eventHistFiller::fillEvent(const EVENT::LCEvent* event){
	try {
		PFOCollection = event->getCollection(PFOCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|eventHistFiller]\tCan't find collection: " << PFOCollectionName << endl;
		return -1;
	}

	if (config::vm.count("debug"))
		cout << "[INFO]	eventHistFiller::fillEvent: " << event->getEventNumber() << endl;

	double truthTheta = -666;  
	double truthPhi = -666;  
	double cosTruthTheta = -666; 
	double truthEnergy = -666;
	bool reclusteringIsDone = false;

	EVENT::MCParticle* genPart = truthCondition::instance()->getGunParticle();
	const double *partMom = genPart->getMomentum();
	TVector3 v1;
	v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
	double partTheta = 180.*v1.Theta()/TMath::Pi();
	double cosPartTheta = TMath::Cos(TMath::Pi()*partTheta/180.);
	truthTheta = partTheta;
	truthPhi = 180.*v1.Phi()/TMath::Pi();
	cosTruthTheta = cosPartTheta;

	getHistFromMap("truthParticle_isDecayedInTracker")->Fill(genPart->isDecayedInTracker());
	nSelectecTruthParticles++;
	getHistFromMap("nTruthPartsVsCosTheta")->Fill(cosPartTheta);
	getHistFromMap("nTruthPartsVsTheta")->Fill(partTheta);
	truthEnergy = genPart->getEnergy();
	getHistFromMap("nTruthPartsVsEnergy")->Fill(truthEnergy);

	vector<EVENT::ReconstructedParticle*> recoPFOs = getObjVecFromCollection(PFOCollection);
	getHistFromMap("nPFOs")->Fill(recoPFOs.size());
	if (recoPFOs.size()==0) 
		return 0; // no reco PFOs

	map <string, unsigned int> pfoCounter;
	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++)
		pfoCounter[it->second] = 0;

	bool efficiencyHistWasAlreadyFilled = false;
	for (int i=0; i<recoPFOs.size(); i++){

		int pfoType = abs(recoPFOs[i]->getType());
		const double *partMom = recoPFOs[i]->getMomentum();
		TVector3 v1;
		v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
		double partTheta = 180.*v1.Theta()/TMath::Pi();
		double partPhi = 180.*v1.Phi()/TMath::Pi();
		double cosPartTheta = TMath::Cos(TMath::Pi()*partTheta/180.);
		double pfoE = recoPFOs[i]->getEnergy();
		double dPhi = TMath::Pi()*(partPhi-truthPhi)/180.0;
		double dTheta = TMath::Pi()*(partTheta-truthTheta)/180.0;

		getHistFromMap("PFOType")->Fill(pfoType);
		getHistFromMap("nPFOsVsCosTheta_all")->Fill(cosTruthTheta);
		getHistFromMap("nPFOsVsTheta_all")->Fill(truthTheta);
		getHistFromMap("totalEnergyVsTheta")->Fill(truthTheta,recoPFOs[i]->getEnergy());

		for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
			if (abs(pfoType)==it->first) {
				getHistFromMap("nPFOsVsCosTheta_"+it->second)->Fill(cosTruthTheta);
				getHistFromMap("nPFOsVsTheta_"+it->second)->Fill(truthTheta);
				getHistFromMap("thetaResolution_"+it->second)->Fill(dTheta);
				getHistFromMap("phiResolution_"+it->second)->Fill(dPhi);
				getHistFromMap("energyResolution_"+it->second)->Fill(pfoE);
				pfoCounter[it->second]++;
			}
		}

		if (efficiencyHistWasAlreadyFilled==false && i==(recoPFOs.size()-1)){ 
			bool partOfTypeIsPresent = false;
			for (int kkk=0; kkk<recoPFOs.size(); kkk++){
				if (abs(recoPFOs[kkk]->getType())==abs(genPart->getPDG()))
					partOfTypeIsPresent = true;
			}
			if (partOfTypeIsPresent==false){
				getHistFromMap("efficiencyVsCosThetaFailType")->Fill(cosTruthTheta);
				for (int kkk=0; kkk<recoPFOs.size(); kkk++){
					for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
						if (abs(recoPFOs[kkk]->getType())==it->first) {
							getHistFromMap("nPFOsVsCosThetaFailType_"+it->second)->Fill(cosTruthTheta);
						}
					}
				}
			}
		}

		if (pfoType!=abs(genPart->getPDG()) || efficiencyHistWasAlreadyFilled==true)
			continue;

		// cout << "PFOmergeMap["<<mergeTag<<"].size(): " << PFOmergeMap[mergeTag].size() << endl;
		if (PFOmergeMap[mergeTag].size()>0) {
			for (int j=0; j<recoPFOs.size(); j++){
				if (i==j) 
					continue;
				for (int iMerge=0; iMerge<PFOmergeMap[mergeTag].size(); iMerge++){
					PFOMergeSettings mergeSettings = PFOmergeMap[mergeTag][iMerge];
					if (abs(recoPFOs[j]->getType())!=mergeSettings.pfoTypeToMerge)
						continue;
					const double *partMom = recoPFOs[i]->getMomentum();
					TVector3 v1;
					v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
					double partTheta_i = 180.*v1.Theta()/TMath::Pi();
					double partPhi_i = 180.*v1.Phi()/TMath::Pi();

					partMom = recoPFOs[j]->getMomentum();
					v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
					double partTheta_j = 180.*v1.Theta()/TMath::Pi();
					double partPhi_j = 180.*v1.Phi()/TMath::Pi();

					bool passTheta = abs(partTheta_i-partTheta_j)<mergeSettings.thetaCone;
					bool passPhi = abs(partPhi_i-partPhi_j)<mergeSettings.phiCone;
					if (mergeSettings.phiConeMomentumDep!=0.0){
						passPhi = abs(partPhi_i-partPhi_j)<(mergeSettings.phiConeMomentumDep/truthEnergy + mergeSettings.phiCone);
						if (config::vm.count("debug")){
							cout << "[CHECK] phiCone: " << mergeSettings.phiConeMomentumDep/truthEnergy + mergeSettings.phiCone << endl;
							cout << "[CHECK] phiCone1: " << mergeSettings.phiConeMomentumDep/truthEnergy << endl;
							cout << "[CHECK] phiCone2: " << mergeSettings.phiCone << endl;
						}
					}

					if (config::vm.count("debug"))
						cout << "[INFO] passTheta: " << passTheta << "; passPhi: " << passPhi << "; E: " << recoPFOs[j]->getEnergy() << " GeV; thetaCone: " << mergeSettings.thetaCone << "; dTheta: " << abs(partTheta_i-partTheta_j) << "; phiCone: " << mergeSettings.phiCone << "; dPhi: " << abs(partPhi_i-partPhi_j) << endl;
					if (passTheta && passPhi) {
						pfoE+=recoPFOs[j]->getEnergy();
					}
				}
			}
		}
		if (abs(pfoE-truthEnergy)<0.75*sqrt(truthEnergy)) {
			// WARNING OLD angular requirements:
			// if (recoPFOs.size()==1 && (abs(dPhi)>0.004 || abs(dTheta)>0.002)) {
			// WARNING NEW angular requirements:
			if (recoPFOs.size()==1 && (abs(dPhi)>0.02 || abs(dTheta)>0.01)) {
				getHistFromMap("efficiencyVsCosThetaFailAngularMatching")->Fill(cosTruthTheta);
				continue;
			}
			getHistFromMap("efficiencyVsTheta")->Fill(truthTheta);
			getHistFromMap("efficiencyVsCosTheta")->Fill(cosTruthTheta);
			getHistFromMap("efficiencyVsEnergy")->Fill(truthEnergy);
			// getHistFromMap("matchedEnergyVsTheta"]->Fill(truthTheta,pfoE);
			efficiencyHistWasAlreadyFilled=true;
			if (config::vm.count("debug")){
				cout << "[INFO]	eventHistFiller::fillEvent(" << event->getEventNumber() << ") eff filled with energy: " << pfoE << " GeV" << endl;
			}
			getHistFromMap("mergedEnergy")->Fill(pfoE);
		}
		else
			getHistFromMap("efficiencyVsCosThetaFailEnergyMatching")->Fill(cosTruthTheta);
	}

	// std::map<unsigned int, std::string> pfoTypeIntStringMap = {{11,"Electron"}, {13,"Muon"},{22,"Photon"},{211,"Pion"},{2112,"NeutralHadron"}};
	if (pfoCounter["Electron"]==1 && pfoCounter["Photon"]==0 && pfoCounter["Pion"]==0 && pfoCounter["NeutralHadron"]==0 && pfoCounter["Muon"]==0)
		getHistFromMap("efficiencyVsCosThetaCat1")->Fill(cosTruthTheta);
	else if (pfoCounter["Electron"]==1 && pfoCounter["Photon"]==1 && pfoCounter["Pion"]==0 && pfoCounter["NeutralHadron"]==0 && pfoCounter["Muon"]==0)
		getHistFromMap("efficiencyVsCosThetaCat2")->Fill(cosTruthTheta);
	else if (pfoCounter["Electron"]==1 && pfoCounter["Photon"]>=2 && pfoCounter["Pion"]==0 && pfoCounter["NeutralHadron"]==0 && pfoCounter["Muon"]==0)
		getHistFromMap("efficiencyVsCosThetaCat3")->Fill(cosTruthTheta);
	else if (pfoCounter["Electron"]==0 && pfoCounter["Photon"]==0 && pfoCounter["Pion"]==1 && pfoCounter["NeutralHadron"]==0 && pfoCounter["Muon"]==0)
		getHistFromMap("efficiencyVsCosThetaCat4")->Fill(cosTruthTheta);
	else if (pfoCounter["Electron"]==0 && pfoCounter["Photon"]==1 && pfoCounter["Pion"]==1 && pfoCounter["NeutralHadron"]==0 && pfoCounter["Muon"]==0)
		getHistFromMap("efficiencyVsCosThetaCat5")->Fill(cosTruthTheta);
	else if (pfoCounter["Electron"]==0 && pfoCounter["Photon"]>=2 && pfoCounter["Pion"]==1 && pfoCounter["NeutralHadron"]==0 && pfoCounter["Muon"]==0)
		getHistFromMap("efficiencyVsCosThetaCat6")->Fill(cosTruthTheta);
	else if (pfoCounter["Electron"]==0 && pfoCounter["Photon"]==1 && pfoCounter["Pion"]==0 && pfoCounter["NeutralHadron"]==0 && pfoCounter["Muon"]==0)
		getHistFromMap("efficiencyVsCosThetaCat7")->Fill(cosTruthTheta);
	else if (pfoCounter["Electron"]==0 && pfoCounter["Photon"]>=2 && pfoCounter["Pion"]==0 && pfoCounter["NeutralHadron"]==0 && pfoCounter["Muon"]==0)
		getHistFromMap("efficiencyVsCosThetaCat8")->Fill(cosTruthTheta);
	else if (pfoCounter["NeutralHadron"]>=1)
		getHistFromMap("efficiencyVsCosThetaCat9")->Fill(cosTruthTheta);
		


}


int eventHistFiller::writeToFile(TFile* outFile){

	if (config::vm.count("debug"))
		cout << "[INFO]	eventHistFiller::writeToFile(" << outFile->GetName() << ")" << endl;

	getHistFromMap("totalEnergyVsTheta")->Divide(getHistFromMap("nTruthPartsVsTheta"));
	getHistFromMap("nTruthPartsVsTheta")->Sumw2();
	getHistFromMap("nTruthPartsVsCosTheta")->Sumw2();
	getHistFromMap("nTruthPartsVsEnergy")->Sumw2();

	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
		getHistFromMap("nPFOsVsCosTheta_"+it->second)->Sumw2();
		getHistFromMap("nPFOsVsCosTheta_"+it->second)->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
		getHistFromMap("nPFOsVsTheta_"+it->second)->Sumw2();
		getHistFromMap("nPFOsVsTheta_"+it->second)->Divide(getHistFromMap("nTruthPartsVsTheta"));
	}

	getHistFromMap("nPFOsVsTheta_all")->Sumw2();
	getHistFromMap("nPFOsVsTheta_all")->Divide(getHistFromMap("nTruthPartsVsTheta"));
	getHistFromMap("nPFOsVsCosTheta_all")->Sumw2();
	getHistFromMap("nPFOsVsCosTheta_all")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	// getHistFromMap("totalEnergyVsTheta")->Sumw2();

	getHistFromMap("efficiencyVsTheta")->Sumw2();
	getHistFromMap("efficiencyVsTheta")->Divide(getHistFromMap("nTruthPartsVsTheta"));
	getHistFromMap("efficiencyVsCosTheta")->Sumw2();
	getHistFromMap("efficiencyVsCosTheta")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	getHistFromMap("efficiencyVsEnergy")->Sumw2();
	getHistFromMap("efficiencyVsEnergy")->Divide(getHistFromMap("nTruthPartsVsEnergy"));

	getHistFromMap("efficiencyVsCosThetaFailType")->Sumw2();
	getHistFromMap("efficiencyVsCosThetaFailType")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	getHistFromMap("efficiencyVsCosThetaFailEnergyMatching")->Sumw2();
	getHistFromMap("efficiencyVsCosThetaFailEnergyMatching")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	getHistFromMap("efficiencyVsCosThetaFailAngularMatching")->Sumw2();
	getHistFromMap("efficiencyVsCosThetaFailAngularMatching")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	
	getHistFromMap("efficiencyVsCosThetaSum")->Add(getHistFromMap("efficiencyVsCosTheta"));
	getHistFromMap("efficiencyVsCosThetaSum")->Add(getHistFromMap("efficiencyVsCosThetaFailType"));
	getHistFromMap("efficiencyVsCosThetaSum")->Add(getHistFromMap("efficiencyVsCosThetaFailEnergyMatching"));
	getHistFromMap("efficiencyVsCosThetaSum")->Add(getHistFromMap("efficiencyVsCosThetaFailAngularMatching"));

	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
		getHistFromMap("nPFOsVsCosThetaFailType_"+it->second)->Sumw2();
		getHistFromMap("nPFOsVsCosThetaFailType_"+it->second)->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
		getHistFromMap("nPFOsVsCosThetaFailType_"+it->second)->Divide(getHistFromMap("efficiencyVsCosThetaFailType"));
	}


	for (int iCount=1; iCount<=9; iCount++){
		getHistFromMap("efficiencyVsCosThetaCat"+DoubToStr(iCount))->Sumw2();
		getHistFromMap("efficiencyVsCosThetaCat"+DoubToStr(iCount))->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
		getHistFromMap("efficiencyVsCosThetaCat0")->Add(getHistFromMap("efficiencyVsCosThetaCat"+DoubToStr(iCount)));
	}

	if (!outFile->IsOpen()){
		cout << "[ERROR|writeToFile]\tno output file is found!" << endl;
		return -1;
	}
	outFile->cd();
	TDirectory *mainDir = outFile->mkdir(outDirName.c_str());
	mainDir->cd();

	map<string,unsigned int> prefixCounter;
	map<string,string> namePrefixMap;
	map<string,bool> isPrefixSubdirCreated;
	map<string,string> nameWithoutPrefixMap;
	for(auto const &it : histMap) {
		string histName = it.first;
		vector<string> tmpStrVec = GetSplittedWords(histName,"_");
		if (tmpStrVec.size()<2) 
			continue;
		string prefix = "";
		for (int i=0; i<tmpStrVec.size()-1; i++){
			if (i==tmpStrVec.size()-2)
				prefix += tmpStrVec[i];
			else
				prefix += tmpStrVec[i] + "_";
		}
		nameWithoutPrefixMap[histName] = tmpStrVec[tmpStrVec.size()-1];
		prefixCounter[prefix] += 1;
		isPrefixSubdirCreated[prefix] = false;
		namePrefixMap[histName] = prefix;
	}
	

	for(auto const &it : histMap) {
		string histName = it.first;
		string prefix = namePrefixMap[histName];
		if (prefixCounter[prefix]<2){
			mainDir->cd();
			it.second->Write();
		}
		else{
			if (isPrefixSubdirCreated[prefix]==false){
				mainDir->mkdir(prefix.c_str());
				isPrefixSubdirCreated[prefix]=true;
			}
			mainDir->cd(prefix.c_str());
			it.second->SetName(nameWithoutPrefixMap[histName].c_str());
			it.second->Write();
			mainDir->cd();
		}
	}
	outFile->cd();
	return 0;


	// objectFill::writeToFile(outFile);
	// for (auto it2 = inVec.begin(); it2 != inVec.end(); it2++) {
	//         bool isAllowed = false;
	//         for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
	//                 if (it->first==(*it2))
	//                         isAllowed = true;
	//         }
	//         if (isAllowed == false) {
	//                 cout << "[ERROR]\teventHistFiller::checkPfoType\tnot allowed PFO type" << endl;
	//                 exit(1);
	//         }
	// }
}
/*===========================================================================*/
/*===============================[ function implementations ]===============================*/
/*===========================================================================*/

// void eventHistFiller::get90(const TH1 *const pTH1F){
//     static const float FLOAT_MAX(std::numeric_limits<float>::max());
//
//     double fractionEventsToIgnorePerSide = 0.05;
//
//     if (NULL == pTH1F)
//         return;
//
//     if (5 > pTH1F->GetEntries())
//     {
//         std::cout << pTH1F->GetName() << " (" << pTH1F->GetEntries() << " entries) - skipped" << std::endl;
//         return;
//     }
//
//     // Calculate raw properties of distribution
//     float sum = 0., total = 0.;
//     double sx = 0., sxx = 0.;
//     const unsigned int nbins(pTH1F->GetNbinsX());
//
//     for (unsigned int i = 0; i <= nbins; ++i)
//     {
//         const float binx(pTH1F->GetBinLowEdge(i) + (0.5 * pTH1F->GetBinWidth(i)));
//         const float yi(pTH1F->GetBinContent(i));
//         sx += yi * binx;
//         sxx += yi * binx * binx;
//         total += yi;
//     }
//
//     const float rawMean(sx / total);
//     const float rawMeanSquared(sxx / total);
//     const float rawRms(std::sqrt(rawMeanSquared - rawMean * rawMean));
//
//     sum = 0.;
//     unsigned int is0 = 0;
//
//     for (unsigned int i = 0; (i <= nbins) && (sum < total * fractionEventsToIgnorePerSide ); ++i)
//     {
//         sum += pTH1F->GetBinContent(i);
//         is0 = i;
//     }
//
//     // Calculate truncated properties
//     float rmsmin(FLOAT_MAX), sigma(FLOAT_MAX), sigmasigma(FLOAT_MAX), frac(FLOAT_MAX), efrac(FLOAT_MAX), mean(FLOAT_MAX), low(FLOAT_MAX), rms(FLOAT_MAX);
//     float high(0.f);
//
//     for (unsigned int istart = 0; istart <= is0; ++istart)
//     {
//         double sumn = 0.;
//         double csum = 0.;
//         double sumx = 0.;
//         double sumxx = 0.;
//         unsigned int iend = 0;
//
//         for (unsigned int i = istart; (i <= nbins) && (csum < (1.0 - fractionEventsToIgnorePerSide) * total); ++i)
//         {
//             const float binx(pTH1F->GetBinLowEdge(i) + (0.5 * pTH1F->GetBinWidth(i)));
//             const float yi(pTH1F->GetBinContent(i));
//             csum += yi;
//
//             if (sumn < (1.0 - fractionEventsToIgnorePerSide) * total)
//             {
//                 sumn += yi;
//                 sumx += yi * binx;
//                 sumxx+= yi * binx * binx;
//                 iend = i;
//             }
//         }
//
//         const float localMean(sumx / sumn);
//         const float localMeanSquared(sumxx / sumn);
//         const float localRms(std::sqrt(localMeanSquared - localMean * localMean));
//
//         if (localRms < rmsmin)
//         {
//             mean = localMean;
//             rms = localRms;
//             low = pTH1F->GetBinLowEdge(istart);
//             high = pTH1F->GetBinLowEdge(iend);
//             rmsmin = localRms;
//
//
//             sigma = rms;
//             sigmasigma = sigma / std::sqrt(total);
//         }
//     }
//
//     std::cout << pTH1F->GetName() << " (" << pTH1F->GetEntries() << " entries), rawrms: " << rawRms << ", rms90: " << rmsmin
//               << " (" << low << "-" << high << "), mean: " << mean << ", sigma: " << sigma << "+-" << sigmasigma << endl << endl;
//
// }
		

void eventHistFiller::setClusterMerging(string _mergeTag){
	if (!IsInVector<string>(_mergeTag,effType))
		cout << "[eventHistFiller::setClusterMerging]\tERROR mergeTag <" << _mergeTag << "> is not supported!!!" << endl;
	else
		mergeTag = _mergeTag;
}
