#include "eventHistFiller.h"

vector<string> typeOfCutFails = {"FailType","FailAngularMatching","FailEnergyMatching"};

int eventHistFiller::init(){

	if (config::vm.count("debug")){
		cout << "[INFO]\teventHistFiller::init(); outDirName: " << outDirName << endl;
	}

	createTH1I("nPFOs","Number of PFOs in event; Number of PFOs; Counts",5,0,5);
	createTH1I("PFOType","PFO particle type; Type; Counts",2200,0,2200); // max part.type = 2112 (neutron)
	
	createTH1D("nPFOsVsTheta_all","nPFOs vs Theta; Theta; Counts per Event",180*2,0,180);
	createTH1D("nPFOsVsCosTheta_all","nPFOs vs cos(#theta); cos(#theta); Counts per Event",15*2,-1,1);

	createTH1D("totalEnergyVsTheta","Sum of PFOs Energy vs Theta; Theta; Energy [GeV]",180*2,0,180);
	createTH1D("matchedEnergyVsTheta","Sum of matched PFOs Energy vs Theta; Theta; Energy [GeV]",180*2,0,180);

	for (auto const &it : typeOfCutFails)


	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
		createTH1D(("nPFOsVsTheta_"+it->second).c_str(),("n"+it->second+"s vs Theta; Theta; Counts per Event").c_str(),180*2,0,180);
		createTH1D(("nPFOsVsCosTheta_"+it->second).c_str(),("n"+it->second+"s vs cos(#theta); cos(#theta); Counts per Event").c_str(),15*2,-1,1);
		createTH1D(("nPFOsVsCosThetaFailType_"+it->second).c_str(),("n"+it->second+"s vs cos(#theta); cos(#theta); Counts per Event").c_str(),15*2,-1,1);

	}
	
	createTH1D("nTruthPartsVsTheta","nTruthParts vs Theta; Theta; Counts per Event",180*2,0,180);
	createTH1D("nTruthPartsVsCosTheta","nTruthParts vs cos(#theta); cos(#theta); Counts per Event",15*2,-1,1);
	createTH1D("nTruthPartsVsEnergy","nTruthParts vs Energy ;Energy [GeV]; Counts per Event",100,0.5,100.5); 

	createTH1D("efficiencyVsTheta","efficiency vs Theta; Theta; Counts per Event",180*2,0,180);
	createTH1D("efficiencyVsCosTheta","efficiency vs cos(#theta); cos(#theta); Counts per Event",15*2,-1,1);
	createTH1D("efficiencyVsEnergy","efficiency vs Energy; Energy [GeV]; Counts per Event",100,0.5,100.5);
	createTH1D("efficiencyVsEnergy_onlyType","efficiency vs Energy; Energy [GeV]; Counts per Event",100,0.5,100.5);


	createTH1D("efficiencyVsCosThetaFailType_all","efficiency vs cos(#theta); cos(#theta); Counts per Event",15*2,-1,1);
	createTH1D("efficiencyVsCosThetaFailType_onlyPion","efficiency vs cos(#theta); cos(#theta); Counts per Event",15*2,-1,1);
	createTH1D("efficiencyVsCosThetaFailType_onlyMuon","efficiency vs cos(#theta); cos(#theta); Counts per Event",15*2,-1,1);
	createTH1D("efficiencyVsCosThetaFailType_onlyElectron","efficiency vs cos(#theta); cos(#theta); Counts per Event",15*2,-1,1);
	createTH1D("efficiencyVsCosThetaFailType_noChargedParts","efficiency vs cos(#theta); cos(#theta); Counts per Event",15*2,-1,1);
	createTH1D("efficiencyVsCosThetaFailType_chargePartsOfTwoOrMoreTypes","efficiency vs cos(#theta); cos(#theta); Counts per Event",15*2,-1,1);
	createTH1D("efficiencyVsCosThetaFailAngularMatching","efficiency vs cos(#theta); cos(#theta); Counts per Event",15*2,-1,1);
	createTH1D("efficiencyVsCosThetaFailEnergyMatching","efficiency vs cos(#theta); cos(#theta); Counts per Event",15*2,-1,1);

	createTH1D("efficiencyVsCosThetaSum","efficiency vs cos(#theta); cos(#theta); Counts per Event",15*2,-1,1);

	for (int iCount=0; iCount<=9; iCount++){
		createTH1D(("efficiencyVsCosThetaCat"+DoubToStr(iCount)).c_str(),"efficiency vs cos(#theta); cos(#theta); Counts per Event",15*2,-1,1);
	}

	createTH1I("truthParticle_isDecayedInTracker","isDecayedInTracker; isDecayedInTracker flag; Counts",2,-0.5,1.5);
	createTH1D("truthParticle_convRadius","convRadius; convRadius; Counts",250,0,2500);

	// cout << "debug1" << endl;

	createTH1D("PFO_passed_eff_E","Candidate Energy After Reclustering; E [GeV]; Counts",1250,0,125);
	createTH1D("PFO_passed_eff_Pt","Candidate pT After Reclustering; p_T [GeV]; Counts",1250,0,125);
	createTH1D("candidateEnergyBeforeRecl","Candidate Energy Before Reclustering; E [GeV]; Counts",1250,0,125);
	createTH1D("totalRecoEnergy","Total Reconstructed energy; E [GeV]; Counts",1250,0,125);

	createTH1D("PFO_passed_eff_dTheta","#Delta #theta; Theta [rad]; Counts",1000,-0.025,0.025);
	createTH1D("PFO_passed_eff_dPhi","#Delta #phi; Phi [rad]; Counts",40000,-0.5,0.5);
	createTH1D("PFO_passed_eff_dPt","#(Delta pt)/pt; dPt/Pt; Counts",300,-0.15,0.15);
	createTH1D("PFO_passed_eff_dE","#(Delta E)/E; dE/E; Counts",200,-1.0,1.0);
	

	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++){
		createTH1D(("phiResolution_"+it->second).c_str(),(it->second+" Phi resolution; dPhi [rad]; Counts").c_str(),20000,-0.2,0.2);
		createTH1D(("thetaResolution_"+it->second).c_str(),(it->second+" Theta resolution; dTheta [rad]; Counts").c_str(),2000,-0.05,0.05);
		createTH1D(("energyResolution_"+it->second).c_str(),(it->second+" Energy resolution; E [GeV]; Counts").c_str(),625,0,125);
		createTH1D(("energyResolution2_"+it->second).c_str(),(it->second+" Energy resolution; E [GeV]; Counts").c_str(),625,0,125);
		createTH1D(("energyResolution3_"+it->second).c_str(),(it->second+" Energy resolution; E [GeV]; Counts").c_str(),625,0,125);
		createTH1D(("energyResolution4_"+it->second).c_str(),(it->second+" Energy resolution; E [GeV]; Counts").c_str(),625,0,125);
	}

	// cout << "debug2" << endl;

	vector<PFOMergeSettings> tmpVec;
	PFOMergeSettings tmpPFOMergeSettings;

	PFOmergeMap["nominal"] = tmpVec;
	tmpVec.clear();

	tmpPFOMergeSettings.pfoTypeToMerge = 22;
	tmpPFOMergeSettings.thetaCone = 0.01*TMath::RadToDeg();
	tmpPFOMergeSettings.phiCone = 0.2*TMath::RadToDeg();
	tmpVec.push_back(tmpPFOMergeSettings);
	PFOmergeMap["photonMerge"] = tmpVec;

	tmpPFOMergeSettings.pfoTypeToMerge = 2112;
	tmpPFOMergeSettings.thetaCone = 0.01*TMath::RadToDeg();
	tmpPFOMergeSettings.phiCone = 0.2*TMath::RadToDeg();
	tmpVec.push_back(tmpPFOMergeSettings);
	PFOmergeMap["photonAndNeutralMerge"] = tmpVec;
	tmpVec.clear();

	tmpPFOMergeSettings.pfoTypeToMerge = 22;
	tmpPFOMergeSettings.thetaCone = 0.01*TMath::RadToDeg();
	tmpPFOMergeSettings.phiCone = 0.2*TMath::RadToDeg();
	tmpVec.push_back(tmpPFOMergeSettings);
	tmpPFOMergeSettings.pfoTypeToMerge = 2112;
	tmpPFOMergeSettings.thetaCone = 0.035*TMath::RadToDeg();
	tmpPFOMergeSettings.phiCone = 0.2*TMath::RadToDeg();
	tmpVec.push_back(tmpPFOMergeSettings);
	PFOmergeMap["photonAndNeutralLooseMerge"] = tmpVec;
	tmpVec.clear();

	return 0;

}


int eventHistFiller::writeToFile(TFile* outFile){

	if (config::vm.count("debug"))
		cout << "[INFO]	eventHistFiller::writeToFile(" << outFile->GetName() << ")" << endl;

	getHistFromMap("nTruthPartsVsTheta")->Sumw2();
	getHistFromMap("nTruthPartsVsCosTheta")->Sumw2();
	getHistFromMap("nTruthPartsVsEnergy")->Sumw2();

	getHistFromMap("totalEnergyVsTheta")->Sumw2();
	getHistFromMap("totalEnergyVsTheta")->Divide(getHistFromMap("nTruthPartsVsTheta"));

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

	createTEff("efficiencyVsTheta","nTruthPartsVsTheta");
	createTEff("efficiencyVsCosTheta","nTruthPartsVsCosTheta");
	createTEff("efficiencyVsEnergy","nTruthPartsVsEnergy");
	createTEff("efficiencyVsEnergy_onlyType","nTruthPartsVsEnergy");

	getHistFromMap("efficiencyVsTheta")->Sumw2();
	getHistFromMap("efficiencyVsTheta")->Divide(getHistFromMap("nTruthPartsVsTheta"));
	getHistFromMap("efficiencyVsCosTheta")->Sumw2();
	getHistFromMap("efficiencyVsCosTheta")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	getHistFromMap("efficiencyVsEnergy")->Sumw2();
	getHistFromMap("efficiencyVsEnergy")->Divide(getHistFromMap("nTruthPartsVsEnergy"));
	getHistFromMap("efficiencyVsEnergy_onlyType")->Sumw2();
	getHistFromMap("efficiencyVsEnergy_onlyType")->Divide(getHistFromMap("nTruthPartsVsEnergy"));

	// getHistFromMap("efficiencyVsCosThetaFailType_all")->Sumw2();
	// getHistFromMap("efficiencyVsCosThetaFailType_all")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	// getHistFromMap("efficiencyVsCosThetaFailType_noChargedParts")->Sumw2();
	// getHistFromMap("efficiencyVsCosThetaFailType_noChargedParts")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	// getHistFromMap("efficiencyVsCosThetaFailType_pion")->Sumw2();
	// getHistFromMap("efficiencyVsCosThetaFailType_pion")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	// getHistFromMap("efficiencyVsCosThetaFailType_muon")->Sumw2();
	// getHistFromMap("efficiencyVsCosThetaFailType_muon")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	// getHistFromMap("efficiencyVsCosThetaFailEnergyMatching")->Sumw2();
	// getHistFromMap("efficiencyVsCosThetaFailEnergyMatching")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	// getHistFromMap("efficiencyVsCosThetaFailAngularMatching")->Sumw2();
	// getHistFromMap("efficiencyVsCosThetaFailAngularMatching")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	//
	// getHistFromMap("efficiencyVsCosThetaSum")->Add(getHistFromMap("efficiencyVsCosTheta"));
	// getHistFromMap("efficiencyVsCosThetaSum")->Add(getHistFromMap("efficiencyVsCosThetaFailType_all"));
	// getHistFromMap("efficiencyVsCosThetaSum")->Add(getHistFromMap("efficiencyVsCosThetaFailEnergyMatching"));
	// getHistFromMap("efficiencyVsCosThetaSum")->Add(getHistFromMap("efficiencyVsCosThetaFailAngularMatching"));

	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
		getHistFromMap("nPFOsVsCosThetaFailType_"+it->second)->Sumw2();
	}

	getHistFromMap("efficiencyVsCosThetaFailType_all")->Sumw2();
	getHistFromMap("efficiencyVsCosThetaFailType_all")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	getHistFromMap("efficiencyVsCosThetaFailType_noChargedParts")->Sumw2();
	getHistFromMap("efficiencyVsCosThetaFailType_noChargedParts")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	getHistFromMap("efficiencyVsCosThetaFailType_onlyPion")->Sumw2();
	getHistFromMap("efficiencyVsCosThetaFailType_onlyPion")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	getHistFromMap("efficiencyVsCosThetaFailType_onlyElectron")->Sumw2();
	getHistFromMap("efficiencyVsCosThetaFailType_onlyElectron")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	getHistFromMap("efficiencyVsCosThetaFailType_onlyMuon")->Sumw2();
	getHistFromMap("efficiencyVsCosThetaFailType_onlyMuon")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	getHistFromMap("efficiencyVsCosThetaFailEnergyMatching")->Sumw2();
	getHistFromMap("efficiencyVsCosThetaFailEnergyMatching")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	getHistFromMap("efficiencyVsCosThetaFailAngularMatching")->Sumw2();
	getHistFromMap("efficiencyVsCosThetaFailAngularMatching")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
	
	getHistFromMap("efficiencyVsCosThetaSum")->Add(getHistFromMap("efficiencyVsCosTheta"));
	getHistFromMap("efficiencyVsCosThetaSum")->Add(getHistFromMap("efficiencyVsCosThetaFailType_all"));
	getHistFromMap("efficiencyVsCosThetaSum")->Add(getHistFromMap("efficiencyVsCosThetaFailEnergyMatching"));
	getHistFromMap("efficiencyVsCosThetaSum")->Add(getHistFromMap("efficiencyVsCosThetaFailAngularMatching"));

	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
		getHistFromMap("nPFOsVsCosThetaFailType_"+it->second)->Sumw2();
		getHistFromMap("nPFOsVsCosThetaFailType_"+it->second)->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
		getHistFromMap("nPFOsVsCosThetaFailType_"+it->second)->Divide(getHistFromMap("efficiencyVsCosThetaFailType_all"));
	}


	for (int iCount=1; iCount<=4; iCount++){
		getHistFromMap("efficiencyVsCosThetaCat"+DoubToStr(iCount))->Sumw2();
		getHistFromMap("efficiencyVsCosThetaCat"+DoubToStr(iCount))->Divide(getHistFromMap("nTruthPartsVsCosTheta"));
		getHistFromMap("efficiencyVsCosThetaCat0")->Add(getHistFromMap("efficiencyVsCosThetaCat"+DoubToStr(iCount)));
	}
	getHistFromMap("efficiencyVsCosThetaCat5")->Sumw2();
	getHistFromMap("efficiencyVsCosThetaCat5")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));

	if (!outFile->IsOpen()){
		cout << "[ERROR|writeToFile]\tno output file is found!" << endl;
		return -1;
	}

	objectFill::writeToFile(outFile);
	return 0;

}

