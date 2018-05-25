/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/src/jetPfoStudy.cpp
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	05th Dec 2017
 * 	Last Update:	05th Dec 2017
 */

#include "jetPfoStudy.h"

/*===========================================================================*/
/*===============================[ function implementations ]================*/
/*===========================================================================*/

int jetPfoStudy::init(){
	debugFlag = false;

	if (config::vm.count("debug"))
		cout << "[INFO]\tjetPfoStudy::init()" << endl;

	TH1* tmpHist;
	tmpHist = new TH1I("nPFOs","Number of PFOs in event; Number of PFOs; Counts",100,0,100);
	histMap["nPFOs"] = tmpHist;
	tmpHist= new TH1I("PFOType","PFO particle type; Type; Counts",2200,0,2200); // max part.type = 2112 (neutron)
	histMap["PFOType"] = tmpHist;
	
	tmpHist = new TH1D("nPFOsVsTheta_all","nPFOs vs Theta; Theta; Counts per Event",180,0,180);
	histMap["nPFOsVsTheta_all"] = tmpHist;
	tmpHist = new TH1D("nPFOsVsCosTheta_all","nPFOs vs cos(#theta); cos(#theta); Counts per Event",180,-1,1);
	histMap["nPFOsVsCosTheta_all"] = tmpHist;

	tmpHist = new TH1D("totalEnergyVsTheta","Sum of PFOs Energy vs Theta; Theta; Energy [GeV]",180,0,180);
	histMap["totalEnergyVsTheta"] = tmpHist;
	tmpHist = new TH1D("totalEnergyVsCosTheta","Sum of PFOs Energy vs CosTheta; Theta; Energy [GeV]",180,-1,1);
	histMap["totalEnergyVsCosTheta"] = tmpHist;


	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
		tmpHist = new TH1D(("nPFOsVsTheta_"+it->second).c_str(),("n"+it->second+"s vs Theta; Theta; Counts per Event").c_str(),180,0,180);
		histMap["nPFOsVsTheta_"+it->second] = tmpHist;
		tmpHist = new TH1D(("nPFOsVsCosTheta_"+it->second).c_str(),("n"+it->second+"s vs cos(#theta); cos(#theta); Counts per Event").c_str(),180,-1,1);
		histMap["nPFOsVsCosTheta_"+it->second] = tmpHist;

	}
	
	tmpHist = new TH1D("nTruthPartsVsTheta","nTruthParts vs Theta; Theta; Counts per Event",180,0,180);
	histMap["nTruthPartsVsTheta"] = tmpHist;
	tmpHist = new TH1D("nTruthPartsVsCosTheta","nTruthParts vs cos(#theta); cos(#theta); Counts per Event",180,-1,1);
	histMap["nTruthPartsVsCosTheta"] = tmpHist;

	tmpHist = new TH1D("totalRecoEnergy","Total Reconstructed energy; E [GeV]; Counts",1250,0,125);
	histMap["totalRecoEnergy"] = tmpHist;

	return 0;
}

int jetPfoStudy::fillEvent(const EVENT::LCEvent* event){
	mcQuarkVector.clear();
	jetTheta = std::numeric_limits<double>::max();
	jetCosTheta = std::numeric_limits<double>::max();

	if (debugFlag)
		cout << "before reading collections\n";

	try {
		if (debugFlag)
			cout << "Try to read <" << MCCollectionName << "> collection...\n";
		EVENT::LCCollection* mcCollection = event->getCollection(MCCollectionName);
		if (debugFlag){
			cout << "Collection pointer: " << mcCollection << endl;
			cout << "Try to access it...\n";
			cout << "nElements: " << mcCollection->getNumberOfElements() << endl;
		}
		for (unsigned int i = 0, nElements = mcCollection->getNumberOfElements(); i < nElements; ++i)
		{
			if (debugFlag)
				cout << "Try to access element: " << i << endl;
			const EVENT::MCParticle *pMCParticle = dynamic_cast<EVENT::MCParticle*>(mcCollection->getElementAt(i));

			if (NULL == pMCParticle)
				throw EVENT::Exception("Collection type mismatch");

			const int absPdgCode(std::abs(pMCParticle->getPDG()));

			if ((absPdgCode >= 1) && (absPdgCode <= 6) && pMCParticle->getParents().empty())
				mcQuarkVector.push_back(pMCParticle);
		}
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|jetPfoStudy]\tCan't find collection: " << MCCollectionName << endl;
		return -1;
	}

	try {
		if (debugFlag)
			cout << "Try to read <" << PFOCollectionName << "> collection\n";
		PFOCollection = event->getCollection(PFOCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|jetPfoStudy]\tCan't find collection: " << PFOCollectionName << endl;
		return -1;
	}

	if (!mcQuarkVector.empty()){
		int m_qPdg = std::abs(mcQuarkVector[0]->getPDG());
		float energyTot(0.f);
		float costTot(0.f);

		for (unsigned int i = 0; i < mcQuarkVector.size(); ++i)
		{
		    const float px(mcQuarkVector[i]->getMomentum()[0]);
		    const float py(mcQuarkVector[i]->getMomentum()[1]);
		    const float pz(mcQuarkVector[i]->getMomentum()[2]);
		    const float energy(mcQuarkVector[i]->getEnergy());
		    const float p(std::sqrt(px * px + py * py + pz * pz));
		    const float cost(std::fabs(pz) / p);
		    energyTot += energy;
		    costTot += cost * energy;
		}

		jetCosTheta = costTot / energyTot;
		jetTheta = acos(jetCosTheta);
	}

	if (config::vm.count("debug"))
		cout << "[INFO]	jetPfoStudy::fillEvent: " << event->getEventNumber() << endl;

	const double truthTheta = jetTheta*TMath::RadToDeg();
	const double cosTruthTheta = jetCosTheta;

	getHistFromMap("nTruthPartsVsCosTheta")->Fill(cosTruthTheta);
	getHistFromMap("nTruthPartsVsTheta")->Fill(truthTheta);

	vector<EVENT::ReconstructedParticle*> recoPFOs = getObjVecFromCollection<EVENT::ReconstructedParticle*>(PFOCollection);
	getHistFromMap("nPFOs")->Fill(recoPFOs.size());

	map <string, unsigned int> pfoCounter;
	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++)
		pfoCounter[it->second] = 0;

	double totalRecoEnergy = 0.0;

	for (int i=0; i<recoPFOs.size(); i++){

		const int pfoType = abs(recoPFOs[i]->getType());
		const double *partMom = recoPFOs[i]->getMomentum();
		TVector3 vPartMom(partMom[0],partMom[1],partMom[2]);
		const double partTheta = vPartMom.Theta()*TMath::RadToDeg();
		const double partPhi = vPartMom.Phi()*TMath::RadToDeg();
		const double partPt = vPartMom.Pt();
		const double cosPartTheta = TMath::Cos(partTheta*TMath::DegToRad());
		const double partEnergy = recoPFOs[i]->getEnergy();

		getHistFromMap("PFOType")->Fill(pfoType);
		getHistFromMap("nPFOsVsCosTheta_all")->Fill(cosTruthTheta);
		getHistFromMap("nPFOsVsTheta_all")->Fill(truthTheta);
		getHistFromMap("totalEnergyVsTheta")->Fill(truthTheta,recoPFOs[i]->getEnergy());
		getHistFromMap("totalEnergyVsCosTheta")->Fill(cosTruthTheta,recoPFOs[i]->getEnergy());
		totalRecoEnergy += recoPFOs[i]->getEnergy();

		for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
			if (abs(pfoType)==it->first) {
				getHistFromMap("nPFOsVsCosTheta_"+it->second)->Fill(cosTruthTheta);
				getHistFromMap("nPFOsVsTheta_"+it->second)->Fill(truthTheta);
				pfoCounter[it->second]++;
			}
		}
	}

	getHistFromMap("totalRecoEnergy")->Fill(totalRecoEnergy);

}


int jetPfoStudy::writeToFile(TFile* outFile){

	if (config::vm.count("debug"))
		cout << "[INFO]	jetPfoStudy::writeToFile(" << outFile->GetName() << ")" << endl;

	getHistFromMap("nTruthPartsVsTheta")->Sumw2();
	getHistFromMap("nTruthPartsVsCosTheta")->Sumw2();

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
	getHistFromMap("totalEnergyVsTheta")->Sumw2();
	getHistFromMap("totalEnergyVsTheta")->Divide(getHistFromMap("nTruthPartsVsTheta"));
	getHistFromMap("totalEnergyVsCosTheta")->Sumw2();
	getHistFromMap("totalEnergyVsCosTheta")->Divide(getHistFromMap("nTruthPartsVsCosTheta"));


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
}

