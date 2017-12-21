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
// void eventHistFiller::initHistStructs(){
//         // all histogram NEED to have unique names!!!
//         const unsigned int nBinsPerGeV = 10;
//         singleParticleHistStructMap["Energy"] = histStruct("Particle Energy; Energy [GeV]; Counts",250*nBinsPerGeV,0.0,250.0);
//         singleParticleHistStructMap["phiResolution"] = histStruct("Phi resolution; dPhi [rad]; Counts",20000,-0.2,0.2);
//         singleParticleHistStructMap["thetaResolution"] = histStruct("Theta resolution; dTheta [rad]; Counts",400,-0.01,0.01);
// }

int eventHistFiller::init(){

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

	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
		tmpHist = new TH1D(("nPFOsVsTheta_"+it->second).c_str(),("n"+it->second+"s vs Theta; Theta; Counts per Event").c_str(),180*2,0,180);
		histMap["nPFOsVsTheta_"+it->second] = tmpHist;
		tmpHist = new TH1D(("nPFOsVsCosTheta_"+it->second).c_str(),("n"+it->second+"s vs Cos(#Theta); Cos(#Theta); Counts per Event").c_str(),180*2,-1,1);
		histMap["nPFOsVsCosTheta_"+it->second] = tmpHist;
		tmpHist = new TH1D(("phiResolution_"+it->second).c_str(),(it->second+" Phi resolution; dPhi [rad]; Counts").c_str(),20000,-0.2,0.2);
		histMap["phiResolution_"+it->second] = tmpHist;
		tmpHist = new TH1D(("thetaResolution_"+it->second).c_str(),(it->second+" Theta resolution; dTheta [rad]; Counts").c_str(),400,-0.01,0.01);
		histMap["thetaResolution_"+it->second] = tmpHist;
		tmpHist = new TH1D(("energyResolution_"+it->second).c_str(),(it->second+" Energy resolution; E [GeV]; Counts").c_str(),1250,0,125);
		histMap["energyResolution_"+it->second] = tmpHist;
	}
	
	tmpHist = new TH1D("nTruthPartsVsTheta","nTruthParts vs Theta; Theta; Counts per Event",180*2,0,180);
	histMap["nTruthPartsVsTheta"] = tmpHist;
	tmpHist = new TH1D("nTruthPartsVsCosTheta","nTruthParts vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["nTruthPartsVsCosTheta"] = tmpHist;
	tmpHist = new TH1D("nTruthPartsVsEnergy","nTruthParts vs Energy ;Energy [GeV]; Counts per Event",100,0.5,100.5); 
	histMap["nTruthPartsVsEnergy"] = tmpHist;

	tmpHist = new TH1I("truthParticle_isDecayedInTracker","Sim. Status; Sim. Status; Counts",35,0,35);
	histMap["truthParticle_isDecayedInTracker"] = tmpHist;

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

	histMap["truthParticle_isDecayedInTracker"]->Fill(genPart->isDecayedInTracker());
	nSelectecTruthParticles++;
	histMap["nTruthPartsVsCosTheta"]->Fill(cosPartTheta);
	histMap["nTruthPartsVsTheta"]->Fill(partTheta);
	truthEnergy = genPart->getEnergy();
	histMap["nTruthPartsVsEnergy"]->Fill(truthEnergy);

	vector<EVENT::ReconstructedParticle*> recoPFOs = getObjVecFromCollection(PFOCollection);
	histMap["nPFOs"]->Fill(recoPFOs.size());
	if (recoPFOs.size()==0) 
		return 0; // no reco PFOs

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

		histMap["PFOType"]->Fill(pfoType);
		histMap["nPFOsVsCosTheta_all"]->Fill(cosTruthTheta);
		histMap["nPFOsVsTheta_all"]->Fill(truthTheta);
		histMap["totalEnergyVsTheta"]->Fill(truthTheta,recoPFOs[i]->getEnergy());

		for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
			if (abs(pfoType)==it->first) {
				histMap["nPFOsVsCosTheta_"+it->second]->Fill(cosTruthTheta);
				histMap["nPFOsVsTheta_"+it->second]->Fill(truthTheta);
				histMap["thetaResolution_"+it->second]->Fill(dTheta);
				histMap["phiResolution_"+it->second]->Fill(dPhi);
				histMap["energyResolution_"+it->second]->Fill(pfoE);
			}
		}

	}
}


int eventHistFiller::writeToFile(TFile* outFile){

	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
		histMap["nPFOsVsCosTheta_"+it->second]->Divide(histMap["nTruthPartsVsCosTheta"]);
		histMap["nPFOsVsTheta_"+it->second]->Divide(histMap["nTruthPartsVsTheta"]);
	}
	histMap["nPFOsVsTheta_all"]->Divide(histMap["nTruthPartsVsTheta"]);
	histMap["nPFOsVsCosTheta_all"]->Divide(histMap["nTruthPartsVsCosTheta"]);
	
	histMap["totalEnergyVsTheta"]->Divide(histMap["nTruthPartsVsTheta"]);

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
	// return 0;

}

/*===========================================================================*/
/*===============================[ function implementations ]===============================*/
/*===========================================================================*/


