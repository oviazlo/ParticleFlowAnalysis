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

	cout << "[INFO]\teventHistFiller::init()" << endl;

	TH1* tmpHist;
	tmpHist = new TH1I("nPFOs","Number of PFOs in event; Number of PFOs; Counts",5,0,5);
	histMap["nPFOs"] = tmpHist;
	tmpHist= new TH1I("PFOType","PFO particle type; Type; Counts",2200,0,2200); // max part.type = 2112 (neutron)
	histMap["PFOType"] = tmpHist;
	
	tmpHist = new TH1D("nPFOsVsTheta","nPFOs vs Theta; Theta; Counts per Event",180*2,0,180);
	histMap["nPFOsVsTheta"] = tmpHist;
	tmpHist = new TH1D("nPFOsVsCosTheta","nPFOs vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["nPFOsVsCosTheta"] = tmpHist;
	
	tmpHist = new TH1D("nElectronsVsTheta","nElectrons vs Theta; Theta; Counts per Event",180*2,0,180);
	histMap["nElectronsVsTheta"] = tmpHist;
	tmpHist = new TH1D("nElectronsVsCosTheta","nElectrons vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["nElectronsVsCosTheta"] = tmpHist;
	
	tmpHist = new TH1D("nMuonsVsTheta","nMuons vs Theta; Theta; Counts per Event",180*2,0,180);
	histMap["nMuonsVsTheta"] = tmpHist;
	tmpHist = new TH1D("nMuonsVsCosTheta","nMuons vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["nMuonsVsCosTheta"] = tmpHist;

	tmpHist = new TH1D("nPhotonsVsTheta","nPhotons vs Theta; Theta; Counts per Event",180*2,0,180);
	histMap["nPhotonsVsTheta"] = tmpHist;
	tmpHist = new TH1D("nPhotonsVsCosTheta","nPhotons vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["nPhotonsVsCosTheta"] = tmpHist;

	tmpHist = new TH1D("nNeutralHadronsVsTheta","nNeutrals vs Theta; Theta; Counts per Event",180*2,0,180);
	histMap["nNeutralHadronsVsTheta"] = tmpHist;
	tmpHist = new TH1D("nNeutralHadronsVsCosTheta","nNeutrals vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["nNeutralHadronsVsCosTheta"] = tmpHist;

	tmpHist = new TH1D("nPionsVsTheta","nPions vs Theta; Theta; Counts per Event",180*2,0,180);
	histMap["nPionsVsTheta"] = tmpHist;
	tmpHist = new TH1D("nPionsVsCosTheta","nPions vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["nPionsVsCosTheta"] = tmpHist;
	
	tmpHist = new TH1D("nTruthPartsVsTheta","nTruthParts vs Theta; Theta; Counts per Event",180*2,0,180);
	histMap["nTruthPartsVsTheta"] = tmpHist;
	tmpHist = new TH1D("nTruthPartsVsCosTheta","nTruthParts vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["nTruthPartsVsCosTheta"] = tmpHist;

	tmpHist = new TH1D("efficiencyVsTheta","efficiency vs Theta; Theta; Counts per Event",180*2,0,180);
	histMap["efficiencyVsTheta"] = tmpHist;
	tmpHist = new TH1D("efficiencyVsCosTheta","efficiency vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["efficiencyVsCosTheta"] = tmpHist;

	tmpHist = new TH1D("phiResolution_photon","Photon Phi resolution; dPhi [rad]; Counts",400,-0.2,0.2);
	histMap["phiResolution_photon"] = tmpHist;
	tmpHist = new TH1D("thetaResolution_photon","Photon Theta resolution; dTheta [rad]; Counts",400,-0.01,0.01);
	histMap["thetaResolution_photon"] = tmpHist;

	tmpHist = new TH1D("phiResolution_neutral","Photon Phi resolution; dPhi [rad]; Counts",400,-0.2,0.2);
	histMap["phiResolution_neutral"] = tmpHist;
	tmpHist = new TH1D("thetaResolution_neutral","Neutral Theta resolution; dTheta [rad]; Counts",400,-0.01,0.01);
	histMap["thetaResolution_neutral"] = tmpHist;
	tmpHist = new TH1I("truthParticle_isDecayedInTracker","Sim. Status; Sim. Status; Counts",35,0,35);
	histMap["truthParticle_isDecayedInTracker"] = tmpHist;

	nSelectecTruthParticles = 0;
	return 0;

}

int eventHistFiller::fillEvent(const EVENT::LCEvent* event){
	// read collections
	try {
		PFOCollection = event->getCollection(PFOCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|eventHistFiller]\tCan't find collection: " << PFOCollectionName << endl;
		return -1;
	}
	try {
		MCTruthCollection = event->getCollection(MCTruthCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|eventHistFiller]\tCan't find collection: " << MCTruthCollectionName << endl;
		return -1;
	}

	// find primary generated MC particle which is the only findable in the event by definition [particle gun]
	

	double truthTheta = -666;  
	double truthPhi = -666;  
	double cosTruthTheta = -666; 
	double truthEnergy = -666;
	bool reclusteringIsDone = false;

	EVENT::MCParticle* genPart = NULL;
	int nElements = MCTruthCollection->getNumberOfElements();
	for(int j=0; j < nElements; j++) {
		auto part = dynamic_cast<EVENT::MCParticle*>(MCTruthCollection->getElementAt(j));
		if (part->getGeneratorStatus()==1){
			const double *partMom = part->getMomentum();
			TVector3 v1;
			v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
			double partTheta = 180.*v1.Theta()/TMath::Pi();
			double cosPartTheta = TMath::Cos(TMath::Pi()*partTheta/180.);
			truthTheta = partTheta;
			truthPhi = 180.*v1.Phi()/TMath::Pi();
			cosTruthTheta = cosPartTheta;
			if (partTheta<8 || partTheta>172) 
				return 0;
			if (part->isDecayedInTracker() && flagDiscardConvertion==true)
				return 0;
			if ((!part->isDecayedInTracker()) && flagSelectConvertion==true)
				return 0;
			histMap["truthParticle_isDecayedInTracker"]->Fill(part->isDecayedInTracker());
			genPart = part;
			nSelectecTruthParticles++;
			histMap["nTruthPartsVsCosTheta"]->Fill(cosPartTheta);
			histMap["nTruthPartsVsTheta"]->Fill(partTheta);
			truthEnergy = part->getEnergy();
			break;
		}
	}

	vector<EVENT::ReconstructedParticle*> recoPFOs = getObjVecFromCollection(PFOCollection);
	histMap["nPFOs"]->Fill(recoPFOs.size());
	if (recoPFOs.size()==0) 
		return 0; // no reco PFOs
	/// 11,-11,13,-13,211,-211,22,2112

	for (int i=0; i<recoPFOs.size(); i++){
		int pfoType = recoPFOs[i]->getType();
		const double *partMom = recoPFOs[i]->getMomentum();
		TVector3 v1;
		v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
		double partTheta = 180.*v1.Theta()/TMath::Pi();
		double partPhi = 180.*v1.Phi()/TMath::Pi();
		double cosPartTheta = TMath::Cos(TMath::Pi()*partTheta/180.);

		// cout << "partTheta: " << partTheta << "; truthTheta: " << truthTheta << endl;

		// WARNING TODO FIXME: wrong when calculate dTheta!!!
		partTheta = truthTheta;
		cosPartTheta = cosTruthTheta;

		histMap["PFOType"]->Fill(pfoType);
		histMap["nPFOsVsCosTheta"]->Fill(cosPartTheta);
		histMap["nPFOsVsTheta"]->Fill(partTheta);


		for (auto it = pfoTypeIntStringMap.begin(); it != pfoTypeIntStringMap.end(); it++) {
			if (abs(pfoType)==it->first) {
				histMap["n"+it->second+"sVsCosTheta"]->Fill(cosPartTheta);
				histMap["n"+it->second+"sVsTheta"]->Fill(partTheta);
			}
		}
		switch(abs(pfoType)){
			case 11 : histMap["nElectronsVsCosTheta"]->Fill(cosPartTheta);
				  histMap["nElectronsVsTheta"]->Fill(partTheta);
				  break;
			case 13 : histMap["nMuonsVsCosTheta"]->Fill(cosPartTheta);
				  histMap["nMuonsVsTheta"]->Fill(partTheta);
				  break;
			case 211: histMap["nPionsVsCosTheta"]->Fill(cosPartTheta);
				  histMap["nPionsVsTheta"]->Fill(partTheta);
				  break;
			case 22 : histMap["nPhotonsVsCosTheta"]->Fill(cosPartTheta);
				  histMap["nPhotonsVsTheta"]->Fill(partTheta);
				  // histMap["thetaResolution_photon"]->Fill((partTheta-truthTheta)*TMath::Pi()/180.0);
				  histMap["thetaResolution_photon"]->Fill((partTheta-truthTheta));
				  histMap["phiResolution_photon"]->Fill((partPhi-truthPhi)*TMath::Pi()/180.0);
				  break;
			case 2112:histMap["nNeutralHadronsVsCosTheta"]->Fill(cosPartTheta);
				  histMap["nNeutralHadronsVsTheta"]->Fill(partTheta);
				  // histMap["thetaResolution_neutral"]->Fill((partTheta-truthTheta)*TMath::Pi()/180.0);
				  histMap["thetaResolution_neutral"]->Fill((partTheta-truthTheta));
				  histMap["phiResolution_neutral"]->Fill((partPhi-truthPhi)*TMath::Pi()/180.0);
				  break;
		}
		
		if (abs(pfoType)==22 && reclusteringIsDone==false) {
			double pfoE = recoPFOs[i]->getEnergy();
			reclusteringIsDone = true;
			if (photonReclustering==true || neutralReclustering==true) {
				for (int j=0; j<recoPFOs.size(); j++){
					if (i==j) 
						continue;
					if (!((recoPFOs[j]->getType()==22) || (recoPFOs[j]->getType()==2112)))
							continue;
					if (recoPFOs[j]->getType()==22 && photonReclustering==false)
						continue;
					if (recoPFOs[j]->getType()==2112 && neutralReclustering==false)
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
					if (abs(partTheta_i-partTheta_j)<0.6 && abs(partPhi_i-partPhi_j)<10.0) {
						pfoE+=recoPFOs[j]->getEnergy();
					}
					// if (recoPFOs[j]->getType()==22){
					//         TODO
					//         ???
					// }
				}
			}
			if (abs(pfoE-truthEnergy)<0.75*sqrt(truthEnergy)) {
				histMap["efficiencyVsTheta"]->Fill(partTheta);
				histMap["efficiencyVsCosTheta"]->Fill(cosPartTheta);
			}
		}
		
	}
}


int eventHistFiller::writeToFile(TFile* outFile){

	for (auto it = pfoTypeIntStringMap.begin(); it != pfoTypeIntStringMap.end(); it++) {
		if (abs(pfoType)==it->first) {
			histMap["n"+it->second+"sVsCosTheta"]->Divide(histMap["nTruthPartsVsCosTheta"]);
			histMap["n"+it->second+"sVsTheta"]->Divide(histMap["nTruthPartsVsTheta"]);
		}
	}
	histMap["nPFOsVsCosTheta"]->Divide(histMap["nTruthPartsVsCosTheta"]);
	histMap["nPFOsVsTheta"]->Divide(histMap["nTruthPartsVsTheta"]);
	histMap["nElectronsVsCosTheta"]->Divide(histMap["nTruthPartsVsCosTheta"]);
	histMap["nElectronsVsTheta"]->Divide(histMap["nTruthPartsVsTheta"]);
	histMap["nMuonsVsCosTheta"]->Divide(histMap["nTruthPartsVsCosTheta"]);
	histMap["nMuonsVsTheta"]->Divide(histMap["nTruthPartsVsTheta"]);
	histMap["nPionsVsCosTheta"]->Divide(histMap["nTruthPartsVsCosTheta"]);
	histMap["nPionsVsTheta"]->Divide(histMap["nTruthPartsVsTheta"]);
	histMap["nPhotonsVsCosTheta"]->Divide(histMap["nTruthPartsVsCosTheta"]);
	histMap["nPhotonsVsTheta"]->Divide(histMap["nTruthPartsVsTheta"]);
	histMap["nNeutralHadronsVsCosTheta"]->Divide(histMap["nTruthPartsVsCosTheta"]);
	histMap["nNeutralHadronsVsTheta"]->Divide(histMap["nTruthPartsVsTheta"]);
	
	histMap["efficiencyVsCosTheta"]->Divide(histMap["nTruthPartsVsCosTheta"]);
	histMap["efficiencyVsTheta"]->Divide(histMap["nTruthPartsVsTheta"]);

	objectFill::writeToFile(outFile);
	return 0;

}
/*===========================================================================*/
/*===============================[ function implementations ]===============================*/
/*===========================================================================*/

