#include <particleFill.h>

#include <IMPL/ReconstructedParticleImpl.h>

#define nBinsPerGeV 1
#define nExampleClusterHists 10
#define maxEnergyRatioExampleClusterHist 2.0

void particleFill::initHistStructs(){
	// all histogram NEED to have unique names!!!
	singleParticleHistStructMap["Energy"] = histStruct("Particle Energy; Energy [GeV]; Counts",250*nBinsPerGeV,0.0,250.0);
	singleParticleHistStructMap["Phi"] = histStruct("Particle Phi; Phi; Counts",360,-180,180);
	singleParticleHistStructMap["Theta"] = histStruct("Particle Theta; Theta; Counts",180*2,0,180);
	singleParticleHistStructMap["CosTheta"] = histStruct("Particle Cos(#Theta); Cos(#Theta); Counts",180*2,-1,1);
	singleParticleHistStructMap["Type"] = histStruct("Particle Type; Type; Counts",2200,0,2200,"TH1I");

	twoParticleCorrelationHistStructMap["dPhi"] = histStruct("dPhi; Phi; Counts",360,-180,180);
	twoParticleCorrelationHistStructMap["dR"] = histStruct("dR; #Delta R; Counts",200,0,2);
	twoParticleCorrelationHistStructMap["dE"] = histStruct("E_{1}-E_{2}; E_{1}-E_{2}; Counts",100*nBinsPerGeV,0.0,100.0);
	twoParticleCorrelationHistStructMap["sumE"] = histStruct("E_{1}+E_{2}; E_{1}+E_{2} [GeV] ; Counts",250*nBinsPerGeV,0.0,250.0);
	twoParticleCorrelationHistStructMap["TwoParticleType"] = histStruct("Particle Type; Type; Counts",2200,0,2200,"TH1I");

	singleRecoParticleClustersHistStructMap["ECALtoHCALRatio"] = histStruct("ECAL/HCAL Energy Ratio; ECAL/HCAL ratio; Counts",250,0,100);
	singleRecoParticleClustersHistStructMap["ECALEnergy"] = histStruct("ECAL Energy; ECAL Energy [GeV]; Counts",125*nBinsPerGeV,0.0,125.0);
	singleRecoParticleClustersHistStructMap["HCALEnergy"] = histStruct("HCAL Energy; HCAL Energy [GeV]; Counts",125*nBinsPerGeV,0.0,125.0);
	singleRecoParticleClustersHistStructMap["TestCorrectedEnergy"] = histStruct("CAL Energy; CAL Energy [GeV]; Counts",125*nBinsPerGeV,0.0,125.0);
	singleRecoParticleClustersHistStructMap["clusterPositionZR"] = histStruct("Cluster Position; Z [mm]; R [mm]",800,-4000,4000,"TH2D",400,0,4000);
	singleRecoParticleClustersHistStructMap["clusterPositionMoreECALEnergy"] = histStruct("Cluster Position; Z [mm]; R [mm]",800,-4000,4000,"TH2D",400,0,4000);
	singleRecoParticleClustersHistStructMap["clusterPositionMoreHCALEnergy"] = histStruct("Cluster Position; Z [mm]; R [mm]",800,-4000,4000,"TH2D",400,0,4000);
	for (unsigned int i; i<nExampleClusterHists; i++)
		singleRecoParticleClustersHistStructMap["hitPositionsOneEventExample"+DoubToStr(i)] = histStruct("Calo hits position; Z [mm]; R [mm]",800,-4000,4000,"TH2D",400,0,4000);
	

}

int particleFill::init(){

	initHistStructs();
	for (int i=0; i<nExampleClusterHists; i++)
		boolVecDefaultFalse.push_back(false);

	if (IsInWord("MCParticle",collectionName)){ 
		createSingleParticleHists("truthParticle_");
	}
	else if(IsInWord("SiTracks",collectionName)){
		// add track-related hists. TODO rewrite it with functions
		TH1I* tmpHist = new TH1I("tracking_nTracks","Number of tracks; # of tracks; Counts", 5, 0, 5);
		histMap["tracking_nTracks"] = tmpHist;
		TH1D* tmpHist2 = new TH1D("tracking_allTracksMomentum","Track momentum; Momentum [GeV]; Counts",250*nBinsPerGeV,0.0,250);
		histMap["tracking_allTracksMomentum"] = tmpHist2;
		tmpHist2 = new TH1D("tracking_1st_trackMomentum","Track momentum; Momentum [GeV]; Counts",250*nBinsPerGeV,0.0,250);
		histMap["tracking_1st_trackMomentum"] = tmpHist2;
		tmpHist2 = new TH1D("tracking_2nd_trackMomentum","Track momentum; Momentum [GeV]; Counts",250*nBinsPerGeV,0.0,250);
		histMap["tracking_2nd_trackMomentum"] = tmpHist2;
		tmpHist2 = new TH1D("tracking_3rd_trackMomentum","Track momentum; Momentum [GeV]; Counts",250*nBinsPerGeV,0.0,250);
		histMap["tracking_3rd_trackMomentum"] = tmpHist2;
		tmpHist2 = new TH1D("tracking_1st_trackPt","Track Pt; Pt [GeV]; Counts",250*nBinsPerGeV,0.0,250);
		histMap["tracking_1st_trackPt"] = tmpHist2;
	}
	else{ // Reconstructed particles
		createSingleParticleHists("allPartsOfType_");
		createSingleParticleHists("signalPFO_");
		createSingleParticleHists("signalPFO_truthParts_");
		createSingleParticleHists("signalPFO_noAdditionalPFOs_");
		createSingleParticleHists("signalPFO_noAdditionalPFOs_truthParts_");
		createSingleParticleHists("signalPFO_thereAreAdditionalPFOs_");
		createSingleParticleHists("firstEnergeticPartOfWrongType_");
		createSingleParticleHists("secondEnergeticPart_");

		createSingleParticleHists("PFOofType_noAdditionalPFOs_");
		createSingleParticleHists("PFOofType_thereAreAdditionalPFOs_");
		createSingleParticleHists("mostEnergeticPFO_noPFOofType_");
		// truth hists:
		createSingleParticleHists("truthParticle_onlyOneRecoPFO_");
		createSingleParticleHists("truthParticle_twoOrMoreRecoPFO_");

		createTwoParticleCorrelationHists("signalPFO_secondEnergeticPFO_");
		createTwoParticleCorrelationHists("signalPFO_secondEnergeticPhoton_");

		createSingleRecoParticleClustersHists("signalPFO_clusterInfo_");
		createSingleRecoParticleClustersHists("signalPFO_noAdditionalPFOs_clusterInfo_");

		TH1I* h_setPartType = new TH1I("setPFOPartType","Set PFO particle type to use; Type; Counts",2200,0,2200); // max part.type = 2112 (neutron)
		histMap["setPFOPartType"] = h_setPartType;
		for (auto i=0; i<partTypeToSelect.size(); i++)
			histMap["setPFOPartType"]->Fill(partTypeToSelect[i]);

		if (partTypeToSelect.size()==0)
			cout << "[ERROR|particleFill]\tNo PFO type particle are set up for collection: " << collectionName << "" << endl;
	}



	for(auto const &iMapElement : histMap) {
		iMapElement.second->AddDirectory(kFALSE);
	}

	dPhiMergeValue = 0;
	debugFlag = false;
	collection = NULL;

	return 0;
}



int particleFill::fillParticleCorrelations (const EVENT::ReconstructedParticle* inPart1, const EVENT::ReconstructedParticle* inPart2, const string prefix){
	// TODO implement dE, dPhi between properly reconstructed and second energetic PFOs
	const double *partMom1 = inPart1->getMomentum();
	TVector3 v1, v2;
	v1.SetXYZ(partMom1[0],partMom1[1],partMom1[2]);
	const double *partMom2 = inPart2->getMomentum();
	v2.SetXYZ(partMom2[0],partMom2[1],partMom2[2]);
	double dPhi = v1.DeltaPhi(v2)*180./TMath::Pi();
	double dR = v1.DeltaR(v2);

	double dE = inPart1->getEnergy() - inPart2->getEnergy();
	double sumE = inPart1->getEnergy() + inPart2->getEnergy();

	fillHist(dPhi, "dPhi", prefix);
	fillHist(dR, "dR", prefix);
	fillHist(dE, "dE", prefix);
	fillHist(sumE, "sumE", prefix);

	fillHist(inPart1->getType(), "TwoParticleType", prefix);
	fillHist(inPart2->getType(), "TwoParticleType", prefix);

	return 0;
}

vector<EVENT::MCParticle*> particleFill::getTruthMCParticlesFromCollection(const EVENT::LCEvent* event, const string truthCollectionName){
	// FIXME hardcoding of the collection name
	vector<EVENT::MCParticle*> outPartVec;
	EVENT::LCCollection *truthCollection;
	try {
		truthCollection = event->getCollection(truthCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|particleFill]\tCan't find collection: " << truthCollectionName << endl;
		return outPartVec;
	}

	const int nElements = truthCollection->getNumberOfElements();
	// const int nElements = collection->getNumberOfElements();
	for(int j=0; j < nElements; j++) {
		auto part = dynamic_cast<EVENT::MCParticle*>(truthCollection->getElementAt(j));
		outPartVec.push_back(part);
	}
	return outPartVec;
}

int particleFill::fillEvent(const EVENT::LCEvent* event){
	try {
		collection = event->getCollection(collectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|particleFill]\tCan't find collection: " << collectionName << endl;
		return -1;
	}

	auto truthParts = getTruthMCParticlesFromCollection(event);
	// for (int iTr=0; iTr<truthParts.size(); iTr++){
	//         cout << iTr << ": " << truthParts[iTr] << endl;
	// }

	if( collection ) {
		string collectionType = collection->getTypeName();
		const int nElements = collection->getNumberOfElements();
		if (nElements<=0)
			return 0;
		// unsigned int countE_1 = 0;
		// unsigned int countE_2 = 0;
		// double E_1 = 0.0;
		// double E_2 = 0.0;
		
		double E_1 = 0.0;
		double E_2 = 0.0;

		intMap["id_firstEnergeticPFO"] = -1;
		intMap["id_secondEnergeticPFO"] = -1;
		intMap["id_mostEnergeticPFOofType"] = -1;
		intMap["n_PFOofType"] = 0;
		intMap["n_PFO"] = nElements;

		
      		for(int j=0; j < nElements; j++) {
			if (collectionType=="MCParticle"){
				auto part = dynamic_cast<EVENT::MCParticle*>(collection->getElementAt(j));
				if (debugFlag) dumpTruthPart(part, j);
				// dumpTruthPart(part, j);
				if (part->isCreatedInSimulation()!=0) continue;
				if (part->getGeneratorStatus()!=1) continue;
				fillPart(part);
			}
			else if (collectionType=="Track"){
				auto track = dynamic_cast<EVENT::Track*>(collection->getElementAt(j));
				double omtrack = track->getOmega();
				double tLtrack = track->getTanLambda();
				double recoThetaRad = 3.141592/2.0 - atan(tLtrack);
				double m_magneticField = 2.0;
				double recoPt = 0.3*m_magneticField/(fabs(omtrack)*1000.);
				double momentum = recoPt/sin(recoThetaRad);
				if (j==0){
					histMap["tracking_nTracks"]->Fill(nElements);
					histMap["tracking_1st_trackMomentum"]->Fill(momentum);
					histMap["tracking_1st_trackPt"]->Fill(recoPt);
				}
				if (j==1)
					histMap["tracking_2nd_trackMomentum"]->Fill(momentum);
				if (j==2)
					histMap["tracking_3rd_trackMomentum"]->Fill(momentum);
				histMap["tracking_allTracksMomentum"]->Fill(momentum);
			}
			else if (collectionType=="ReconstructedParticle"){
				auto part = dynamic_cast<EVENT::ReconstructedParticle*>(collection->getElementAt(j));
				if (debugFlag) dumpReconstructedPart(part, j);
				if(find(partTypeToSelect.begin(), partTypeToSelect.end(), part->getType()) != partTypeToSelect.end()){
					fillPart(part,"allPartsOfType_");
				}


				clasiffyPFO(part);

				// if (j==0){
				//         E_1 = part->getEnergy();
				// }
				// else{
				//         double E = part->getEnergy();
				//         if (E>E_1){
				//                 intMap["id_secondEnergeticPFO"] = intMap["id_firstEnergeticPFO"];
				//                 E_2 = E_1;
				//                 E_1 = E;
				//                 intMap["id_firstEnergeticPFO"] = j;
				//         }
				//         if (E>E_2 && E<E_1){
				//                 E_2 = E;
				//                 intMap["id_secondEnergeticPFO"] = j;
				//         }
				// }
				// cout << j << ":\tE1: " << E_1 << "; E2: " << E_2 << endl;
			}
			else{
				cout << "[ERROR|fillEvent]\tuknown particle collection: " << collectionType << endl;
			}
		}
		if (collectionType=="ReconstructedParticle"){
			auto part = dynamic_cast<IMPL::ReconstructedParticleImpl*>(collection->getElementAt(countE_1));
			
			// FIXME WARNING testing merging of clusters...
			if (dPhiMergeValue > 0.0){
				auto tmpPart = dynamic_cast<IMPL::ReconstructedParticleImpl*>(collection->getElementAt(countE_1));
				part = new IMPL::ReconstructedParticleImpl();
				part->setEnergy(tmpPart->getEnergy());
				part->setMomentum(tmpPart->getMomentum());
				part->setType(tmpPart->getType());

				for(int j=0; j < nElements; j++) {
					if (j==countE_1) continue;
					TVector3 v1, v2;
					const double *partMom1 = part->getMomentum();
					v1.SetXYZ(partMom1[0],partMom1[1],partMom1[2]);

					auto part2 = dynamic_cast<EVENT::ReconstructedParticle*>(collection->getElementAt(j));
					const double *partMom2 = part2->getMomentum();
					v2.SetXYZ(partMom2[0],partMom2[1],partMom2[2]);
					double dPhi = v1.DeltaPhi(v2)*180./TMath::Pi();
					if (abs(dPhi)<2.0)
						part->setEnergy(part->getEnergy()+part2->getEnergy());
				}
			}
			// FIXME
			if(find(partTypeToSelect.begin(), partTypeToSelect.end(), part->getType()) != partTypeToSelect.end()){

				fillPart(part,"signalPFO_");
				for (int iTr=0; iTr<truthParts.size(); iTr++)
					fillPart(truthParts[iTr],"signalPFO_truthParts_");
				fillClusterInfo(part,"signalPFO_clusterInfo_");
				if (nElements>=2 && countE_1!=countE_2){
					auto part2 = dynamic_cast<EVENT::ReconstructedParticle*>(collection->getElementAt(countE_2));
					fillPart(part2,"secondEnergeticPart_");
					fillParticleCorrelations(part,part2,"signalPFO_secondEnergeticPFO_");
					if (part2->getType()==22)
						fillParticleCorrelations(part,part2,"signalPFO_secondEnergeticPhoton_");
				}
				// if (nElements>=2 && countE_1==countE_2){
				//         cout << "countE: " << countE_2 << "; E1: " << E_1 << "; E2: " << E_2 << endl;
				// }


				if (nElements>=2){
					fillPart(part,"signalPFO_thereAreAdditionalPFOs_");
					fillPart(truthParts[0],"truthParticle_twoOrMoreRecoPFO_");
				}
				else{
					fillPart(part,"signalPFO_noAdditionalPFOs_");
					for (int iTr=0; iTr<truthParts.size(); iTr++)
						fillPart(truthParts[iTr],"signalPFO_noAdditionalPFOs_truthParts_");
					fillPart(truthParts[0],"truthParticle_onlyOneRecoPFO_");
					fillClusterInfo(part,"signalPFO_noAdditionalPFOs_clusterInfo_");
					
				}


			}
			else
				fillPart(part,"firstEnergeticPartOfWrongType_");
				

		}
	}
	return 0;
}

int particleFill::fillPart (const EVENT::MCParticle* inPart, const string prefix) {
	const double *partMom = inPart->getMomentum();
	TVector3 v1;
	v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
	double partPhi = 180.*v1.Phi()/TMath::Pi();
	double partTheta = 180.*v1.Theta()/TMath::Pi();
	double partEnergy = inPart->getEnergy();
	fillHist(partEnergy, "Energy", prefix);
	fillHist(partPhi, "Phi", prefix);
	fillHist(partTheta, "Theta", prefix);
	fillHist(TMath::Cos(TMath::Pi()*partTheta/180.), "CosTheta", prefix);
	fillHist(inPart->getPDG(), "Type", prefix);
	return 0;
}

int particleFill::fillPart (const EVENT::ReconstructedParticle* inPart, const string prefix){
	const double *partMom = inPart->getMomentum();
	TVector3 v1;
	v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
	double partPhi = 180.*v1.Phi()/TMath::Pi();
	double partTheta = 180.*v1.Theta()/TMath::Pi();
	double partEnergy = inPart->getEnergy();
	int partType= inPart->getType();
	fillHist(partEnergy, "Energy", prefix);
	fillHist(partPhi, "Phi", prefix);
	fillHist(partTheta, "Theta", prefix);
	fillHist(TMath::Cos(TMath::Pi()*partTheta/180.), "CosTheta", prefix);
	fillHist(partType, "Type", prefix);
	return 0;
}

int particleFill::fillClusterInfo (const EVENT::ReconstructedParticle* inPart, const string prefix){
	vector<EVENT::Cluster*> clusterVec = inPart->getClusters();
	for (int i=0; i<clusterVec.size(); i++){
		const float *pos = clusterVec[i]->getPosition();
		fillHist(pos[2],sqrt(pos[1]*pos[1]+pos[0]*pos[0]), "clusterPositionZR", prefix);
		// "ecal", "hcal", "yoke", "lcal", "lhcal", "bcal"
		vector<float> subdetEnergies = clusterVec[i]->getSubdetectorEnergies();
		double ratio = subdetEnergies[0]/subdetEnergies[1];
		fillHist(subdetEnergies[0], "ECALEnergy", prefix);
		fillHist(subdetEnergies[1], "HCALEnergy", prefix);
		double testCorrectedEnergy = 1.1191707392*(subdetEnergies[0]/1.02425625854) + 1.18038308419*(subdetEnergies[1]/1.02425625854);
		fillHist(testCorrectedEnergy, "TestCorrectedEnergy", prefix);
		fillHist(ratio, "ECALtoHCALRatio", prefix);
		if (ratio>1)
			fillHist(pos[2],sqrt(pos[1]*pos[1]+pos[0]*pos[0]), "clusterPositionMoreECALEnergy", prefix);
		else
			fillHist(pos[2],sqrt(pos[1]*pos[1]+pos[0]*pos[0]), "clusterPositionMoreHCALEnergy", prefix);

		// EXAMPLE CLUSTER HISTS
		for (unsigned int iExample = 0; iExample<nExampleClusterHists; iExample++){
			if (ratio > (maxEnergyRatioExampleClusterHist/nExampleClusterHists*iExample) && ratio < (maxEnergyRatioExampleClusterHist/nExampleClusterHists*(iExample+1)) && boolVecDefaultFalse[iExample]==false){
				boolVecDefaultFalse[iExample]=true;
				histMap[prefix+"hitPositionsOneEventExample"+DoubToStr(iExample)]->SetTitle(("E: " + DoubToStr(clusterVec[i]->getEnergy()) + " GeV, ECAL to HCAL ratio: " + DoubToStr(ratio)).c_str());
				vector<EVENT::CalorimeterHit*> caloHits = clusterVec[i]->getCalorimeterHits();
				for (int j=0; j<caloHits.size(); j++){
					const float *hitPos = caloHits[j]->getPosition();
					fillHist(hitPos[2],sqrt(hitPos[1]*hitPos[1]+hitPos[0]*hitPos[0]), "hitPositionsOneEventExample"+DoubToStr(iExample), prefix);
				}

			}
		}

	}

	return 0;
}

void particleFill::dumpTruthPart(const EVENT::MCParticle* part, const int counter){
	const double *partMom = part->getMomentum();
	TVector3 v1;
	v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
	double partPhi = 180.*v1.Phi()/TMath::Pi();
	double partTheta = 180.*v1.Theta()/TMath::Pi();
	bool inTracker = part->isDecayedInTracker();
	bool inCal = part->isDecayedInCalorimeter();
	int genStatus = part->getGeneratorStatus();
	int pdgId = part->getPDG();
	cout << "t" << counter << ": pdg: " << pdgId << ": E: " << ((int)(100*part->getEnergy()))/100.0 << " GeV; theta: " << partTheta << "; phi: " << partPhi << "; inTracker: " << inTracker << "; inCal: " << inCal << "; genStatus: " << genStatus << "; pointer: " << part << endl;

}

void particleFill::dumpReconstructedPart(const EVENT::ReconstructedParticle* part, const int counter){
	const double *partMom = part->getMomentum();
	TVector3 v1;
	v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
	double partPhi = 180.*v1.Phi()/TMath::Pi();
	double partTheta = 180.*v1.Theta()/TMath::Pi();
	int partType= part->getType();
	cout << "r" << counter << ": E: " << ((int)(100*part->getEnergy()))/100.0 << " GeV; theta: " << partTheta << "; phi: " << partPhi << "; partType: " << partType << endl;
}

int particleFill::fillHist(const double inVal, const string baseString, const string prefix){
	try {
		if (histMap[prefix+baseString]==NULL){
			cout << "[FATAL]\t[particleFill] Can't find histogram with name: " << prefix+baseString << endl;
			throw 1;
		}
	}
	catch (const std::exception& e) {
		cout << "[ERROR]\t[particleFill] Exception was caught, during accessing element: <" << prefix << baseString << "> from histMap:"
			  << endl << e.what() << endl;
	}

	return histMap[prefix+baseString]->Fill(inVal);
}

int particleFill::fillHist(const double inVal1, const double inVal2, const string baseString, const string prefix){
	try {
		if (histMap[prefix+baseString]==NULL){
			cout << "[FATAL]\t[particleFill] Can't find histogram with name: " << prefix+baseString << endl;
			throw 1;
		}
	}
	catch (const std::exception& e) {
		cout << "[ERROR]\t[particleFill] Exception was caught, during accessing element: <" << prefix << baseString << "> from histMap:"
			  << endl << e.what() << endl;
	}

	return histMap[prefix+baseString]->Fill(inVal1,inVal2);
}


double particleFill::getMeanEnergy(){
	return histMap["allPartsOfType_Energy"]->GetMean();
}

double particleFill::getMeanTheta(){
	return histMap["allPartsOfType_Theta"]->GetMean();
}

void particleFill::createHistsFromMap(const map<string,histStruct> inHistStructMap, const string prefix){
	for(auto const &ent1 : inHistStructMap){
		TH1* tmpHist;
		if (ent1.second.histType=="TH1D")
			tmpHist = new TH1D((prefix+ent1.first).c_str(),ent1.second.title.c_str(),ent1.second.nBins,ent1.second.xLow,ent1.second.xHigh);
		if (ent1.second.histType=="TH1I")
			tmpHist = new TH1I((prefix+ent1.first).c_str(),ent1.second.title.c_str(),ent1.second.nBins,ent1.second.xLow,ent1.second.xHigh);
		if (ent1.second.histType=="TH2D")
			tmpHist = new TH2D((prefix+ent1.first).c_str(),ent1.second.title.c_str(),ent1.second.nBins,ent1.second.xLow,ent1.second.xHigh,ent1.second.ynBins,ent1.second.yLow,ent1.second.yHigh);
		histMap[prefix+ent1.first] = tmpHist;
	}
	
}



void particleFill::createSingleRecoParticleClustersHists(const string prefix){
	createHistsFromMap(singleRecoParticleClustersHistStructMap,prefix);
}
void particleFill::createSingleParticleHists(const string prefix){
	createHistsFromMap(singleParticleHistStructMap,prefix);
}
void particleFill::createTwoParticleCorrelationHists(const string prefix){
	createHistsFromMap(twoParticleCorrelationHistStructMap,prefix);
}


int particleFill::writeToFile(TFile* outFile){
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
		if (tmpStrVec.size()<2) continue;
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
void particleFill::clasiffyPFO(EVENT::ReconstructedParticle* inPFO){
	
	int PFOtype = inPFO->getType();
	double PFOenergy = inPFO->getEnergy();

	if (intMap["id_firstEnergeticPFO"] == -1){
		intMap["id_firstEnergeticPFO"] = 0;
		if(find(partTypeToSelect.begin(), partTypeToSelect.end(), inPFO->getType()) != partTypeToSelect.end()){
			intMap["id_mostEnergeticPFOofType"] = 0;
			intMap["n_PFOofType"] = 1;
		}
	}
	else{
		if (PFOenergy>PFOtypeAndEnergyVec.at(intMap["id_firstEnergeticPFO"])){
			intMap["id_secondEnergeticPFO"] = intMap["id_firstEnergeticPFO"];
			intMap["id_firstEnergeticPFO"] = PFOtypeAndEnergyVec.size();
		}
		if(find(partTypeToSelect.begin(), partTypeToSelect.end(), inPFO->getType()) != partTypeToSelect.end()){
			if (PFOenergy>PFOtypeAndEnergyVec.at(intMap["mostEnergeticPFOofType"]))
				intMap["mostEnergeticPFOofType"] = PFOtypeAndEnergyVec.size();
			intMap["n_PFOofType"]++;
		}
	}
	
	PFOtypeAndEnergyVec.push_back( pair(PFOtype, PFOenergy) );

}
