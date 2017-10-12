#include <particleFill.h>

#include <IMPL/ReconstructedParticleImpl.h>

#define nBinsPerGeV 10

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

	// singleRecoParticleClustersHistStructMap["clusterPositionZY"] = histStruct("Cluster Position; Z [mm]; Y [mm]",800,-4000,4000,"TH2D",800,-4000,4000);
	singleRecoParticleClustersHistStructMap["clusterPositionZR"] = histStruct("Cluster Position; Z [mm]; R [mm]",800,-4000,4000,"TH2D",400,0,4000);
	// singleRecoParticleClustersHistStructMap["clusterPositionXY"] = histStruct("Cluster Position; X [mm]; Y [mm]",800,-4000,4000,"TH2D",800,-4000,4000);

}

int particleFill::init(){

	initHistStructs();

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
		createSingleParticleHists("firstEnergeticPartOfType_");
		createSingleParticleHists("firstEnergeticPartOfType_noAdditionalPFOs_");
		createSingleParticleHists("firstEnergeticPartOfType_thereAreAdditionalPFOs_");
		createSingleParticleHists("firstEnergeticPartOfWrongType_");
		createSingleParticleHists("secondEnergeticPart_");
		// truth hists:
		createSingleParticleHists("truthParticle_onlyOneRecoPFO_");
		createSingleParticleHists("truthParticle_twoOrMoreRecoPFO_");

		createTwoParticleCorrelationHists("signalPFO_secondEnergeticPFO_");
		createTwoParticleCorrelationHists("signalPFO_secondEnergeticPhoton_");

		createSingleRecoParticleClustersHists("test_");

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

	return 0;
}



int particleFill::fillParticleCorrelations (EVENT::ReconstructedParticle* inPart1, EVENT::ReconstructedParticle* inPart2, string prefix){
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

vector<EVENT::MCParticle*> particleFill::getTruthMCParticlesFromCollection(EVENT::LCEvent* event, string truthCollectionName){
	// FIXME hardcoding of the collection name
	vector<EVENT::MCParticle*> outPartVec;
	EVENT::LCCollection *truthCollection;
	try {
		truthCollection = event->getCollection(truthCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|particleFill]\tCan't find collection: " << truthCollectionName << endl;
		return outPartVec;
	}

	const int nElements = collection->getNumberOfElements();
	for(int j=0; j < nElements; j++) {
		auto part = static_cast<EVENT::MCParticle*>(truthCollection->getElementAt(j));
		outPartVec.push_back(part);
	}
	return outPartVec;
}

int particleFill::fillEvent(EVENT::LCEvent* event){
	try {
		collection = event->getCollection(collectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|particleFill]\tCan't find collection: " << collectionName << endl;
		return -1;
	}

	auto truthParts = getTruthMCParticlesFromCollection(event);

	if( collection ) {
		string collectionType = collection->getTypeName();
		const int nElements = collection->getNumberOfElements();
		if (nElements<=0)
			return 0;
		unsigned int countE_1 = 0;
		unsigned int countE_2 = 0;
		double E_1 = 0.0;
		double E_2 = 0.0;

      		for(int j=0; j < nElements; j++) {
			if (collectionType=="MCParticle"){
				auto part = static_cast<EVENT::MCParticle*>(collection->getElementAt(j));
				if (debugFlag) dumpTruthPart(part, j);
				if (part->isCreatedInSimulation()!=0) continue;
				if (part->getGeneratorStatus()!=1) continue;
				fillPart(part);
			}
			else if (collectionType=="Track"){
				auto track = static_cast<EVENT::Track*>(collection->getElementAt(j));
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
				auto part = static_cast<EVENT::ReconstructedParticle*>(collection->getElementAt(j));
				if (debugFlag) dumpReconstructedPart(part, j);
				if(find(partTypeToSelect.begin(), partTypeToSelect.end(), part->getType()) != partTypeToSelect.end()){
					fillPart(part,"allPartsOfType_");
				}
				if (j==0){
					E_1 = part->getEnergy();	
				}
				else{
					double E = part->getEnergy();
					if (E>E_1){
						countE_2 = countE_1;
						E_2 = E_1;
						E_1 = E;
						countE_1 = j;
					}
					if (E>E_2 && E<E_1){
						E_2 = E;
						countE_2 = j;
					}
				}
				// cout << j << ":\tE1: " << E_1 << "; E2: " << E_2 << endl;
			}
			else{
				cout << "[ERROR|fillEvent]\tuknown particle collection: " << collectionType << endl;
			}
		}
		if (collectionType=="ReconstructedParticle"){
			auto part = static_cast<IMPL::ReconstructedParticleImpl*>(collection->getElementAt(countE_1));
			
			// FIXME WARNING testing merging of clusters...
			if (dPhiMergeValue > 0.0){
				auto tmpPart = static_cast<IMPL::ReconstructedParticleImpl*>(collection->getElementAt(countE_1));
				part = new IMPL::ReconstructedParticleImpl();
				part->setEnergy(tmpPart->getEnergy());
				part->setMomentum(tmpPart->getMomentum());
				part->setType(tmpPart->getType());

				for(int j=0; j < nElements; j++) {
					if (j==countE_1) continue;
					TVector3 v1, v2;
					const double *partMom1 = part->getMomentum();
					v1.SetXYZ(partMom1[0],partMom1[1],partMom1[2]);

					auto part2 = static_cast<EVENT::ReconstructedParticle*>(collection->getElementAt(j));
					const double *partMom2 = part2->getMomentum();
					v2.SetXYZ(partMom2[0],partMom2[1],partMom2[2]);
					double dPhi = v1.DeltaPhi(v2)*180./TMath::Pi();
					if (abs(dPhi)<2.0)
						part->setEnergy(part->getEnergy()+part2->getEnergy());
				}
			}
			// FIXME
			if(find(partTypeToSelect.begin(), partTypeToSelect.end(), part->getType()) != partTypeToSelect.end()){

				fillPart(part,"firstEnergeticPartOfType_");
				fillClusterInfo(part,"test_");
				if (nElements>=2 && countE_1!=countE_2){
					auto part2 = static_cast<EVENT::ReconstructedParticle*>(collection->getElementAt(countE_2));
					fillPart(part2,"secondEnergeticPart_");
					fillParticleCorrelations(part,part2,"signalPFO_secondEnergeticPFO_");
					if (part2->getType()==22)
						fillParticleCorrelations(part,part2,"signalPFO_secondEnergeticPhoton_");
				}
				// if (nElements>=2 && countE_1==countE_2){
				//         cout << "countE: " << countE_2 << "; E1: " << E_1 << "; E2: " << E_2 << endl;
				// }


				if (nElements>=2){
					fillPart(part,"firstEnergeticPartOfType_thereAreAdditionalPFOs_");
					fillPart(truthParts[0],"truthParticle_twoOrMoreRecoPFO_");
				}
				else{
					fillPart(part,"firstEnergeticPartOfType_noAdditionalPFOs_");
					fillPart(truthParts[0],"truthParticle_onlyOneRecoPFO_");
				}


			}
			else
				fillPart(part,"firstEnergeticPartOfWrongType_");
				

		}
	}
	return 0;
}

int particleFill::fillPart (EVENT::MCParticle* inPart, string prefix) {
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

int particleFill::fillPart (EVENT::ReconstructedParticle* inPart, string prefix){
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

int particleFill::fillClusterInfo (EVENT::ReconstructedParticle* inPart, string prefix){
	vector<EVENT::Cluster*> clusterVec = inPart->getClusters();
	for (int i=0; i<clusterVec.size(); i++){
		 const float *pos = clusterVec[i]->getPosition();
		 histMap[prefix+"clusterPositionZR"]->Fill(pos[2],sqrt(pos[1]*pos[1]+pos[0]*pos[0]));
	}

	return 0;
}

void particleFill::dumpTruthPart(EVENT::MCParticle* part, int counter){
	const double *partMom = part->getMomentum();
	TVector3 v1;
	v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
	double partPhi = 180.*v1.Phi()/TMath::Pi();
	double partTheta = 180.*v1.Theta()/TMath::Pi();
	bool inTracker = part->isDecayedInTracker();
	bool inCal = part->isDecayedInCalorimeter();
	int genStatus = part->getGeneratorStatus();
	cout << "t" << counter << ": E: " << ((int)(100*part->getEnergy()))/100.0 << " GeV; theta: " << partTheta << "; phi: " << partPhi << "; inTracker: " << inTracker << "; inCal: " << inCal << "; genStatus: " << genStatus << endl;

}

void particleFill::dumpReconstructedPart(EVENT::ReconstructedParticle* part, int counter){
	const double *partMom = part->getMomentum();
	TVector3 v1;
	v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
	double partPhi = 180.*v1.Phi()/TMath::Pi();
	double partTheta = 180.*v1.Theta()/TMath::Pi();
	int partType= part->getType();
	cout << "r" << counter << ": E: " << ((int)(100*part->getEnergy()))/100.0 << " GeV; theta: " << partTheta << "; phi: " << partPhi << "; partType: " << partType << endl;
}

int particleFill::fillHist(double inVal, string baseString, string prefix){
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

double particleFill::getMeanEnergy(){
	return histMap["allPartsOfType_Energy"]->GetMean();
}

double particleFill::getMeanTheta(){
	return histMap["allPartsOfType_Theta"]->GetMean();
}

void particleFill::createHistsFromMap(map<string,histStruct> inHistStructMap, string prefix){
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



void particleFill::createSingleRecoParticleClustersHists(string prefix){
	createHistsFromMap(singleRecoParticleClustersHistStructMap,prefix);
}
void particleFill::createSingleParticleHists(string prefix){
	createHistsFromMap(singleParticleHistStructMap,prefix);
}
void particleFill::createTwoParticleCorrelationHists(string prefix){
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
		for (int i=0; i<tmpStrVec.size()-1; i++)
			prefix += tmpStrVec[i];
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

