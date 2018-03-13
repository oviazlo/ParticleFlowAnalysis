/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/src/electronStudy.cpp
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	05th Dec 2017
 * 	Last Update:	05th Dec 2017
 */

#include "electronStudy.h"

/*===========================================================================*/
/*===============================[ function implementations ]================*/
/*===========================================================================*/

int electronStudy::init(){

	if (config::vm.count("debug"))
		cout << "[INFO]\telectronStudy::init()" << endl;

	// TODO hardcoded
	// trackCollectionName = "SiTracks_Refitted";
	trackCollectionName = "SiTracks";
	ecalBarrelCollectionName = "ECalBarrelCollection";

	TH1* tmpHist;
	tmpHist = new TH1I("nTracks","Number of Tracks; Number of tracks; Counts",5,0,5);
	histMap["nTracks"] = tmpHist;

	map <string, unsigned int> tmpMap;
	string tmpString = "";

	tmpMap.clear();
	tmpString = "OneElectron";
	tmpMap["Electron"] = 1;
	tmpMap["Total"] = 1;
	categoryMap[tmpString] = tmpMap;
	tmpMap.clear();

	tmpMap.clear();
	tmpString = "OnePion";
	tmpMap["Pion"] = 1;
	tmpMap["Total"] = 1;
	categoryMap[tmpString] = tmpMap;
	tmpMap.clear();

	tmpMap.clear();
	tmpString = "OnePhoton";
	tmpMap["Photon"] = 1;
	tmpMap["Total"] = 1;
	categoryMap[tmpString] = tmpMap;
	tmpMap.clear();

	tmpMap.clear();
	tmpString = "OneNeutralHadron";
	tmpMap["NeutralHadron"] = 1;
	tmpMap["Total"] = 1;
	categoryMap[tmpString] = tmpMap;
	tmpMap.clear();

	tmpMap.clear();
	tmpString = "TwoPhotons";
	tmpMap["Photon"] = 2;
	tmpMap["Total"] = 2;
	categoryMap[tmpString] = tmpMap;
	tmpMap.clear();

	tmpMap.clear();
	tmpString = "OneElectronOnePhoton";
	tmpMap["Photon"] = 1;
	tmpMap["Electron"] = 1;
	tmpMap["Total"] = 2;
	categoryMap[tmpString] = tmpMap;
	tmpMap.clear();

	tmpMap.clear();
	tmpString = "OneElectronTwoPhotons";
	tmpMap["Photon"] = 2;
	tmpMap["Electron"] = 1;
	tmpMap["Total"] = 3;
	categoryMap[tmpString] = tmpMap;
	tmpMap.clear();

	tmpMap.clear();
	tmpString = "OneElectronPlusMore";
	tmpMap["Electron"] = 1;
	categoryMap[tmpString] = tmpMap;
	tmpMap.clear();

	tmpMap.clear();
	tmpString = "NoElectron";
	tmpMap["Electron"] = 0;
	categoryMap[tmpString] = tmpMap;
	tmpMap.clear();

	tmpMap.clear();
	tmpString = "NoElectronNoPion";
	tmpMap["Electron"] = 0;
	tmpMap["Pion"] = 0;
	categoryMap[tmpString] = tmpMap;
	tmpMap.clear();

	for (auto it = categoryMap.begin(); it != categoryMap.end(); it++){
		string prefix = it->first;

		// HISTS with type + energyIndex postfix
		for (auto iType = config::pfoTypeIntStringMap.begin(); iType != config::pfoTypeIntStringMap.end(); iType++) {
			for (int ii=0; ii<20; ii++){
				string postFix = iType->second + DoubToStr(ii);

				tmpString = prefix + "_" + "dPhi" + "-" + postFix;
				tmpHist = new TH1I(tmpString.c_str(),"#phi_{reco}-#phi_{truth}; Phi [rad]; Counts",400,-0.5,0.5);
				histMap[tmpString] = tmpHist;
				
				tmpString = prefix + "_" + "dTheta" + "-" + postFix;
				tmpHist = new TH1I(tmpString.c_str(),"#Theta_{reco}-#Theta_{truth}; Theta [rad]; Counts",1000,-0.025,0.025);
				histMap[tmpString] = tmpHist;

				tmpString = prefix + "_" + "dPhiTrack" + "-" + postFix;
				tmpHist = new TH1I(tmpString.c_str(),"#phi_{pfo}-#phi_{track}; Phi [rad]; Counts",1200,-TMath::Pi(),TMath::Pi());
				histMap[tmpString] = tmpHist;
				
				tmpString = prefix + "_" + "dThetaTrack" + "-" + postFix;
				tmpHist = new TH1I(tmpString.c_str(),"#Theta_{pfo}-#Theta_{track}; Theta [rad]; Counts",1200,-TMath::Pi(),TMath::Pi());
				histMap[tmpString] = tmpHist;

				tmpString = prefix + "_" + "E" + "-" + postFix;
				tmpHist = new TH1I(tmpString.c_str(),"PFO energy; E [GeV]; Counts",5000,0,125);
				histMap[tmpString] = tmpHist;

				tmpString = prefix + "_" + "nClusters" + "-" + postFix;
				tmpHist = new TH1I(tmpString.c_str(),"Number of clusters per PFO; Number of clusters; Counts",5,0,5);
				histMap[tmpString] = tmpHist;

				tmpString = prefix + "_" + "trackClusterParallelDistance" + "-" + postFix;
				tmpHist = new TH1I(tmpString.c_str(),"Track-Cluster Parallel Distance; Distance [mm]; Counts",1000,0.,250);
				histMap[tmpString] = tmpHist;

				tmpString = prefix + "_" + "trackClusterPerpendicularDistance" + "-" + postFix;
				tmpHist = new TH1I(tmpString.c_str(),"Track-Cluster Perpendicular Distance; Distance [mm]; Counts",1000,0.,250);
				histMap[tmpString] = tmpHist;

				tmpString = prefix + "_" + "trackClusterPerpendicularDistanceWithCut" + "-" + postFix;
				tmpHist = new TH1I(tmpString.c_str(),"Track-Cluster Perpendicular Distance with Cut; Distance [mm]; Counts",1000,0.,250);
				histMap[tmpString] = tmpHist;
				
				tmpString = prefix + "_" + "trackClusterCosAngle" + "-" + postFix;
				tmpHist = new TH1I(tmpString.c_str(),"Track-Cluster Cosine Angle; Cos(Opening Angle); Counts",1000,-1.,1.);
				histMap[tmpString] = tmpHist;

				tmpString = prefix + "_" + "trackZ0" + "-" + postFix;
				tmpHist = new TH1I(tmpString.c_str(),"Track Z0; Z0 [mm]; Counts",800,-100,100);
				histMap[tmpString] = tmpHist;
				
				tmpString = prefix + "_" + "trackD0" + "-" + postFix;
				tmpHist = new TH1I(tmpString.c_str(),"Track D0; D0 [mm]; Counts",800,-100,100);
				histMap[tmpString] = tmpHist;

				tmpString = prefix + "_" + "trackR0" + "-" + postFix;
				tmpHist = new TH1I(tmpString.c_str(),"Radius of the innermost hit that has been used in the track fit; R0 [mm]; Counts",400,0,100);
				histMap[tmpString] = tmpHist;
			}
		}

		// OTHER HISTS
		tmpString = prefix + "_" + "nTracks";
		tmpHist = new TH1I(tmpString.c_str(),"Number of Tracks; Number of tracks; Counts",5,0,5);
		histMap[tmpString] = tmpHist;

		tmpString = prefix + "_" + "nTracks";
		tmpHist = new TH1I(tmpString.c_str(),"Number of Tracks; Number of tracks; Counts",5,0,5);
		histMap[tmpString] = tmpHist;

		tmpString = prefix + "_" + "truthCosTheta";
		tmpHist = new TH1I(tmpString.c_str(),"Truth Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
		histMap[tmpString] = tmpHist;

		tmpString = prefix + "_" + "truthPhi";
		tmpHist = new TH1I(tmpString.c_str(),"Truth #Phi; #Phi; Counts per Event",180*2,-180,180);
		histMap[tmpString] = tmpHist;
	}




	// else if (pfoCounter["Electron"]==1 && pfoCounter["Total"]>1){
	//
	// }
	// else if (pfoCounter["Electron"]==0 && pfoCounter["Photon"]>2){
	//
	// }


	return 0;
}

int electronStudy::fillEvent(const EVENT::LCEvent* event){
	try {
		PFOCollection = event->getCollection(PFOCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|electronStudy]\tCan't find collection: " << PFOCollectionName << endl;
		return -1;
	}
	try {
		trackCollection = event->getCollection(trackCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|electronStudy]\tCan't find collection: " << trackCollectionName << endl;
		return -1;
	}
	try {
		ecalBarrelCollection = event->getCollection(ecalBarrelCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|electronStudy]\tCan't find collection: " << ecalBarrelCollectionName << endl;
		return -1;
	}
	_encoder = new UTIL::BitField64(ecalBarrelCollection->getParameters().getStringVal( EVENT::LCIO::CellIDEncoding ));

	// cout << "Event " << event->getEventNumber() << endl;
	// cout << "NumberOfElements in " << PFOCollectionName << " collection: " << PFOCollection->getNumberOfElements() << endl;

	// init pfoCounter
	pfoCounter.clear();
	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++)
		pfoCounter[it->second] = 0;

	if (config::vm.count("debug"))
		cout << "[INFO]	electronStudy::fillEvent: " << event->getEventNumber() << endl;

	EVENT::MCParticle* genPart = truthCondition::instance()->getGunParticle();
	const double *partMom = genPart->getMomentum();
	TVector3 v1;
	v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
	double truthTheta = 180.*v1.Theta()/TMath::Pi();
	double truthPhi = 180.*v1.Phi()/TMath::Pi();
	double cosTruthTheta = TMath::Cos(TMath::Pi()*truthTheta/180.);
	double truthEnergy = genPart->getEnergy();
	double truthPt = v1.Pt();

	recoPFOs = getObjVecFromCollection<EVENT::ReconstructedParticle*>(PFOCollection);
	if (recoPFOs.size()==0){
		cout << "[WARNING]\t Event " << event->getEventNumber() << " has no PFOs!!!" << endl;
		return 0; // no reco PFOs
	}

	pfoIdSortedByEnergyAndType.clear();
	for (int i=0; i<recoPFOs.size(); i++){
		unsigned int pfoType = abs(recoPFOs[i]->getType());
		double pfoEnergy = recoPFOs[i]->getEnergy();
		pfoIdSortedByEnergyAndType.push_back( make_pair( pfoEnergy, make_pair(i, pfoType) ) );
		// double energyRanking
	}
	sort(pfoIdSortedByEnergyAndType.begin(),pfoIdSortedByEnergyAndType.end());
	reverse(pfoIdSortedByEnergyAndType.begin(),pfoIdSortedByEnergyAndType.end());
	// cout << "nPFOs: " << pfoIdSortedByEnergyAndType.size() << endl;
	performEnergyPfoTypeRanking();

	vector<EVENT::Track*> tracks = getObjVecFromCollection<EVENT::Track*>(trackCollection);
	getHistFromMap("nTracks")->Fill(tracks.size());
	fillHistsPerCategory("nTracks",tracks.size(),-1);
	fillHistsPerCategory("truthCosTheta",cosTruthTheta,-1);
	fillHistsPerCategory("truthPhi",truthPhi,-1);

	for (int i=0; i<recoPFOs.size(); i++){
		int pfoType = abs(recoPFOs[i]->getType());
		const double *partMom = recoPFOs[i]->getMomentum();
		TVector3 v1;
		v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
		double partTheta = 180.*v1.Theta()/TMath::Pi();
		double partPhi = 180.*v1.Phi()/TMath::Pi();
		double partPt = v1.Pt();
		double cosPartTheta = TMath::Cos(TMath::Pi()*partTheta/180.);
		double pfoE = recoPFOs[i]->getEnergy();
		double dPhi = TMath::Pi()*(partPhi-truthPhi)/180.0;
		double dTheta = TMath::Pi()*(partTheta-truthTheta)/180.0;

		fillHistsPerCategory("dPhi",dPhi,i);
		fillHistsPerCategory("dTheta",dTheta,i);
		fillHistsPerCategory("E",pfoE,i);
		fillHistsPerCategory("nClusters",recoPFOs[i]->getClusters().size(),i);


		if (config::vm.count("debug")) cout << endl << "access cluster" << endl;
		auto clusterVec = recoPFOs[i]->getClusters();
		if (config::vm.count("debug")) cout << "number of clusters: " << clusterVec.size() << "; partType: " << pfoType << endl;
		if (clusterVec.size()>0){
			auto cluster = clusterVec[0];
			if (config::vm.count("debug")) cout << "access caloHit vector" << endl;
			auto tmpVec = cluster->getCalorimeterHits();
			if (config::vm.count("debug")) cout << "size of caloHit vector: " << tmpVec.size() << endl;
			if (config::vm.count("debug")) cout << "continue further" << endl;
			if (config::vm.count("debug")) cout << endl;
		}
		// cout << "Size of caloHits vector: " << tmpVec.size() << endl;
		
		// cout << "layerNumber: " << getLayerNumber(tmpVec[0]) << endl;
		// pandora::Cluster* myCluster = new pandora::Cluster(recoPFOs[i]->getClusters()[0]);
		// for (auto it = tmpVec.begin(); it!=tmpVec.end(); it++){
		//         myCluster->AddCaloHit(it);
		// }
		
		if (tracks.size()==1){
			auto trackStateMomentum = getTrackStateMomentum(tracks[0]);
			fillHistsPerCategory("dPhiTrack",v1.Phi()-trackStateMomentum->Phi(),i);
			fillHistsPerCategory("dThetaTrack",v1.Theta()-trackStateMomentum->Theta(),i);

			auto trackClusterDistanceMap = getTrackClusterDistance(recoPFOs[i],tracks[0]);

			fillHistsPerCategory("trackClusterPerpendicularDistance",trackClusterDistanceMap["trackClusterPerpendicularDistance"],i);
			fillHistsPerCategory("trackClusterPerpendicularDistanceWithCut",trackClusterDistanceMap["trackClusterPerpendicularDistanceWithCut"],i);
			fillHistsPerCategory("trackClusterParallelDistance",trackClusterDistanceMap["trackClusterParallelDistance"],i);
			fillHistsPerCategory("trackClusterCosAngle",trackClusterDistanceMap["trackClusterCosAngle"],i);

			fillHistsPerCategory("trackZ0",tracks[0]->getZ0(),i);
			fillHistsPerCategory("trackD0",tracks[0]->getD0(),i);
			fillHistsPerCategory("trackR0",tracks[0]->getRadiusOfInnermostHit(),i);
			
			// for (auto it=trackClusterDistanceMap.begin(); it!=trackClusterDistanceMap.end(); it++){
			//         cout << "[Sasha]\t " << (*it).first << ": " << (*it).second << endl;
			// }
			// cout << endl;

		}
		


	}

	return 0;
}

void electronStudy::fillHistsPerCategory(string histNameCore, double fillValue, int pfoId){
	for (auto it = categoryMap.begin(); it != categoryMap.end(); it++){
		bool passRequirement = true;
		for (auto itTmp = it->second.begin(); itTmp != it->second.end(); itTmp++){
			// cout << "itTmp->first = " << itTmp->first << "itTmp->second = " << itTmp->second << "; pfoCounter[itTmp->first] = " << pfoCounter[itTmp->first] << endl;
			if (itTmp->second != pfoCounter[itTmp->first]){
				passRequirement = false;
				break;
			}

		}
		if (passRequirement==false)
			continue;
		string prefix = it->first;
		
		string tmpString = prefix + "_" + histNameCore;
		if (pfoId>=0)
			tmpString = tmpString + "-" + pfoIdEnergyTypeMap[pfoId];
		getHistFromMap(tmpString)->Fill(fillValue);
	}
}





int electronStudy::writeToFile(TFile* outFile){

	if (config::vm.count("debug"))
		cout << "[INFO]	electronStudy::writeToFile(" << outFile->GetName() << ")" << endl;

	// getHistFromMap("totalEnergyVsTheta")->Divide(getHistFromMap("nTruthPartsVsTheta"));
	// getHistFromMap("nTruthPartsVsTheta")->Sumw2();
	

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
		if (it.second->GetEntries()==0)
			continue;
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

void electronStudy::performEnergyPfoTypeRanking(){
	for (auto iType = config::pfoTypeIntStringMap.begin(); iType != config::pfoTypeIntStringMap.end(); iType++) {
		unsigned int pdgId = iType->first;
		string typeString = iType->second;
		for (int iCount=0; iCount<pfoIdSortedByEnergyAndType.size(); iCount++){
			pair<unsigned int, unsigned int> idTypePair = pfoIdSortedByEnergyAndType[iCount].second;
			unsigned int pfoId = idTypePair.first;
			unsigned int pfoType = idTypePair.second;

			if (pdgId == pfoType){
				pfoIdEnergyTypeMap[ pfoId ] = typeString + DoubToStr(pfoCounter[typeString]+1);
				pfoCounter[typeString]++;
				pfoCounter["Total"]++;
				// cout << "ID: " << pfoId << "; Type: " << pfoType << "; Energy: " << pfoIdSortedByEnergyAndType[iCount].first << "; label: " << pfoIdEnergyTypeMap[ pfoId ] << endl;
			}
		}
	}
}

TVector3* electronStudy::getTrackStateMomentum(EVENT::Track *inTrack){
	double bField = 2.0;
	auto pTrackState = inTrack->getTrackState (EVENT::TrackState::AtCalorimeter);
	const double pt(bField * 2.99792e-4 / std::fabs(pTrackState->getOmega()));

	const double px(pt * std::cos(pTrackState->getPhi()));
	const double py(pt * std::sin(pTrackState->getPhi()));
	const double pz(pt * pTrackState->getTanLambda());

	TVector3* outVec = new TVector3(px,py,pz);
	return outVec;
}


TVector3* electronStudy::getTrackStatePosition(EVENT::Track *inTrack){
	auto pTrackState = inTrack->getTrackState (EVENT::TrackState::AtCalorimeter);
	const double xs(pTrackState->getReferencePoint()[0]);
	const double ys(pTrackState->getReferencePoint()[1]);
	const double zs(pTrackState->getReferencePoint()[2]);
	TVector3* outVec = new TVector3(xs,ys,zs);
	return outVec;
}

int electronStudy::getLayerNumber(EVENT::CalorimeterHit* calHit){
	lcio::long64 cellId = long( calHit->getCellID0() & 0xffffffff ) | ( long( calHit->getCellID1() ) << 32 );
	_encoder->setValue(cellId);
	return (*_encoder)["layer"].value();
}


pandora::TrackState* electronStudy::getPandoraTrackState(const EVENT::TrackState *const pTrackState)
{
	if (!pTrackState)
		throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

	double m_bField = 2.0;
	const double pt(m_bField * 2.99792e-4 / std::fabs(pTrackState->getOmega()));

	const double px(pt * std::cos(pTrackState->getPhi()));
	const double py(pt * std::sin(pTrackState->getPhi()));
	const double pz(pt * pTrackState->getTanLambda());

	const double xs(pTrackState->getReferencePoint()[0]);
	const double ys(pTrackState->getReferencePoint()[1]);
	const double zs(pTrackState->getReferencePoint()[2]);

	// cout << "[DEBUG|electronStudy::CopyTrackState]\t" << " xs: " << xs  << " ys: " << ys  << " zs: " << zs  << " px: " << px  << " py: " << py  << " pz: " << pz << endl;

	return new pandora::TrackState(xs, ys, zs, px, py, pz);
}


map <string, double> electronStudy::getTrackClusterDistance(const EVENT::ReconstructedParticle* const inPart, const EVENT::Track* const inTrack){
	map<string,double> outMap;
	auto clusterVec = inPart->getClusters();
	float minDistanceSquared(std::numeric_limits<float>::max());
	float minDistanceSquaredWithCut(std::numeric_limits<float>::max());
	float minTrackClusterCosAngle(std::numeric_limits<float>::max());
	pandora::TrackState* pTrackState = getPandoraTrackState(inTrack->getTrackState (EVENT::TrackState::AtCalorimeter));
	if (clusterVec.size()>0){
		auto cluster = clusterVec[0];
		auto caloHitVec = cluster->getCalorimeterHits();
		pandora::CartesianVector* tmpVec = new pandora::CartesianVector(pTrackState->GetPosition());
		const pandora::CartesianVector &trackPosition(pTrackState->GetPosition());
		const pandora::CartesianVector trackDirection(pTrackState->GetMomentum().GetUnitVector());

		const pandora::CartesianVector* clusterDirection = new pandora::CartesianVector(inPart->getMomentum()[0],inPart->getMomentum()[1],inPart->getMomentum()[2]);
		double trackClusterCosAngle = trackDirection.GetCosOpeningAngle(clusterDirection->GetUnitVector());

		double minParallelDistance = std::numeric_limits<double>::max();

		EVENT::CalorimeterHitVec::const_iterator savedHitIter;

		unsigned int hitCounter = 0;
		unsigned int previousHitLayer = 0;
		for (auto hitIter=caloHitVec.begin(); hitIter!=caloHitVec.end(); hitIter++){
			unsigned int hitLayer = getLayerNumber(*hitIter);
			// if (hitLayer>maxCaloSearchLayer)
			//         continue;
			auto caloHitPositionVec = (*hitIter)->getPosition();
			// cout << "caloHitID: " << hitCounter << "\tlayer: " << hitLayer << endl;
			const pandora::CartesianVector* caloHitPosition = new pandora::CartesianVector(caloHitPositionVec[0],caloHitPositionVec[1],caloHitPositionVec[2]);
			// const pandora::CartesianVector &caloHitPosition((*hitIter)->getPosition()[0],(*hitIter)->getPosition()[1],(*hitIter)->getPosition()[2]);
			const pandora::CartesianVector positionDifference(*caloHitPosition - trackPosition);
			double parallelDistance = std::fabs(trackDirection.GetDotProduct(positionDifference));
			const float perpendicularDistanceSquared((trackDirection.GetCrossProduct(positionDifference)).GetMagnitudeSquared());

			// cout << "HitID: " << hitCounter << "\tlayer: " << hitLayer << "\tpos[x,y,z]: " << caloHitPositionVec[0] << "\t" << caloHitPositionVec[1] << "\t" << caloHitPositionVec[2] << "\tr: " << sqrt(caloHitPositionVec[0]*caloHitPositionVec[0]+caloHitPositionVec[1]*caloHitPositionVec[1]) << "\ttrackClusterDistance: " << sqrt(perpendicularDistanceSquared) <<  endl;
			hitCounter++;
			// TODO remove 2 lines below and uncomment the same lines above!!!
			// if (hitLayer>maxCaloSearchLayer)
			//         continue;

			if (hitLayer>maxCaloSearchLayer)
				break;
			if (hitLayer<previousHitLayer)
				break;

			previousHitLayer = hitLayer;
			if (minParallelDistance>parallelDistance)
				minParallelDistance = parallelDistance;
			if (perpendicularDistanceSquared < minDistanceSquared)
				minDistanceSquared = perpendicularDistanceSquared;
			if (perpendicularDistanceSquared < minDistanceSquaredWithCut && parallelDistance<100.0 && trackClusterCosAngle>0.0){
				minDistanceSquaredWithCut = perpendicularDistanceSquared;
				savedHitIter = hitIter;
			}
		}

		outMap["trackClusterPerpendicularDistance"] = std::sqrt(minDistanceSquared);
		outMap["trackClusterPerpendicularDistanceWithCut"] = std::sqrt(minDistanceSquaredWithCut);
		outMap["trackClusterParallelDistance"] = minParallelDistance;
		outMap["trackClusterCosAngle"] = trackClusterCosAngle;

		// cout << "[Sasha]\tnCaloHitsInCluster: " << caloHitVec.size() << endl;
		// cout << "[Sasha]\ttrackDirection: " << trackDirection << endl;
		// cout << "[Sasha]\ttrackPosition: " << trackPosition << endl;
		// auto caloHitPositionVec = (*savedHitIter)->getPosition();
		// const pandora::CartesianVector* caloHitPosition = new pandora::CartesianVector(caloHitPositionVec[0],caloHitPositionVec[1],caloHitPositionVec[2]);
		// cout << "[Sasha]\tcaloHitPosition: " << (*caloHitPosition) << endl;
		// const pandora::CartesianVector positionDifference(*caloHitPosition - trackPosition);
		// cout << "[Sasha]\tpositionDifference: " << positionDifference << endl;
		// cout << "[Sasha]\tcaloHitLayer: " << getLayerNumber(*savedHitIter) << endl;



	}
	// cout << "1" << endl;
	return outMap;
}

 
