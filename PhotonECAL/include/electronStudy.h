/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/include/electronStudy.h
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	05th Dec 2017
 * 	Last Update:	05th Dec 2017
 */
#ifndef ELECTRONSTUDY_H
#define ELECTRONSTUDY_H

#include <objectFill.h>
#include <energyFill.h>
#include <globalConfig.h>
#include "truthCondition.h"
#include "Objects/CartesianVector.h"
#include "Objects/TrackState.h"
#include "EVENT/TrackState.h"
// #include "Objects/Cluster.h"

// struct PFOMergeSettings{
//         unsigned int pfoTypeToMerge;
//         double thetaCone;
//         double phiCone;
//         double phiConeMomentumDep = 0.0;
// };

class electronStudy : public objectFill{
	public:
		electronStudy(string _outDirName, string _PFOCollectionName = "PandoraPFOs") : objectFill(_outDirName) {PFOCollectionName = _PFOCollectionName;}
		~electronStudy(){}

		int init();
		int fillEvent(const EVENT::LCEvent*);
		int writeToFile(TFile* outFile);

	private:
		unsigned int nSelectecTruthParticles;
		string PFOCollectionName;
		string trackCollectionName;
		string ecalBarrelCollectionName;
		EVENT::LCCollection *PFOCollection;
		EVENT::LCCollection *trackCollection;
		EVENT::LCCollection *ecalBarrelCollection;

		UTIL::BitField64* _encoder;

		map <string,map <string, unsigned int> > categoryMap;
		// if pfoId<0 - dont apply prefix
		void fillHistsPerCategory(string histNameCore, double fillValue, int pfoId);
		map <string, unsigned int> pfoCounter;
		map <unsigned int, unsigned int> energyRanking;
		// (energy, pfoID, pfoType)
		vector<pair<double, pair<unsigned int, unsigned int> > > pfoIdSortedByEnergyAndType; 
		map<unsigned int, string> pfoIdEnergyTypeMap;
		vector<EVENT::ReconstructedParticle*> recoPFOs;

		void performEnergyPfoTypeRanking();

		TVector3* getTrackStateMomentum(EVENT::Track *inTrack);
		TVector3* getTrackStatePosition(EVENT::Track *inTrack);

		// void CopyTrackState(const EVENT::TrackState *const pTrackState, pandora::TrackState &inputTrackState) const;
		pandora::TrackState* getPandoraTrackState(const EVENT::TrackState *const pTrackState);

		int getLayerNumber(EVENT::CalorimeterHit* calHit);
		unsigned int maxCaloSearchLayer = 9;

		map <string, double> getTrackClusterDistance(const EVENT::ReconstructedParticle* const inPart, const EVENT::Track* const inTrack);

};
#endif // ELECTRONSTUDY_H
