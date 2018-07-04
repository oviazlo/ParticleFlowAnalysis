/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/include/muonEffInJets.h
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	05th Dec 2017
 * 	Last Update:	05th Dec 2017
 */
#ifndef MUONEFFINJETS_H
#define MUONEFFINJETS_H

#include <objectFill.h>
#include <globalConfig.h>
#include "truthZWCondition.h"

#include <cmath>
#include <gsl/gsl_sf_gamma.h>
#include <vector>

#include <marlin/Global.h>
#include <GeometryUtil.h>

#include <DD4hep/DD4hepUnits.h>
#include <DD4hep/DetType.h>
#include <DD4hep/DetectorSelector.h>
#include <DD4hep/Detector.h>
#include <DDRec/DetectorData.h>

#include "EVENT/ReconstructedParticle.h"
#include "EVENT/MCParticle.h"
#include "EVENT/Cluster.h"
#include <lcio.h>

#include <EVENT/CalorimeterHit.h>
#include "CalorimeterHitType.h"

#include <ClusterShapes.h>

struct yokeHitsStruct{
	unsigned int nHits;	
	unsigned int clusterLayerSpan;
	unsigned int nLayers;
};

struct truthPartParameters{
	double truthTheta;
	double cosTruthTheta;
	double truthPhi;
	double truthPt;
	double truthEnergy;
	double vertexR;
	double vertexZ;
};

// partParameters extractTruthParticleParameters(MCParticle* genPart){
//         truthPartParameters tmpPars;
//         TVector3 vTruthMom(genPart->getMomentum());
//         tmpPars.truthTheta = vTruthMom.Theta()*TMath::RadToDeg();
//         tmpPars.cosTruthTheta = TMath::Cos(truthTheta*TMath::DegToRad());
//         tmpPars.truthPhi = vTruthMom.Phi()*TMath::RadToDeg();
//         tmpPars.truthPt = vTruthMom.Pt();
//         tmpPars.truthEnergy = genPart->getEnergy();
//
//         const double *vertexPos = genPart->getVertex();
//         tmpPars.vertexR = sqrt(vertexPos[0]*vertexPos[0]+vertexPos[1]*vertexPos[1]);
//         tmpPars.vertexZ = vertexPos[2];
//         return tmpPars;
// }

class muonEffInJets : public objectFill{
	public:
		muonEffInJets(string _outDirName, string _PFOCollectionName = "PandoraPFO") : objectFill(_outDirName) {PFOCollectionName = _PFOCollectionName;}
		~muonEffInJets(){}

		int init();
		int fillEvent(const EVENT::LCEvent*);
		int writeToFile(TFile* outFile);

	private:
		string PFOCollectionName;
		EVENT::LCCollection *PFOCollection = nullptr;

		EVENT::MCParticle* genPart = nullptr;

		unsigned int nPFOs = 0;

		int fillMuonClusterInfo();
		EVENT::FloatVec getClusterShape(EVENT::Cluster* pCluster);
		yokeHitsStruct getNumberOfYokeHits(EVENT::Cluster* pCluster); 

		vector<EVENT::ReconstructedParticle*> matchedRecoPartsOfType;
		vector<EVENT::ReconstructedParticle*> nonMatchedRecoPartsOfType;
		vector<EVENT::MCParticle*> matchedTruthPartsOfType;
		vector<EVENT::MCParticle*> nonMatchedTruthPartsOfType;

		void findMatchedRecoParts();
		void fillTruthInfo();
		void fillTruthInfoHelper(string histNamePrefix, vector<EVENT::MCParticle*> inTruthPartVec);

};
#endif 



