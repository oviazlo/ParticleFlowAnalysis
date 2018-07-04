/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/include/muonClusterFiller.h
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	05th Dec 2017
 * 	Last Update:	05th Dec 2017
 */
#ifndef MUONCLUSTERFILLER_H
#define MUONCLUSTERFILLER_H

#include <objectFill.h>
#include <globalConfig.h>
#include "truthCondition.h"

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

class muonClusterFiller : public objectFill{
	public:
		muonClusterFiller(string _outDirName, string _PFOCollectionName = "PandoraPFO") : objectFill(_outDirName) {PFOCollectionName = _PFOCollectionName;}
		~muonClusterFiller(){}

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
};
#endif 



