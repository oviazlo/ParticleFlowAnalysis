/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/include/trackerTiming.h
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	05th Dec 2017
 * 	Last Update:	05th Dec 2017
 */
#ifndef TRACKERTIMING_H
#define TRACKERTIMING_H

#include <objectFill.h>
#include <globalConfig.h>
#include "truthCondition.h"

#include <cmath>
#include <vector>

#include "EVENT/ReconstructedParticle.h"
#include <lcio.h>

#include <UTIL/CellIDDecoder.h>
#include <EVENT/TrackerHit.h>

class trackerTiming : public objectFill{
	public:
		trackerTiming(string _outDirName, string _PFOCollectionName = "PandoraPFO") : objectFill(_outDirName) {PFOCollectionName = _PFOCollectionName;}
		~trackerTiming(){}

		int init();
		int fillEvent(const EVENT::LCEvent*);
		int writeToFile(TFile* outFile);

	private:
		string PFOCollectionName;
		EVENT::LCCollection *PFOCollection = nullptr;
		EVENT::LCCollection *trackCollection = nullptr;

		EVENT::MCParticle* genPart = nullptr;

		unsigned int nPFOs = 0;

		int fillTrackerTimingInfo();
		inline float getSimpleTimeOfFlight(float x, float y, float z);
		int getLayer(EVENT::TrackerHit* hit, UTIL::BitField64 &encoder);
		inline float getHitRadius(float x, float y);
};
#endif 



