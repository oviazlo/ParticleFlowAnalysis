#include "trackerTiming.h"
#include <UTIL/LCTrackerConf.h>

int trackerTiming::fillEvent(const EVENT::LCEvent* event){
	try {
		PFOCollection = event->getCollection(PFOCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|trackerTiming]\tCan't find collection: " << PFOCollectionName << endl;
		return -1;
	
	}

	string trackCollectionName = "SiTracks_Refitted";
	// trackCollectionName = "SiTracks";
	try {
		trackCollection = event->getCollection(trackCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|electronStudy]\tCan't find collection: " << trackCollectionName << endl;
		return -1;
	}

	nPFOs = PFOCollection->getNumberOfElements();

	fillTrackerTimingInfo();
	return 0;

}


int trackerTiming::fillTrackerTimingInfo(){
	if (nPFOs!=1)
		return 1;

	UTIL::BitField64 encoder( lcio::LCTrackerCellID::encoding_string() ) ;
	encoder.reset();

	vector<EVENT::Track*> tracks = getObjVecFromCollection<EVENT::Track*>(trackCollection);
	getHistFromMap("nTracks")->Fill(tracks.size());

	if (tracks.size()!=1)
		return 1;
	auto track = tracks[0];
	auto trackHits = track->getTrackerHits();
	// cout << "\nnHits: " << trackHits.size() << endl;
	for (auto iHit: trackHits){
		const double *hitPos = iHit->getPosition();
		double hitRadius = getHitRadius(hitPos[0],hitPos[1]);
		unsigned int layer = getLayer(iHit,encoder);
		// cout << "radius: " << hitRadius << "; layer: " << layer << endl;
		if (hitRadius>2100)
			getHistFromMap("lastBarrelLayerTiming")->Fill(iHit->getTime());
	}
}


inline float trackerTiming::getSimpleTimeOfFlight(float x, float y, float z){
	return std::sqrt((x * x) + (y * y) + (z * z))/299.792458;
}

inline float trackerTiming::getHitRadius(float x, float y){
	return std::sqrt((x * x) + (y * y));
}

int trackerTiming::getLayer(EVENT::TrackerHit* hit, UTIL::BitField64 &encoder){
  const int celId = hit->getCellID0() ;
  encoder.setValue(celId) ;
  int layer = encoder[lcio::LCTrackerCellID::layer()].value();
  return layer;
}
