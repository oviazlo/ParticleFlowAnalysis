#ifndef energyFill_H
#define energyFill_H

#include <objectFill.h>

#include <UTIL/CellIDEncoder.h>

class energyFill : public objectFill{

	public:
		energyFill(string _outDirName) : objectFill(_outDirName) {}
		~energyFill(){}
		int init();
		int fillEvent(const EVENT::LCEvent*);
		void setCollectionName(const string _collectionName){collectionName = _collectionName;}
		int writeToFile(TFile* outFile);

	private:
		int fillEnergy(const double energy);
		int fillNHits(const int nHits);
		int fillMaxLayer(const int maxLayer);
		int createHists();
		EVENT::LCCollection *collection;
		string collectionName;
		
		//Initialize CellID encoder
		// UTIL::BitField64 m_encoder;
};


#endif
