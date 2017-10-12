#ifndef energyFill_H
#define energyFill_H

#include <objectFill.h>

#include <UTIL/CellIDEncoder.h>

class energyFill : public objectFill{

	public:
		energyFill(string _outDirName) : objectFill(_outDirName) {}
		~energyFill(){}
		int init();
		int fillEvent(EVENT::LCEvent*);
		void setCollectionName(string _collectionName){collectionName = _collectionName;}
		int writeToFile(TFile* outFile);

	private:
		int fillEnergy(double energy);
		int fillNHits(int nHits);
		int fillMaxLayer(int maxLayer);
		int createHists();
		EVENT::LCCollection *collection;
		string collectionName;
		
		//Initialize CellID encoder
		// UTIL::BitField64 m_encoder;
};


#endif
