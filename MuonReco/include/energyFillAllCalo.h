#ifndef energyFillAllCalo_H
#define energyFillAllCalo_H

#include <objectFill.h>

#include <UTIL/CellIDEncoder.h>

class energyFillAllCalo : public objectFill{

	public:
		energyFillAllCalo(string _outDirName) : objectFill(_outDirName) {}
		~energyFillAllCalo(){}
		int init();
		int fillEvent(const EVENT::LCEvent*);
		// void setCollectionName(const string _collectionName){collectionName = _collectionName;}
		int writeToFile(TFile* outFile);

	private:
		EVENT::LCCollection *collection;
		string collectionName;
};


#endif
