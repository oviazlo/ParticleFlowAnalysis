/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/include/eventHistFiller.h
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	05th Dec 2017
 * 	Last Update:	05th Dec 2017
 */
#ifndef EVENTHISTFILLER_H
#define EVENTHISTFILLER_H

#include <objectFill.h>
#include <globalConfig.h>

class eventHistFiller : public objectFill{
	public:
		eventHistFiller(string _outDirName) : objectFill(_outDirName) {flagDiscardConvertion = false; flagSelectConvertion = false; photonReclustering = false; neutralReclustering = false;}
		~eventHistFiller(){}

		int init();
		int fillEvent(const EVENT::LCEvent*);
		int writeToFile(TFile* outFile);
		void setPFOCollection(const string _collectionName){PFOCollectionName = _collectionName;}
		void setMCTruthCollection(const string _collectionName){MCTruthCollectionName = _collectionName;}
		void setDiscardConvertions(bool inFlag){flagDiscardConvertion = inFlag;}
		void setSelectConvertions(bool inFlag){flagSelectConvertion = inFlag;}
		void setPhotonReclustering(bool inFlag){photonReclustering = inFlag;}
		void setNeutralReclustering(bool inFlag){neutralReclustering = inFlag;}

	private:
		bool flagDiscardConvertion;
		bool flagSelectConvertion;
		bool photonReclustering;
		bool neutralReclustering;
		unsigned int nSelectecTruthParticles;
		string PFOCollectionName;
		string MCTruthCollectionName;
		EVENT::LCCollection *PFOCollection;
		EVENT::LCCollection *MCTruthCollection;
};
#endif // EVENTHISTFILLER_H
