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
#include "truthCondition.h"

class eventHistFiller : public objectFill{
	public:
		eventHistFiller(string _outDirName, unsigned int _pfoTypeToSelect) 
			: objectFill(_outDirName) 
		{flagDiscardConvertion = false; flagSelectConvertion = false; pfoTypeToSelect.push_back(_pfoTypeToSelect); thetaMergingCut = 0.6;}
		eventHistFiller(string _outDirName, vector<unsigned int> _pfoTypeToSelect) 
			: objectFill(_outDirName) 
		{flagDiscardConvertion = false; flagSelectConvertion = false; pfoTypeToSelect.insert(pfoTypeToSelect.begin(),_pfoTypeToSelect.begin(),_pfoTypeToSelect.end()); thetaMergingCut = 0.6;}
		~eventHistFiller(){}

		int init();
		int fillEvent(const EVENT::LCEvent*);
		int writeToFile(TFile* outFile);
		void setPFOCollection(const string _collectionName){PFOCollectionName = _collectionName;}
		void setMCTruthCollection(const string _collectionName){MCTruthCollectionName = _collectionName;}
		void setDiscardConvertions(bool inFlag){flagDiscardConvertion = inFlag;}
		void setSelectConvertions(bool inFlag){flagSelectConvertion = inFlag;}
		// void setPhotonReclustering(bool inFlag){photonReclustering = inFlag;}
		// void setNeutralReclustering(bool inFlag){neutralReclustering = inFlag;}
		void setMergePfoType(vector<unsigned int> inVec);
		void setMergePfoType(unsigned int inVal);
		void setThetaMergingCut(double inCut){thetaMergingCut = inCut;}

	private:
		void get90(const TH1 *const pTH1F);
		vector<unsigned int> pfoTypeToSelect;
		vector<unsigned int> pfoTypeToMerge; 
		bool flagDiscardConvertion;
		bool flagSelectConvertion;
		// bool photonReclustering;
		// bool neutralReclustering;
		unsigned int nSelectecTruthParticles;
		string PFOCollectionName;
		string MCTruthCollectionName;
		EVENT::LCCollection *PFOCollection;
		EVENT::LCCollection *MCTruthCollection;
		int checkPfoType(vector<unsigned int> inVec);
		double thetaMergingCut;
};
#endif // EVENTHISTFILLER_H
