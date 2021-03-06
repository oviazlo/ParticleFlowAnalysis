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

#include <cmath>
// #include <gsl/gsl_sf_gamma.h>
#include <vector>

// #include <marlin/Global.h>
// #include <GeometryUtil.h>

// #include <DD4hep/DD4hepUnits.h>
// #include <DD4hep/DetType.h>
// #include <DD4hep/DetectorSelector.h>
// #include <DD4hep/Detector.h>
// #include <DDRec/DetectorData.h>

#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Cluster.h"
#include <lcio.h>

// #include <EVENT/CalorimeterHit.h>
// #include "CalorimeterHitType.h"
//
// #include <ClusterShapes.h>

struct PFOMergeSettings{
	unsigned int pfoTypeToMerge;
	double thetaCone;
	double phiCone;
};

// struct mergedEnergyContainer{
//         double mergedEnergy;
//         mergedEnergyContainer(): mergedEnergy(0.0){}
//         mergedEnergyContainer(double partE): mergedEnergy(partE){}
//         vector<PFOMergeSettings> > PFOmergeMap;
//         bool tryToMergeParticle(EVENT::ReconstructedParticle* inPart){
//
//         }
// };

class eventHistFiller : public objectFill{
	public:
		eventHistFiller(string _outDirName, string _PFOCollectionName = "PandoraPFO") : objectFill(_outDirName) {PFOCollectionName = _PFOCollectionName; mergeTag = "nominal"; applyAngularMatching = true; applyEnergyMatching = true;}
		~eventHistFiller(){}

		int init();
		int fillEvent(const EVENT::LCEvent*);
		int writeToFile(TFile* outFile);

		void setClusterMerging(string _mergeTag);
		void SetApplyAngularMatching(bool _applyAngularMatching){applyAngularMatching = _applyAngularMatching;}
		void SetApplyEnergyMatching(bool _applyEnergyMatching){applyEnergyMatching = _applyEnergyMatching;}
		void SetUseCaloInfoForEnergyMerging(bool _useCaloInfoDuringEnergyMerging){useCaloInfoDuringEnergyMerging = _useCaloInfoDuringEnergyMerging;}
		void SetUseCaloCutInsteadMomentum(bool _useCaloCutInsteadMomentum){useCaloCutInsteadMomentum = _useCaloCutInsteadMomentum;}

	private:
		string PFOCollectionName;
		EVENT::LCCollection *PFOCollection = nullptr;
		int checkPfoType(vector<unsigned int> inVec);
		// map<string,histStruct> singleParticleHistStructMap;
		// vector<string> effType = {"nominal","photonMerge","photonAndNeutralMerge","photonAndNeutralLooseMerge","bremRecovery"};
		vector<string> effType = {"nominal","photonMerge","photonAndNeutralMerge","photonAndNeutralLooseMerge","photonMergeMomentumDep"};
		// vector<string> effType = {"nominal","photonMerge"};
		map<string,vector<PFOMergeSettings> > PFOmergeMap;
		string mergeTag;

		bool applyAngularMatching = true;
		bool applyEnergyMatching = true;
		bool useCaloInfoDuringEnergyMerging = false;
		bool useCaloCutInsteadMomentum = false;

		EVENT::MCParticle* genPart = nullptr;
		IMPL::ReconstructedParticleImpl* partCandidate = nullptr;
		unsigned int idOfPartCandidate = std::numeric_limits<unsigned int>::max();

		double truthTheta = std::numeric_limits<double>::max();
		double cosTruthTheta = std::numeric_limits<double>::max();
		double truthPhi = std::numeric_limits<double>::max();
		double truthPt = std::numeric_limits<double>::max();
		double truthEnergy = std::numeric_limits<double>::max();
		int truthType = std::numeric_limits<int>::max();

		map <string, unsigned int> pfoCounter;
		unsigned int nPFOs = 0;

		double totalRecoEnergy = 0.0;

		void fillOtherHists();
		void angularAndEnergyMatching();
		void fillEventsFailedSelection();
		vector <unsigned int> mergeClusters();
		void fillPfoCounterAndFillGeneralPfoInfo(unsigned int partId);
		void initTruthInfoAndFillIt();


		void createTH1I(string histName, string histTitle, unsigned int nBins, double leftRange, double rightRange){
			// string finalHistName = outDirName+"-"+histName;
			string finalHistName = histName;
			delete gROOT->FindObject(finalHistName.c_str());
			TH1I* tmpHist = new TH1I(finalHistName.c_str(),histTitle.c_str(),nBins,leftRange,rightRange);
			tmpHist->SetDirectory(0);
			histMap[finalHistName] = tmpHist;
		}
		void createTH1D(string histName, string histTitle, unsigned int nBins, double leftRange, double rightRange){
			string finalHistName = histName;
			delete gROOT->FindObject(finalHistName.c_str());
			TH1D* tmpHist = new TH1D(finalHistName.c_str(),histTitle.c_str(),nBins,leftRange,rightRange);
			tmpHist->SetDirectory(0);
			histMap[finalHistName] = tmpHist;
		}

		TH1* getHistFromMap(string _histID){
			string histID = _histID;
			if (histMap[histID]==NULL)
				cout << "[ERROR]\teventHistFiller::getHistFromMap(" << histID << ") no hist in the histMap with name <" << histID << ">" << endl;
			return histMap[histID];
		}

};
#endif // EVENTHISTFILLER_H



