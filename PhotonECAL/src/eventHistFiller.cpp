/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/src/eventHistFiller.cpp
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	05th Dec 2017
 * 	Last Update:	05th Dec 2017
 */

#include "eventHistFiller.h"

/*===========================================================================*/
/*===============================[ function implementations ]================*/
/*===========================================================================*/

int eventHistFiller::init(){

	cout << "[INFO]\teventHistFiller::init()" << endl;

	TH1* tmpHist;
	tmpHist = new TH1I("nPFOs","Number of PFOs in event; Number of PFOs; Counts",5,0,5);
	histMap["nPFOs"] = tmpHist;
	tmpHist= new TH1I("PFOType","PFO particle type; Type; Counts",2200,0,2200); // max part.type = 2112 (neutron)
	histMap["PFOType"] = tmpHist;
	
	tmpHist = new TH1D("nPFOsVsTheta","nPFOs vs Theta; Theta; Counts per Event",180*2,0,180);
	histMap["nPFOsVsTheta"] = tmpHist;
	tmpHist = new TH1D("nPFOsVsCosTheta","nPFOs vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["nPFOsVsCosTheta"] = tmpHist;

	tmpHist = new TH1D("totalEnergyVsTheta","Sum of PFOs Energy vs Theta; Theta; Energy [GeV]",180*2,0,180);
	histMap["totalEnergyVsTheta"] = tmpHist;
	tmpHist = new TH1D("matchedEnergyVsTheta","Sum of matched PFOs Energy vs Theta; Theta; Energy [GeV]",180*2,0,180);
	histMap["matchedEnergyVsTheta"] = tmpHist;


	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
		tmpHist = new TH1D(("n"+it->second+"sVsCosTheta").c_str(),("n"+it->second+"s vs Theta; Theta; Counts per Event").c_str(),180*2,0,180);
		histMap["n"+it->second+"sVsTheta"] = tmpHist;
		tmpHist = new TH1D(("n"+it->second+"sVsCosTheta").c_str(),("n"+it->second+"s vs Cos(#Theta); Cos(#Theta); Counts per Event").c_str(),180*2,-1,1);
		histMap["n"+it->second+"sVsCosTheta"] = tmpHist;

	}
	
	tmpHist = new TH1D("nTruthPartsVsTheta","nTruthParts vs Theta; Theta; Counts per Event",180*2,0,180);
	histMap["nTruthPartsVsTheta"] = tmpHist;
	tmpHist = new TH1D("nTruthPartsVsCosTheta","nTruthParts vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["nTruthPartsVsCosTheta"] = tmpHist;
	tmpHist = new TH1D("nTruthPartsVsEnergy","nTruthParts vs Energy ;Energy [GeV]; Counts per Event",100,0.5,100.5); 
	histMap["nTruthPartsVsEnergy"] = tmpHist;

	tmpHist = new TH1D("efficiencyVsTheta","efficiency vs Theta; Theta; Counts per Event",180*2,0,180);
	histMap["efficiencyVsTheta"] = tmpHist;
	tmpHist = new TH1D("efficiencyVsCosTheta","efficiency vs Cos(#Theta); Cos(#Theta); Counts per Event",180*2,-1,1);
	histMap["efficiencyVsCosTheta"] = tmpHist;
	tmpHist = new TH1D("efficiencyVsEnergy","efficiency vs Energy; Energy [GeV]; Counts per Event",100,0.5,100.5);
	histMap["efficiencyVsEnergy"] = tmpHist;

	// cout << "debug1" << endl;

	for (auto it = pfoTypeToSelect.begin(); it != pfoTypeToSelect.end(); it++){
		tmpHist = new TH1D(("phiResolution_"+config::pfoTypeIntStringMap[*it]).c_str(),(config::pfoTypeIntStringMap[*it]+" Phi resolution; dPhi [rad]; Counts").c_str(),20000,-0.2,0.2);
		histMap["phiResolution_"+config::pfoTypeIntStringMap[*it]] = tmpHist;
		tmpHist = new TH1D(("thetaResolution_"+config::pfoTypeIntStringMap[*it]).c_str(),(config::pfoTypeIntStringMap[*it]+" Theta resolution; dTheta [rad]; Counts").c_str(),400,-0.01,0.01);
		histMap["thetaResolution_"+config::pfoTypeIntStringMap[*it]] = tmpHist;
		tmpHist = new TH1D(("energyResolution_"+config::pfoTypeIntStringMap[*it]).c_str(),(config::pfoTypeIntStringMap[*it]+" Energy resolution; E [GeV]; Counts").c_str(),1250,0,125);
		histMap["energyResolution_"+config::pfoTypeIntStringMap[*it]] = tmpHist;
	}

	// cout << "debug2" << endl;

	tmpHist = new TH1I("truthParticle_isDecayedInTracker","Sim. Status; Sim. Status; Counts",35,0,35);
	histMap["truthParticle_isDecayedInTracker"] = tmpHist;

	nSelectecTruthParticles = 0;
	return 0;

}

int eventHistFiller::fillEvent(const EVENT::LCEvent* event){
	// read collections
	try {
		PFOCollection = event->getCollection(PFOCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|eventHistFiller]\tCan't find collection: " << PFOCollectionName << endl;
		return -1;
	}

	// find primary generated MC particle which is the only findable in the event by definition [particle gun]
	

	double truthTheta = -666;  
	double truthPhi = -666;  
	double cosTruthTheta = -666; 
	double truthEnergy = -666;
	bool reclusteringIsDone = false;

	EVENT::MCParticle* genPart = truthCondition::instance()->getGunParticle();
	const double *partMom = genPart->getMomentum();
	TVector3 v1;
	v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
	double partTheta = 180.*v1.Theta()/TMath::Pi();
	double cosPartTheta = TMath::Cos(TMath::Pi()*partTheta/180.);
	truthTheta = partTheta;
	truthPhi = 180.*v1.Phi()/TMath::Pi();
	cosTruthTheta = cosPartTheta;
	if (partTheta<8 || partTheta>172) 
		return 0;
	if (genPart->isDecayedInTracker() && flagDiscardConvertion==true)
		return 0;
	if ((!genPart->isDecayedInTracker()) && flagSelectConvertion==true)
		return 0;
	histMap["truthParticle_isDecayedInTracker"]->Fill(genPart->isDecayedInTracker());
	nSelectecTruthParticles++;
	histMap["nTruthPartsVsCosTheta"]->Fill(cosPartTheta);
	histMap["nTruthPartsVsTheta"]->Fill(partTheta);
	truthEnergy = genPart->getEnergy();
	histMap["nTruthPartsVsEnergy"]->Fill(truthEnergy);

	vector<EVENT::ReconstructedParticle*> recoPFOs = getObjVecFromCollection(PFOCollection);
	histMap["nPFOs"]->Fill(recoPFOs.size());
	if (recoPFOs.size()==0) 
		return 0; // no reco PFOs

	for (int i=0; i<recoPFOs.size(); i++){
		int pfoType = abs(recoPFOs[i]->getType());
		const double *partMom = recoPFOs[i]->getMomentum();
		TVector3 v1;
		v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
		double partTheta = 180.*v1.Theta()/TMath::Pi();
		double partPhi = 180.*v1.Phi()/TMath::Pi();
		double cosPartTheta = TMath::Cos(TMath::Pi()*partTheta/180.);

		// cout << "partTheta: " << partTheta << "; truthTheta: " << truthTheta << endl;

		// WARNING TODO FIXME: wrong when calculate dTheta!!!
		// partTheta = truthTheta;
		// cosPartTheta = cosTruthTheta;

		histMap["PFOType"]->Fill(pfoType);
		histMap["nPFOsVsCosTheta"]->Fill(cosTruthTheta);
		histMap["nPFOsVsTheta"]->Fill(truthTheta);

		histMap["totalEnergyVsTheta"]->Fill(truthTheta,recoPFOs[i]->getEnergy());

		for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
			if (abs(pfoType)==it->first) {
				histMap["n"+it->second+"sVsCosTheta"]->Fill(cosTruthTheta);
				histMap["n"+it->second+"sVsTheta"]->Fill(truthTheta);
			}
		}

		// for (auto it = pfoTypeToMerge.begin(); it != pfoTypeToMerge.end(); it++) {
		//         double dPhi = TMath::Pi()*(partPhi-truthPhi)/180.0;
		//         double dTheta = TMath::Pi()*(partTheta-truthTheta)/180.0;
		//         histMap["thetaResolution_"+config::pfoTypeIntStringMap[*it]]->Fill(dTheta);
		//         histMap["phiResolution_"+config::pfoTypeIntStringMap[*it]]->Fill(dPhi);
		// }
	
		// cout << "debug3" << endl;

		if (IsInVector<unsigned int>(abs(pfoType),pfoTypeToSelect) && reclusteringIsDone==false) {
			double pfoE = recoPFOs[i]->getEnergy();
			reclusteringIsDone = true;

			// cout << "debug31" << endl;
			if (pfoTypeToMerge.size()>0) {
				for (int j=0; j<recoPFOs.size(); j++){
					if (i==j) 
						continue;
					if (IsInVector<unsigned int>(recoPFOs[j]->getType(),pfoTypeToMerge)==false)
						continue;
					const double *partMom = recoPFOs[i]->getMomentum();
					TVector3 v1;
					v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
					double partTheta_i = 180.*v1.Theta()/TMath::Pi();
					double partPhi_i = 180.*v1.Phi()/TMath::Pi();

					partMom = recoPFOs[j]->getMomentum();
					v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
					double partTheta_j = 180.*v1.Theta()/TMath::Pi();
					double partPhi_j = 180.*v1.Phi()/TMath::Pi();
					// cout << "partTheta_i-partTheta_j: " << partTheta_i-partTheta_j << ", partPhi_i-partPhi_j: " << partPhi_i-partPhi_j << endl;
					if (abs(partTheta_i-partTheta_j)<thetaMergingCut && abs(partPhi_i-partPhi_j)<10.0) {
					// if (abs(partTheta_i-partTheta_j)<2.0 && abs(partPhi_i-partPhi_j)<10.0) {
						pfoE+=recoPFOs[j]->getEnergy();
					}
				}
			}

			double dPhi = TMath::Pi()*(partPhi-truthPhi)/180.0;
			double dTheta = TMath::Pi()*(partTheta-truthTheta)/180.0;

			if (abs(pfoE-truthEnergy)<0.75*sqrt(truthEnergy)) {
				if (recoPFOs.size()==1 && (abs(dPhi)>0.004 || abs(dTheta)>0.002)) {
					continue;
				}
				histMap["efficiencyVsTheta"]->Fill(truthTheta);
				histMap["efficiencyVsCosTheta"]->Fill(cosTruthTheta);
				histMap["efficiencyVsEnergy"]->Fill(truthEnergy);
				histMap["matchedEnergyVsTheta"]->Fill(truthTheta,pfoE);
			}
			
			histMap["thetaResolution_"+config::pfoTypeIntStringMap[pfoType]]->Fill(dTheta);
			histMap["phiResolution_"+config::pfoTypeIntStringMap[pfoType]]->Fill(dPhi);
			histMap["energyResolution_"+config::pfoTypeIntStringMap[pfoType]]->Fill(pfoE);

		}
	}
}


int eventHistFiller::writeToFile(TFile* outFile){

	for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
		histMap["n"+it->second+"sVsCosTheta"]->Divide(histMap["nTruthPartsVsCosTheta"]);
		histMap["n"+it->second+"sVsTheta"]->Divide(histMap["nTruthPartsVsTheta"]);
	}
	histMap["nPFOsVsTheta"]->Divide(histMap["nTruthPartsVsTheta"]);
	histMap["nPFOsVsCosTheta"]->Divide(histMap["nTruthPartsVsCosTheta"]);
	
	histMap["efficiencyVsCosTheta"]->Divide(histMap["nTruthPartsVsCosTheta"]);
	histMap["efficiencyVsTheta"]->Divide(histMap["nTruthPartsVsTheta"]);
	histMap["efficiencyVsEnergy"]->Divide(histMap["nTruthPartsVsEnergy"]);
	histMap["totalEnergyVsTheta"]->Divide(histMap["nTruthPartsVsTheta"]);
	histMap["matchedEnergyVsTheta"]->Divide(histMap["nTruthPartsVsTheta"]);

	// for (auto it = pfoTypeToMerge.begin(); it != pfoTypeToMerge.end(); it++) {
	//         get90(histMap["thetaResolution_"+config::pfoTypeIntStringMap[*it]]);
	//         get90(histMap["phiResolution_"+config::pfoTypeIntStringMap[*it]]);
	// }

	objectFill::writeToFile(outFile);
	return 0;

}

void eventHistFiller::setMergePfoType(vector<unsigned int> inVec){
	pfoTypeToMerge.clear();
	pfoTypeToMerge.insert(pfoTypeToMerge.begin(),inVec.begin(),inVec.end());
	checkPfoType(pfoTypeToMerge);
}
void eventHistFiller::setMergePfoType(unsigned int inVal){
	pfoTypeToMerge.clear();
	pfoTypeToMerge.push_back(inVal);
	checkPfoType(pfoTypeToMerge);
}
int eventHistFiller::checkPfoType(vector<unsigned int> inVec){
	for (auto it2 = inVec.begin(); it2 != inVec.end(); it2++) {
		bool isAllowed = false;
		for (auto it = config::pfoTypeIntStringMap.begin(); it != config::pfoTypeIntStringMap.end(); it++) {
			if (it->first==(*it2))
				isAllowed = true;
		}
		if (isAllowed == false) {
			cout << "[ERROR]\teventHistFiller::checkPfoType\tnot allowed PFO type" << endl;
			exit(1);
		}
	}
}
/*===========================================================================*/
/*===============================[ function implementations ]===============================*/
/*===========================================================================*/

void eventHistFiller::get90(const TH1 *const pTH1F){
    static const float FLOAT_MAX(std::numeric_limits<float>::max());

    double fractionEventsToIgnorePerSide = 0.05;

    if (NULL == pTH1F)
        return;

    if (5 > pTH1F->GetEntries())
    {
        std::cout << pTH1F->GetName() << " (" << pTH1F->GetEntries() << " entries) - skipped" << std::endl;
        return;
    }

    // Calculate raw properties of distribution
    float sum = 0., total = 0.;
    double sx = 0., sxx = 0.;
    const unsigned int nbins(pTH1F->GetNbinsX());

    for (unsigned int i = 0; i <= nbins; ++i)
    {
        const float binx(pTH1F->GetBinLowEdge(i) + (0.5 * pTH1F->GetBinWidth(i)));
        const float yi(pTH1F->GetBinContent(i));
        sx += yi * binx;
        sxx += yi * binx * binx;
        total += yi;
    }

    const float rawMean(sx / total);
    const float rawMeanSquared(sxx / total);
    const float rawRms(std::sqrt(rawMeanSquared - rawMean * rawMean));

    sum = 0.;
    unsigned int is0 = 0;

    for (unsigned int i = 0; (i <= nbins) && (sum < total * fractionEventsToIgnorePerSide ); ++i)
    {
        sum += pTH1F->GetBinContent(i);
        is0 = i;
    }

    // Calculate truncated properties
    float rmsmin(FLOAT_MAX), sigma(FLOAT_MAX), sigmasigma(FLOAT_MAX), frac(FLOAT_MAX), efrac(FLOAT_MAX), mean(FLOAT_MAX), low(FLOAT_MAX), rms(FLOAT_MAX);
    float high(0.f);

    for (unsigned int istart = 0; istart <= is0; ++istart)
    {
        double sumn = 0.;
        double csum = 0.;
        double sumx = 0.;
        double sumxx = 0.;
        unsigned int iend = 0;

        for (unsigned int i = istart; (i <= nbins) && (csum < (1.0 - fractionEventsToIgnorePerSide) * total); ++i)
        {
            const float binx(pTH1F->GetBinLowEdge(i) + (0.5 * pTH1F->GetBinWidth(i)));
            const float yi(pTH1F->GetBinContent(i));
            csum += yi;

            if (sumn < (1.0 - fractionEventsToIgnorePerSide) * total)
            {
                sumn += yi;
                sumx += yi * binx;
                sumxx+= yi * binx * binx;
                iend = i;
            }
        }

        const float localMean(sumx / sumn);
        const float localMeanSquared(sumxx / sumn);
        const float localRms(std::sqrt(localMeanSquared - localMean * localMean));

        if (localRms < rmsmin)
        {
            mean = localMean;
            rms = localRms;
            low = pTH1F->GetBinLowEdge(istart);
            high = pTH1F->GetBinLowEdge(iend);
            rmsmin = localRms;

            
	    sigma = rms;
	    sigmasigma = sigma / std::sqrt(total);
        }
    }

    std::cout << pTH1F->GetName() << " (" << pTH1F->GetEntries() << " entries), rawrms: " << rawRms << ", rms90: " << rmsmin
    	  << " (" << low << "-" << high << "), mean: " << mean << ", sigma: " << sigma << "+-" << sigmasigma << endl << endl;

}
		

