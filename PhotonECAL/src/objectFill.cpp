#include <objectFill.h>

objectFill::objectFill(string _outDirName){
	outDirName = _outDirName;
}

objectFill::~objectFill(){
	for(auto const &it : histMap) {
		it.second->Delete();
	}
}

int objectFill::writeToFile(TFile* outFile){
	if (!outFile->IsOpen()){
		cout << "[ERROR|writeToFile]\tno output file is found!" << endl;
		return -1;
	}
	outFile->cd();
	TDirectory *dir = outFile->mkdir(outDirName.c_str());
	dir->cd();
	for(auto const &it : histMap) {
		it.second->Write();
	}
	outFile->cd();
	return 0;
}

// vector<EVENT::ReconstructedParticle*> objectFill::getObjVecFromCollection(const EVENT::LCCollection* inCollection){
//         int nElements = inCollection->getNumberOfElements();
//         vector<EVENT::ReconstructedParticle*> outVec;
//         for(int j=0; j < nElements; j++) {
//                 auto part = dynamic_cast<EVENT::ReconstructedParticle*>(inCollection->getElementAt(j));
//                 outVec.push_back(part);
//         }
//         return outVec;
// }

double objectFill::get_dPhi(double phi_reco, double phi_truth){
	double dPhi = phi_reco - phi_truth;
	if (fabs(dPhi)<=180.0)
		return dPhi;
	if (dPhi>0.0)
		return dPhi - 360.0;
	else
		return 360.0 + dPhi;
}

void objectFill::createHistsFromMap(const map<string,histStruct> inHistStructMap, const string prefix){
	for(auto const &ent1 : inHistStructMap){
		TH1* tmpHist;
		if (ent1.second.histType=="TH1D")
			tmpHist = new TH1D((prefix+ent1.first).c_str(),ent1.second.title.c_str(),ent1.second.nBins,ent1.second.xLow,ent1.second.xHigh);
		if (ent1.second.histType=="TH1I")
			tmpHist = new TH1I((prefix+ent1.first).c_str(),ent1.second.title.c_str(),ent1.second.nBins,ent1.second.xLow,ent1.second.xHigh);
		if (ent1.second.histType=="TH2D")
			tmpHist = new TH2D((prefix+ent1.first).c_str(),ent1.second.title.c_str(),ent1.second.nBins,ent1.second.xLow,ent1.second.xHigh,ent1.second.ynBins,ent1.second.yLow,ent1.second.yHigh);
		histMap[prefix+ent1.first] = tmpHist;
	}
	
}

void objectFill::DeleteHists(){
	for(auto const &mapElement : histMap)
		delete mapElement.second;
}

TH1* objectFill::getHistFromMap(string histID){
	if (histMap[histID]==NULL)
		cout << "[ERROR]\tobjectFill::getHistFromMap(" << histID << ") no hist in the histMap with name <" << histID << ">" << endl;
	return histMap[histID];
}

IMPL::ReconstructedParticleImpl* objectFill::CopyReconstructedParticle (const EVENT::ReconstructedParticle* const pfo_orig ) {
	// copy this in an ugly fashion to be modifiable - a versatile copy constructor would be much better!
	IMPL::ReconstructedParticleImpl* pfo = new IMPL::ReconstructedParticleImpl();
	pfo->setMomentum(pfo_orig->getMomentum());
	pfo->setEnergy(pfo_orig->getEnergy());
	pfo->setType(pfo_orig->getType());
	pfo->setCovMatrix(pfo_orig->getCovMatrix());
	pfo->setMass(pfo_orig->getMass());
	pfo->setCharge(pfo_orig->getCharge());
	pfo->setParticleIDUsed(pfo_orig->getParticleIDUsed());
	pfo->setGoodnessOfPID(pfo_orig->getGoodnessOfPID());
	pfo->setStartVertex(pfo_orig->getStartVertex());
	return pfo;
}
