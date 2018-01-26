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

vector<EVENT::ReconstructedParticle*> objectFill::getObjVecFromCollection(const EVENT::LCCollection* inCollection){
	int nElements = inCollection->getNumberOfElements();
	vector<EVENT::ReconstructedParticle*> outVec;
	for(int j=0; j < nElements; j++) {
		auto part = dynamic_cast<EVENT::ReconstructedParticle*>(inCollection->getElementAt(j));
		outVec.push_back(part);
	}
	return outVec;
}

// template <class T> vector<T> objectFill::getObjVecFromCollection(EVENT::LCCollection* inCollection){
//         int nElements = inCollection->getNumberOfElements();
//         vector<T> outVec;
//         for(int j=0; j < nElements; j++) {
//                 auto part = dynamic_cast<T>(inCollection->getElementAt(j));
//                 outVec.push_back(part);
//         }
//         return outVec;

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

TH1* objectFill::getHistFromMap(string histID){
	if (histMap[histID]==NULL)
		cout << "[ERROR]\tobjectFill::getHistFromMap(" << histID << ") no hist in the histMap with name <" << histID << ">" << endl;
	return histMap[histID];
}
