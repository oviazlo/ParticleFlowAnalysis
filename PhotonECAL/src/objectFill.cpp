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

vector<EVENT::ReconstructedParticle*> objectFill::getObjVecFromCollection(EVENT::LCCollection* inCollection){
	int nElements = inCollection->getNumberOfElements();
	vector<EVENT::ReconstructedParticle*> outVec;
	for(int j=0; j < nElements; j++) {
		auto part = static_cast<EVENT::ReconstructedParticle*>(inCollection->getElementAt(j));
		outVec.push_back(part);
	}
	return outVec;
}

// template <class T> vector<T> objectFill::getObjVecFromCollection(EVENT::LCCollection* inCollection){
//         int nElements = inCollection->getNumberOfElements();
//         vector<T> outVec;
//         for(int j=0; j < nElements; j++) {
//                 auto part = static_cast<T>(inCollection->getElementAt(j));
//                 outVec.push_back(part);
//         }
//         return outVec;
// }
