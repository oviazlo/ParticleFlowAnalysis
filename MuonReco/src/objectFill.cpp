#include <objectFill.h>
#include <TGraphAsymmErrors.h>

objectFill::objectFill(string _outDirName){
	outDirName = _outDirName;
}

objectFill::~objectFill(){
	for(auto const &it : histMap) {
		it.second->Delete();
	}
	// for(auto const &it : tEffMap) {
	//         it.second->~TEfficiency();
	// }
}

int objectFill::writeToFile(TFile* outFile){
	if (!outFile->IsOpen()){
		cout << "[ERROR|writeToFile]\tno output file is found!" << endl;
		return -1;
	}
	outFile->cd();
	TDirectory *mainDir = outFile->mkdir(outDirName.c_str());
	mainDir->cd();

	mainDir->mkdir("tEff");
	mainDir->cd("tEff");
	for(auto const &it : tEffMap){
		// cout << "tEff name: " << it.second->GetName();
		TGraphAsymmErrors *tmpGr = it.second->CreateGraph();
		tmpGr->SetName(it.second->GetName());
		tmpGr->Write();
	}
	mainDir->cd();

	map<string,unsigned int> prefixCounter;
	map<string,string> namePrefixMap;
	map<string,bool> isPrefixSubdirCreated;
	map<string,string> nameWithoutPrefixMap;
	for(auto const &it : histMap) {
		string histName = it.first;
		vector<string> tmpStrVec = GetSplittedWords(histName,"_");
		if (tmpStrVec.size()<2) 
			continue;
		string prefix = "";
		for (int i=0; i<tmpStrVec.size()-1; i++){
			if (i==tmpStrVec.size()-2)
				prefix += tmpStrVec[i];
			else
				prefix += tmpStrVec[i] + "_";
		}
		nameWithoutPrefixMap[histName] = tmpStrVec[tmpStrVec.size()-1];
		prefixCounter[prefix] += 1;
		isPrefixSubdirCreated[prefix] = false;
		namePrefixMap[histName] = prefix;
	}
	

	for(auto const &it : histMap) {
		string histName = it.first;
		string prefix = namePrefixMap[histName];
		if (prefixCounter[prefix]<2){
			mainDir->cd();
			it.second->Write();
		}
		else{
			if (isPrefixSubdirCreated[prefix]==false){
				mainDir->mkdir(prefix.c_str());
				isPrefixSubdirCreated[prefix]=true;
			}
			mainDir->cd(prefix.c_str());
			it.second->SetName(nameWithoutPrefixMap[histName].c_str());
			it.second->Write();
			mainDir->cd();
		}
	}
	outFile->cd();
	return 0;
	// if (!outFile->IsOpen()){
	//         cout << "[ERROR|writeToFile]\tno output file is found!" << endl;
	//         return -1;
	// }
	// outFile->cd();
	// TDirectory *dir = outFile->mkdir(outDirName.c_str());
	// dir->cd();
	// for(auto const &it : histMap) {
	//         it.second->Write();
	// }
	// outFile->cd();
	// return 0;
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
	for(auto &mapElement : histMap){
		delete mapElement.second;
		histMap.erase(mapElement.first);    
	}
	for(auto &mapElement : tEffMap){
		delete mapElement.second;
		histMap.erase(mapElement.first);    
	}
}

TH1* objectFill::getHistFromMap(string histID){
	if (histMap[histID]==NULL)
		cout << "[ERROR]\t" + className + "::getHistFromMap(" << histID << ") no hist in the histMap with name <" << histID << ">" << endl;
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

void objectFill::createTH1I(string histName, string histTitle, unsigned int nBins, double leftRange, double rightRange){
	// string finalHistName = outDirName+"-"+histName;
	string finalHistName = histName;
	delete gROOT->FindObject(finalHistName.c_str());
	TH1I* tmpHist = new TH1I(finalHistName.c_str(),histTitle.c_str(),nBins,leftRange,rightRange);
	tmpHist->SetDirectory(0);
	histMap[finalHistName] = tmpHist;
}
void objectFill::createTH1D(string histName, string histTitle, unsigned int nBins, double leftRange, double rightRange){
	string finalHistName = histName;
	delete gROOT->FindObject(finalHistName.c_str());
	TH1D* tmpHist = new TH1D(finalHistName.c_str(),histTitle.c_str(),nBins,leftRange,rightRange);
	tmpHist->SetDirectory(0);
	histMap[finalHistName] = tmpHist;
}
