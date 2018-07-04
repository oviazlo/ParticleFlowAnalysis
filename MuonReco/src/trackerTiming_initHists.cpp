#include "trackerTiming.h"

int trackerTiming::init(){

	className = "trackerTiming";
	if (config::vm.count("debug")){
		cout << "[INFO]"+className+"::init(); outDirName: " << outDirName << endl;
	}

	createTH1D("lastBarrelLayerTiming","Time information from the outermost OT layer; Time [ns]; Counts",10000,0,10);
	createTH1I("nTracks","Number of tracks per event, nPFOs==1; nTracks; Counts",10,0,10);

	return 0;
}


int trackerTiming::writeToFile(TFile* outFile){

	if (config::vm.count("debug"))
		cout << "[INFO]"+className+"::writeToFile(" << outFile->GetName() << ")" << endl;

	objectFill::writeToFile(outFile);
	return 0;


}

