#include "muonClusterFiller.h"

int muonClusterFiller::init(){

	className = "muonClusterFiller";
	if (config::vm.count("debug")){
		cout << "[INFO]"+className+"::init(); outDirName: " << outDirName << endl;
	}

	createTH1D("cluster_xt90","radius where 90\% of the cluster energy exists; R [mm]; Counts",5000,0,500);
	createTH1D("cluster_depth","depth of the cluster; L [mm]; Counts",5000,0,5000);
	createTH1D("cluster_RhitMean","mean of the radius of the hits wrt cog; <R_{hit}> [mm]; Counts",5000,0,5000);
	createTH1D("cluster_RhitRMS","RMS of the radius of the hits wrt cog; RMS(R_{hit}) [mm]; Counts",5000,0,5000);
	createTH1D("cluster_nYokeHits","Number of yoke hits in cluster; nHits; Counts",100,0,100);
	createTH1D("cluster_nLayers","Number of yoke layers in cluster; nLayers; Counts",10,0,10);
	createTH1D("cluster_clusterLayerSpan","Number of yoke hits in cluster; nLayers; Counts",10,0,10);

	return 0;
}


int muonClusterFiller::writeToFile(TFile* outFile){

	if (config::vm.count("debug"))
		cout << "[INFO]"+className+"::writeToFile(" << outFile->GetName() << ")" << endl;

	objectFill::writeToFile(outFile);
	return 0;


}

