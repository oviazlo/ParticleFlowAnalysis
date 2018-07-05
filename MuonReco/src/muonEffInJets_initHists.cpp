#include "muonEffInJets.h"

int muonEffInJets::init(){

	className = "muonEffInJets";
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

	vector<string> truthPartDirPrefix = {"truthPartAll","truthPartMatched","truthPartNotMatched","eff_truthPartMatched"};
	for (auto iDir: truthPartDirPrefix){
		createTH1D(iDir+"_E","Truth Energy; Energy [GeV]; Counts",50,0,500);
		createTH1D(iDir+"_pt","Truth Pt; pT [GeV]; Counts",15,0,150);
		createTH1D(iDir+"_theta","Truth Theta; Theta [#circ]; Counts",45,0,180);
		createTH1D(iDir+"_cosTheta","Truth cos(Theta); cos(Theta); Counts",50,-1,1);
		createTH1D(iDir+"_phi","Truth Phi; Phi [#circ]; Counts",45,-180,180);
		createTH1D(iDir+"_vertexR","Truth Vertex R; R [mm]; Counts",100,0,50);
		createTH1D(iDir+"_vertexZ","Truth Vertex Z; Z [mm]; Counts",100,0,50);
	}

	return 0;
}


int muonEffInJets::writeToFile(TFile* outFile){

	if (config::vm.count("debug"))
		cout << "[INFO]"+className+"::writeToFile(" << outFile->GetName() << ")" << endl;

	string histNamePrefix = "eff_truthPartMatched";
	string histNameDenominator = "truthPartAll";

	getHistFromMap(histNamePrefix + "_E")->Sumw2();
	getHistFromMap(histNamePrefix + "_pt")->Sumw2();
	getHistFromMap(histNamePrefix + "_theta")->Sumw2();
	getHistFromMap(histNamePrefix + "_cosTheta")->Sumw2();
	getHistFromMap(histNamePrefix + "_phi")->Sumw2();
	getHistFromMap(histNamePrefix + "_vertexR")->Sumw2();
	getHistFromMap(histNamePrefix + "_vertexZ")->Sumw2();

	createTEff(histNamePrefix + "_E",histNameDenominator + "_E");
	createTEff(histNamePrefix + "_pt",histNameDenominator + "_pt");
	createTEff(histNamePrefix + "_theta",histNameDenominator + "_theta");
	createTEff(histNamePrefix + "_cosTheta",histNameDenominator + "_cosTheta");
	createTEff(histNamePrefix + "_phi",histNameDenominator + "_phi");
	createTEff(histNamePrefix + "_vertexR",histNameDenominator + "_vertexR");
	createTEff(histNamePrefix + "_vertexZ",histNameDenominator + "_vertexZ");

	getHistFromMap(histNamePrefix + "_E")->Divide(getHistFromMap(histNameDenominator + "_E"));
	getHistFromMap(histNamePrefix + "_pt")->Divide(getHistFromMap(histNameDenominator + "_pt"));
	getHistFromMap(histNamePrefix + "_theta")->Divide(getHistFromMap(histNameDenominator + "_theta"));
	getHistFromMap(histNamePrefix + "_cosTheta")->Divide(getHistFromMap(histNameDenominator + "_cosTheta"));
	getHistFromMap(histNamePrefix + "_phi")->Divide(getHistFromMap(histNameDenominator + "_phi"));
	getHistFromMap(histNamePrefix + "_vertexR")->Divide(getHistFromMap(histNameDenominator + "_vertexR"));
	getHistFromMap(histNamePrefix + "_vertexZ")->Divide(getHistFromMap(histNameDenominator + "_vertexZ"));

	objectFill::writeToFile(outFile);
	return 0;


}

