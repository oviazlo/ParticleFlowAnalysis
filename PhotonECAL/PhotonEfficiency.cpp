//custom libs
#include <particleFill.h>
#include <energyFill.h>
#include <photonEffCalculator.h>

using namespace std;
int main (int argn, char* argv[]) {
	// Read collections
	std::vector<std::string> m_fileNames;

	if (argn < 2) {
		std::cout << "[WARNING]\tNo input arguments. Exit!" << std::endl;
		return 0;
	}
	else{
		for (int i = 1; i < argn; ++i) {
			vector<string> tmpStrVec = getFilesMatchingPattern(argv[i]);
			m_fileNames.insert(m_fileNames.end(), tmpStrVec.begin(), tmpStrVec.end());
		}
	}

	cout << "[INFO]\tNumber of input files to be used: " << m_fileNames.size() << " files" << endl;

	bool printVerbose = false;

	// Open Files
	auto m_reader( IOIMPL::LCFactory::getInstance()->createLCReader());

	try{
		m_reader->open( m_fileNames );
	} catch (IO::IOException &e)  {
		std::cerr << "Error opening files: " << e.what() << std::endl;
		return 1;
	}

	map<string, objectFill*> objFillMap;
	vector<string> particleFillCollections = {"MCParticlesSkimmed","LooseSelectedPandoraPFOs","TightSelectedPandoraPFOs","SelectedPandoraPFOs"};
	vector<string> energyFillCollections = {"ECALBarrel","ECALEndcap"/*, "ECalBarrelCollection", "ECalEndcapCollection"*/};
	for (int i; i<particleFillCollections.size(); i++){
		particleFill* tmpPartFill = new particleFill(particleFillCollections[i]);
		tmpPartFill->setCollectionName(particleFillCollections[i]);
		objFillMap[particleFillCollections[i]] = tmpPartFill;
	}
	for (int i; i<energyFillCollections.size(); i++){
		energyFill* tmpEnergyFill = new energyFill(energyFillCollections[i]);
		tmpEnergyFill->setCollectionName(energyFillCollections[i]); 
		objFillMap[energyFillCollections[i]] = tmpEnergyFill;
	}

	photonEffCalculator* effCalculator = new photonEffCalculator("photonEfficiency");
	effCalculator->setPFOCollection("LooseSelectedPandoraPFOs");
	effCalculator->setMCTruthCollection("MCParticlesSkimmed");
	objFillMap["photonEfficiency"] = effCalculator;

	for(auto const &mapElement : objFillMap){
		mapElement.second->init();
	}

	// static_cast<particleFill*>(objFillMap["MCParticlesSkimmed"])->setDebugFlag(true);
	// static_cast<particleFill*>(objFillMap["SelectedPandoraPFOs"])->setDebugFlag(true);

	EVENT::LCEvent *event = NULL;
	int eventCounter = 0;
	while ( ( event = m_reader->readNextEvent() ) ) {
		// if (eventCounter>100) break;
		if (printVerbose) cout << endl << "Event " << eventCounter << ":" << endl;
		eventCounter++;
		
		for(auto const &mapElement : objFillMap){
			mapElement.second->fillEvent(event);
		} 
	}

	double scale = 0.1; 
	string meanEnergy = DoubToStr( floor(static_cast<particleFill*>(objFillMap["MCParticlesSkimmed"])->getMeanEnergy() / scale + 0.5)*scale );
	string meanTheta = DoubToStr( floor(static_cast<particleFill*>(objFillMap["MCParticlesSkimmed"])->getMeanTheta() / scale + 0.5)*scale );
	TFile *outFile = new TFile(("ECAL_photonGun_E"+meanEnergy+"_theta"+meanTheta+".root").c_str(), "RECREATE");

	for(auto const &mapElement : objFillMap){
		mapElement.second->writeToFile(outFile);
	} 

	outFile->Close();

	for(auto const &mapElement : objFillMap){
		delete mapElement.second;
	} 
	
}

