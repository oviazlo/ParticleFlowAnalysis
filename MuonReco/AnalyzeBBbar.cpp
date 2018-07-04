//custom libs
#include <muonEffInJets.h>
#include <trackerTiming.h>
#include <boostServiceFunctions.h>
#include <serviceFunctions.h>
#include "truthZWCondition.h"
// #include <globalConfig.h>

using namespace std;
using namespace config; 
 
// COLLECTIONS TO USE
vector<string> particleFillCollections = {"MCParticle","PandoraPFOs"};
vector<string> energyFillCollections = {};
vector<string> additionalCollections = {"SiTracks_Refitted","LCRelation","RecoMCTruthLink","PandoraClusters"};

int main (int argn, char* argv[]) {

	map<string, objectFill*> objFillMap;
	objectFill* objFill = NULL;

	objFill = new muonEffInJets("muonEffInJets", "PandoraPFOs");
	objFillMap["muonEffInJets"] = objFill;

	for(auto const &mapElement : objFillMap){
		cout << "Init mapElement: " << mapElement.first << endl;
		mapElement.second->init();
	}

	po::options_description desc("Options");
	desc.add_options()
		("help,h", "Example:\n ./MuonReconMisID -f \"/ssd/viazlo/data/FCCee_o5_v04_ILCSoft-2017-07-27_gcc62_photons_cosTheta_v1_files/FCCee_o5_v04_ILCSoft-2017-07-27_gcc62_photons_cosTheta_v1_E10_*\" -n 10 --energy 9 11 --theta 50 60 --phi 0 10")
		("filesTemplate,f", po::value<string>()->required(), "file template")
		("nFiles,n", po::value<unsigned int>(), "Set up limit on number of files to read")
		("energy", po::value<vector<double> >()->multitoken(), "To specify energy ranges. \nFormat: 10 500")
		("theta", po::value<vector<double> >()->multitoken(), "To specify theta ranges. \nFormat: 0 45")
		("phi", po::value<vector<double> >()->multitoken(), "To specify phi ranges. \nFormat: 0 90")
		("minE", po::value<double>(), "minimum energy")
		("maxE", po::value<double>(), "maximum energy")
		("minTh", po::value<double>(), "minimum theta")
		("maxTh", po::value<double>(), "maximum theta")
		("minPhi", po::value<double>(), "minimum phi")
		("maxPhi", po::value<double>(), "maximum phi")
		("motherPDG", po::value<unsigned int>(), "Set mother PDG")
		("partTypeToSelect", po::value<unsigned int>(), "Set particle type to select")
		("debug,d", "debug flag")
		;

	/// get global input arguments
	const size_t returnedMessage = parseOptionsWithBoost(vm,argn,argv, desc);
	if (returnedMessage!=SUCCESS) 
		std::exit(returnedMessage);

	// Read collections
	std::vector<std::string> m_fileNames = getFilesMatchingPattern(vm["filesTemplate"].as<string>());
	if (m_fileNames.size()==0){
		cout << "[ERROR]\t[StudyElectronPerformance] No input files found..." << endl;
		return 0;
	}


	if (vm.count("nFiles"))
		if (vm["nFiles"].as<unsigned int>()<m_fileNames.size())
			m_fileNames.resize(vm["nFiles"].as<unsigned int>());

	cout << endl << "[INFO]\tNumber of input files to be used: " << m_fileNames.size() << " files" << endl;

	// Open Files
	auto m_reader( IOIMPL::LCFactory::getInstance()->createLCReader());
	try{
		m_reader->open( m_fileNames );
	} catch (IO::IOException &e)  {
		std::cerr << "Error opening files: " << e.what() << std::endl;
		return 1;
	}

	if (vm.count("debug")){
		cout << "First file to be read: " << m_fileNames[0] << endl;
		cout << "Number of events to be read: " << m_reader->getNumberOfEvents() << endl;
	}

	vector<string> collectionsToRead = {};
	if (!vm.count("accessCaloHitInfo")) {
		energyFillCollections = {};
	}
	else{
		energyFillCollections = {"ECALBarrel","ECALEndcap", "HCALBarrel","HCALEndcap","ECalBarrelCollection", "ECalEndcapCollection","HCalBarrelCollection", "HCalEndcapCollection"};
	}
	collectionsToRead.insert(collectionsToRead.end(),energyFillCollections.begin(),energyFillCollections.end());
	collectionsToRead.insert(collectionsToRead.end(),particleFillCollections.begin(),particleFillCollections.end());
	collectionsToRead.insert(collectionsToRead.end(),additionalCollections.begin(),additionalCollections.end());

	cout << endl << "Collections to be read:" << endl;
	for (int kk=0; kk<collectionsToRead.size(); kk++){
		cout << "- " << collectionsToRead[kk] << endl;
	}
	cout << endl;
	// m_reader->setReadCollectionNames(collectionsToRead);
 
        
	vector<double> energyRange = {0.0,std::numeric_limits<double>::max()};
	vector<double> thetaRange = {-180.0,180.0};
	vector<double> phiRange = {-360.0,360.0};

	if (vm.count("energy"))
		energyRange = vm["energy"].as<vector<double> >();
	if (vm.count("theta"))
		thetaRange = vm["theta"].as<vector<double> >();
	if (vm.count("phi"))
		phiRange = vm["phi"].as<vector<double> >();
	if (energyRange.size()!=2 || thetaRange.size()!=2 || phiRange.size()!=2){
		cout << "[ERROR]\tWrong input range of energy/theta/phi!!!" << endl;
		return 1;
	}

	if (vm.count("minE") && vm.count("maxE")){
		energyRange = {};
		energyRange.push_back(vm["minE"].as<double>());
		energyRange.push_back(vm["maxE"].as<double>());
	}

	if (vm.count("minTh") && vm.count("maxTh")){
		thetaRange = {};
		thetaRange.push_back(vm["minTh"].as<double>());
		thetaRange.push_back(vm["maxTh"].as<double>());
	}
	if (vm.count("minPhi") && vm.count("maxPhi")){
		phiRange = {};
		phiRange.push_back(vm["minPhi"].as<double>());
		phiRange.push_back(vm["maxPhi"].as<double>());
	}

	// LOOP OVER EVENTS
	if (vm.count("debug")) cout << "Reading first event..." << endl;
	EVENT::LCEvent *event = m_reader->readNextEvent();
	int eventCounter = 0;
	if (vm.count("debug"))
		cout << "First event pointer: " << event << endl;

	truthZWCondition::instance()->setMCTruthCollectionName("MCParticlesSkimmed");
	// truthZWCondition::instance()->setMCTruthCollectionName("MCParticle");
	truthZWCondition::instance()->setPartTypeToSelect(13);
	if (vm.count("partTypeToSelect"))
		truthZWCondition::instance()->setPartTypeToSelect(vm["partTypeToSelect"].as<unsigned int>());
	truthZWCondition::instance()->setEnergySelectionCut(1);
	truthZWCondition::instance()->setThetaSelectionCut(18);
	truthZWCondition::instance()->setMotherPDG(24);
	if (vm.count("motherPDG"))
		truthZWCondition::instance()->setMotherPDG(vm["motherPDG"].as<unsigned int>());
	if (vm.count("debug"))
		truthZWCondition::instance()->setDebugFlag(true);

	unsigned int nEvents = m_reader->getNumberOfEvents() * m_fileNames.size();
	cout << "[INFO]\nNumber of events to process: " << nEvents << endl;

	while ( event != NULL ) {
		if (eventCounter % 1000 == 0)
			cout << "[INFO]\tProcess " << eventCounter << " events; Done: " << round(100*eventCounter/nEvents) << "%" << endl;
		bool passKinematicCuts = true;
		if (vm.count("debug")) 
			cout << endl << "[INFO]\t *****EVENT: " << event->getEventNumber() << " *****" <<endl;
		truthZWCondition::instance()->processEvent(event);
		eventCounter++;

		vector<EVENT::MCParticle*> genParts = truthZWCondition::instance()->getGenParticles();

		// const double *partMom = part->getMomentum();
		// TVector3 v1;
		// v1.SetXYZ(partMom[0],partMom[1],partMom[2]);
		// double partTheta = 180.*v1.Theta()/TMath::Pi();
		// double partMomMag = v1.Mag();
		// double partPhi = 180.*v1.Phi()/TMath::Pi();
                //
		// if ((partMomMag<energyRange[0]) || (partMomMag>energyRange[1]))
		//         passKinematicCuts = false;
		// if ((partTheta<thetaRange[0]) || (partTheta>thetaRange[1]))
		//         passKinematicCuts = false;
		// if ((partPhi<phiRange[0]) || (partPhi>phiRange[1]))
		//         passKinematicCuts = false;

		if (passKinematicCuts==true){
			for(auto const &mapElement : objFillMap){
				if (genParts.size()>0)
					mapElement.second->fillEvent(event);
			}
		}

		event = m_reader->readNextEvent();
	}

	string outFileName = "muonRecoMisID";
	if (energyRange[1]!=std::numeric_limits<double>::max())
		outFileName += "_E" + DoubToStr(energyRange[0]) + "_" + DoubToStr(energyRange[1]);
	outFileName += ".root";
	TFile *outFile = new TFile(outFileName.c_str(), "RECREATE");
	for(auto const &mapElement : objFillMap){
		mapElement.second->writeToFile(outFile);
	}
	outFile->Close();

	for(auto const &mapElement : objFillMap){
		mapElement.second->DeleteHists();
		delete mapElement.second;
	} 
}


