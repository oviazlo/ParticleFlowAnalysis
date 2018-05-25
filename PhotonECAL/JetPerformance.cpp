//custom libs
#include <jetPfoStudy.h>
#include <truthParticleSelector.h>
#include <boostServiceFunctions.h>
// #include "truthCondition.h"
// #include <globalConfig.h>


using namespace std;
using namespace config; 

vector<string> collectionsToRead = {};
// vector<string> particleFillCollections = {"PandoraPFOs","LooseSelectedPandoraPFOs","SelectedPandoraPFOs","TightSelectedPandoraPFOs"};
vector<string> particleFillCollections = {"PandoraPFOs"};
vector<string> energyFillCollections = {};
vector<string> additionalCollections = {"MCParticle"};

int main (int argn, char* argv[]) {

	po::options_description desc("Options");
	desc.add_options()
		("help,h", "Example:\n ./JetPerformance -f \"/ssd/viazlo/data/FCCee_o5_v04_ILCSoft-2017-07-27_gcc62_photons_cosTheta_v1_files/FCCee_o5_v04_ILCSoft-2017-07-27_gcc62_photons_cosTheta_v1_E10_*\" -n 10")
		("filesTemplate,f", po::value<string>()->required(), "file template")
		("nFiles,n", po::value<unsigned int>(), "Set up limit on number of files to read")
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

	collectionsToRead.insert(collectionsToRead.end(),energyFillCollections.begin(),energyFillCollections.end());
	collectionsToRead.insert(collectionsToRead.end(),particleFillCollections.begin(),particleFillCollections.end());
	collectionsToRead.insert(collectionsToRead.end(),additionalCollections.begin(),additionalCollections.end());

	cout << endl << "Collections to be read:" << endl;
	for (int kk=0; kk<collectionsToRead.size(); kk++){
		cout << "- " << collectionsToRead[kk] << endl;
	}
	cout << endl;
	// m_reader->setReadCollectionNames(collectionsToRead);
        
	vector<jetPfoStudy*> selectors;
	for (auto colIt=particleFillCollections.begin(); colIt!=particleFillCollections.end(); colIt++){
		jetPfoStudy *sel = new jetPfoStudy("jetStudy",additionalCollections[0],*colIt);
		// jetPfoStudy *sel = new jetPfoStudy("jetStudy_"+*colIt,additionalCollections[0],*colIt);
		sel->init();
		if (vm.count("debug")) 
			sel->setDebugFlag(true);
		selectors.push_back(sel);
	}


	TFile *outFile = new TFile("jetStudy.root", "RECREATE");
	outFile->Close();

	// LOOP OVER EVENTS
	if (vm.count("debug")) cout << "Reading first event..." << endl;
	EVENT::LCEvent *event = m_reader->readNextEvent();
	// EVENT::LCEvent *event = m_reader->readEvent(0,1);
	int eventCounter = 0;
	if (vm.count("debug"))
		cout << "First event pointer: " << event << endl;

	while ( event != NULL ) {
		if (vm.count("debug")) 
			cout << endl << "Event " << eventCounter << ":" << endl;
	
		for(auto i=0; i<selectors.size(); i++){
			selectors[i]->fillEvent(event);
		}
		event = m_reader->readNextEvent();
		eventCounter++;
		if (eventCounter%10==0)
			cout << "Event processed: " << eventCounter << ":" << endl;
	}

	for(auto i=0; i<selectors.size(); i++){
		TFile *outFile = new TFile("jetStudy.root", "UPDATE");
		selectors[i]->writeToFile(outFile);
		outFile->Close();
	}

}

