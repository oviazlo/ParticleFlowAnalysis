//custom libs
#include <particleFill.h>
#include <energyFill.h>
#include <photonEffCalculator.h>
#include <truthParticleSelector.h>
#include <boostServiceFunctions.h>

using namespace std;

// COLLECTIONS TO USE
// vector<string> particleFillCollections = {"SiTracks","MCParticlesSkimmed","PandoraPFOs","PandoraPFOs","PandoraPFOs","PandoraPFOs","PandoraPFOs","PandoraPFOs","PandoraPFOs","PandoraPFOs","PandoraPFOs","PandoraPFOs","PandoraPFOs"};
vector<string> particleFillCollections = {"MCParticlesSkimmed","PandoraPFOs","PandoraPFOs"};
// 11:		electron
// 13:		muon
// 22:		photon
// 211: 	pi+
// 2112:	neutron
// vector<vector<int> > PFOPartTypes = {{},{},{11},{-11},{13},{-13},{22},{-211},{211},{2112},{11,22},{11,-11,13,-13,211,-211},{22,2112}};
vector<vector<int> > PFOPartTypes = {{},{22},{2112}};
// int efficiencyPFOType = 11;
// vector<string> energyFillCollections = {"ECALBarrel","ECALEndcap"[>, "ECalBarrelCollection", "ECalEndcapCollection"<]};
// vector<string> energyFillCollections = {"ECALBarrel","ECALEndcap","ECalBarrelCollection", "ECalEndcapCollection", "HCALBarrel","HCALEndcap","HCalBarrelCollection", "HCalEndcapCollection"};
vector<string> energyFillCollections = {"ECALBarrel","ECALEndcap", "HCALBarrel","HCALEndcap"};
po::variables_map vm;

int main (int argn, char* argv[]) {

	po::options_description desc("Options");
	desc.add_options()
		("help,h", "Example:\n ./AnalyzeIsotropPhotonSample -f \"/ssd/viazlo/data/FCCee_o5_v04_ILCSoft-2017-07-27_gcc62_photons_cosTheta_v1_files/FCCee_o5_v04_ILCSoft-2017-07-27_gcc62_photons_cosTheta_v1_E10_*\" -n 10 --energy 9 11 --theta 50 60 70 80 90 --phi 0 10 20")
		("filesTemplate,f", po::value<string>()->required(), "file template")
		("nFiles,n", po::value<unsigned int>(), "Set up limit on number of files to read")
		("energy", po::value<vector<double> >()->multitoken(), "To specify energy ranges. \nFormat: 10 50 100")
		("theta", po::value<vector<double> >()->multitoken(), "To specify theta ranges. \nFormat: 0 45 90")
		("phi", po::value<vector<double> >()->multitoken(), "To specify phi ranges. \nFormat: 0 90 180")
		("minE", po::value<double>(), "minimum energy")
		("maxE", po::value<double>(), "maximum energy")
		("minTh", po::value<double>(), "minimum theta")
		("maxTh", po::value<double>(), "maximum theta")
		("minPhi", po::value<double>(), "minimum phi")
		("maxPhi", po::value<double>(), "maximum phi")
		("nE", po::value<unsigned int>(), "number of energy bins")
		("nTh", po::value<unsigned int>(), "number of theta bins")
		("nPhi", po::value<unsigned int>(), "number of phi bins")
		("effPfoType", po::value<int>(), "PFO type to use for efficiency calculation")
		("noFSR", "discard events with FSR (only one truth particle allowed)")
		("dPhiMerge", po::value<double>(), "dPhi value in degrees to merge clusters")
		;

	/// get global input arguments
	const size_t returnedMessage = parseOptionsWithBoost(vm,argn,argv, desc);
	if (returnedMessage!=SUCCESS) std::exit(returnedMessage);

	// Read collections
	std::vector<std::string> m_fileNames = getFilesMatchingPattern(vm["filesTemplate"].as<string>());
	if (m_fileNames.size()==0){
		cout << "[ERROR]\t[AnalyzeIsotropPhotonSample] No input files found..." << endl;
		return 0;
	}


	if (vm.count("nFiles"))
		if (vm["nFiles"].as<unsigned int>()<m_fileNames.size())
			m_fileNames.resize(vm["nFiles"].as<unsigned int>());

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

	vector<double> energyRanges = {9.9,10.1};
	vector<double> thetaRanges = {-180.0,180.0};
	vector<double> phiRanges = {-360.0,360.0};

	if (vm.count("energy"))
		energyRanges = vm["energy"].as<vector<double> >();
	if (vm.count("theta"))
		thetaRanges = vm["theta"].as<vector<double> >();
	if (vm.count("phi"))
		phiRanges = vm["phi"].as<vector<double> >();

	if (vm.count("minE") && vm.count("maxE") && vm.count("nE")){
		energyRanges = {};
		for (auto i=0; i<vm["nE"].as<unsigned int>()+1; i++){
			double iE = vm["minE"].as<double>() + i*( vm["maxE"].as<double>() - vm["minE"].as<double>() )/vm["nE"].as<unsigned int>();
			energyRanges.push_back(iE);
		}
	}
	
	if (vm.count("minTh") && vm.count("maxTh") && vm.count("nTh")){
		thetaRanges = {};
		for (auto i=0; i<vm["nTh"].as<unsigned int>()+1; i++){
			double iTh = vm["minTh"].as<double>() + i*( vm["maxTh"].as<double>() - vm["minTh"].as<double>() )/vm["nTh"].as<unsigned int>();
			thetaRanges.push_back(iTh);
		}
	}
	if (vm.count("minPhi") && vm.count("maxPhi") && vm.count("nPhi")){
		phiRanges = {};
		for (auto i=0; i<vm["nPhi"].as<unsigned int>()+1; i++){
			double iPhi = vm["minPhi"].as<double>() + i*( vm["maxPhi"].as<double>() - vm["minPhi"].as<double>() )/vm["nPhi"].as<unsigned int>();
			phiRanges.push_back(iPhi);
		}
	}

	vector<truthParticleSelector*> selectors;
	for(auto iE=0; iE<energyRanges.size()-1;iE++){
	for(auto iTheta=0; iTheta<thetaRanges.size()-1;iTheta++){
	for(auto iPhi=0; iPhi<phiRanges.size()-1;iPhi++){
		
		truthParticleSelector *sel1 = new truthParticleSelector("MCParticlesSkimmed");
		sel1->setEnergyRange(energyRanges[iE],energyRanges[iE+1]);
		sel1->setThetaRange(thetaRanges[iTheta],thetaRanges[iTheta+1]);
		sel1->setPhiRange(phiRanges[iPhi],phiRanges[iPhi+1]);
		sel1->setEfficiencyCollection("PandoraPFOs");
		if (vm.count("effPfoType"))
			sel1->setEfficiencyPFOType(vm["effPfoType"].as<int>());
		if (vm.count("dPhiMerge"))
			sel1->setDPhiMergeValue(vm["dPhiMerge"].as<double>());
		sel1->setParticleFillCollections(particleFillCollections);
		sel1->setPFOTypes(PFOPartTypes);
		sel1->setEnergyFillCollections(energyFillCollections);
		sel1->init();
		if (vm.count("noFSR"))
			sel1->setDiscardFSREvents(true);
		selectors.push_back(sel1);

	}
	}
	}


	// LOOP OVER EVENTS
	EVENT::LCEvent *event = NULL;
	int eventCounter = 0;
	while ( ( event = m_reader->readNextEvent() ) ) {
		// if (eventCounter>100) break;
		if (printVerbose) cout << endl << "Event " << eventCounter << ":" << endl;
		eventCounter++;
	
		for(auto i=0; i<selectors.size(); i++){
			selectors[i]->selectEvent(event);
		}
	}

	for(auto i=0; i<selectors.size(); i++){
		TFile *outFile = new TFile(("particleGun_"+ selectors[i]->getPostFixString() +".root").c_str(), "RECREATE");
		selectors[i]->writeToFile(outFile); 
		outFile->Close();
	}

	// SAVE OUTPUT HISTOGRAMS
	// double scale = 0.1;
	// string meanEnergy = DoubToStr( floor(static_cast<particleFill*>(objFillMap["MCParticlesSkimmed"])->getMeanEnergy() / scale + 0.5)*scale );
	// string meanTheta = DoubToStr( floor(static_cast<particleFill*>(objFillMap["MCParticlesSkimmed"])->getMeanTheta() / scale + 0.5)*scale );
	// TFile *outFile = new TFile(("ECAL_photonGun_E"+meanEnergy+"_theta"+meanTheta+".root").c_str(), "RECREATE");
        //
	// for(auto const &mapElement : objFillMap){
	//         mapElement.second->writeToFile(outFile);
	// }

	// outFile->Close();

}

