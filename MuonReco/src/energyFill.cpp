#include <energyFill.h>
int energyFill::init(){
	TH1D *h_ECALTotalEnergy = new TH1D("CAL_Energy","energy; Total Energy [GeV]; Counts",250,0.0,125);
	TH1I *h_ECALNHits = new TH1I("CAL_Nhits","Number of hits; Number of ECAL hits; Counts",3000,0,3000);
	TH1I *h_ECALMaxLayer = new TH1I("CAL_maxLayer","Max fired layer; max. layer; Counts",45,0,45);
	TH1D *h_ECALEnergyPerLayers = new TH1D("CAL_EnergyPerLayers","Energy per layer; Layer number; Deposited energy [???]",50,0,50);
	TH1D *h_ECALAverageEnergyPerLayers = new TH1D("CAL_AverageEnergyPerLayers","Average Energy per layer; Layer number; Deposited energy fraction",50,0,50);

	histMap["CAL_Energy"] = h_ECALTotalEnergy;
	histMap["CAL_Nhits"] = h_ECALNHits;
	histMap["CAL_maxLayer"] = h_ECALMaxLayer;
	histMap["CAL_EnergyPerLayers"] = h_ECALEnergyPerLayers;
	histMap["CAL_AverageEnergyPerLayers"] = h_ECALAverageEnergyPerLayers;

	createTH1D("CAL_hitTiming","Hit Timing; Time [ns]; Counts",5000,0,500);
	createTH1D("CAL_energyVsHitTiming","Energy vs. Hit Timing; Time [ns]; Energy [GeV / 1 ns]",500,0,500);
	createTH1D("CAL_energyVsRadius","Hit Energy vs Radius; Radius [mm]; Hit Energy [GeV / 50 mm]",80,0,4000);
	createTH1D("CAL_energyVsZ","Hit Energy vs Z; Z [mm]; Hit Energy [GeV / 51 mm]",160,-4088,4088);
	createTH1D("CAL_energyVsZ2","Hit Energy vs Z; Z [mm]; Hit Energy [GeV / 51 mm]",160,-4090.5,4085.5);
	createTH1D("CAL_energyVsZ3","Hit Energy vs Z; Z [mm]; Hit Energy [GeV / 51 mm]",26,-3900,3900);
	createTH1D("CAL_energyVsTheta","Hit Energy vs Theta; Theta [deg]; Hit Energy [GeV / 5 deg]",36,0,180);
	createTH1D("CAL_energyVsCosTheta","Hit Energy vs cosTheta; cos(Theta); Hit Energy [GeV / 0.02]",40,-1,1);
	createTH1D("CAL_energySpectrum","Hit Energy Spectrum; Enegy [GeV]; Counts",10000,0,1);

	for(auto const &iMapElement : histMap) {
		iMapElement.second->AddDirectory(kFALSE);
	}
	//Initialize CellID encoder
	// m_encoder = new UTIL::BitField64(lcio::LCTrackerCellID::encoding_string());

	return 0;
}


int energyFill::fillEvent(const EVENT::LCEvent* event){

	try {
		collection = event->getCollection(collectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|energyFill]\tCan't find collection: " << collectionName << endl;
		return -1;
	}

	if( collection ) {
		string collectionType = collection->getTypeName();
		const int nElements = collection->getNumberOfElements();
		double energyCount = 0.0;
		double totalEnergyDeposited = 0.0;
		int maxLayer = 0;
		for(int j=0; j < nElements; j++){
			auto calHit = dynamic_cast<EVENT::CalorimeterHit*>(collection->getElementAt(j));
			totalEnergyDeposited += calHit->getEnergy();
			// int nContributions = calHit->getNMCContributions();
			// for (int iCont = 0; iCont < nContributions;++iCont){
			//         energyCount = calHit->getEnergyCont(iCont);
			//         totalEnergyDeposited += energyCount;
			// }

			// cout << endl << "totalEnergy: " << calHit->getEnergy() << endl;

			UTIL::BitField64 _encoder(collection->getParameters().getStringVal( EVENT::LCIO::CellIDEncoding ));
			lcio::long64 cellId = long( calHit->getCellID0() & 0xffffffff ) | ( long( calHit->getCellID1() ) << 32 );
			_encoder.setValue(cellId);
			int layer=_encoder["layer"].value();
			if (layer>maxLayer)
				maxLayer = layer;
			histMap["CAL_EnergyPerLayers"]->Fill(layer,calHit->getEnergy());
			histMap["CAL_AverageEnergyPerLayers"]->Fill(layer,calHit->getEnergy());
			auto caloRawHit = dynamic_cast<EVENT::SimCalorimeterHit*>(calHit->getRawHit());
			const unsigned int n = caloRawHit->getNMCContributions(); //number of subhits of this SimHit
			for (unsigned int i_t=0; i_t<n; i_t++){
				float timei = caloRawHit->getTimeCont(i_t); //absolute hit timing of current subhit
				getHistFromMap("CAL_hitTiming")->Fill(timei);
				getHistFromMap("CAL_energyVsHitTiming")->Fill(timei,caloRawHit->getEnergyCont(i_t));
				// cout << "deposit_" << i_t << ": " << caloRawHit->getEnergyCont(i_t) << endl;
			}
			auto hitPos = calHit->getPosition();
			double hitRadius = sqrt((hitPos[0] * hitPos[0]) + (hitPos[1] * hitPos[1]));
			double hitZ = hitPos[2];
			double hitTheta = TVector3(hitPos[0],hitPos[1],hitPos[2]).Theta();
			double hitCosTheta = cos(hitTheta);
			getHistFromMap("CAL_energyVsRadius")->Fill(hitRadius,calHit->getEnergy());
			getHistFromMap("CAL_energyVsZ")->Fill(hitZ,calHit->getEnergy());
			getHistFromMap("CAL_energyVsZ2")->Fill(hitZ,calHit->getEnergy());
			getHistFromMap("CAL_energyVsZ3")->Fill(hitZ,calHit->getEnergy());
			getHistFromMap("CAL_energyVsTheta")->Fill(hitTheta*TMath::RadToDeg(),calHit->getEnergy());
			getHistFromMap("CAL_energyVsCosTheta")->Fill(hitCosTheta,calHit->getEnergy());
			getHistFromMap("CAL_energySpectrum")->Fill(calHit->getEnergy());

		}
		fillMaxLayer(maxLayer);
		fillEnergy(totalEnergyDeposited);
		fillNHits(nElements);
	}
	return 0;
}

int energyFill::fillEnergy(const double energy){
	return histMap["CAL_Energy"]->Fill(energy);
}
int energyFill::fillNHits(const int nHits){
	return histMap["CAL_Nhits"]->Fill(nHits);
}
int energyFill::fillMaxLayer(const int maxLayer){
	return histMap["CAL_maxLayer"]->Fill(maxLayer);
}

int energyFill::writeToFile(TFile* outFile){
	histMap["CAL_AverageEnergyPerLayers"]->Scale(1.0/histMap["CAL_Energy"]->GetEntries()/histMap["CAL_Energy"]->GetMean());

	unsigned int nEvents = getHistFromMap("CAL_Energy")->GetEntries();
	getHistFromMap("CAL_energyVsRadius")->Scale(1.0/nEvents);
	getHistFromMap("CAL_energyVsZ")->Scale(1.0/nEvents);
	getHistFromMap("CAL_energyVsZ2")->Scale(1.0/nEvents);
	getHistFromMap("CAL_energyVsZ3")->Scale(1.0/nEvents);
	getHistFromMap("CAL_energyVsTheta")->Scale(1.0/nEvents);
	getHistFromMap("CAL_energyVsCosTheta")->Scale(1.0/nEvents);
	getHistFromMap("CAL_energyVsHitTiming")->Scale(1.0/nEvents);

	objectFill::writeToFile(outFile);
}
