#include <energyFill.h>

int energyFill::init(){
	TH1D *h_ECALTotalEnergy = new TH1D("CAL_Energy","energy; Total Energy [GeV]; Counts",250,0.0,125);
	TH1I *h_ECALNHits = new TH1I("CAL_Nhits","Number of hits; Number of ECAL hits; Counts",3000,0,3000);
	TH1I *h_ECALMaxLayer = new TH1I("CAL_maxLayer","Max fired layer; max. layer; Counts",45,0,45);
	TH1D *h_ECALEnergyPerLayers = new TH1D("CAL_EnergyPerLayers","Energy per layer; Layer number; Deposited energy [MeV]",50,0,50);
	TH1D *h_ECALAverageEnergyPerLayers = new TH1D("CAL_AverageEnergyPerLayers","Average Energy per layer; Layer number; Deposited energy fraction",50,0,50);

	histMap["CAL_Energy"] = h_ECALTotalEnergy;
	histMap["CAL_Nhits"] = h_ECALNHits;
	histMap["CAL_maxLayer"] = h_ECALMaxLayer;
	histMap["CAL_EnergyPerLayers"] = h_ECALEnergyPerLayers;
	histMap["CAL_AverageEnergyPerLayers"] = h_ECALAverageEnergyPerLayers;

	for(auto const &iMapElement : histMap) {
		iMapElement.second->AddDirectory(kFALSE);
	}
	//Initialize CellID encoder
	// m_encoder = new UTIL::BitField64(lcio::LCTrackerCellID::encoding_string());

	return 0;
}


int energyFill::fillEvent(EVENT::LCEvent* event){

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
			auto calHit = static_cast<EVENT::SimCalorimeterHit*>(collection->getElementAt(j));
			totalEnergyDeposited += calHit->getEnergy();
			// int nContributions = calHit->getNMCContributions();
			// for (int iCont = 0; iCont < nContributions;++iCont){
			//         energyCount = calHit->getEnergyCont(iCont);
			//         totalEnergyDeposited += energyCount;
			// }

			UTIL::BitField64 _encoder(collection->getParameters().getStringVal( EVENT::LCIO::CellIDEncoding ));
			lcio::long64 cellId = long( calHit->getCellID0() & 0xffffffff ) | ( long( calHit->getCellID1() ) << 32 );
			_encoder.setValue(cellId);
			int layer=_encoder["layer"].value();
			if (layer>maxLayer)
				maxLayer = layer;
			histMap["CAL_EnergyPerLayers"]->Fill(layer,calHit->getEnergy());
			histMap["CAL_AverageEnergyPerLayers"]->Fill(layer,calHit->getEnergy());

		}
		fillMaxLayer(maxLayer);
		fillEnergy(totalEnergyDeposited);
		fillNHits(nElements);
	}
	return 0;
}

int energyFill::fillEnergy(double energy){
	return histMap["CAL_Energy"]->Fill(energy);
}
int energyFill::fillNHits(int nHits){
	return histMap["CAL_Nhits"]->Fill(nHits);
}
int energyFill::fillMaxLayer(int maxLayer){
	return histMap["CAL_maxLayer"]->Fill(maxLayer);
}

int energyFill::writeToFile(TFile* outFile){
	histMap["CAL_AverageEnergyPerLayers"]->Scale(1.0/histMap["CAL_Energy"]->GetEntries()/histMap["CAL_Energy"]->GetMean());
	objectFill::writeToFile(outFile);
}
