#include <energyFillAllCalo.h>
int energyFillAllCalo::init(){

	createTH1D("CALratio_EcalToHcal_BarrelAndEndcaps","Ecal/Hcal Energy Ratio; Ecal/Hcal; Counts",500,0,50);
	createTH1D("CALratio_HcalToEcal_BarrelAndEndcaps","Hcal/Ecal Energy Ratio; Hcal/Ecal; Counts",500,0,50);
	createTH1D("CALratio_EcalToHcal_BarrelOnly","Ecal/Hcal Barrel Energy Ratio; Barrel Ecal/Hcal; Counts",500,0,50);
	createTH1D("CALratio_HcalToEcal_BarrelOnly","Hcal/Ecal Barrel Energy Ratio; Barrel Hcal/Ecal; Counts",500,0,50);
	createTH1D("CALratio_EcalToHcal_EndcapOnly","Ecal/Hcal Endcap Energy Ratio; Endcap Ecal/Hcal; Counts",500,0,50);
	createTH1D("CALratio_HcalToEcal_EndcapOnly","Hcal/Ecal Endcap Energy Ratio; Endcap Hcal/Ecal; Counts",500,0,50);

	for(auto const &iMapElement : histMap) {
		iMapElement.second->AddDirectory(kFALSE);
	}

	return 0;
}

int energyFillAllCalo::fillEvent(const EVENT::LCEvent* event){

	vector<string> energyFillCollections = {"ECALBarrel","ECALEndcap", "HCALBarrel","HCALEndcap"};
	map<string, double> totalEnergyMap;

	for (auto collectionName: energyFillCollections){
		try {
			collection = event->getCollection(collectionName);
		} catch (EVENT::DataNotAvailableException &e) {
			cout << "[ERROR|energyFillAllCalo]\tCan't find collection: " << collectionName << endl;
			return -1;
		}

		if( collection ) {
			const int nElements = collection->getNumberOfElements();
			double  totalEnergyDeposited = 0.0;
			for(int j=0; j < nElements; j++){
				auto calHit = dynamic_cast<EVENT::CalorimeterHit*>(collection->getElementAt(j));
				totalEnergyDeposited += calHit->getEnergy();
			}
			totalEnergyMap[collectionName] = totalEnergyDeposited;
		}
	}

	getHistFromMap("CALratio_EcalToHcal_BarrelAndEndcaps")->Fill((totalEnergyMap["ECALBarrel"] + totalEnergyMap["ECALEndcap"]) / (totalEnergyMap["HCALBarrel"] + totalEnergyMap["HCALEndcap"]));
	getHistFromMap("CALratio_HcalToEcal_BarrelAndEndcaps")->Fill((totalEnergyMap["HCALBarrel"] + totalEnergyMap["HCALEndcap"]) / (totalEnergyMap["ECALBarrel"] + totalEnergyMap["ECALEndcap"]));

	getHistFromMap("CALratio_EcalToHcal_BarrelOnly")->Fill(totalEnergyMap["ECALBarrel"]/totalEnergyMap["HCALBarrel"]);
	getHistFromMap("CALratio_HcalToEcal_BarrelOnly")->Fill(totalEnergyMap["HCALBarrel"]/totalEnergyMap["ECALBarrel"]);
	getHistFromMap("CALratio_EcalToHcal_EndcapOnly")->Fill(totalEnergyMap["ECALEndcap"]/totalEnergyMap["HCALEndcap"]);
	getHistFromMap("CALratio_HcalToEcal_EndcapOnly")->Fill(totalEnergyMap["HCALEndcap"]/totalEnergyMap["ECALEndcap"]);

	return 0;
}

int energyFillAllCalo::writeToFile(TFile* outFile){

	// unsigned int nEvents = getHistFromMap("CAL_Energy")->GetEntries();

	objectFill::writeToFile(outFile);
}
