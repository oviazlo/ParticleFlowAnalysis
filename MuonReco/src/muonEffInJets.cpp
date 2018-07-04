#include "muonEffInJets.h"

int muonEffInJets::fillEvent(const EVENT::LCEvent* event){
	try {
		PFOCollection = event->getCollection(PFOCollectionName);
	} catch (EVENT::DataNotAvailableException &e) {
		cout << "[ERROR|muonEffInJets]\tCan't find collection: " << PFOCollectionName << endl;
		return -1;
	
	}
	nPFOs = PFOCollection->getNumberOfElements();
	// cout << "nPFOs: " << nPFOs << endl;

	matchedRecoPartsOfType.clear();
	nonMatchedRecoPartsOfType.clear();
	matchedTruthPartsOfType.clear();
	nonMatchedTruthPartsOfType.clear();

	findMatchedRecoParts();
	fillTruthInfo();
	fillMuonClusterInfo();
	return 0;

}






void muonEffInJets::fillTruthInfo(){
	fillTruthInfoHelper("truthPartAll",matchedTruthPartsOfType);
	fillTruthInfoHelper("truthPartAll",nonMatchedTruthPartsOfType);
	fillTruthInfoHelper("truthPartMatched",matchedTruthPartsOfType);
	fillTruthInfoHelper("eff_truthPartMatched",matchedTruthPartsOfType);
	fillTruthInfoHelper("truthPartNotMatched",nonMatchedTruthPartsOfType);
}

void muonEffInJets::fillTruthInfoHelper(string histNamePrefix, vector<EVENT::MCParticle*> inTruthPartVec){
	for (auto genPart: inTruthPartVec){
		TVector3 vTruthMom(genPart->getMomentum());
		double truthTheta = vTruthMom.Theta()*TMath::RadToDeg();
		double cosTruthTheta = TMath::Cos(truthTheta*TMath::DegToRad());
		double truthPhi = vTruthMom.Phi()*TMath::RadToDeg();
		double truthPt = vTruthMom.Pt();
		double truthEnergy = genPart->getEnergy();

		const double *vertexPos = genPart->getVertex();
		double vertexR = sqrt(vertexPos[0]*vertexPos[0]+vertexPos[1]*vertexPos[1]);
		double vertexZ = vertexPos[2];

		getHistFromMap(histNamePrefix + "_E")->Fill(truthEnergy);
		getHistFromMap(histNamePrefix + "_pt")->Fill(truthPt);
		getHistFromMap(histNamePrefix + "_theta")->Fill(truthTheta);
		getHistFromMap(histNamePrefix + "_cosTheta")->Fill(cosTruthTheta);
		getHistFromMap(histNamePrefix + "_phi")->Fill(truthPhi);
		getHistFromMap(histNamePrefix + "_vertexR")->Fill(vertexR);
		getHistFromMap(histNamePrefix + "_vertexZ")->Fill(vertexZ);
	}
}

void muonEffInJets::findMatchedRecoParts(){
	nonMatchedTruthPartsOfType = truthZWCondition::instance()->getGenParticles();
	for (int iPart=0; iPart<PFOCollection->getNumberOfElements(); iPart++){
		EVENT::ReconstructedParticle* recoPart = static_cast<IMPL::ReconstructedParticleImpl*>(PFOCollection->getElementAt(iPart));
		if (abs(recoPart->getType())!=truthZWCondition::instance()->getPartTypeToSelect())
			continue;
		TVector3 vPartMom(recoPart->getMomentum());
		int nTruthPartToLoop = nonMatchedTruthPartsOfType.size();
		unsigned int counter = 0;
		while(nTruthPartToLoop>0){
		// for (auto genPart: nonMatchedTruthPartsOfType){
			auto genPart = nonMatchedTruthPartsOfType[counter];
			TVector3 vTruthMom(genPart->getMomentum());
			double dAngle = TMath::RadToDeg()*vPartMom.Angle(vTruthMom);
			double dPt = abs(vPartMom.Pt() - vTruthMom.Pt());

			if (dAngle<1.0 && dPt < 0.05*vTruthMom.Pt()){
			// if (dAngle<1.0){
				matchedRecoPartsOfType.push_back(recoPart);
				matchedTruthPartsOfType.push_back(genPart);
				nonMatchedTruthPartsOfType.erase(nonMatchedTruthPartsOfType.begin()+counter);
				counter--;
			}
			else{
				nonMatchedRecoPartsOfType.push_back(recoPart);
			}
			counter++;
			nTruthPartToLoop--;
		}
	}
}



int muonEffInJets::fillMuonClusterInfo(){

	auto allRecoPartsOfType(matchedRecoPartsOfType);
	allRecoPartsOfType.insert(allRecoPartsOfType.end(),nonMatchedRecoPartsOfType.begin(),nonMatchedRecoPartsOfType.end());

	for (auto recoPart : allRecoPartsOfType){
		auto clusters = recoPart->getClusters();
		// cout << "[MuonReco]\tnClusters: " << clusters.size() << endl;
		
		auto cluster = clusters[0];
		// auto clusterShape = getClusterShape(cluster);
		// if (clusterShape.size()>3){
		//         cout << "xt90: " << clusterShape[0] << "; depth of the cluster: " << clusterShape[1] << "; RhitMean: " << clusterShape[2] << "; RhitRMS: " << clusterShape[3] << endl;
		//         getHistFromMap("cluster_xt90")->Fill(clusterShape[0]);
		//         getHistFromMap("cluster_depth")->Fill(clusterShape[1]);
		//         getHistFromMap("cluster_RhitMean")->Fill(clusterShape[2]);
		//         getHistFromMap("cluster_RhitRMS")->Fill(clusterShape[3]);
		//         getHistFromMap("cluster_nYokeHits")->Fill(clusterShape[4]);
		// }
		yokeHitsStruct yokeHitCount = getNumberOfYokeHits(cluster);
		getHistFromMap("cluster_nYokeHits")->Fill(yokeHitCount.nHits);
		getHistFromMap("cluster_nLayers")->Fill(yokeHitCount.nLayers);
		getHistFromMap("cluster_clusterLayerSpan")->Fill(yokeHitCount.clusterLayerSpan);
	}
}

yokeHitsStruct muonEffInJets::getNumberOfYokeHits(EVENT::Cluster* pCluster) { 
	unsigned int nYokeHits = 0;
        unsigned int nHitsInCluster(pCluster->getCalorimeterHits().size());
	std::set<unsigned int> hitLayer;
	for (unsigned int iHit = 0; iHit < nHitsInCluster; ++iHit)
	{
		EVENT::CalorimeterHit *pCalorimeterHit = (EVENT::CalorimeterHit*)(pCluster->getCalorimeterHits()[iHit]);
		if (CHT(pCalorimeterHit->getType()).caloID()==CHT::yoke){
			hitLayer.insert(CHT(pCalorimeterHit->getType()).layer());
			nYokeHits++;
		}
	}
	yokeHitsStruct yokeHitCount;
	yokeHitCount.clusterLayerSpan = *hitLayer.end() - *hitLayer.begin();
	yokeHitCount.nLayers = hitLayer.size();
	yokeHitCount.nHits = nYokeHits;
	// cout << "nYokeHits: " << nYokeHits << "; nLayers: " << nLayers << "; clusterLayerSpan: " << clusterLayerSpan << endl;
	return yokeHitCount;
}

EVENT::FloatVec muonEffInJets::getClusterShape(EVENT::Cluster* pCluster) { 
  
      EVENT::FloatVec shapes;
  const unsigned int ecal_Index(0) ;
  const unsigned int hcal_Index(1) ;
  const unsigned int yoke_Index(2) ;
  const unsigned int lcal_Index(3) ;
  const unsigned int lhcal_Index(4);
  const unsigned int bcal_Index(5) ;

  bool useOnlyYokeHits = false;

      // const unsigned int nHitsInCluster(pCluster->getCalorimeterHits().size());
      unsigned int nHitsInCluster(pCluster->getCalorimeterHits().size());
      // cout << "nAllHits: " << nHitsInCluster << "; ";
      
      if (useOnlyYokeHits==true){
	      vector<unsigned int> idOfYokeHits;
	      for (unsigned int iHit = 0; iHit < nHitsInCluster; ++iHit)
		{
		  EVENT::CalorimeterHit *pCalorimeterHit = (EVENT::CalorimeterHit*)(pCluster->getCalorimeterHits()[iHit]);
		  if (CHT(pCalorimeterHit->getType()).caloID()==CHT::yoke){
			idOfYokeHits.push_back(iHit);
		  }
		}

	      nHitsInCluster = idOfYokeHits.size();
	      // cout << "; nYokeHits: " << nHitsInCluster << endl;
      }
      if (nHitsInCluster==0)
	      return shapes;

      float clusterEnergy(0.);
      float *pHitE = new float[nHitsInCluster];
      float *pHitX = new float[nHitsInCluster];
      float *pHitY = new float[nHitsInCluster];
      float *pHitZ = new float[nHitsInCluster];
      int *typ = new int[nHitsInCluster];

      unsigned int nYokeHits = 0;

      for (unsigned int iHit = 0; iHit < nHitsInCluster; ++iHit)
      // for (auto iHit: idOfYokeHits)
	{
	  EVENT::CalorimeterHit *pCalorimeterHit = (EVENT::CalorimeterHit*)(pCluster->getCalorimeterHits()[iHit]);

	  const float caloHitEnergy(pCalorimeterHit->getEnergy());
	  
	  pHitE[iHit] = caloHitEnergy;
	  pHitX[iHit] = pCalorimeterHit->getPosition()[0];
	  pHitY[iHit] = pCalorimeterHit->getPosition()[1];
	  pHitZ[iHit] = pCalorimeterHit->getPosition()[2];
	  clusterEnergy += caloHitEnergy;

	  switch (CHT(pCalorimeterHit->getType()).caloID())
	    {
	    case CHT::ecal:  typ[iHit]=ecal_Index; break;
	    case CHT::hcal:  typ[iHit]=hcal_Index; break;
	    case CHT::yoke:  typ[iHit]=yoke_Index; nYokeHits++; break;
	    case CHT::lcal:  typ[iHit]=lcal_Index; break;
	    case CHT::lhcal: typ[iHit]=lhcal_Index; break;
	    case CHT::bcal:  typ[iHit]=bcal_Index; break;
	    default: cout << "[DEBUG]\tno subdetector found for hit with type: " << pCalorimeterHit->getType() << std::endl;
	    }
	}
    
      ClusterShapes* pClusterShapes = new ClusterShapes(nHitsInCluster, pHitE, pHitX, pHitY, pHitZ);
      pClusterShapes->setHitTypes(typ);   //set hit types

      //here is cluster shape study - cluster transverse & longitudinal information
      //define variables
      float chi2,a,b,c,d,xl0,CoG[3],xStart[3]; //for fitting parameters
      float X0[2]={0,0};  //in mm. //this is the exact value of tangsten and iron
      float Rm[2]={0,0};  //in mm. need to change to estimate correctly times 2
      float _X01 = 3.50;
      float _X02 = 17.57;
      float _Rm1 = 9.0;
      float _Rm2 = 17.19;
      X0[0]=_X01;
      X0[1]=_X02;
      Rm[0]=_Rm1;
      Rm[1]=_Rm2;

      //get barrel detector surfce
      //Get ECal Barrel extension by type, ignore plugs and rings
      // const dd4hep::rec::LayeredCalorimeterData * eCalBarrelExtension= MarlinUtil::getLayeredCalorimeterData( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::BARREL),
													      // ( dd4hep::DetType::AUXILIARY  |  dd4hep::DetType::FORWARD ) );
      // float ecalrad = eCalBarrelExtension->extent[0]/dd4hep::mm;

      //get endcap detector surfce
      //Get ECal Endcap extension by type, ignore plugs and rings
      // const dd4hep::rec::LayeredCalorimeterData * eCalEndcapExtension= MarlinUtil::getLayeredCalorimeterData( ( dd4hep::DetType::CALORIMETER | dd4hep::DetType::ELECTROMAGNETIC | dd4hep::DetType::ENDCAP),
													      // ( dd4hep::DetType::AUXILIARY |  dd4hep::DetType::FORWARD  ) );
      // float plugz = eCalEndcapExtension->extent[2]/dd4hep::mm;
      
      //looking for the hit which corresponds to the nearest hit from IP in the direction of the center of gravity
      int index_xStart=0;
      float lCoG=0.0,tmpcos=0.0,tmpsin=0.0,detsurface=0.0;
      CoG[0]=pClusterShapes->getEigenVecInertia()[0];
      CoG[1]=pClusterShapes->getEigenVecInertia()[1];
      CoG[2]=pClusterShapes->getEigenVecInertia()[2];
      //CoG2[0]=pCluster->getPosition()[0];
      //CoG2[1]=pCluster->getPosition()[1];
      //CoG2[2]=pCluster->getPosition()[2];
      
      lCoG=sqrt(CoG[0]*CoG[0]+CoG[1]*CoG[1]+CoG[2]*CoG[2]);
      tmpcos=CoG[2]/lCoG;
      tmpsin=sqrt(CoG[0]*CoG[0]+CoG[1]*CoG[1])/lCoG;
      pClusterShapes->fit3DProfile(chi2,a,b,c,d,xl0,xStart,index_xStart,X0,Rm);  //is this good??
      float lxstart=sqrt(xStart[0]*xStart[0]+xStart[1]*xStart[1]);
      //calculate detector surface
      // if(fabs(xStart[2])<plugz){   //if in the barrel
	// detsurface=(lxstart-ecalrad)/tmpsin;
      // }else{  //if in plug
	// detsurface=(fabs(xStart[2])-plugz)/fabs(tmpcos);
      // }
      
      //float maxed=a*pow(b/c,b)*exp(-b);   //for simple fitting
      float maxed = a*c*gsl_sf_gammainv(b)*pow(b-1,b-1)*exp(-b+1);  //for advanced multiply with fabs(d) to avoid NaN
      float maxlength_pho=(1.0*std::log(clusterEnergy/(X0[0] * 0.021/Rm[0]))-0.5);  //this definition, +0.5 if gamma

      TVector3 clusdirection(CoG[0]-xStart[0],CoG[1]-xStart[1],CoG[2]-xStart[2]);


      //these variables are fit based variables
      // shapes.push_back(chi2);
      // shapes.push_back(maxed);
      // // shapes.push_back(((b-1.0)*X0[0]/c+xl0+detsurface)/(2.0*X0[0]));
      // shapes.push_back(1/fabs(d));
      // shapes.push_back(maxlength_pho);
      // shapes.push_back(((b-1.0)/c)/maxlength_pho);
      // shapes.push_back(Rm[0]*2.0);
      // // shapes.push_back(detsurface);
      // shapes.push_back(xl0);
      // shapes.push_back(a);
      // shapes.push_back(b);
      // shapes.push_back(c);
      // shapes.push_back(d);
      // //these variables are detector based variables
      // shapes.push_back(pClusterShapes->getEmax(xStart,index_xStart,X0,Rm));
      // shapes.push_back(pClusterShapes->getsmax(xStart,index_xStart,X0,Rm));
      shapes.push_back(pClusterShapes->getxt90(xStart,index_xStart,X0,Rm));
      // shapes.push_back(pClusterShapes->getxl20(xStart,index_xStart,X0,Rm));
      
      //20150708 add variables by Hale
      shapes.push_back(clusdirection.Mag());  // depth of the cluster
      shapes.push_back(pClusterShapes->getRhitMean(xStart,index_xStart,X0,Rm));  // mean of the radius of the hits wrt cog
      shapes.push_back(pClusterShapes->getRhitRMS(xStart,index_xStart,X0,Rm));   // RMS of the radius of the hits wrt cog
      shapes.push_back(nYokeHits);

      delete pClusterShapes;
      delete[] pHitE; delete[] pHitX; delete[] pHitY; delete[] pHitZ;

      return shapes;
}

