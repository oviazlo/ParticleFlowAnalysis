INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR}/include )

ADD_EXECUTABLE ( MuonRecoMisID MuonRecoMisID.cpp src/objectFill.cpp src/particleFill.cpp src/globalConfig.cpp src/eventHistFiller.cpp src/eventHistFiller_initHists.cpp src/eventHistFiller_mainFunctions.cpp src/truthCondition.cpp src/muonClusterFiller.cpp src/muonClusterFiller_initHists.cpp)
ADD_EXECUTABLE ( TimingInfo TimingInfo.cpp src/objectFill.cpp src/globalConfig.cpp src/eventHistFiller.cpp src/eventHistFiller_initHists.cpp src/eventHistFiller_mainFunctions.cpp src/truthCondition.cpp src/trackerTiming.cpp src/trackerTiming_initHists.cpp)
ADD_EXECUTABLE ( AnalyzeBBbar AnalyzeBBbar.cpp src/objectFill.cpp src/globalConfig.cpp src/truthZWCondition.cpp src/muonEffInJets.cpp src/muonEffInJets_initHists.cpp)
ADD_EXECUTABLE ( OverlayAnalysis OverlayAnalysis.cpp src/objectFill.cpp src/energyFillAllCalo.cpp src/energyFill.cpp src/globalConfig.cpp src/eventHistFiller.cpp src/eventHistFiller_initHists.cpp src/eventHistFiller_mainFunctions.cpp src/truthCondition.cpp )

TARGET_LINK_LIBRARIES ( MuonRecoMisID ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${Boost_LIBRARIES} ${PandoraSDK_LIBRARIES} MyServiceFunctions BoostServiceFunctions )
TARGET_LINK_LIBRARIES ( TimingInfo ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${Boost_LIBRARIES} ${PandoraSDK_LIBRARIES} MyServiceFunctions BoostServiceFunctions )
TARGET_LINK_LIBRARIES ( AnalyzeBBbar ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${Boost_LIBRARIES} ${PandoraSDK_LIBRARIES} MyServiceFunctions BoostServiceFunctions )
TARGET_LINK_LIBRARIES ( OverlayAnalysis ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${Boost_LIBRARIES} ${PandoraSDK_LIBRARIES} MyServiceFunctions BoostServiceFunctions )
