INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR}/include )

# ADD_EXECUTABLE ( PhotonECAL PhotonECAL.cpp src/objectFill.cpp src/particleFill.cpp src/energyFill.cpp src/globalConfig.cpp)
# ADD_EXECUTABLE ( PhotonEfficiency PhotonEfficiency.cpp src/objectFill.cpp src/particleFill.cpp src/energyFill.cpp src/photonEffCalculator.cpp src/globalConfig.cpp)
# ADD_EXECUTABLE ( AnalyzeIsotropPhotonSample AnalyzeIsotropPhotonSample.cpp src/objectFill.cpp src/particleFill.cpp src/energyFill.cpp src/photonEffCalculator.cpp src/truthParticleSelector.cpp src/globalConfig.cpp src/eventHistFiller.cpp src/truthCondition.cpp src/electronStudy.cpp)
ADD_EXECUTABLE ( StudyElectronPerformance StudyElectronPerformance.cpp src/objectFill.cpp src/particleFill.cpp src/energyFill.cpp src/truthParticleSelector.cpp src/globalConfig.cpp src/eventHistFiller.cpp src/eventHistFiller_initHists.cpp src/eventHistFiller_mainFunctions.cpp src/truthCondition.cpp src/electronStudy.cpp)
# ADD_EXECUTABLE ( testRead testRead.cpp src/truthCondition.cpp src/globalConfig.cpp)
ADD_EXECUTABLE ( JetPerformance JetPerformance.cpp src/objectFill.cpp src/globalConfig.cpp src/jetPfoStudy.cpp)

# TARGET_LINK_LIBRARIES ( PhotonECAL ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${Boost_LIBRARIES} MyServiceFunctions BoostServiceFunctions)
# TARGET_LINK_LIBRARIES ( PhotonEfficiency ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${Boost_LIBRARIES} MyServiceFunctions BoostServiceFunctions)
# TARGET_LINK_LIBRARIES ( AnalyzeIsotropPhotonSample ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${Boost_LIBRARIES} MyServiceFunctions BoostServiceFunctions )
TARGET_LINK_LIBRARIES ( StudyElectronPerformance ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${Boost_LIBRARIES} ${PandoraSDK_LIBRARIES} MyServiceFunctions BoostServiceFunctions )
# TARGET_LINK_LIBRARIES ( testRead ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${Boost_LIBRARIES} MyServiceFunctions BoostServiceFunctions )
TARGET_LINK_LIBRARIES ( JetPerformance ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${Boost_LIBRARIES} ${PandoraSDK_LIBRARIES} MyServiceFunctions BoostServiceFunctions )

