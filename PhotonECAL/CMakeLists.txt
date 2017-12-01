INCLUDE_DIRECTORIES ( ${CMAKE_CURRENT_SOURCE_DIR}/include )

# ADD_EXECUTABLE ( PhotonECAL PhotonECAL.cpp src/objectFill.cpp src/particleFill.cpp src/energyFill.cpp src/globalConfig.cpp)
# ADD_EXECUTABLE ( PhotonEfficiency PhotonEfficiency.cpp src/objectFill.cpp src/particleFill.cpp src/energyFill.cpp src/photonEffCalculator.cpp src/globalConfig.cpp)
ADD_EXECUTABLE ( AnalyzeIsotropPhotonSample AnalyzeIsotropPhotonSample.cpp src/objectFill.cpp src/particleFill.cpp src/energyFill.cpp src/photonEffCalculator.cpp src/truthParticleSelector.cpp src/globalConfig.cpp)

# TARGET_LINK_LIBRARIES ( PhotonECAL ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${Boost_LIBRARIES} MyServiceFunctions BoostServiceFunctions)
# TARGET_LINK_LIBRARIES ( PhotonEfficiency ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${Boost_LIBRARIES} MyServiceFunctions BoostServiceFunctions)
TARGET_LINK_LIBRARIES ( AnalyzeIsotropPhotonSample ${ROOT_LIBRARIES} ${LCIO_LIBRARIES} ${Boost_LIBRARIES} MyServiceFunctions BoostServiceFunctions )
