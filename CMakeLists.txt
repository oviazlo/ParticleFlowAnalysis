CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

PROJECT( ECalAnalysis )

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

SET( PROJECT_VERSION_MAJOR 0 )
SET( PROJECT_VERSION_MINOR 1 )
#FIXME remove line below
SET( MarlinUtil_DIR "/cvmfs/clicdp.cern.ch/iLCSoft/builds/2019-04-17/x86_64-slc6-gcc62-opt/MarlinUtil/HEAD/")
#FIXME remove line below
SET( PandoraSDK_DIR "/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/modStandartPackages/PandoraPFANew_modTrackClusterAssCuts")

# FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )
# INCLUDE( ilcsoft_default_settings )

set(Boost_INCLUDE_DIR "/cvmfs/clicdp.cern.ch/software/Boost/1.62.0/x86_64-slc6-gcc62-opt/include/")
set(Boost_LIBRARY_DIR "/cvmfs/clicdp.cern.ch/software/Boost/1.62.0/x86_64-slc6-gcc62-opt/lib/")

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

INCLUDE_DIRECTORIES ( ${BOOST_ROOT}/include )
INCLUDE_DIRECTORIES ( ${PROJECT_SOURCE_DIR}/include )
INCLUDE_DIRECTORIES ( ${PROJECT_SOURCE_DIR}/../MyServiceFunctions/include )
INCLUDE_DIRECTORIES ( ${PROJECT_SOURCE_DIR}/../BoostServiceFunctions/include )

FIND_PACKAGE( MarlinUtil 1.4 REQUIRED )
FIND_PACKAGE( Marlin 1.0 REQUIRED )
INCLUDE_DIRECTORIES( SYSTEM ${Marlin_INCLUDE_DIRS} )
LINK_LIBRARIES( ${Marlin_LIBRARIES} )
ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )

#FIXME uncomment line below
# FIND_PACKAGE( GSL REQUIRED )

find_package( DD4hep REQUIRED COMPONENTS DDRec )
set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}  ${DD4hep_ROOT}/cmake )
include( DD4hep )
include(DD4hep_XML_setup)
INCLUDE_DIRECTORIES( SYSTEM ${DD4hep_INCLUDE_DIRS} )
LINK_LIBRARIES( ${DD4hep_LIBRARIES} ${DD4hep_COMPONENT_LIBRARIES} )

#FIXME uncomment line below
# FOREACH( pkg LCIO PandoraSDK DD4hep Marlin MarlinUtil GSL)
#FIXME remove line below
FOREACH( pkg LCIO PandoraSDK DD4hep Marlin MarlinUtil)
    FIND_PACKAGE (${pkg} REQUIRED )
    IF( ${pkg}_FOUND )
	INCLUDE_DIRECTORIES( SYSTEM ${${pkg}_INCLUDE_DIRS} )
	LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
	ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
	# ADD_DEFINITIONS( -DUSE_${pkg} )
    ENDIF()
ENDFOREACH()


LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/../MyServiceFunctions/build/)
LINK_DIRECTORIES(${PROJECT_SOURCE_DIR}/../BoostServiceFunctions/build/)

SET( CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR} )

# SET(GCC_COVERAGE_COMPILE_FLAGS "-O3 -Wall -Wextra -Wshadow -ansi -Wno-long-long -Wuninitialized -fopenmp -flto -Weffc++ --std=c++11 -ggdb" )
# SET(GCC_COVERAGE_LINK_FLAGS    "-flto -O3 -ggdb" )
# SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}" )
# SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}" )

# LIST ( APPEND CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wshadow -ansi -Wno-long-long -Wuninitialized -fopenmp -flto -Weffc++ --std=c++11 -ggdb" )
LIST ( APPEND CMAKE_CXX_FLAGS "-O3 -ansi -fopenmp -flto --std=c++11 -ggdb" )
LIST ( APPEND CMAKE_LD_FLAGS "-flto -O3 -ggdb" )


find_package ( ROOT REQUIRED )
include(${ROOT_USE_FILE})

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
# find_package(Boost 1.45.0 COMPONENTS program_options)
find_package(Boost 1.62.0 COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})

ADD_SUBDIRECTORY ( ${PROJECT_SOURCE_DIR}/PhotonECAL )
ADD_SUBDIRECTORY ( ${PROJECT_SOURCE_DIR}/MuonReco )

MESSAGE (STATUS "*** Build Type ${CMAKE_BUILD_TYPE} ***" )
MESSAGE (STATUS "*** Compiler Flags:  ${CMAKE_CXX_FLAGS}"  )
MESSAGE (STATUS "*** Compiler Flags:  ${CMAKE_CXX_FLAGS_DEBUG}"  )

