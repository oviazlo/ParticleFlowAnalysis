CMAKE_MINIMUM_REQUIRED( VERSION 2.8 )

PROJECT( ECalAnalysis )

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

SET( PROJECT_VERSION_MAJOR 0 )
SET( PROJECT_VERSION_MINOR 1 )

#INCLUDE_DIRECTORIES ( ${BOOST_ROOT}/include )
INCLUDE_DIRECTORIES ( ${PROJECT_SOURCE_DIR}/include )
INCLUDE_DIRECTORIES ( ${PROJECT_SOURCE_DIR}/../MyServiceFunctions/include )
INCLUDE_DIRECTORIES ( ${PROJECT_SOURCE_DIR}/../BoostServiceFunctions/include )

FOREACH( pkg LCIO )
    FIND_PACKAGE (${pkg} REQUIRED )
    IF( ${pkg}_FOUND )
	INCLUDE_DIRECTORIES( SYSTEM ${${pkg}_INCLUDE_DIRS} )
#	LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
	ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
	ADD_DEFINITIONS( -DUSE_${pkg} )
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
find_package(Boost 1.45.0 COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})

ADD_SUBDIRECTORY ( ${PROJECT_SOURCE_DIR}/PhotonECAL )

MESSAGE (STATUS "*** Build Type ${CMAKE_BUILD_TYPE} ***" )
MESSAGE (STATUS "*** Compiler Flags:  ${CMAKE_CXX_FLAGS}"  )
MESSAGE (STATUS "*** Compiler Flags:  ${CMAKE_CXX_FLAGS_DEBUG}"  )

