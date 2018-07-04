/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/src/globalConfig.cpp
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	30th Nov 2017
 * 	Last Update:	30th Nov 2017
 */

/*===========================================================================*/
/*===============================[ globalConfig ]===============================*/
/*===========================================================================*/

#include "globalConfig.h"
namespace config
{
	int some_config_int = 123;
	std::string some_config_string = "foo";
	boost::program_options::variables_map vm;
	std::map<unsigned int, std::string> pfoTypeIntStringMap = {{11,"Electron"}, {13,"Muon"},{22,"Photon"},{211,"Pion"},{2112,"NeutralHadron"}};
}

bool config::loadConfigBoostOptions(boost::program_options::variables_map &vm){
	return true;
}

/*===========================================================================*/
/*===============================[ globalConfig ]===============================*/
/*===========================================================================*/



