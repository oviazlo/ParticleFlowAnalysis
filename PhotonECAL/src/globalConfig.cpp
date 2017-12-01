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
}

bool config::loadConfigBoostOptions(boost::program_options::variables_map &vm){
	return true;
}

/*===========================================================================*/
/*===============================[ globalConfig ]===============================*/
/*===========================================================================*/



