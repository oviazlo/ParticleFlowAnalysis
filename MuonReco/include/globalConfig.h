/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/include/globalConfig.h
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	30th Nov 2017
 * 	Last Update:	30th Nov 2017
 */
#ifndef GLOBALCONFIG_H
#define GLOBALCONFIG_H

#include <boost/program_options.hpp>

namespace config
{
	extern std::map<unsigned int, std::string> pfoTypeIntStringMap;
	extern std::string some_config_string;
	extern int some_config_int;

	extern boost::program_options::variables_map vm;
	bool loadConfigBoostOptions(boost::program_options::variables_map &vm);
	// bool loadConfigFile();
}


#endif // GLOBALCONFIG_H
