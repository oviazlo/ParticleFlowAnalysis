/**@file	/afs/cern.ch/work/v/viazlo/analysis/PFAAnalysis/PhotonECAL/include/configStore.h
 * @author	viazlo
 * @version	800
 * @date
 * 	Created:	29th Nov 2017
 * 	Last Update:	29th Nov 2017
 */
#ifndef CONFIGSTORE_H
#define CONFIGSTORE_H

#include <boost/program_options.hpp>
using namespace std;

/*===========================================================================*/
/*======================[ Class configStore ]======================*/
/*===========================================================================*/
class configStore
{
/*-----------------------------[ Construction ]------------------------------*/
    /**@name
     */
    //@{
public:
	
    //@}

/*----------------------------[ Public methods ]-----------------------------*/
    /**@name
     */
    //@{
public:
    //@}
    
	static configStore *instance(){
		if (!s_instance)
			s_instance = new configStore;
		return s_instance;
	}
	static ConfigStore& get(){
		static ConfigStore instance;
		return instance;
	}
	// void parseFile(std::ifstream& inStream);
	void parseBoostOptions(boost::program_options::variables_map &vm);
	template<typename _T>
	_T getValue(std::string key);

/*---------------------------[ Internal methods ]----------------------------*/
    /**@name
     */
    //@{
protected:
    //@}

/*---------------------------------[ Data ]----------------------------------*/
    /**@name
     */
    //@{
private:
	ConfigStore(){};
	ConfigStore(const ConfigStore&);
	ConfigStore& operator=(const ConfigStore&);
	std::map<std::string,std::string> storedConfig;
    //@}

/*-----------------------------[ End of class ]------------------------------*/
};
#endif // CONFIGSTORE_H
