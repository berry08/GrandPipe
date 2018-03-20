#ifndef _GLOBAL_CMD_PARA_H
#define _GLOBAL_CMD_PARA_H

#include <string>
#include <map>
#include <vector>
using namespace::std;

class cmd_para{
public:
	cmd_para();
	cmd_para(string a,string b,string c,string d);
	string bin_dir;
	string output_dir;
	string output_prefix;
	map<string,vector<string> > m_config;
	map<string,string> m_config2;
};
#endif