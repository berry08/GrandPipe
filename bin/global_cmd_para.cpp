#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include "global_cmd_para.h"
#include "gc.h"
using namespace::std;

cmd_para::cmd_para(string a,string b,string c,string d){
	bin_dir=a;
	output_dir=b;
	output_prefix=c;
	ifstream if_config(d.c_str());
	string config_line;
	string step_name;
	while(getline(if_config,config_line)){
		if(config_line.find("#")==0){
			step_name=config_line;
			step_name.erase(0,1);
			continue;
		}
		m_config[step_name].push_back(config_line);
	}
	m_config2.clear();
//	for(map<string,vector<string> >::iterator ix=m_config.begin();ix!=m_config.end();ix++)
//		cout<<ix->first<<"\t"<<join_vector(ix->second,'\t')<<endl;
}
cmd_para::cmd_para(){
	bin_dir="";
	output_dir="";
	output_prefix="";
	m_config.clear();
	m_config2.clear();
}