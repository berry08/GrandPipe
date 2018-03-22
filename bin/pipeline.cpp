#include <string>
#include <fstream>
#include <map>
#include <iostream>
#include <stdlib.h>
#include "pipeline.h"
using namespace::std;
map<string,string> create_m_name_filepath(string bin_dir,string pipe_name){
	map<string,string> m_name_filepath;
	string pipeline_file_dir=bin_dir+"/../pipe/";
	string new_pipe_name;
	for(string::size_type ix=0;ix!=pipe_name.size();ix++){
		new_pipe_name.insert(new_pipe_name.end(),tolower(pipe_name[ix]));
	}
	string pipe_path=pipeline_file_dir+new_pipe_name+".pipe";
	ifstream if_pipe(pipe_path.c_str());
	if(!if_pipe){
		cerr<<"Error:no such pipe file,"<<pipe_path<<endl;
		exit(1);
	}
	if_pipe.close();
	m_name_filepath[pipe_name]=pipe_path;
	return m_name_filepath;
}