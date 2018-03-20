#include <string>
#include <map>
#include "pipeline.h"
using namespace::std;
map<string,string> create_m_name_filepath(string bin_dir){
	map<string,string> m_name_filepath;
	string pipeline_file_dir=bin_dir+"/../pipe/";
	string wes_pipe_path=pipeline_file_dir+"wes.pipe";
	m_name_filepath["WES"]=wes_pipe_path;
	return m_name_filepath;
}