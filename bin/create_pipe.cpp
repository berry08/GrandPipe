#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "PipeGraph.h"
#include "pipeline.h"
#include "gc.h"
using namespace::std;

Graph createPipeGraph(string bin_dir,string pipeline_name,cmd_para aaa){
	map<string,string> m_name_filepath=create_m_name_filepath(bin_dir);
	if(m_name_filepath.find(pipeline_name)==m_name_filepath.end())
		cerr<<"Error:no such pipeline"<<endl;
	string pipeline_file=m_name_filepath[pipeline_name];
	ifstream if_pipe_file;
	if_pipe_file.open(pipeline_file.c_str());
	if(!if_pipe_file)
		cerr<<"Error:no such file,"<<pipeline_name<<endl;
	vector<string> steps;
	multimap<string,string> relations;
	set<string> s_steps;
	string pipe_file_line;
	int step_flag(0),relation_flag(0);
	int relation_num=0;
	while(getline(if_pipe_file,pipe_file_line)){
		if(pipe_file_line.find("#steps")==0){
			step_flag=1;
			relation_flag=0;
			continue;
		}
		if(pipe_file_line.find("#relations")==0){
			relation_flag=1;
			step_flag=0;
			continue;
		}
		if(step_flag==1){
			steps.push_back(pipe_file_line);
			s_steps.insert(pipe_file_line);
		}
		if(relation_flag==1){
			vector<string> line_eles;
			line_split(pipe_file_line,'=',line_eles);
			if(line_eles.size()!=2)
				cerr<<"Error,relation format error"<<endl;
			relations.insert(make_pair(line_eles[0],line_eles[1]));
			relation_num++;
		}
	}
	for(multimap<string,string>::iterator ix=relations.begin();ix!=relations.end();ix++){
		if(s_steps.find(ix->first)==s_steps.end() || s_steps.find(ix->second)==s_steps.end())
			cerr<<"Error:at least one step is unrecognized,"<<ix->first<<","<<ix->second<<endl;
	}
	Graph pipe_graph(steps,relations,relation_num,aaa);
	return pipe_graph;
}