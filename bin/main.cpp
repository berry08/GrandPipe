#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <sys/types.h>
#include <unistd.h>
#include "gc.h"
#include "PipeGraph.h"
#include "create_pipe.h"
#include "check_graph_exists.h"
#include "pipeline.h"
#include "global_cmd_para.h"
#include "postprocess.h"
Graph pipe_graph;
int main(int argc,char* argv[]){
	if(argc!=5 && argc!=6){
		cerr<<endl;
		cerr<<"Program: GrandPipe (Tools for process pipelines)"<<endl;
		cerr<<"Version: 1.0"<<endl;
		cerr<<"Usage:   GrandPipe  <pipeline name>  <config file path>  <output directory>  <output prefix>  <other option>\n"<<endl;
		return -1;
	}

	string pipeline_name(argv[1]);
	string config_file_path(argv[2]);
	string output_dir_name(argv[3]);
	string output_prefix(argv[4]);
	string exe_path(argv[0]);
	string bin_dir=exe_path;
	if(exe_path.find("/")==string::npos)
		bin_dir=get_absolute_path(exe_path);
	bin_dir.erase(bin_dir.begin()+bin_dir.find_last_of("/"),bin_dir.end());
	cmd_para env_parameter(bin_dir,output_dir_name,output_prefix,config_file_path);

	if(!check_pipeGraph(output_dir_name,output_prefix))
		pipe_graph=createPipeGraph(bin_dir,pipeline_name,env_parameter);
	//print_cmd_para(bin_dir,output_dir_name,output_prefix);
	//pipe_graph.DFSTraverse2();
	pipe_graph.Traverse();
	post_process(env_parameter);
/*	while(!check_done())
		run_uncomplished();
*/	cout<<get_local_time()<<"\tAll analysis done!"<<endl;
	return 0;
}