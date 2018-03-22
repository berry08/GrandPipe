#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <stdlib.h>
#include <unistd.h>
#include "run.h"
#include "add_done_flag.h"
#include "wes_run.h"
#include "wgs_run.h"
#include "gc.h"
using namespace::std;

//int run_container(string fid,string root_dir,string output_dir,string output_prefix){

void run_container(string fid,cmd_para &env_para){
	/*string dir=env_para.output_dir;
	string cmd="mkdir -p "+dir;
	if(system(cmd.c_str())==-1)
		cerr<<"system error"<<endl;
	string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
	string out_test=whole_out+".test.out";
	ofstream of_out(out_test.c_str(),ios::app);
	if(!of_out){
		cerr<<"Error:cannot write to such file,"<<out_test<<endl;
		exit(1);
	}
	if(fid=="a0"){
		of_out<<"here is a0."<<endl;
	}
	if(fid=="a1"){
		of_out<<"here is a1."<<endl;
	}
	if(fid=="a2"){
		sleep(10);
		of_out<<"here is a2."<<endl;
		
	}
	if(fid=="a3"){
		sleep(40);
		of_out<<"here is a3."<<endl;
		
	}
	if(fid=="a4"){
		of_out<<"here is a4."<<endl;
	}
	if(fid=="a5"){
		of_out<<"here is a5."<<endl;
	}
	if(fid=="a6"){
		of_out<<"here is a6."<<endl;
	}
	if(fid=="a7"){
		of_out<<"here is a7."<<endl;
	}
	if(fid=="a8"){
		sleep(20);
		of_out<<"here is a8."<<endl;
		
	}
	if(fid=="a9"){
		of_out<<"here is a9."<<endl;
	}
	if(fid=="a10"){
		of_out<<"here is a10."<<endl;
	}
	if(fid=="a11"){
		sleep(10);
		of_out<<"here is a11."<<endl;
		
	}
	if(fid=="a12"){
		of_out<<"here is a12."<<endl;
	}
	if(fid=="a13"){
		of_out<<"here is a13."<<endl;
	}
	if(fid=="a14"){
		of_out<<"here is a14."<<endl;
	}
	if(fid=="a15"){
		sleep(10);
		of_out<<"here is a15."<<endl;
		
	}
	of_out.close();
	*/
	if(env_para.pipe_name=="wes"){
		run_wes_container(fid,env_para);		
	}
	if(env_para.pipe_name=="wgs"){
		run_wgs_container(fid,env_para);
	}

}
