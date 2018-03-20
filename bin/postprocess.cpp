#include <string>
#include <stdlib.h>
#include <iostream>
#include "global_cmd_para.h"
#include "gc.h"
using namespace::std;
void post_process(cmd_para env_parameter){
	string toberm_file1=env_parameter.output_dir+"/trim.log";
	if(file_exist_and_not_empty(toberm_file1)>0){
		string rmcmd1="rm "+toberm_file1;
		if(system(rmcmd1.c_str())==-1){
			cerr<<"Error:cannot remove such file,"<<toberm_file1;
			exit(1);
		}
	}
}