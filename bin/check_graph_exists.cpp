#include <iostream>
#include <string>
#include <fstream>
#include "gc.h"
using namespace::std;

int check_pipeGraph(string output_dir_name,string output_prefix){
	string check_file=output_dir_name+"/"+output_prefix+".pipeline.sh";
	if(file_exist_and_not_empty(check_file)){
		return 1;
	}else{
		return 0;
	}
}