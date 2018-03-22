#include <iostream>
#include <fstream>
#include <string>
#include "global_cmd_para.h"
using namespace::std;
void add_done_flag(string fid,cmd_para &env_para){
	ofstream of_done_log;
	string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
	string done_log=whole_out+".total.log";
	of_done_log.open(done_log.c_str(),ios::app);
	of_done_log<<fid<<" done!"<<endl;
	of_done_log.close();
}