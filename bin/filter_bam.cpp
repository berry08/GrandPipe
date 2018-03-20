#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fstream>
#include "gc.h"
using namespace::std;

int main(int argc,char* argv[]){
	if(argc!=3){
		cerr<<"wrong arguments!\n";
		cerr<<"Usage: "<<argv[0]<<"  <bam file>  <out sam file>"<<endl;
		cerr<<"If <out sam file> equal \"_STD_\", the result will output to stdout"<<endl;
		exit(1);
	}
	string find_exe_path=argv[0];
	find_exe_path="which "+find_exe_path;

	FILE *fp = popen(find_exe_path.c_str(), "r");
	char exe_whole_path[200];
	string whole_path;
	if(fgets(exe_whole_path,sizeof(exe_whole_path),fp)!=NULL){
		whole_path=exe_whole_path;
	}
	whole_path.erase(whole_path.size()-1);
	string tmp_whole_path=whole_path;
	string bin_dir;
	int iter=0;
	fclose(fp);
	while(1){
		if(tmp_whole_path[tmp_whole_path.size()-1]=='/'){
			iter++;
			tmp_whole_path.erase(tmp_whole_path.size()-1);
			if(iter==1){
				bin_dir=tmp_whole_path;
				break;
			}
		}else{
			tmp_whole_path.erase(tmp_whole_path.size()-1);
		}
	}
	string s_cmd=bin_dir;
	s_cmd=s_cmd+"/samtools view -h ";
	string bamfile=argv[1];
	s_cmd=s_cmd+bamfile;
	FILE* read_bam=popen(s_cmd.c_str(),"r");
	string out_s(argv[2]);
	ofstream out_sam;
	if(out_s.find("_STD_")==string::npos){
		out_sam.open(argv[2]);
	}
	char line_info[10000];
//	cout<<cmd<<endl;
	while(fgets(line_info,10000,read_bam)!=NULL){
		line_info[strlen(line_info)-1]=0;
		if(line_info[0]=='@'){
			if(out_s.find("_STD_")==string::npos){
				out_sam<<line_info<<endl;
			}else{
				cout<<line_info<<endl;
			}
			continue;
		}
		string line(line_info);
		vector<string> line_ele;
		line_split(line_info,'\t',line_ele);
		int flag=atoi(line_ele[1].c_str());
		int mapq=atoi(line_ele[4].c_str());
		if((flag&4)!=4 && mapq>=20){
			if(out_s.find("_STD_")==string::npos){
				out_sam<<line_info<<endl;
			}else{
				cout<<line_info<<endl;
			}
		}
	}
	pclose(read_bam);
	out_sam.close();
}
