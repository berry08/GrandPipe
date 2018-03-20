#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <unistd.h>
#include <sys/stat.h>
#include "gc.h"
#include "callVCFfromBam.h"
using namespace::std;

string link_dir_file(string dir,string file);
string get_absolute_path(string path);
void get_mt_bam(string raw_bam_path,string mt_bam_path,string root_path);
string call_mt_snp(string mt_bam_path);
int check_bam_ok(string bin_path,string bam_file);
string root_dir="";
string bin_dir="";
int main(int argc,char* argv[]){
	if(argc!=2){
		cerr<<"Arguments error!"<<endl;
		cerr<<"This program is used to add MT analysis based on wes, including callVCF and annotation"<<endl;
		cerr<<"Usage:"<<endl;
		cerr<<"\t"<<argv[0]<<"\t<sample name>"<<endl;
		cerr<<"\tThis program should be run at the grandmule directory"<<endl;
		cerr<<"\tThe annotation output will be generated at the grandanno directory with \"addMT\" keywords"<<endl;
		exit(1);
	}
	char cwd[200];
	if(getcwd(cwd,200)==NULL){
		cerr<<"ERROR:get cwd error"<<endl;
		exit(1);
	}
	string wd(cwd);
	string sample_name(argv[1]);
	string bam_path=wd+"/"+sample_name+"_result/"+sample_name+".sort.rmdup.filter.bam";
	string program_path(argv[0]);
	string total_path=program_path;
	if(program_path.find("/")!=0){
		total_path=get_absolute_path(program_path);
	}
	string::size_type last_slash=total_path.find_last_of("/");
	bin_dir=total_path.substr(0,last_slash);
	root_dir=total_path.substr(0,bin_dir.find_last_of("/"));
	int flag_check_bam=check_bam_ok(bin_dir,bam_path);
	if(flag_check_bam<0){
		cerr<<"ERROR:bam file not exists or malformed"<<endl;
		exit(1);
	}
	string mt_bam_path=wd+"/"+sample_name+"_result/"+sample_name+".sort.rmdup.filter.mt.bam";
	int flag_check_mt_bam=check_bam_ok(bin_dir,mt_bam_path);
	if(flag_check_mt_bam<0){
		get_mt_bam(bam_path,mt_bam_path,bin_dir);
	}
	string vcf_file=wd+"/"+sample_name+"_result/"+sample_name+".mt.vcf";
	callVCF(root_dir,sample_name,mt_bam_path,vcf_file);
	string anno_exe=root_dir+"/../grandanno/bin/mt_anno";
	string anno_config=root_dir+"/../grandanno/bin/mt.config.ini";
	string anno_out_prefix=wd+"/../grandanno/"+sample_name+".addMT";
	string anno_cmd=anno_exe+" "+vcf_file+" "+anno_out_prefix+" "+anno_config;
/*	if(system(anno_cmd.c_str())==-1){
		cerr<<"ERROR:anno MT error"<<endl;
		exit(1);
	}
*/	return 0;
}
int check_bam_ok(string bin_path,string bam_file){
	string check_cmd=bin_path+"/samtools quickcheck "+bam_file+" && echo 'all ok' || echo 'fail!'";
	FILE* s_check=popen(check_cmd.c_str(),"r");
	char line[100];
	while(fgets(line,100,s_check)!=NULL){
		string s_line(line);
		if(s_line.find("ok")!=string::npos){
			return 1;
		}
		if(s_line.find("fail")!=string::npos){
			return -1;
		}
	}
	cerr<<"ERROR:bam check error"<<endl;
	exit(1);
}
void get_mt_bam(string bam_path,string mt_bam_path,string bin_path){
	string samtools_path=bin_path+"/samtools";
	string s_cmd=samtools_path+" view -h "+bam_path;
	string mt_sam_path=mt_bam_path;
	mt_sam_path[mt_sam_path.size()-3]='s';
	ofstream of_mt_sam(mt_sam_path.c_str());
	FILE* s_view=popen(s_cmd.c_str(),"r");
	char line[10000];
	int iter=0;
//	string start_time=get_local_time();
//	cout<<start_time<<endl;
	while(fgets(line,10000,s_view)!=NULL){
		line[strlen(line)-1]=0;
/*		iter++;
		if(iter%1000000==0){
			cout<<"iter:"<<iter<<"\t"<<line<<endl;
			cout<<get_local_time()<<endl;
		}
*/		string s_line(line);
		if(s_line.find("@")==0){
			if(s_line.find("SN:")!=string::npos){
				if(s_line.find("SN:MT")!=string::npos){
					of_mt_sam<<s_line<<endl;
				}
			}else{
				of_mt_sam<<s_line<<endl;
			}
		}else{
/*			vector<string> line_eles;
			line_split(s_line,'\t',line_eles);
			if(line_eles[2]=="MT"){
				of_mt_sam<<s_line<<endl;
			}
*/
			int tab_num=0;
			string chr="";
			for(string::size_type ix=0;ix!=s_line.size();ix++){
				if(s_line[ix]=='\t'){
					tab_num++;
					continue;
				}
				if(tab_num==2){
					chr+=s_line[ix];
				}
				if(tab_num>2){
					break;
				}
			}
			if(chr=="MT"){
				of_mt_sam<<s_line<<endl;
			}
		}
	}
	string samtobam_cmd=samtools_path+" view -b "+mt_sam_path+" > "+mt_bam_path;
	pclose(s_view);
	of_mt_sam.close();
	if(system(samtobam_cmd.c_str())==-1){
		cerr<<"ERROR:sam transfer to bam error"<<endl;
		exit(1);
	}
}
string get_absolute_path(string path){
	char* _path=getenv("PATH");
	string s_path(_path);
	vector<string> path_eles;
	line_split(s_path,':',path_eles);
	struct stat buf;
	for(vector<string>::iterator ix=path_eles.begin();ix!=path_eles.end();ix++){
		string process_path=link_dir_file(*ix,path);
		if(lstat(process_path.c_str(),&buf)<0){
			continue;
		}
		if(S_ISREG(buf.st_mode) && S_IXUSR&(buf.st_mode)){
			return process_path;
		}
	}
	cerr<<"no such file: "<<path<<endl;
	exit(1);
}
string link_dir_file(string dir,string file){
	string::size_type index;
	// remove null char in the start or end;
	dir.erase(0,dir.find_first_not_of(" "));  
	dir.erase(dir.find_last_not_of(" ")+1);
/*	while((index=dir.find(" ",index))!=string::npos){       //remove all null char
		dir.erase(index,1);
	}
*/
	if(dir[dir.size()-1]!='/'){
		dir=dir+"/";
	}
	if(file.find("/")!=string::npos){
		string::size_type ix=file.find("/");
		file.erase(0,ix+1);
	}
	string out=dir+file;
	return out;
}
