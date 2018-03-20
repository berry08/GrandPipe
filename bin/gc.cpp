//#include<iostream>
//#include<algorithm>
#include <string>
#include <vector>
#include <sstream>
#include <cctype>
#include <set>
#include <time.h>
#include <sstream>
#include <unistd.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <sys/stat.h>
#include <algorithm>
#include <fstream>
//#include<stdio.h>

using namespace::std;

//char separator='\t';
void remove_space(string &a){
	string b;
	for(string::size_type ix=0;ix!=a.size();ix++){
		if(!isspace(a[ix])){
			b+=a[ix];
		}
	}
	a=b;
}
int file_exist_and_not_empty(string file_name){
	int return_value=1;
	ifstream if_input;
	if_input.open(file_name.c_str());
	if(!if_input){
		return_value=0;
		return return_value;
	}
	if(if_input.peek()==EOF){
		return_value=0;
		return return_value;
	}
	return return_value;
}
void uniq_vector(vector<string> &a){
	sort(a.begin(),a.end());
	vector<string>::iterator end_iter=unique(a.begin(),a.end());
	a.erase(end_iter,a.end());
}
string get_local_time(){
	time_t a;
	time(&a);
	struct tm* b=localtime(&a);
	int cur_year=b->tm_year+1900;
	int cur_mon=b->tm_mon+1;
	int cur_day=b->tm_mday;
	ostringstream out_time;
	out_time<<cur_year<<"-"<<cur_mon<<"-"<<cur_day<<"  "<<b->tm_hour<<":"<<b->tm_min<<":"<<b->tm_sec<<endl;
	string out=out_time.str();
	out.erase(out.end()-1);
	return out;
}
void line_split(string line_info,char sep,vector<string> &elements){
	elements.clear();
	string element;
	for(string::size_type ix=0;ix!=line_info.size();ix++){
		if(line_info[ix]!=sep){
			element+=line_info[ix];
		}else{
			elements.push_back(element);
			element="";
		}
	}
	elements.push_back(element);
}
void line_split(string line_info,vector<string> &elements){
	elements.clear();
	string element;
	for(string::size_type ix=0;ix!=line_info.size();ix++){
		if(!isspace(line_info[ix])){
			element+=line_info[ix];
		}else{
			elements.push_back(element);
			element="";
		}
	}
	elements.push_back(element);
}
void line_split(string line_info,char sep,set<string> &elements){
	elements.clear();
	string element;
	for(string::size_type ix=0;ix!=line_info.size();ix++){
		if(line_info[ix]!=sep){
			element+=line_info[ix];
		}else{
			elements.insert(element);
			element="";
		}
	}
	elements.insert(element);
}
int count_gc(string a){
	string::size_type t_index=0,t_index2;
	int gc_counts=0;
	while((t_index2=a.find_first_of("GC",t_index))!=string::npos){
		t_index=t_index2+1;
		gc_counts++;
	}
	return gc_counts;
}
void int2string(int &a,string &b){
        stringstream ss; 
        ss<<a;
        b=ss.str();
}
void float2string(float &a,string &b){
		stringstream ss;
		ss<<a;
		b=ss.str();
}
void double2string(double &a,string &b){
	stringstream ss;
	ss<<a;
	b=ss.str();
}
void string2double(string &a,double &b){
	stringstream ss;
	ss<<a;
	ss>>b;
}
void string2upper(string &a){
	for(string::iterator ix=a.begin();ix!=a.end();ix++){
		*ix=toupper(*ix);
	}
}
string join_vector(vector<string> a,char sep){
	string result;
	for(vector<string>::iterator ix=a.begin();ix!=a.end();ix++){
		if(ix!=a.end()-1){
			result+=*ix+sep;
		}else{
			result+=*ix;
		}
	}
	return result;
}
string join_vector(vector<int> a,char sep){
	string result;
	string tmp;
	for(vector<int>::iterator ix=a.begin();ix!=a.end();ix++){
		int2string(*ix,tmp);
		if(ix!=a.end()-1){
			result+=tmp+sep;
		}else{
			result+=tmp;
		}
	}
	return result;
}
string join_vector(vector<double> a,char sep){
	string result;
	for(vector<double>::iterator ix=a.begin();ix!=a.end();ix++){
		string tmp_a;
		double2string(*ix,tmp_a);
		if(ix!=a.end()-1){
			result+=tmp_a+sep;
		}else{
			result+=tmp_a;
		}
	}
	return result;
}
string join_vector(set<string> a,char sep){
	string result;
	for(set<string>::iterator ix=a.begin();ix!=a.end();ix++){
		result+=*ix+sep;
	}
	result.erase(result.size()-1,1);
	return result;
}
string join_vector(vector<string> a,string sep){
	string result;
	for(vector<string>::iterator ix=a.begin();ix!=a.end();ix++){
		if(ix!=a.end()-1){
			result+=*ix+sep;
		}else{
			result+=*ix;
		}
	}
	return result;
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
int check_bam_ok(string bin_path,string bam_file){
	ifstream if_pre_check(bam_file.c_str());
	if(!if_pre_check)
		return 0;
	if_pre_check.close();
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
int check_bai_ok(string bin_path,string bam_file){
	string check_cmd=bin_path+"/samtools tview -d T "+bam_file+" 2>&1";
	FILE* s_check=popen(check_cmd.c_str(),"r");
	char line[1000];
	while(fgets(line,1000,s_check)!=NULL){
		string s_line(line);
		if(s_line.find("Cannot read index for")!=string::npos){
			return 0;
		}
	}
	return 1;
}
string get_exe_path(string path){
	if(path.find("/")!=string::npos){
		return path.substr(0,path.find_last_of("/")+1);
	}else{
		string abs_path=get_absolute_path(path);
		return path.substr(0,abs_path.find_last_of("/")+1);
	}
}
void chomp_space(string& a,string type){
	if(type=="head"){
		int s_num=0;
		for(string::size_type ix=0;ix!=a.size();ix++){
			if(!isspace(a[ix])){
				break;
			}else{
				s_num++;
			}
		}
		a.erase(0,s_num);
	}else if(type=="tail"){
		string::size_type s_iter=0;;
		int s_num=0;
		for(string::size_type ix=a.size()-1;ix>=0;ix++){
			if(!isspace(a[ix])){
				s_iter=ix;
				break;
			}else{
				s_num++;
			}
		}
		a.erase(s_iter,s_num);
	}else{
		int s_num=0;
		for(string::size_type ix=0;ix!=a.size();ix++){
			if(!isspace(a[ix])){
				break;
			}else{
				s_num++;
			}
		}
		a.erase(0,s_num);
		string::size_type s_iter=0;;
		int s_num2=0;
		for(string::size_type ix=a.size()-1;ix>=0;ix++){
			if(!isspace(a[ix])){
				s_iter=ix;
				break;
			}else{
				s_num++;
			}
		}
		a.erase(s_iter,s_num2);
	}
}