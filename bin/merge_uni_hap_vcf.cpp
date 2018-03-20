#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <map>
#include <algorithm>
#include "gc.h"
using namespace::std;
int main(int argc,char* argv[]){
	ifstream uni_vcf,hap_vcf;
	ofstream merged_vcf;
	uni_vcf.open(argv[1]);
	hap_vcf.open(argv[2]);
	merged_vcf.open(argv[3]);
	string uni_line,hap_line;
	map<string,string> hap_info;
	vector<string> output_vcf_line,output_head;
	set<string> hap_head;
	while(getline(hap_vcf,hap_line)){
		if(hap_line.find("#")==0){
			if(hap_line.find("##INFO")==0){
				hap_head.insert(hap_line);
				continue;
			}else{
				continue;
			}
		}
		vector<string> line_eles;
		line_split(hap_line,'\t',line_eles);
		string key="";
		key=line_eles[0]+"\t"+line_eles[1]+"\t"+line_eles[3]+"\t"+line_eles[4];
		hap_info[key]=hap_line;
		output_vcf_line.push_back(hap_line);
	}
	int flag=0,flag2=0;
	while(getline(uni_vcf,uni_line)){
		if(uni_line.find("#")==0){
			if(uni_line.find("##INFO")==0){
				hap_head.insert(uni_line);
				flag=1;
			}else{
				if(flag==1 && flag2!=1){
					for(set<string>::iterator ix=hap_head.begin();ix!=hap_head.end();ix++){
						merged_vcf<<*ix<<endl;
						flag2=1;
					}
				}
				merged_vcf<<uni_line<<endl;
			}
		}else{
			vector<string> line_eles2;
			line_split(uni_line,'\t',line_eles2);
			string key="";
			key=line_eles2[0]+"\t"+line_eles2[1]+"\t"+line_eles2[3]+"\t"+line_eles2[4];
			if(hap_info.find(key)==hap_info.end()){
				output_vcf_line.push_back(uni_line);
			}
		}
	}
	multimap<double,string> value_line;
	vector<double> v_value;
//	vector<string> sorted_output_vcf_line;
	for(vector<string>::iterator ix=output_vcf_line.begin();ix!=output_vcf_line.end();ix++){
		vector<string> line_eles;
		line_eles.clear();
		line_split(*ix,'\t',line_eles);
		double chr=0;
		if(line_eles[0].find("X")!=string::npos){
			chr=23.0;
		}else if(line_eles[0].find("Y")!=string::npos){
			chr=24.0;
		}else if(line_eles[0]=="M" || line_eles[0]=="MT"){
			chr=25.0;
		}else{
			chr=atof(line_eles[0].c_str());
		}
		double value=chr*1000000000+atof(line_eles[1].c_str());
		value_line.insert(make_pair<double,string>(value,(*ix)));
		v_value.push_back(value);
	}
	unique(v_value.begin(),v_value.end());
	sort(v_value.begin(),v_value.end());
	double tmp=0;
	for(vector<double>::iterator ix=v_value.begin();ix!=v_value.end();ix++){
		if(*ix==tmp){
			continue;
		}
		tmp=*ix;
		multimap<double,string>::iterator lower=value_line.lower_bound(*ix);
		multimap<double,string>::iterator upper=value_line.upper_bound(*ix);
		if(lower==upper){
			merged_vcf<<lower->second<<endl;
		}else{
			while(lower!=upper){
				merged_vcf<<lower->second<<endl;
				++lower;
			}
		}
	}
	return 0;
}