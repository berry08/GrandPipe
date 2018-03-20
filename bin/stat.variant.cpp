#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include "/home/gongchun/bin/gc.h"
using namespace::std;

int main(int argc,char* argv[]){
	if(argc!=4){
		cerr<<argv[0]<<"  <anno.result>  <output_result>  <type>\n";
		cerr<<"type:\"snp\" or \"indel\"\n";
		exit(0);
	}
	char* anno_vcf=argv[1];
	ifstream input_file;
	input_file.open(anno_vcf);
	ofstream output_file;
	output_file.open(argv[2]);
	output_file<<anno_vcf<<"_"<<argv[3]<<"\n";
	char* c_type=argv[3];
	string type(c_type);
	vector<string> line_ele;
	string line_info;
	int total_variant=0,dbsnp_number=0,novel_number=0,hom_number=0,het_number=0,ti_number=0,tv_number=0,ti_dbsnp=0,tv_dbsnp=0,ti_novel=0,tv_novel=0;
	int ag(0),ac(0),at(0),ga(0),gc(0),gt(0),ca(0),cg(0),ct(0),ta(0),tg(0),tc(0);
	map<string,int> map_function;
	map<string,int> map_exonic_func,map_exonic_func_stat;
	vector<int> indel_size;
	while(getline(input_file,line_info)){
		if(line_info.find("Chr")==0){
			continue;
		}
		line_info.erase(line_info.end()-1);
		line_split(line_info,'\t',line_ele);
		string ref=line_ele[3];
		string alt=line_ele[4];
		string function_stat=line_ele[5];
		string exonic_func_stat=line_ele[8];
		if(function_stat.find("exonic")==0){
			if(map_exonic_func_stat.count(exonic_func_stat)){
				map_exonic_func_stat[exonic_func_stat]++;
			}else{
				map_exonic_func_stat.insert(make_pair(exonic_func_stat,1));
			}
		}else if(function_stat.find("splicing")!=string::npos){
			if(map_exonic_func_stat.count("splicing")){
				map_exonic_func_stat["splicing"]++;
			}else{
				map_exonic_func_stat.insert(make_pair("splicing",1));
			}
		}else{
		}
		if(ref.size()==1 && ref.find("-")==string::npos && alt.size()==1 && alt.find("-")==string::npos){
			if(type=="snp"){
				total_variant++;
			}else if(type=="indel"){
				continue;
			}else{
				cerr<<"unrecognized type:"<<type<<". The type should be \"snp\" or \"indel\""<<endl;
				exit(0);
			}
		}else{
			if(type=="snp"){
				continue;
			}else if(type=="indel"){
				total_variant++;
			}else{
				cerr<<"unrecognized type:"<<type<<". The type should be \"snp\" or \"indel\""<<endl;
				exit(0);
			}
		}
		string dbsnp=line_ele[30];
		if(type=="indel"){
			goto label_indel;
		}
		{string combine_ref_alt=ref+alt;
		if(combine_ref_alt.find("AG")!=string::npos || combine_ref_alt.find("GA")!=string::npos || combine_ref_alt.find("CT")!=string::npos ||  combine_ref_alt.find("TC")!=string::npos){
			ti_number++;
			if(dbsnp.find("rs")==0){
				ti_dbsnp++;
			}else{
				ti_novel++;
			}
		}else{
			tv_number++;
			if(dbsnp.find("rs")==0){
				tv_dbsnp++;
			}else{
				tv_novel++;
			}
		}
		if(ref=="A" && alt=="G"){
			ag++;
		}else if(ref=="A" && alt=="C"){
			ac++;
		}else if(ref=="A" && alt=="T"){
			at++;
		}else if(ref=="G" && alt=="A"){
			ga++;
		}else if(ref=="G" && alt=="C"){
			gc++;
		}else if(ref=="G" && alt=="T"){
			gt++;
		}else if(ref=="C" && alt=="A"){
			ca++;
		}else if(ref=="C" && alt=="G"){
			cg++;
		}else if(ref=="C" && alt=="T"){
			ct++;
		}else if(ref=="T" && alt=="A"){
			ta++;
		}else if(ref=="T" && alt=="G"){
			tg++;
		}else if(ref=="T" && alt=="C"){
			tc++;
		}else{
			cerr<<"unknown base:\t"<<ref<<"\t"<<alt<<endl;
			exit(0);
		}
		}
		label_indel:
		if(ref=="-"){
			indel_size.push_back(alt.size()-ref.size()+1);
		}else if(alt=="-"){
			indel_size.push_back(alt.size()-ref.size()-1);
		}else{
			if(type=="snp"){
			}else{
				cerr<<"unknown tpye:\t"<<ref<<"\t"<<alt<<endl;
				exit(0);
			}
		}
		string function=line_ele[5];	
		string exonic_func=line_ele[8];
		if(dbsnp.find("rs")!=string::npos){
			dbsnp_number++;
		}else{
			novel_number++;
		}
		string info=line_ele[line_ele.size()-1];
		string genotype;
		for(string::iterator ix=info.begin();ix!=info.end();ix++){
			if(*ix!=':'){
				genotype+=*ix;
			}else{
				break;
			}
		}
		if(genotype=="0/1" || genotype=="1/2"){
			het_number++;
		}else if(genotype=="1/1"){
			hom_number++;
		}else{
			cerr<<"unrecognized genotype:\t";
			cerr<<genotype<<endl;
			cerr<<line_info<<endl;
			//exit(0);

		}
		if(map_function.count(function)){
			map_function[function]++;
		}else{
			map_function.insert(make_pair(function,1));
		}
		if(function.find("exonic")==0){
			if(map_exonic_func.count(exonic_func)){
				map_exonic_func[exonic_func]++;
			}else{
				map_exonic_func.insert(make_pair(exonic_func,1));
			}
		}
	}
	for(map<string,int>::iterator ix=map_exonic_func_stat.begin();ix!=map_exonic_func_stat.end();ix++){
		cout<<ix->first<<":\t"<<ix->second<<endl;
	}
	output_file<<"total "<<type<<":\t"<<total_variant<<"\n";
	output_file<<"dbsnp number:\t"<<dbsnp_number<<"\n";
	output_file<<"novel number:\t"<<novel_number<<"\n";
	output_file<<"hom number:\t"<<hom_number<<"\n";
	output_file<<"het number:\t"<<het_number<<"\n";
	for(map<string,int>::iterator ix=map_function.begin();ix!=map_function.end();ix++){
		output_file<<ix->first<<":\t"<<ix->second<<"\n";
	}
	for(map<string,int>::iterator ix=map_exonic_func.begin();ix!=map_exonic_func.end();ix++){
		output_file<<ix->first<<":\t"<<ix->second<<"\n";
	}
	if(type=="snp"){
		output_file<<"Ti number:\t"<<ti_number<<"\n";
		output_file<<"Tv number:\t"<<tv_number<<"\n";
		output_file<<"Ti/Tv:\t"<<(float)ti_number/tv_number<<"\n";
		output_file<<"dbsnp Ti/Tv:\t"<<(float)ti_dbsnp/tv_dbsnp<<"\n";
		output_file<<"novel Ti/Tv:\t"<<(float)ti_novel/tv_novel<<endl;
		output_file<<"AG:\t"<<ag<<"\n";
		output_file<<"AC:\t"<<ac<<"\n";
		output_file<<"AT:\t"<<at<<"\n";
		output_file<<"GA:\t"<<ga<<"\n";
		output_file<<"GC:\t"<<gc<<"\n";
		output_file<<"GT:\t"<<gt<<"\n";
		output_file<<"CA:\t"<<ca<<"\n";
		output_file<<"CG:\t"<<cg<<"\n";
		output_file<<"CT:\t"<<ct<<"\n";
		output_file<<"TA:\t"<<ta<<"\n";
		output_file<<"TG:\t"<<tg<<"\n";
		output_file<<"TC:\t"<<tc<<"\n";
	}else if(type=="indel"){
		ofstream output_file2;
		char anno_vcf_indel_size[100];
		strcpy(anno_vcf_indel_size,anno_vcf);
		strcat(anno_vcf_indel_size,".indel_size");
		output_file2.open(anno_vcf_indel_size);
		for(vector<int>::iterator ix=indel_size.begin();ix!=indel_size.end();ix++){
			output_file2<<*ix<<endl;
		}
		output_file2.close();
		ofstream output_rscript;
		char rscript[100];
		strcpy(rscript,anno_vcf);
		strcat(rscript,".R");
		output_rscript.open(rscript);
		output_rscript<<"jpeg(file=\""<<anno_vcf<<".jpg\",width=600,height=400);\n";
		output_rscript<<"a=read.table(\""<<anno_vcf_indel_size<<"\");\n";
		output_rscript<<"hist(a[,1],col=\"green\",xlim=c(-20,20),xlab=\"length\",ylab=\"count\",main=\"Insertion and deletion length distribution\",breaks=seq(-100,100,1));\n";
		output_rscript<<"dev.off();\n";
	}else{
		cerr<<"unrecognized type:\t"<<type<<endl;
		exit(0);
	}
	return 0;
}
