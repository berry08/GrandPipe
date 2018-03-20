#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <stdlib.h>
#include <math.h>
#include "gc.h"
using namespace::std;
double sigmod(double number);
string gender_decide(string bedcov_file){
	string gender;
	set<string> y_gene;
	y_gene.insert("AMELY");
	y_gene.insert("DDX3Y");
	y_gene.insert("EIF1AY");
	y_gene.insert("NLGN4Y");
	y_gene.insert("PCDH11Y");
	y_gene.insert("RPS4Y1");
	y_gene.insert("RPS4Y2");
	y_gene.insert("SRY");
	y_gene.insert("TBL1Y");
	y_gene.insert("TGIF2LY");
	y_gene.insert("USP9Y");
	y_gene.insert("ZFY");
	map<string,double> gene_line_count,gene_value;
	ifstream input;
	input.open(bedcov_file.c_str());
	if(!input){
		cerr<<"Error:can not open file\t"<<bedcov_file<<endl;
		exit(1);
	}
	string line;
	double region_size,total_count;
	while(getline(input,line)){
		vector<string> line_eles;
		line_eles.clear();
		line_split(line,'\t',line_eles);
		double start,end,count;
		string2double(line_eles[1],start);
		string2double(line_eles[2],end);
		string2double(line_eles[line_eles.size()-1],count);
		region_size+=end-start;
		total_count+=count;
	}
	input.close();
	double mean_depth=total_count/region_size;
	input.open(bedcov_file.c_str());
	while(getline(input,line)){
		vector<string> line_eles;
		line_eles.clear();
		line_split(line,'\t',line_eles);
		double start,end;
		string2double(line_eles[1],start);
		string2double(line_eles[2],end);
		double value;
		if(y_gene.find(line_eles[3])!=y_gene.end()){
			gene_line_count[line_eles[3]]++;
			string2double(line_eles[4],value);
			gene_value[line_eles[3]]+=value/(end-start);
		}
	}
	input.close();
	for(set<string>::iterator ix=y_gene.begin();ix!=y_gene.end();ix++){
		if(gene_line_count.find(*ix)==gene_line_count.end() || gene_value.find(*ix)==gene_value.end()){
			cerr<<"Error:no such gene"<<endl;
			exit(1);
		}
		gene_value[*ix]=gene_value[*ix]/gene_line_count[*ix];
		gene_value[*ix]=gene_value[*ix]/mean_depth;
	}
	//double weight[13]={-3.72765,0.989338,0.985697,1.0002,0.962339,0.993416,0.992819,0.997403,0.985218,0.992781,0.992141,0.992109,0.987863};
	//medexome weight	
	//double weight[13]={-15.1821,2.77203,3.5266,2.89889,3.16594,3.12798,3.22509,2.93238,2.24359,2.88721,2.57312,3.598,3.10347};
	//idtexome weigth
	double weight[13]={-13.3879,3.88845,2.77065,2.36206,4.34374,3.16145,3.15498,2.99938,4.78625,3.88096,8.84362,2.60803,3.04213};
	if(gene_value.size()!=12){
		cerr<<"Error:gene number error\t"<<gene_value.size()<<endl;
		exit(1);
	}
	int i=0;
	double result=weight[0];
	for(set<string>::iterator ix=y_gene.begin();ix!=y_gene.end();ix++){
		i++;
		result+=gene_value[*ix]*weight[i];
	}
	double final_value=sigmod(result);
	if(final_value>0.5){
		gender="Male";
	}else if(final_value<0.5){
		gender="Female";
	}else{
		gender="Unkonwn_gender";
	}
	string s_value;
	double2string(final_value,s_value);
	gender=gender+"("+s_value+")"; 
	return gender;
}
double sigmod(double number){
	double y_value=exp(number)/(1+exp(number));
	return y_value;
}
