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

int main(int argc,char* argv[]){
	if(argc!=3){
		cerr<<"Arguments error!"<<endl;
		cerr<<"Usage:\n\t"<<argv[0]<<"  <transcript coverage>  <final stat file>"<<endl;
		exit(1);
	}
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
	ifstream if_cov(argv[1]);
	ifstream if_final(argv[2]);
	string final_stat_line,ts_cov_line;
	while(getline(if_final,final_stat_line)){
		if(final_stat_line.find("Gender(no control)")==0)
			return 0;
		if(final_stat_line.find("depth")==0){
			vector<string> line_eles;
			line_split(final_stat_line,'\t',line_eles);
			double depth;
			string2double(line_eles[1],depth);
			if(depth<50)
				cerr<<"Error:depth is only "<<depth<<". It may make the gengder auto-detection not available"<<endl;
		}
	}
	if_final.close();
	int y_gene_to_detect_num=12;
	int y_gene_pass_num(0);
	while(getline(if_cov,ts_cov_line)){
		vector<string> line_eles;
		line_split(ts_cov_line,'\t',line_eles);
		if(y_gene.find(line_eles[0])!=y_gene.end()){
			double cov_10x(0);
			string2double(line_eles[7],cov_10x);
			if(cov_10x>0.8)
				y_gene_pass_num++;
		}
	}
	if_cov.close();
	if(y_gene_pass_num>y_gene_to_detect_num/5){
		gender="male";
	}else{
		gender="female";
	}
	ofstream of_final(argv[2],ofstream::app);
	of_final<<"Gender(no control):\t"<<gender<<endl;
	of_final.close();
	return 0;
}
