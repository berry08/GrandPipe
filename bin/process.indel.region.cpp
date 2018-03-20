#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "gc.h"
#include <boost/lexical_cast.hpp>
using namespace::std;
using namespace::boost;
int main(int argc,char* argv[]){
	if(argc!=3){
		cerr<<"Arguments wrong!"<<endl;
		cerr<<argv[0]<<"\t<input vcf file>\t<output bed file>"<<endl;
		exit(1);
	}
	ifstream raw_vcf_file;
	ofstream out_bed_file;
	out_bed_file.open(argv[2]);
	raw_vcf_file.open(argv[1]);
	string vcf_line;
	while(getline(raw_vcf_file,vcf_line)){
		int start(0),end(0);
		int expand=0;
		if(vcf_line.find("#")==0){
			continue;
		}
		vector<string> line_eles;
		line_split(vcf_line,'\t',line_eles);
		if(line_eles[3].size()==1 && line_eles[4].size()==1){
			continue;
		}
		if(line_eles[3].find(",")!=string::npos){
			cerr<<"error in vcf file: ref should not contain comma"<<endl;
			exit(1);
		}else if(line_eles[4].find(",")!=string::npos){
			vector<string> alts;
			line_split(line_eles[4],',',alts);
			for(vector<string>::iterator ix=alts.begin();ix!=alts.end();ix++){
				if((*ix).size()>expand){
					expand=(*ix).size();
				}
			}
			if(line_eles[3].size()>expand){
				expand=line_eles[3].size();
			}
		}else{
			expand=line_eles[3].size()>line_eles[4].size()?line_eles[3].size():line_eles[4].size();
		}
		string chr=line_eles[0];
		start=atoi(line_eles[1].c_str())-2;
		end=atoi(line_eles[1].c_str())+expand+2;
		out_bed_file<<chr<<"\t"<<start<<"\t"<<end<<endl;
	}
	return 0;
}