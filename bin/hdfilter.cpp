#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "gc.h"
#include <boost/lexical_cast.hpp>
using namespace::std;
using namespace::boost;

float snp_qdcutoff=2.0;
float snp_fscutoff=60.0;
float snp_mqcutoff=40.0;
float snp_mqranksumcutoff=-12.5;
float snp_readposranksum=-8.0;
float indel_qdcutoff=2.0;
float indel_fscutoff=200.0;
float indel_readposranksum=-20.0;
int main(int argc,char* argv[]){
	ifstream input_vcf;
	ofstream anno_vcf;
	input_vcf.open(argv[1]);
	anno_vcf.open(argv[2]);
	if(argc!=3){
		cerr<<"wrong arguments!"<<endl;
		cerr<<argv[0]<<"\tinput_vcf  filtered_vcf"<<endl;
		exit(1);
	}
	vector<string> line_eles;
	string line_info;
	while(getline(input_vcf,line_info)){
		if(line_info.find("#")==0){
			anno_vcf<<line_info<<endl;
		}else{
			float qd(3),fs(30),mq(50),mqranksum(0),readposranksum(0);
			line_split(line_info,'\t',line_eles);
			vector<string> para_value;
			line_split(line_eles[7],';',para_value);
			vector<string> flag;
			for(vector<string>::iterator ix=para_value.begin();ix!=para_value.end();ix++){
				vector<string> p_v;
				line_split(*ix,'=',p_v);
				if(p_v[0]=="QD"){
					qd=lexical_cast<float>(p_v[1]);
				}
				if(p_v[0]=="FS"){
					fs=lexical_cast<float>(p_v[1]);
				}
				if(p_v[0]=="MQ"){
					mq=lexical_cast<float>(p_v[1]);
				}
				if(p_v[0]=="MQRankSum"){
					mqranksum=lexical_cast<float>(p_v[1]);
				}
				if(p_v[0]=="ReadPosRankSum"){
					readposranksum=lexical_cast<float>(p_v[1]);
				}
			}
			if(line_eles[3].size()==1 && line_eles[4].size()==1){
				//snp filter
				if(qd<snp_qdcutoff){
					flag.push_back("QDfilter");
				}
				if(fs>snp_fscutoff){
					flag.push_back("FSfilter");
				}
				if(mq<snp_mqcutoff){
					flag.push_back("MQfilter");
				}
				if(mqranksum<snp_mqranksumcutoff){
					flag.push_back("MQRSfilter");
				}
				if(readposranksum<snp_readposranksum){
					flag.push_back("RPRSfilter");
				}
				if(flag.size()==0){
					flag.push_back("PASS");
				}
			}else{
				//indel filter
				if(qd<indel_qdcutoff){
					flag.push_back("QDfilter");
				}
				if(fs>indel_fscutoff){
					flag.push_back("FSfilter");
				}
				if(readposranksum<indel_readposranksum){
					flag.push_back("RPRSfilter");
				}
				if(flag.size()==0){
					flag.push_back("PASS");
				}
			}
			string filter_flag=join_vector(flag,",");
			line_eles[6]=filter_flag;
			string new_line_info=join_vector(line_eles,"\t");
			anno_vcf<<new_line_info<<endl;
		}
	}
	return 0;
}