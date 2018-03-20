#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>
#include "gc.h"
using namespace::std;

int main(int argc,char* argv[]){
	if(argc!=4){
		cerr<<argv[0]<<"  <bed_file>  <samtools_depth_file>  <stat_output>\n";
		exit(0);
	}
	ifstream i_bed,i_depth;
	i_bed.open(argv[1]);
	i_depth.open(argv[2]);
/*	
	string bam_file=argv[2];
	bam_file.erase(bam_file.end()-6,bam_file.end());
	//cout<<bam_file<<endl;
	string root_dir=argv[0];
	root_dir.erase(root_dir.end()-14,root_dir.end());
	//string smn_bed=root_dir+"database/smn.medexome.bed";
	//cout<<root_dir<<endl;
	//cout<<smn_bed<<endl;
	//string samtools_path=root_dir+"bin/samtools";
	//string bedcov_cmd=samtools_path+" bedcov "+smn_bed+" "+bam_file+" >"+bam_file+".smn.bedcov";
	//system(bedcov_cmd.c_str());
	string final_stat_file;
	for(string::iterator ix=bam_file.begin();ix!=bam_file.end();ix++){
		if(*ix!='/'){
			final_stat_file+=*ix;
		}else{
			break;
		}
	}
	string sample_name=final_stat_file;
	sample_name.erase(sample_name.end()-7,sample_name.end());
	final_stat_file=final_stat_file+"/"+sample_name+".final.stat";
	ifstream final_stat_input;
	final_stat_input.open(final_stat_file.c_str());
	if(!final_stat_input){
		cerr<<"Error:unable to open file,"<<final_stat_file<<endl;
		exit(1);
	}
	string final_stat_file_line;
	string sample_depth;
	while(getline(final_stat_input,final_stat_file_line)){
		if(final_stat_file_line.find("depth")==0){
			vector<string> depth_eles;
			line_split(final_stat_file_line,'\t',depth_eles);
			sample_depth=depth_eles[1];
			break;
		}
	}
	string awk_cmd="awk '{a=$5/(($3-$2)*"+sample_depth+")}{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5\"\\t\"a\"\\t\"}' "+bam_file+".smn.bedcov >"+bam_file+".smn.bedcov2";
	system(awk_cmd.c_str());
	//cout<<awk_cmd<<endl;
	string smn_database=root_dir+"database/smn.medexome.database";
	string tmp_cmd="paste "+smn_database+" "+bam_file+".smn.bedcov2 "+"|awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$11\"\\t\"$4\"\\t\"$5\"\\t\"$6}' >"+bam_file+".smn.result";
	//cout<<tmp_cmd<<endl;
	system(tmp_cmd.c_str());
	string rm_string="rm "+bam_file+".smn.bedcov "+bam_file+".smn.bedcov2";
	system(rm_string.c_str());
*/
	//exit(1);
	string bed_info,depth_info;
	vector<string> bed_ele,depth_ele;
	int bed_size;
	int level0=0,level1=0,level2=0,level3=0,level4=0,level5=0,level6=0;
	while(getline(i_bed,bed_info)){
		line_split(bed_info,'\t',bed_ele);
		bed_size+=atoi(bed_ele[2].c_str())-atoi(bed_ele[1].c_str());
	}
	i_bed.close();
	while(getline(i_depth,depth_info)){
		line_split(depth_info,'\t',depth_ele);
		int pos_depth;
		pos_depth=(atoi(depth_ele[2].c_str())+9)/10;
		if(pos_depth>5){
			pos_depth=6;
		}
		switch(pos_depth){
			case 0:level0++;break;
			case 1:level1++;break;
			case 2:level2++;break;
			case 3:level3++;break;
			case 4:level4++;break;
			case 5:level5++;break;
			case 6:level6++;break;
			default:cout<<"error:\t"<<pos_depth<<endl;break;
		}
	}
	i_depth.close();
	level0=bed_size-level1-level2-level3-level4-level5-level6;
//	float level1_ratio,level2_ratio,level3_ratio,level4_ratio,level5_ratio,level6_ratio,level0_ratio;
	ofstream o_result;
	o_result.open(argv[3]);
	o_result<<"bed size:\t"<<bed_size<<endl;
	o_result<<"1-10X:\t"<<level1<<"\t"<<(float)level1/bed_size<<endl;
	o_result<<"11-20X:\t"<<level2<<"\t"<<(float)level2/bed_size<<endl;
	o_result<<"21-30X:\t"<<level3<<"\t"<<(float)level3/bed_size<<endl;
	o_result<<"31-40X:\t"<<level4<<"\t"<<(float)level4/bed_size<<endl;
	o_result<<"41-50X:\t"<<level5<<"\t"<<(float)level5/bed_size<<endl;
	o_result<<">50X:\t"<<level6<<"\t"<<(float)level6/bed_size<<endl;
	o_result<<"0X:\t"<<level0<<"\t"<<(float)level0/bed_size<<endl;
}