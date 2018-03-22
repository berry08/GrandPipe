#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <stdlib.h>
#include <unistd.h>
#include "global_cmd_para.h"
#include "check.h"
#include "gc.h"
using namespace::std;

//int run_container(string fid,string root_dir,string output_dir,string output_prefix){
int check_container(string fid,cmd_para &env_para){
	/*string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
	string out_test=whole_out+".test.out";
	if(fid=="a0"){
		string dir=env_para.output_dir;
		string cmd="mkdir -p "+dir;
		string file_path=whole_out+".test.out";
		string cmd2="touch "+file_path;
		if(system(cmd.c_str())==-1)
			cerr<<"system error"<<endl;
		if(system(cmd2.c_str())==-1)
			cerr<<"system error"<<endl;
		return 1;
	}
	ifstream if_out(out_test.c_str());
	if(!if_out){
		cerr<<"Error:cannot open such file,"<<out_test<<endl;
		exit(1);
	}
	string tmp_line;
	int value=0;
	string q_str=fid+".";
	while(getline(if_out,tmp_line)){
		if(tmp_line.find(q_str)!=string::npos){
			value=1;
			break;
		}
	}
	return value;
	*/
	if(fid=="mkdir_p"){
		vector<string> v_raw_fq=env_para.m_config["fastq_file"];
		if(v_raw_fq.size()!=2){
			cerr<<"Error,now only support PE fastq"<<endl;
			exit(1);
		}
		string fq1,fq2;
		int merge1(0),merge2(0);
		for(vector<string>::iterator ix=v_raw_fq.begin();ix!=v_raw_fq.end();ix++){
			vector<string> line_eles;
			line_split(*ix,'=',line_eles);
			if((*ix).find("fq1")==0){
				fq1=line_eles[1];
			}else if((*ix).find("fq2")==0){
				fq2=line_eles[1];
			}else{
				cerr<<"Error,config format error"<<endl;
				exit(1);
			}
		}
		if(fq1.find(",")!=string::npos || fq2.find(",")!=string::npos){
			vector<string> v_fq1s,v_fq2s;
			line_split(fq1,',',v_fq1s);
			line_split(fq2,',',v_fq2s);
			if(v_fq1s.size()!=v_fq2s.size())
				cout<<"Warning! The files number of fq1 and fq2 are different!"<<endl;
			string merge_fq1_cmd=join_vector(v_fq1s,' ');
			string merge_fq2_cmd=join_vector(v_fq2s,' ');
			string merge_fq1,merge_fq2;
			if(v_fq1s[0].find(".gz")==v_fq1s[0].size()-3){
				merge_fq1_cmd="cat "+merge_fq1_cmd+" > "+env_para.output_dir+"/"+env_para.output_prefix+".merged.r1.fq.gz";
				merge_fq2_cmd="cat "+merge_fq2_cmd+" > "+env_para.output_dir+"/"+env_para.output_prefix+".merged.r2.fq.gz";
				merge_fq1=env_para.output_dir+"/"+env_para.output_prefix+".merged.r1.fq.gz";
				merge_fq2=env_para.output_dir+"/"+env_para.output_prefix+".merged.r2.fq.gz";
			}else{
				merge_fq1_cmd="cat "+merge_fq1_cmd+" > "+env_para.output_dir+"/"+env_para.output_prefix+".merged.r1.fq";
				merge_fq2_cmd="cat "+merge_fq2_cmd+" > "+env_para.output_dir+"/"+env_para.output_prefix+".merged.r2.fq";
				merge_fq1=env_para.output_dir+"/"+env_para.output_prefix+".merged.r1.fq";
				merge_fq2=env_para.output_dir+"/"+env_para.output_prefix+".merged.r2.fq";
			}
			ifstream if_m1,if_m2;
			if_m1.open(merge_fq1.c_str());
			if_m2.open(merge_fq2.c_str());
			if(if_m1 && if_m2){
			}else{
				return 0;
			}
			env_para.m_config2["fastq1"]=merge_fq1;
			env_para.m_config2["fastq2"]=merge_fq2;
		}else{
			env_para.m_config2["fastq1"]=fq1;
			env_para.m_config2["fastq2"]=fq2;
		}
		string bed_path;
		if(env_para.m_config.find("bed_file")!=env_para.m_config.end()){
			for(vector<string>::iterator ix=env_para.m_config["bed_file"].begin();ix!=env_para.m_config["bed_file"].end();ix++){
				if((*ix).find("bed")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					bed_path=k_v[1];
				}
			}
		}
		env_para.m_config2["bed"]=bed_path;
		string ref_fasta;
		if(env_para.m_config.find("reference")!=env_para.m_config.end()){
			for(vector<string>::iterator ix=env_para.m_config["reference"].begin();ix!=env_para.m_config["reference"].end();ix++){
				if((*ix).find("hg19_fasta")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					if(k_v[1].find("/")==string::npos){
						ref_fasta=env_para.bin_dir+"/../database/"+k_v[1];
					}else{
						ref_fasta=k_v[1];
					}
				}
			}
		}
		if(file_exist_and_not_empty(ref_fasta)<=0){
			cerr<<"Error:no such file,"<<ref_fasta<<endl;
			exit(1);
		}
		env_para.m_config2["ref_fasta"]=ref_fasta;
		string cmd="mkdir -p "+env_para.output_dir;
		if(system(cmd.c_str())==-1){
			cerr<<"Error:mkdir error"<<endl;
			exit(1);
		}
		//return 1;
	}
	string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
	string done_log=whole_out+".total.log";
	ifstream if_done_log;
	if_done_log.open(done_log.c_str());
	string check_str=fid+" done!";
	string tmp_line;
	int return_value=0;
	while(getline(if_done_log,tmp_line)){
		if(tmp_line==check_str){
			if_done_log.close();
			return 1;
		}
	}
	if_done_log.close();
	return return_value;
	
	
	/*
	if(fid=="trim"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string log_file=whole_out+".trim.e.log";
		ifstream if_log(log_file.c_str());
		string log_line;
		int flag=0;
		while(getline(if_log,log_line)){
			if(log_line.find("Completed successfully")!=string::npos){
				flag=1;
				break;
			}
		}
		if_log.close();
		return flag;
	}
	if(fid=="fastqc"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string log_file=whole_out+".fastqc.log";
		if(file_exist_and_not_empty(log_file)<=0){
			return 0;
		}
		ifstream if_log(log_file.c_str());
		string log_line;
		int iter=0;
		if(!if_log){
			cerr<<"Error:cannot open such file,"<<log_file<<endl;
			exit(1);
		}
		while(getline(if_log,log_line)){
			if(log_line.find("Analysis complete")!=string::npos){
				iter++;
			}
		}
		if(iter==2){
			return 1;
		}else{
			return 0;
		}
		if_log.close();
		
	}
	if(fid=="bwa"){
		string cmd1,cmd2;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_bam=whole_out+".sort.bam";
		string final_bam_bai=whole_out+".sort.bam.bai";
		if(check_bam_ok(env_para.bin_dir,final_bam)==1 && check_bai_ok(env_para.bin_dir,final_bam)==1){
			return 1;
		}else{
			return 0;
		}
	}
	if(fid=="rmdup"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_bam=whole_out+".sort.rmdup.bam";
		if(check_bam_ok(env_para.bin_dir,final_bam)){
			return 1;
		}else{
			return 0;
		}
	}
	if(fid=="filter_bam"){
		string cmd1,cmd2;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_bam=whole_out+".sort.rmdup.filter.bam";
		if(check_bam_ok(env_para.bin_dir,final_bam)==1 && check_bai_ok(env_para.bin_dir,final_bam)==1){
			return 1;
		}else{
			return 0;
		}
	}
	if(fid=="realign1"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".indelcreate.intervals";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="realign2"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_bam=whole_out+".filter.realigned.bam";
		string final_bed=env_para.output_dir+"/add_mt.bed";
		if(file_exist_and_not_empty(final_bed)>0 && check_bam_ok(env_para.bin_dir,final_bam)==1 && check_bai_ok(env_para.bin_dir,final_bam)==1){
			return 1;
		}else{
			return 0;
		}
	}
	if(fid=="bqsr1"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".recalibration.grp";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="bqsr2"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_bam=whole_out+".filter.realigned.recalibrate.bam";
		if(check_bam_ok(env_para.bin_dir,final_bam)==1 && check_bai_ok(env_para.bin_dir,final_bam)==1){
			return 1;
		}else{
			return 0;
		}
	}
	if(fid=="unifiedgenotyper"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".raw.vcf";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="process_indel_region"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".indel_region.bed";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="haplotype"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".hap.vcf";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="merge_vcf"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".merged.raw.vcf";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="vqsr1"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".tranches";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="vqsr2"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".recalibrated.filtered.vcf";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="hdfilter"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".hdfilter.vcf";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="target_reads"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".on_target_reads.txt";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="target_ratio"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".target_ratio";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="depthofcoverage"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+"_target_coverage.sample_summary";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="samtools_depth1"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".sort.rmdup.filter.bam.depth";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="samtools_depth2"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".sort.rmdup.filter.bam.target.depth";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="bedcov"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".filter.realigned.recalibrate.bam.bedcov";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="final_stat"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".final.stat";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="ts_cov"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".transcript.cov";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="stat_depth"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".depth.stat.result";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="copy_vcf"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".final.vcf";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="flagstat"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".aln.stat";
		return file_exist_and_not_empty(final_file);
	}
	if(fid=="flagstat2"){
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string final_file=whole_out+".flagstat";
		return file_exist_and_not_empty(final_file);
	}
	*/
}
