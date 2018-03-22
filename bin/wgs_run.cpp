#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <stdlib.h>
#include <unistd.h>
#include "wgs_run.h"
#include "add_done_flag.h"
#include "gc.h"
using namespace::std;
void run_wgs_container(string fid,cmd_para &env_para){
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
			merge_fq1=env_para.output_dir+"/"+env_para.output_prefix+".merged.r1.fq.gz";
			merge_fq2=env_para.output_dir+"/"+env_para.output_prefix+".merged.r2.fq.gz";
		}else{
			merge_fq1=env_para.output_dir+"/"+env_para.output_prefix+".merged.r1.fq";
			merge_fq2=env_para.output_dir+"/"+env_para.output_prefix+".merged.r2.fq";
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
	if(fid=="mkdir_p"){
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
			if(system(merge_fq1_cmd.c_str())==-1){
				cerr<<"Error:merge fq error"<<endl;
				exit(1);
			}
			if(system(merge_fq2_cmd.c_str())==-1){
				cerr<<"Error:merge fq error"<<endl;
				exit(1);
			}
		}
	}
	if(fid=="trim"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string jar_path=env_para.bin_dir+"/../software/trimmomatic.jar";
		string fq1=env_para.m_config2["fastq1"];
		string fq2=env_para.m_config2["fastq2"];
		if(env_para.m_config.find(fid)!=env_para.m_config.end()){
			string threads,phred_value,other_para_str;
			vector<string> v_config=env_para.m_config[fid];
			for(vector<string>::iterator ix=v_config.begin();ix!=v_config.end();ix++){
				if((*ix).find("threads")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					threads=k_v[1];
				}else if((*ix).find("phred")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					phred_value=k_v[1];
				}else{
					other_para_str=*ix;
				}
			}
			cmd="java -Xmx1750m -jar "+jar_path+" PE -threads "+threads+" -phred"+phred_value+" -trimlog "+env_para.output_dir+"/trim.log "+fq1+" "+fq2+" "+whole_out+".trim.r1.fq.gz "+whole_out+".unpaired.r1.fq.gz "+whole_out+".trim.r2.fq.gz "+whole_out+".unpaired.r2.fq.gz "+other_para_str+" 2>"+whole_out+".trim.e.log";
		}else{
			cmd="java -Xmx1750m -jar "+jar_path+" PE -threads 12 -phred33 -trimlog "+env_para.output_dir+"/trim.log "+fq1+" "+fq2+" "+whole_out+".trim.r1.fq.gz "+whole_out+".unpaired.r1.fq.gz "+whole_out+".trim.r2.fq.gz "+whole_out+".unpaired.r2.fq.gz LEADING:10 TRAILING:10 SLIDINGWINDOW:10:15 MINLEN:50 2>"+whole_out+".trim.e.log";
		}
		if(system(cmd.c_str())==-1){
			cerr<<"Error:trim error"<<endl;
			exit(1);
		}
		
	}
	if(fid=="fastqc"){
		string fq1_cmd,fq2_cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string fq1=env_para.m_config2["fastq1"];
		string fq2=env_para.m_config2["fastq2"];
		string fastqc_exe=env_para.bin_dir+"/../software/FastQC/fastqc";
		string threads="4"; //default value
		if(env_para.m_config.find("fastqc")!=env_para.m_config.end()){
			for(vector<string>::iterator ix=env_para.m_config["fastqc"].begin();ix!=env_para.m_config["fastqc"].end();ix++){
				if((*ix).find("threads")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					threads=k_v[1];
				}
			}
		}
		fq1_cmd=fastqc_exe+" --format fastq -t "+threads+" "+fq1+" -o "+env_para.output_dir+" 1>"+whole_out+".fastqc.log";
		fq2_cmd=fastqc_exe+" --format fastq -t "+threads+" "+fq2+" -o "+env_para.output_dir+" 1>>"+whole_out+".fastqc.log";
		if(system(fq1_cmd.c_str())==-1){
			cerr<<"Error:fastqc error"<<endl;
			exit(1);
		}
		if(system(fq2_cmd.c_str())==-1){
			cerr<<"Error:fastqc error"<<endl;
			exit(1);
		}
	}
	if(fid=="bwa"){
		string cmd1,cmd2;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		//string ref_fasta;
		string bwa_path=env_para.bin_dir+"/bwa";
		string samtools_path=env_para.bin_dir+"/samtools";
		string pl="illumina";//default
		string mem_threads="8";
		string sort_threads="8";
		
		if(env_para.m_config.find("bwa")!=env_para.m_config.end()){
			for(vector<string>::iterator ix=env_para.m_config["bwa"].begin();ix!=env_para.m_config["bwa"].end();ix++){
				if((*ix).find("mem_threads")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					mem_threads=k_v[1];
				}
				if((*ix).find("sort_threads")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					sort_threads=k_v[1];
				}
				if((*ix).find("PL")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					pl=k_v[1];
				}
			}
		}
		cmd1=bwa_path+" mem -R '@RG\\tID:"+env_para.output_prefix+"\\tSM:"+env_para.output_prefix+"\\tPL:"+pl+"' -t "+mem_threads+" "+env_para.m_config2["ref_fasta"]+" "+env_para.m_config2["fastq1"]+" "+env_para.m_config2["fastq1"]+" | "+samtools_path+" view -b - | "+samtools_path+" sort -@ "+sort_threads+" -O bam -T "+env_para.output_dir+"/tmp123 - -o "+whole_out+".sort.bam";
		if(system(cmd1.c_str())==-1){
			cerr<<"Error:bwa error"<<endl;
			exit(1);
		}
		cmd2=samtools_path+" index "+whole_out+".sort.bam";
		if(system(cmd2.c_str())==-1){
			cerr<<"Error:samtools index error"<<endl;
			exit(1);
		}
	}
	if(fid=="flagstat"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string samtools_path=env_para.bin_dir+"/samtools";
		cmd=samtools_path+" flagstat "+whole_out+".sort.bam >"+whole_out+".aln.stat";
		if(system(cmd.c_str())==-1){
			cerr<<"Error:flagstat error"<<endl;
			exit(1);
		}
	}
	if(fid=="rmdup"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string picard_jar=env_para.bin_dir+"/../software/picard.jar";
		string xms="500m";
		string xmx="1750m";
		if(env_para.m_config.find("picard_rmdup")!=env_para.m_config.end()){
			for(vector<string>::iterator ix=env_para.m_config["picard_rmdup"].begin();ix!=env_para.m_config["picard_rmdup"].end();ix++){
				if((*ix).find("xms")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xms=k_v[1];
				}
				if((*ix).find("xmx")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xmx=k_v[1];
				}
			}
		}
		cmd="java -Xms"+xms+" -Xmx"+xmx+" -jar "+picard_jar+" MarkDuplicates INPUT="+whole_out+".sort.bam OUTPUT="+whole_out+".sort.rmdup.bam METRICS_FILE="+whole_out+".rmdup.info.txt REMOVE_DUPLICATES=true";
		if(system(cmd.c_str())==-1){
			cerr<<"Error:picard rmdup error"<<endl;
			exit(1);
		}
	}
	if(fid=="filter_bam"){
		string cmd1,cmd2;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string samtools_path=env_para.bin_dir+"/samtools";
		string filter_bam_exe=env_para.bin_dir+"/filter_bam";
		cmd1=filter_bam_exe+" "+whole_out+".sort.rmdup.bam _STD_ | "+samtools_path+" view -b - >"+whole_out+".sort.rmdup.filter.bam";
		if(system(cmd1.c_str())==-1){
			cerr<<"Error:filter bam error"<<endl;
			exit(1);
		}
		cmd2=samtools_path+" index "+whole_out+".sort.rmdup.filter.bam";
		if(system(cmd2.c_str())==-1){
			cerr<<"Error:samtools index error"<<endl;
			exit(1);
		}
	}
	if(fid=="flagstat2"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string samtools_path=env_para.bin_dir+"/samtools";
		cmd=samtools_path+" flagstat "+whole_out+".sort.rmdup.filter.bam >"+whole_out+".flagstat";
		if(system(cmd.c_str())==-1){
			cerr<<"Error:flagstat error"<<endl;
			exit(1);
		}
	}
	if(fid=="realign1"){
		
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string gatk_jar=env_para.bin_dir+"/../software/GenomeAnalysisTK.jar";
		string xms="500m";
		string xmx="1750m";
		string threads="8";
		if(env_para.m_config.find("realign1")!=env_para.m_config.end()){
			for(vector<string>::iterator ix=env_para.m_config["realign1"].begin();ix!=env_para.m_config["realign1"].end();ix++){
				if((*ix).find("xms")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xms=k_v[1];
				}
				if((*ix).find("xmx")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xmx=k_v[1];
				}
				if((*ix).find("threads")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					threads=k_v[1];
				}
			}
		}
		cmd="java -Xms"+xms+" -Xmx"+xmx+" -jar "+gatk_jar+" -T RealignerTargetCreator -R "+env_para.m_config2["ref_fasta"]+" -L "+env_para.m_config2["bed"]+" -I "+whole_out+".sort.rmdup.filter.bam -o "+whole_out+".indelcreate.intervals -nt "+threads;
		if(system(cmd.c_str())==-1){
			cerr<<"Error:realign1 error"<<endl;
			exit(1);
		}
	}
	if(fid=="realign2"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string gatk_jar=env_para.bin_dir+"/../software/GenomeAnalysisTK.jar";
		string xms="500m";
		string xmx="1750m";
		string other_para_str;
		if(env_para.m_config.find("realign2")!=env_para.m_config.end()){
			for(vector<string>::iterator ix=env_para.m_config["realign2"].begin();ix!=env_para.m_config["realign2"].end();ix++){
				if((*ix).find("xms")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xms=k_v[1];
				}
				if((*ix).find("xmx")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xmx=k_v[1];
				}
				if((*ix).find("=")==string::npos){
					other_para_str=*ix;
				}
			}
		}
		string process_bed_cmd="cp "+env_para.m_config2["bed"]+" "+env_para.output_dir+"/add_mt.bed";
		if(system(process_bed_cmd.c_str())==-1){
			cerr<<"Error:copy error"<<endl;
			exit(1);
		}
		string new_bed_path=env_para.output_dir+"/add_mt.bed";
		ofstream of_bed(new_bed_path.c_str(),ios::app);
		of_bed<<"MT\t1\t16569"<<endl;
		of_bed.close();
		env_para.m_config2["bed2"]=new_bed_path;
		cmd="java -Xms"+xms+" -Xmx"+xmx+" -jar "+gatk_jar+" -T IndelRealigner -R "+env_para.m_config2["ref_fasta"]+" -targetIntervals "+whole_out+".indelcreate.intervals -L "+env_para.m_config2["bed2"]+" -I "+whole_out+".sort.rmdup.filter.bam -o "+whole_out+".filter.realigned.bam "+other_para_str;
		if(system(cmd.c_str())==-1){
			cerr<<"Error:realign2 error"<<endl;
			exit(1);
		}
	}
	if(fid=="bqsr1"){
		string cmd;
		string gatk_jar=env_para.bin_dir+"/../software/GenomeAnalysisTK.jar";
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string xms="500m";
		string xmx="1750m";
		string threads="8";
		string other_para_str;
		if(env_para.m_config.find("bqsr1")!=env_para.m_config.end()){
			for(vector<string>::iterator ix=env_para.m_config["bqsr1"].begin();ix!=env_para.m_config["bqsr1"].end();ix++){
				if((*ix).find("xms")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xms=k_v[1];
				}
				if((*ix).find("xmx")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xmx=k_v[1];
				}
				if((*ix).find("threads")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					threads=k_v[1];
				}
				if((*ix).find("=")==string::npos){
					other_para_str=*ix;
				}
			}
		}
		string new_bed_path=env_para.output_dir+"/add_mt.bed";
		env_para.m_config2["bed2"]=new_bed_path;
		cmd="java -Xms"+xms+" -Xmx"+xmx+" -jar "+gatk_jar+" -T BaseRecalibrator -R "+env_para.m_config2["ref_fasta"]+" -I "+whole_out+".filter.realigned.bam -knownSites "+env_para.bin_dir+"/../database/dbsnp_hg19_138.vcf -o "+whole_out+".recalibration.grp -L "+env_para.m_config2["bed2"]+" -nct "+threads+" "+other_para_str;
		cout<<cmd<<endl;
		if(system(cmd.c_str())==-1){
			cerr<<"Error:bqsr1 error"<<endl;
			exit(1);
		}
	}
	if(fid=="bqsr2"){
		string cmd;
		string gatk_jar=env_para.bin_dir+"/../software/GenomeAnalysisTK.jar";
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string xms="500m";
		string xmx="1750m";
		string threads="8";
		if(env_para.m_config.find("bqsr2")!=env_para.m_config.end()){
			for(vector<string>::iterator ix=env_para.m_config["bqsr2"].begin();ix!=env_para.m_config["bqsr2"].end();ix++){
				if((*ix).find("xms")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xms=k_v[1];
				}
				if((*ix).find("xmx")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xmx=k_v[1];
				}
				if((*ix).find("threads")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					threads=k_v[1];
				}
			}
		}
		cmd="java -Xms"+xms+" -Xmx"+xmx+" -jar "+gatk_jar+" -T PrintReads -R "+env_para.m_config2["ref_fasta"]+" -I "+whole_out+".filter.realigned.bam -BQSR "+whole_out+".recalibration.grp -o "+whole_out+".filter.realigned.recalibrate.bam -nct "+threads;
		cout<<cmd<<endl;
		if(system(cmd.c_str())==-1){
			cerr<<"Error:bqsr2 error"<<endl;
			exit(1);
		}
	}
	if(fid=="unifiedgenotyper"){
		string cmd;
		string gatk_jar=env_para.bin_dir+"/../software/GenomeAnalysisTK.jar";
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string xms="500m";
		string xmx="1750m";
		string other_para_str;
		if(env_para.m_config.find("unifiedgenotyper")!=env_para.m_config.end()){
			for(vector<string>::iterator ix=env_para.m_config["unifiedgenotyper"].begin();ix!=env_para.m_config["unifiedgenotyper"].end();ix++){
				if((*ix).find("xms")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xms=k_v[1];
				}
				if((*ix).find("xmx")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xmx=k_v[1];
				}
				if((*ix).find("=")==string::npos){
					other_para_str=*ix;
				}
			}
		}
		cmd="java -Xms"+xms+" -Xmx"+xmx+" -jar "+gatk_jar+" -T UnifiedGenotyper -R "+env_para.m_config2["ref_fasta"]+" -I "+whole_out+".filter.realigned.recalibrate.bam -o "+whole_out+".raw.vcf "+other_para_str;
		if(system(cmd.c_str())==-1){
			cerr<<"Error:call snv error"<<endl;
			exit(1);
		}
	}
	if(fid=="process_indel_region"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string exe_path=env_para.bin_dir+"/process.indel.region";
		cmd=exe_path+" "+whole_out+".raw.vcf "+whole_out+".indel_region.bed";
		if(system(cmd.c_str())==-1){
			cerr<<"Error:process indel region error"<<endl;
			exit(1);
		}
	}
	if(fid=="haplotype"){
		string cmd;
		string gatk_jar=env_para.bin_dir+"/../software/GenomeAnalysisTK.jar";
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string indel_bed=whole_out+".indel_region.bed";
		string xms="500m";
		string xmx="1750m";
		string threads="8";
		if(env_para.m_config.find("haplotype")!=env_para.m_config.end()){
			for(vector<string>::iterator ix=env_para.m_config["haplotype"].begin();ix!=env_para.m_config["haplotype"].end();ix++){
				if((*ix).find("xms")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xms=k_v[1];
				}
				if((*ix).find("xmx")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xmx=k_v[1];
				}
				if((*ix).find("threads")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					threads=k_v[1];
				}
			}
		}
		cmd="java -Xms"+xms+" -Xmx"+xmx+" -jar "+gatk_jar+" -T HaplotypeCaller -R "+env_para.m_config2["ref_fasta"]+ " -I "+whole_out+".filter.realigned.recalibrate.bam -L "+indel_bed+" -o "+whole_out+".hap.vcf -nct "+threads;
		if(system(cmd.c_str())==-1){
			cerr<<"Error,haplotype error"<<endl;
			exit(1);
		}
	}
	if(fid=="merge_vcf"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string exe_path=env_para.bin_dir+"/merge_uni_hap_vcf";
		cmd=exe_path+" "+whole_out+".raw.vcf "+whole_out+".hap.vcf "+whole_out+".merged.raw.vcf";
		if(system(cmd.c_str())==-1){
			cerr<<"Error:merge vcf error"<<endl;
			exit(1);
		}
	}
	if(fid=="vqsr1"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string gatk_jar=env_para.bin_dir+"/../software/GenomeAnalysisTK.jar";
		string xms="500m";
		string xmx="1750m";
		string threads="8";
		string maxgaussians="4";
		string other_para_str;
		if(env_para.m_config.find("vqsr1")!=env_para.m_config.end()){
			for(vector<string>::iterator ix=env_para.m_config["vqsr1"].begin();ix!=env_para.m_config["vqsr1"].end();ix++){
				if((*ix).find("xms")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xms=k_v[1];
				}
				if((*ix).find("xmx")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xmx=k_v[1];
				}
				if((*ix).find("threads")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					threads=k_v[1];
				}
				if((*ix).find("maxgaussians")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					maxgaussians=k_v[1];
				}
				if((*ix).find("resource:")!=string::npos){
					vector<string> v_resource;
					line_split(*ix,v_resource);
					for(vector<string>::iterator ix=v_resource.begin();ix!=v_resource.end();ix++){
						if((*ix).find("vcf")!=string::npos){
							if((*ix).find("/")==string::npos){
								*ix=env_para.bin_dir+"/../database/"+*ix;
							}
						}
					}
					other_para_str=join_vector(v_resource,' ');
				}
			}
		}
		cmd="java -Xms"+xms+" -Xmx"+xmx+" -jar "+gatk_jar+" -T VariantRecalibrator -R "+env_para.m_config2["ref_fasta"]+" -input "+whole_out+".merged.raw.vcf "+other_para_str+" -recalFile "+whole_out+".vcf.recal -tranchesFile "+whole_out+".tranches -rscriptFile "+whole_out+".plots.R -nt "+threads+" --maxGaussians "+maxgaussians;
		if(system(cmd.c_str())==-1){
			cerr<<"Error:vqsr1 error"<<endl;
			exit(1);
		}
	}
	if(fid=="vqsr2"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string gatk_jar=env_para.bin_dir+"/../software/GenomeAnalysisTK.jar";
		string xms="500m";
		string xmx="1750m";
		string other_para_str;
		if(env_para.m_config.find("vqsr2")!=env_para.m_config.end()){
			for(vector<string>::iterator ix=env_para.m_config["vqsr2"].begin();ix!=env_para.m_config["vqsr2"].end();ix++){
				if((*ix).find("xms")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xms=k_v[1];
				}
				if((*ix).find("xmx")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xmx=k_v[1];
				}
				if((*ix).find("=")==string::npos){
					other_para_str=*ix;
				}
			}
		}
		cmd="java -Xms"+xms+" -Xmx"+xmx+" -jar "+gatk_jar+" -T ApplyRecalibration -R "+env_para.m_config2["ref_fasta"]+" -input "+whole_out+".merged.raw.vcf "+other_para_str+" -tranchesFile "+whole_out+".tranches -recalFile "+whole_out+".vcf.recal -o "+whole_out+".recalibrated.filtered.vcf";
		if(system(cmd.c_str())==-1){
			cerr<<"Error:vqsr2 error"<<endl;
			exit(1);
		}
	}
	if(fid=="hdfilter"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string exe_path=env_para.bin_dir+"/hdfilter";
		cmd=exe_path+" "+whole_out+".merged.raw.vcf "+whole_out+".hdfilter.vcf";
		if(system(cmd.c_str())==-1){
			cerr<<"Error:hdfilter error"<<endl;
			exit(1);
		}
	}
	if(fid=="target_reads"){
		string cmd;
		string exe_path=env_para.bin_dir+"/bedtools";
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		cmd=exe_path+" intersect -bed -u -abam "+whole_out+".sort.rmdup.filter.bam -b "+env_para.m_config2["bed"]+" | wc -l > "+whole_out+".on_target_reads.txt";
		if(system(cmd.c_str())==-1){
			cerr<<"Error:calculate target reads error"<<endl;
			exit(1);
		}
	}
	if(fid=="target_ratio"){
		string cmd;
		string exe_path=env_para.bin_dir+"/cal_on_target_ratio";
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		cmd=exe_path+" "+whole_out+".on_target_reads.txt "+whole_out+".flagstat >"+whole_out+".target_ratio";
		if(system(cmd.c_str())==-1){
			cerr<<"Error:calculate target reads error"<<endl;
			exit(1);
		}
	}
	if(fid=="depthofcoverage"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string gatk_jar=env_para.bin_dir+"/../software/GenomeAnalysisTK.jar";
		string xms="500m";
		string xmx="1750m";
		string other_para_str;
		if(env_para.m_config.find("depthofcoverage")!=env_para.m_config.end()){
			for(vector<string>::iterator ix=env_para.m_config["depthofcoverage"].begin();ix!=env_para.m_config["depthofcoverage"].end();ix++){
				if((*ix).find("xms")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xms=k_v[1];
				}
				if((*ix).find("xmx")==0){
					vector<string> k_v;
					line_split(*ix,'=',k_v);
					if(k_v.size()!=2){
						cerr<<"Error,config format error"<<endl;
						exit(1);
					}
					xmx=k_v[1];
				}
				if((*ix).find("=")==string::npos){
					other_para_str=*ix;
				}
			}
		}
		cmd="java -Xms"+xms+" -Xmx"+xmx+" -jar "+gatk_jar+" -T DepthOfCoverage -R "+env_para.m_config2["ref_fasta"]+" -I "+whole_out+".sort.rmdup.filter.bam -o "+whole_out+"_target_coverage -L "+env_para.m_config2["bed"]+" "+other_para_str;
		if(system(cmd.c_str())==-1){
			cerr<<"Error:depthofcoverage error"<<endl;
			exit(1);
		}
	}
	if(fid=="samtools_depth1"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string exe_path=env_para.bin_dir+"/samtools";
		string refgene_cds_bed=env_para.bin_dir+"/../database/grch37_refgene_cds.sorted.bed";
		cmd=exe_path+" depth -d 300000 -b "+refgene_cds_bed+" "+whole_out+".sort.rmdup.filter.bam >"+whole_out+".sort.rmdup.filter.bam.depth";
		if(system(cmd.c_str())==-1){
			cerr<<"Error:samtools depth1 error"<<endl;
			exit(1);
		}
	}
	if(fid=="samtools_depth2"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string exe_path=env_para.bin_dir+"/samtools";
		cmd=exe_path+" depth -d 300000 -b "+env_para.m_config2["bed"]+" "+whole_out+".sort.rmdup.filter.bam >"+whole_out+".sort.rmdup.filter.bam.target.depth";
		if(system(cmd.c_str())==-1){
			cerr<<"Error:samtools depth2 error"<<endl;
			exit(1);
		}
	}
	if(fid=="bedcov"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string exe_path=env_para.bin_dir+"/samtools";
		string refgene_cds_bed=env_para.bin_dir+"/../database/grch37_refgene_cds.sorted.bed";
		cmd=exe_path+" bedcov "+refgene_cds_bed+" "+whole_out+".sort.rmdup.filter.bam >"+whole_out+".filter.realigned.recalibrate.bam.bedcov";
		if(system(cmd.c_str())==-1){
			cerr<<"Error:samtools bedcov error"<<endl;
			exit(1);
		}
	}
	if(fid=="final_stat"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string exe_path=env_para.bin_dir+"/final_stat";
		cmd=exe_path+" "+whole_out+".trim.e.log "+whole_out+".aln.stat "+whole_out+".rmdup.info.txt "+whole_out+".target_ratio "+whole_out+"_target_coverage.sample_summary "+whole_out+".filter.realigned.recalibrate.bam.bedcov >"+whole_out+".final.stat";
		if(system(cmd.c_str())==-1){
			cerr<<"Error:final stat error"<<endl;
			exit(1);
		}
	}
	if(fid=="ts_cov"){
		string cmd;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		string exe_path=env_para.bin_dir+"/stat.transcript.cov.py";
		cmd="python "+exe_path+" "+whole_out+".sort.rmdup.filter.bam.depth "+env_para.bin_dir+"/../database/grch37.refgene.txt "+env_para.bin_dir+"/../database/gene_id.txt "+whole_out+".transcript.cov";
		if(system(cmd.c_str())==-1){
			cerr<<"Error:transcript coverage calculate error"<<endl;
			exit(1);
		}
	}
	if(fid=="stat_depth"){
		string cmd;
		string exe_path=env_para.bin_dir+"/stat.depth";
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		cmd=exe_path+" "+env_para.m_config2["bed"]+" "+whole_out+".sort.rmdup.filter.bam.target.depth "+whole_out+".final.stat " +whole_out+".depth.stat.result";
		if(system(cmd.c_str())==-1){
			cerr<<"Error:stat depth error"<<endl;
			exit(1);
		}
	}
	if(fid=="copy_vcf"){
		string cmd1,cmd2;
		string whole_out=env_para.output_dir+"/"+env_para.output_prefix;
		cmd1="cp "+whole_out+".hdfilter.vcf "+whole_out+".final.vcf";
		cmd2="cp "+whole_out+".recalibrated.filtered.vcf "+whole_out+".final.vcf";
		if(system(cmd1.c_str())==-1){
			cerr<<"Error:copy hdfilter vcf error"<<endl;
			exit(1);
		}
		if(system(cmd2.c_str())==-1){
			cerr<<"Error:copy vqsr vcf error"<<endl;
			exit(1);
		}
	}
	add_done_flag(fid,env_para);
}