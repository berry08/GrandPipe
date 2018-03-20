#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include "gc.h"
#include <unistd.h>
#include <stdlib.h>
using namespace::std;
void samtools_index(string root_dir,string bam_file);
void callVCF(string root_dir,string sample_name,string bam_file,string vcf_file,...){
	samtools_index(root_dir,bam_file);
	string realign1="java -Xmx1750m -jar "+root_dir+"/software/GenomeAnalysisTK.jar -T RealignerTargetCreator -R "+root_dir+"/database/human_g1k_v37.mt.fasta -I "+bam_file+" -o "+sample_name+"_result/"+sample_name+".mt.indelcreate.intervals";
	cout<<realign1<<endl;
	if(system(realign1.c_str())==-1){
		cerr<<"ERROR:realign error"<<endl;
		exit(1);
	}
	string realign_bam=sample_name+"_result/"+sample_name+".mt.realigned.bam";
	string realign2="java -Xmx1750m -jar "+root_dir+"/software/GenomeAnalysisTK.jar -T IndelRealigner -targetIntervals "+sample_name+"_result/"+sample_name+".mt.indelcreate.intervals -R "+root_dir+"/database/human_g1k_v37.mt.fasta -I "+bam_file+" -o "+realign_bam;
	cout<<realign2<<endl;
	if(system(realign2.c_str())==-1){
		cerr<<"ERROR:realign error"<<endl;
		exit(1);
	}
	samtools_index(root_dir,realign_bam);
	string unifiedgenotyper="java -Xmx1750m -jar "+root_dir+"/software/GenomeAnalysisTK.jar -T UnifiedGenotyper -R "+root_dir+"/database/human_g1k_v37.mt.fasta -I "+realign_bam+" -o "+vcf_file+" -glm BOTH -dt none";
	if(system(unifiedgenotyper.c_str())==-1){
		cerr<<"ERROR:call vcf error"<<endl;
		exit(1);
	}
}
void samtools_index(string root_dir,string bam_file){
	string cmd=root_dir+"/bin/samtools index "+bam_file;
	if(system(cmd.c_str())==-1){
		cerr<<"ERROR:samtools index error "<<bam_file<<endl;
		exit(1);
	}
}
