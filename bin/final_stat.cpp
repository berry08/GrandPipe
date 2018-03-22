#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include "gc.h"
#include "gender_decide.h"
using namespace::std;
void stat_hom_het(string i_vcf_file,int &het,int &hom);
int main(int argc,char* argv[]){
	if(argc!=7){
		cerr<<argv[0]<<"  <qc file>  <aln file>  <dup file>  <target_ratio_file>  <cov_depth_file>  <bedcov file>"<<endl;
		exit(1);
	}
	ifstream qc_file,aln_file,dup_file,target_ratio_file,cov_depth_file;
	ofstream summary_file;
	qc_file.open(argv[1]);
	//ifstream rerun_qc_file;
	string error_log=argv[1];
	if(!qc_file){
		cerr<<"Error:can no open file\t"<<argv[1]<<"\tor\t"<<error_log<<endl;
		exit(1);
	}
	string qc_line;
	while(getline(qc_file,qc_line)){
		if(qc_line.find("Surviving")!=string::npos){
			string::size_type total1=qc_line.find_first_not_of("Input Read Pairs: ");
			string::size_type total2=qc_line.find(" Both Surviving");
			string::size_type survive1=qc_line.find_first_not_of("Both Surviving: ",total2);
			string::size_type survive2=qc_line.find(" (");
			string::size_type ratio1=qc_line.find("(");
			string::size_type ratio2=qc_line.find(")");
			int total_reads=atoi(qc_line.substr(total1,total2-total1).c_str());
			int survive_reads=atoi(qc_line.substr(survive1,survive2-survive1).c_str());
			string survive_rate=qc_line.substr(ratio1+1,ratio2-ratio1-1);
			cout<<"total reads:\t"<<total_reads*2<<endl;
			cout<<"QC passed reads:\t"<<survive_reads*2<<" ("<<survive_rate<<")"<<endl;
		}
	}

	aln_file.open(argv[2]);
	string aln_line;
	while(getline(aln_file,aln_line)){
		string map_rate;
		int reads_l,reads_r;
		if(aln_line.find("mapped (")!=string::npos){
			string::size_type t2=aln_line.find(" +");
			string::size_type u1=aln_line.find_first_not_of(" + ",t2);
			string::size_type u2=aln_line.find(" mapped");
			string::size_type m1=aln_line.find("(");
			string::size_type m2=aln_line.find(" :");
			reads_l=atoi(aln_line.substr(0,t2).c_str());
			reads_r=atoi(aln_line.substr(u1,u2-u1).c_str());
			map_rate=aln_line.substr(m1+1,m2-m1-1);
			cout<<"mapped reads:\t"<<reads_l+reads_r<<" ("<<map_rate<<")"<<endl;
			break;
		}else{
			continue;
		}
	}
	dup_file.open(argv[3]);
	string dup_line;
	while(getline(dup_file,dup_line)){
		if(dup_line.find("PERCENT_DUPLICATION")!=string::npos){
			getline(dup_file,dup_line);
			getline(dup_file,dup_line);
			vector<string> tmp_v;
			line_split(dup_line,'\t',tmp_v);
			cout<<"duplication rate:\t"<<tmp_v[tmp_v.size()-2]<<endl;
		}
	}
	target_ratio_file.open(argv[4]);
	string target_ratio_line;
	while(getline(target_ratio_file,target_ratio_line)){
		if(target_ratio_line.find(".")!=string::npos){
			cout<<"target reads ratio:\t"<<target_ratio_line<<endl;
			break;
		}
	}
	cov_depth_file.open(argv[5]);
	if(!cov_depth_file){
		cerr<<"Error:unable to open file: "<<argv[5]<<endl;
		exit(1);
	}
	string cov_line;
	getline(cov_depth_file,cov_line);
	getline(cov_depth_file,cov_line);
	vector<string> tmp_v2;
	line_split(cov_line,'\t',tmp_v2);
	cout<<"depth:\t"<<tmp_v2[2]<<endl;
	cout<<"%_bases_above_1:\t"<<tmp_v2[6]<<endl;
	cout<<"%_bases_above_5:\t"<<tmp_v2[7]<<endl;
	cout<<"%_bases_above_10:\t"<<tmp_v2[8]<<endl;
	cout<<"%_bases_above_20:\t"<<tmp_v2[9]<<endl;
	cout<<"%_bases_above_30:\t"<<tmp_v2[10]<<endl;
	string bedcov=argv[6];
	string gender=gender_decide(bedcov);
	cout<<"Sample_gender:\t"<<gender<<endl;
	string vcf_file=argv[2];
	string::size_type pos=vcf_file.find(".aln.stat");
	vcf_file.erase(pos,9);
	string vqsr_vcf_file=vcf_file+".recalibrated.filtered.vcf";
	string hd_vcf_file=vcf_file+".hdfilter.vcf";
	ifstream vqsr,hard;
	int vqsr_flag=1;
	vqsr.open(vqsr_vcf_file.c_str());
	if(!vqsr){
		vqsr_flag=0;
	}else{
		if(vqsr.eof()){
			vqsr_flag=0;
		}
	}
	vqsr.close();
	int het_snp(0),hom_snp(0);
	if(vqsr_flag==1){
		stat_hom_het(vqsr_vcf_file,het_snp,hom_snp);
	}else{
		hard.open(hd_vcf_file.c_str());
		if(!hard){
			string final_vcf=vcf_file+".final.vcf";
			stat_hom_het(final_vcf,het_snp,hom_snp);
		}
		stat_hom_het(hd_vcf_file,het_snp,hom_snp);
	}
	cout<<het_snp<<"\t"<<hom_snp<<endl;
	float het_ratio=0;
	if(het_snp!=0 || hom_snp!=0){
		int total_snp=het_snp+hom_snp;
		het_ratio=float(het_snp)/float(total_snp);
	}
	cout<<setprecision(4)<<"Heterozygosis variants ratio:\t"<<het_ratio<<endl;
	string depth_file=vcf_file+".sort.rmdup.filter.bam.depth";
	string boot_dir=argv[0];
	string::size_type pos2=boot_dir.find("bin/final_stat");
	boot_dir.erase(pos2,14);
	string cds_bed_file=boot_dir+"database/grch37_refgene_cds.sorted.txt";
	string stat_cov=boot_dir+"bin/stat.cov.py";
	string out_file=vcf_file+".exon.cov.stat.result";
	string cmd2="python "+stat_cov+" "+cds_bed_file+" "+depth_file+" "+out_file;
	if(system(cmd2.c_str())==-1){
		cerr<<"Error:run cov stat error"<<endl;
		exit(1);
	}
	return 0;
}
void stat_hom_het(string i_vcf_file,int &het,int &hom){
	ifstream vcf_if;
	het=0;
	hom=0;
	vcf_if.open(i_vcf_file.c_str());
	string vcf_line;
	vector<string> line_eles;
	while(getline(vcf_if,vcf_line)){
		if(vcf_line.find("#")==0){
			continue;
		}
		line_split(vcf_line,'\t',line_eles);
		vector<string> last_eles;
		line_split(line_eles[line_eles.size()-1],':',last_eles);
		vector<string> gt_eles;
		line_split(last_eles[0],'/',gt_eles);
		if(last_eles[0]=="./."){
			continue;
		}
		if(gt_eles[0]!=gt_eles[1]){
			het++;
		}else{
			hom++;
		}
	}
	cout<<het<<"\t"<<hom<<endl;
	vcf_if.close();
}
