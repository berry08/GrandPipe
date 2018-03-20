#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include <map>
#include <unistd.h>
#include <stdlib.h>
#include <algorithm>
#include <sys/stat.h>
#include "gc.h"
using namespace::std;

map<int,int> chr_length;



class window{
public:
	int window_size(int end,int start){
		return end-start;
	}
	void print_win(){
		cout<<chr<<"\t"<<start<<"\t"<<end<<endl;
	}
	int chr;
	int start;
	int end;
	bool operator <(const window &w1) const 
	{  
	    if(w1.chr*1000000000+w1.start<chr*1000000000+start)  
	        return true;  
	    else  
	        return false;  
	}  
};
map<window,int> win_stat;
vector<window> gen_window(int window_size);
int init_win_size=1000000;
int main(int argc,char* argv[]){
	if(argc!=3){
		cerr<<"Arguments error!"<<endl;
		cerr<<argv[0]<<"\t<input vcf>\t<output stat>"<<endl;
		exit(1);
	}
	chr_length[1]=249250621;
	chr_length[2]=243199373;
	chr_length[3]=198022430;
	chr_length[4]=191154276;
	chr_length[5]=180915260;
	chr_length[6]=171115067;
	chr_length[7]=159138663;
	chr_length[8]=146364022;
	chr_length[9]=141213431;
	chr_length[10]=135534747;
	chr_length[11]=135006516;
	chr_length[12]=133851895;
	chr_length[13]=115169878;
	chr_length[14]=107349540;
	chr_length[15]=102531392;
	chr_length[16]=90354753;
	chr_length[17]=81195210;
	chr_length[18]=78077248;
	chr_length[19]=59128983;
	chr_length[20]=63025520;
	chr_length[21]=48129895;
	chr_length[22]=51304566;
	chr_length[23]=155270560;
	chr_length[24]=59373566;
	vector<window> v_windows=gen_window(init_win_size);
	ifstream if_vcf;
	if_vcf.open(argv[1]);
	string vcf_line;
	void* file_end=getline(if_vcf,vcf_line);
	vector<window>::iterator ix=v_windows.begin();
	int init_chr=0;
	while(1){
		if(file_end<=0 || ix==v_windows.end()){
			break;
		}
		if(vcf_line.find("#")==0 || vcf_line.find("1/1")!=string::npos || vcf_line.find("PASS")==string::npos){
			file_end=getline(if_vcf,vcf_line);
			continue;
		}
		vector<string> line_eles;
		line_split(vcf_line,'\t',line_eles);
		if(line_eles[0]=="X"){
			line_eles[0]="23";
		}else if(line_eles[0]=="Y"){
			line_eles[0]="24";
		}else if(line_eles[0].find("M")==0){
			line_eles[0]="25";
		}else{
			//line_eles[0]=atoi(line_eles[0].c_str());
		}
		int variant_chr=atoi(line_eles[0].c_str());
		int variant_pos=atoi(line_eles[1].c_str());
		int win_chr=(*ix).chr;
		int win_start=(*ix).start;
		int win_end=(*ix).end;
		if(variant_chr<win_chr){
			file_end=getline(if_vcf,vcf_line);
			continue;
		}else if(variant_chr>win_chr){
			ix++;
			continue;
		}else{
			if(variant_pos>win_end){
				ix++;
				continue;
			}else if(variant_pos<win_start){
				file_end=getline(if_vcf,vcf_line);
				continue;
			}else{
				win_stat[*ix]+=1;
				file_end=getline(if_vcf,vcf_line);
				continue;
			}
		}
	}
	ofstream of_stat;
	of_stat.open(argv[2]);
	for(vector<window>::iterator ix=v_windows.begin();ix!=v_windows.end();ix++){
		if(win_stat.find(*ix)==win_stat.end()){
			win_stat[*ix]=0;
		}
		of_stat<<(*ix).chr<<"\t"<<(*ix).start<<"\t"<<(*ix).end<<"\t"<<win_stat[*ix]<<endl;
	}
/*	string exe_path(argv[0]);
	string ab_path="";
	
	if(exe_path.find("/")!=0){
		ab_path=get_absolute_path(exe_path);
	}else{
		ab_path=exe_path;
	}
	string control_path=ab_path;
	int slash_pos=control_path.find_last_of("/");
	control_path.erase(slash_pos,control_path.size()-slash_pos);
	control_path=control_path+"/../database/idt_homo_stat";
	
	ifstream if_control;
	if_control.open(control_path.c_str());
	vector<string> v_control;
	string control_line;
	while(getline(if_control,control_line)){
		v_control.push_back(control_line);
	}
	if_control.close();
	int iter=0;
	of_stat<<"#Chr\tstart\tend\thet_num_mean\tsd\tsample_het_num\tz_score"<<endl;
	for(vector<window>::iterator ix=v_windows.begin();ix!=v_windows.end();ix++){
		if(win_stat.find(*ix)==win_stat.end()){
			win_stat[*ix]=0;
		}
		vector<string> v_control_line_eles;
		line_split(v_control[iter],'\t',v_control_line_eles);
		float c_mean=atof(v_control_line_eles[3].c_str());
		float c_sd=atof(v_control_line_eles[4].c_str());
		float z_score;
		if(c_sd!=0){
			z_score=(win_stat[*ix]-c_mean)/c_sd;
		}else{
			z_score=0;
		}
		of_stat<<v_control[iter]<<"\t"<<win_stat[*ix]<<"\t"<<z_score<<endl;
		iter++;
	}
*/
	return 0;
}
vector<window> gen_window(int window_size){
	vector<window> v_windows;
	for(int i=1;i<=24;i++){
		int max_length=chr_length[i];
		int j=1;
		for(;j<=max_length;j+=init_win_size){
			//cout<<"here2\t"<<j<<endl;
			window win_1m;
			win_1m.chr=i;
			win_1m.start=j;
			win_1m.end=j+init_win_size-1;
			v_windows.push_back(win_1m);
		}
		if(j!=max_length){
			window win_1m;
			win_1m.chr=i;
			win_1m.start=j;
			win_1m.end=max_length;
			v_windows.push_back(win_1m);
		}
	}
	return v_windows;
}
