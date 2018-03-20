#ifndef GC_H
#define GC_H

#include <string>
#include <vector>
#include <sstream>
#include <cctype>
#include <set>
using namespace::std;

void remove_space(string &a);
int file_exist_and_not_empty(string file_name);
void uniq_vector(vector<string> &a);
string get_local_time();
void line_split(string line_info,char sep,vector<string> &elements);
void line_split(string line_info,char sep,set<string> &elements);
void line_split(string line_info,vector<string> &elements);
int count_gc(string a);
void int2string(int &a,string &b);
void double2string(double &a,string &b);
void float2string(float &a,string &b);
void string2double(string &a,double &b);

void string2upper(string &a);
string join_vector(vector<string> a,char sep);
string join_vector(vector<int> a,char sep);
string join_vector(set<string> a,char sep);
string join_vector(vector<double> a,char sep);
string join_vector(vector<string> a,string seq);
string get_absolute_path(string path);
int check_bam_ok(string bin_path,string bam_file);
int check_bai_ok(string bin_path,string bam_file);
string get_exe_path(string path);
string link_dir_file(string dir,string file);
void chomp_space(string& a,string type);
#endif
