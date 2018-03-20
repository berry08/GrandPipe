#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include "gc.h"
using namespace::std;
using namespace::boost;
int main(int argc,char* argv[]){
	ifstream input_file1,input_file2;
	input_file1.open(argv[1]);
	input_file2.open(argv[2]);
	string line1,line2;
	int on_target=0;
	int total_base=0;
	getline(input_file1,line1);
	on_target=lexical_cast<int>(line1);
	getline(input_file2,line2);
	vector<string> line_eles;
	line_split(line2,' ',line_eles);
	total_base=lexical_cast<int>(line_eles[0]);
	if(total_base==0){
		cerr<<"error:\tno data"<<endl;
		exit(1);
	}
	float ratio=(float)on_target/total_base;
	cout<<ratio<<endl;
	return 0;
}