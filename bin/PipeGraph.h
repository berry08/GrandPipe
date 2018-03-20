#ifndef _PipeGraph_H
#define _PipeGraph_H
#include <string>
#include <vector>
#include <map>
#include <sys/types.h>
#include <set>
#include "run.h"
#include "check.h"
#include "global_cmd_para.h"
using namespace::std;
#define MAXV 100

class anode{
public:
	int targetVex;
	anode *nextarc;
	int weight;
};
class vnode{
public:
	string v_data;
	anode *firstarc;
};

class algraph{
public:
	vnode pipelist[MAXV] ;
	int num_vertex;
	int num_edge;
};

class Graph{
public:
	Graph();
	Graph(vector<string> vertex_names,multimap<string,string> relation,int num_relation,cmd_para aaa);
	//CreateGraph(algraph raw_graph);
	//DestroyGraph(algraph raw_graph);
	//BFSTraverse(algraph raw_graph);
	void Traverse();
	void Tran_sub(int index);
	void DFSTraverse();
	void DFS(int index);
	void DFSTraverse2();
	void DFS_and_run(int index);
	void print_allnodes();
	int get_vnode(string query_str);
	vector<int> get_in_node(int index);
	vector<int> get_out_node(int index);
private:
	algraph raw_graph;
	cmd_para env_para;
};

#endif
