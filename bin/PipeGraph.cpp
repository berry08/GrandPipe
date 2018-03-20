#include <string>
#include <vector>
#include <map>
#include <iostream>
#include "PipeGraph.h"
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>  
#include <sys/wait.h>
#include <set>
#include <fstream>
#include "check.h"
#include "run.h"
#include "gc.h"
using namespace::std;
int visited[MAXV];
/*
#define MAXV 100
int visited[MAXV];
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
	Graph(vector<string> vertex_names,multimap<string,string> relation,int num_relation);
	//CreateGraph(algraph raw_graph);
	//DestroyGraph(algraph raw_graph);
	//BFSTraverse(algraph raw_graph);
	void DFSTraverse();
	void DFS(int index);
	void print_allnodes();
	void print_alledges();
	int get_vnode(string query_str);
private:
	algraph raw_graph;
};
*/
void Graph::print_allnodes(){
	for(int i=0;i<raw_graph.num_vertex;i++){
		cout<<raw_graph.pipelist[i].v_data<<endl;
	}
}
int Graph::get_vnode(string query_str){
	for(int i=0;i<raw_graph.num_vertex;i++){
		if(raw_graph.pipelist[i].v_data==query_str)
			return i;
	}
	return -1;
}

Graph::Graph(){
	raw_graph.num_vertex=0;
	raw_graph.num_edge=0;
	raw_graph.pipelist[0].v_data="0";
	raw_graph.pipelist[0].firstarc=NULL;
	env_para.bin_dir="";
	env_para.output_dir="";
	env_para.output_prefix="";
	env_para.m_config.clear();
	env_para.m_config2.clear();
}
Graph::Graph(vector<string> vertex_names,multimap<string,string> relation,int num_relation,cmd_para aaa){
	raw_graph.num_vertex=vertex_names.size();
	raw_graph.num_edge=num_relation;
	env_para=aaa;
	int i=0;
	for(vector<string>::iterator ix=vertex_names.begin();ix!=vertex_names.end();ix++){
		raw_graph.pipelist[i].v_data=*ix;
		raw_graph.pipelist[i].firstarc=NULL;
		i++;
	}
	for(multimap<string,string>::iterator ix=relation.begin();ix!=relation.end();ix++){
		anode *edge=new anode;
		int head_no=get_vnode(ix->second);
		int tail_no=get_vnode(ix->first);
		if(head_no==-1 || tail_no==-1)
			cerr<<"Error:no such step in pipeline"<<endl;
		edge->targetVex=head_no;
		edge->weight=1;
		edge->nextarc=raw_graph.pipelist[tail_no].firstarc;
		raw_graph.pipelist[tail_no].firstarc=edge;
		//cout<<raw_graph.pipelist[head_no].firstarc<<endl;
	}
}
vector<int> Graph::get_in_node(int index){
	//int iter=0;
	vector<int> in_nodes;
	for(int i=0;i<raw_graph.num_vertex;i++){
		anode* next_p=raw_graph.pipelist[i].firstarc;
		while(next_p!=NULL){
			if(next_p->targetVex==index){
				in_nodes.push_back(i);
			}
			next_p=next_p->nextarc;
		}
	}
	return in_nodes;
}
vector<int> Graph::get_out_node(int index){
	vector<int> out_nodes;
	vnode vertex_=raw_graph.pipelist[index];
	anode* next_p=vertex_.firstarc;
	//int iter=0;
	while(next_p){
		//iter++;
		out_nodes.push_back(next_p->targetVex);
		next_p=next_p->nextarc;
	}
	return out_nodes;
}
void Graph::DFS(int index){
	//cout<<raw_graph.pipelist[index].v_data<<endl;
	visited[index]=1;
	anode *edge=raw_graph.pipelist[index].firstarc;
	//cout<<edge->targetVex<<endl;
	//cout<<edge<<endl;
	while(edge){
		int j=edge->targetVex;
		if(visited[j]==0){
			DFS(j);
		}
		edge=edge->nextarc;
	}
}
void Graph::DFSTraverse(){
	for(int i=0;i<raw_graph.num_vertex;i++)
		visited[i]=0;
	for(int i=0;i<raw_graph.num_vertex;i++){
		if(visited[i]==0)
			DFS(i);
	}
}
void Graph::DFS_and_run(int index){
	//cout<<raw_graph.pipelist[index].v_data<<endl;
	
	string stepid=raw_graph.pipelist[index].v_data;
	//cout<<stepid<<":\t"<<index<<"\tout:"<<join_vector(get_out_node(index),',')<<"\tin:\t"<<join_vector(get_in_node(index),',')<<endl;
	//cout<<stepid<<"\t"<<check_container(stepid,env_para)<<endl;
	if(!(check_container(stepid,env_para))){
		//cout<<"run step:"<<stepid<<endl;
		run_container(stepid,env_para);
	}
	
	visited[index]=1;
	anode *edge=raw_graph.pipelist[index].firstarc;
	while(edge){
		int j=edge->targetVex;
		if(visited[j]==0){
			DFS_and_run(j);
		}
		edge=edge->nextarc;
	}
}
void Graph::DFSTraverse2(){
	for(int i=0;i<raw_graph.num_vertex;i++)
		visited[i]=0;
	for(int i=0;i<raw_graph.num_vertex;i++){
		if(visited[i]==0)
			DFS_and_run(i);
	}
}
void Graph::Tran_sub(int index){
	string stepid=raw_graph.pipelist[index].v_data;
	vector<int> in_nodes=get_in_node(index);
	if(in_nodes.size()>1){
		sleep(5);
		return;
	}
	//cout<<stepid<<":status\t"<<check_container(stepid,env_para)<<endl;
	if(!check_container(stepid,env_para)){
		run_container(stepid,env_para);
	}
	vector<int> out_nodes=get_out_node(index);
	if(out_nodes.size()==0){
		//cout<<"node 0:\t"<<getpid()<<endl;
		sleep(2);
	}else if(out_nodes.size()==1){
		anode *edge=raw_graph.pipelist[index].firstarc;
		//cout<<"node 1\t"<<getpid()<<endl;
		Tran_sub(edge->targetVex);
	}else{
		pid_t pid[out_nodes.size()-1];
		for(int i=0;i<out_nodes.size();i++){
			if((pid[i]=fork())==0){
				//cout<<getpid()<<endl;
				Tran_sub(out_nodes[i]);
				exit(0);
			}else if(pid[i]>0){
				//cout<<"parent pid:\t"<<getpid()<<endl;
			}else{
				cerr<<"fork error"<<endl;
				exit(1);
			}
		}
	}
	
}
void Graph::Traverse(){
	vector<int> zero_innode,one_innode,multi_innode;
	for(int i=0;i<raw_graph.num_vertex;i++){
		visited[i]=0;
		if(get_in_node(i).size()==0){
			zero_innode.push_back(i);
		}else if(get_in_node(i).size()==1){
			one_innode.push_back(i);
		}else{
			multi_innode.push_back(i);
		}
	}
	pid_t top_pid=getpid();
	for(vector<int>::iterator ix=zero_innode.begin();ix!=zero_innode.end();ix++){
//		string stepid=raw_graph.pipelist[*ix].v_data;
		Tran_sub(*ix);
	}
	int not_run=multi_innode.size();
	while(1){
		for(vector<int>::iterator ix=multi_innode.begin();ix!=multi_innode.end();ix++){
			vector<int> innodes=get_in_node(*ix);
			int request_num=innodes.size();
			int ready_num=0;
			string stepid=raw_graph.pipelist[*ix].v_data;
			//cout<<"multi_in_nodes:\t"<<stepid<<":status\t"<<check_container(stepid,env_para)<<endl;
			for(vector<int>::iterator ix2=innodes.begin();ix2!=innodes.end();ix2++){
				string stepid2=raw_graph.pipelist[*ix2].v_data;
				if(check_container(stepid2,env_para)==1){
					ready_num++;
				}
			}

			if(!check_container(stepid,env_para)){
				if(request_num==ready_num){
					//cout<<stepid<<":ready to run"<<endl;
					run_container(stepid,env_para);
					if(check_container(stepid,env_para)>0)
						not_run--;
				}else{
					//cout<<stepid<<":not ready to run"<<endl;
				}
			}else{
				not_run--;
			}
		}
		if(not_run==0)
			break;
		sleep(200);
	}
	int request_num=raw_graph.num_vertex;
	while(1){
		int done_num=0;
		for(int i=0;i<raw_graph.num_vertex;i++){
			string stepid=raw_graph.pipelist[i].v_data;
			//cout<<stepid<<"\t"<<check_container(stepid,env_para)<<endl;
			if(check_container(stepid,env_para)>0){
				done_num++;
			}
		}
		if(request_num==done_num){
			break;
		}
		sleep(200);
	}
}

/*
int main(){
	vector<string> steps;
	abc.push_back("a");
	abc.push_back("b");
	abc.push_back("c");
	abc.push_back("d");
	abc.push_back("e");
	abc.push_back("g");
	abc.push_back("f");
	abc.push_back("h");
	multimap<string,string> relationship;
	def.insert(make_pair("a","b"));
	def.insert(make_pair("a","d"));
	def.insert(make_pair("b","c"));
	//def.insert(make_pair("c","d"));
	def.insert(make_pair("d","e"));
	def.insert(make_pair("c","g"));
	def.insert(make_pair("a","f"));
	def.insert(make_pair("g","h"));

	int 
	Graph new_graph(steps,relationship,6);
	new_graph.DFSTraverse();
	return 0;
}
*/
