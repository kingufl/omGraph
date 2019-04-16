#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <string>
#include <libgen.h> // basename
#include <limits>
#include <cassert>
#include <cmath>


#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "io.hpp"
#include "debruijn_graph_shifted.hpp"
#include "algorithm.hpp"
#include "cosmo-color.hpp"

using namespace std;
using namespace sdsl;

#include <sys/timeb.h>

set<int> rest_nodes;
set<int>::iterator it;
set <int> loopy_nodes;

int getMilliCount(){
  timeb tb;
  ftime(&tb);
  int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
  return nCount;
}


int getMilliSpan(int nTimeStart){
  int nSpan = getMilliCount() - nTimeStart;
  if(nSpan < 0)
    nSpan += 0x100000 * 1000;
  return nSpan;
}

string extension = ".dbg";

int sign(float x){

    if(x>0){
        return 1;
    }
    else if(x<0)
        return -1;
    else
        return 0;
}



void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("Opgraph Copyright (c) Kingshuk Mukherjee 2017", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".packed edge file (output from pack-edges).", true, "", "input_file", cmd);
  TCLAP::ValueArg<std::string> restriction_enzyme_arg("o", "restriction_enzyme",
            "restriction_enzyme", true, "", "restriction_enzyme", cmd);

  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  params.restriction_seq   = restriction_enzyme_arg.getValue();
}

static char base[] = {'?','A','C','G','T'};


const char *const starts[] = {"GCCATACTGCGTCATGTCGCCCTGACGCGC","GCAGGTTCGAATCCTGCACGACCCACCAAT","GCTTAACCTCACAACCCGAAGATGTTTCTT","AAAACCCGCCGAAGCGGGTTTTTACGTAAA","AATCCTGCACGACCCACCAGTTTTAACATC","AGAGTTCCCCGCGCCAGCGGGGATAAACCG","GAATACGTGCGCAACAACCGTCTTCCGGAG"};

void find_restriction_kmers(debruijn_graph<> dbg, string restriction_enzyme)
{
    ofstream ofile;
    ofile.open("restriction_nodes");
    int cnt=0;

    for (size_t i = 0; i < dbg.num_nodes(); i++) {

        if(i%100000==0)
            cout<<i<<" "<<cnt<<endl;
        ssize_t start = i; // place to store start of branch kmer
        std::string start_label(dbg.node_label(start));

        if(start_label.compare(0,restriction_enzyme.length(),restriction_enzyme)==0){
            ofile<<i<<" "<<start_label<<endl;
            cnt++;
        }


    }
}

int loop=0;

int path_len=0;

vector <int> visited;
vector <int> distance_from_start;
vector <int> before;

vector < vector <int> > rmap;

stack <int > dfs_stack;


void read_rmaps(char* filename){

    ifstream file(filename, ifstream::in);
	if (!file.is_open()) {
		cout << "Error opening file" << endl;
		exit(1);
	}

	string str;

	while(getline(file, str))
    {
        istringstream ss(str);
        string map_name;
        ss >> map_name;
        float frag;
        vector <int> tem;
        while(ss>>frag){
            tem.push_back(frag*1000);
        }

        rmap.push_back(tem);

//        cout<<tem<<endl;
    }

//    cout<< "total number of rmaps = "<<rmap.size()<<endl;
    file.close();

}

int check_restrcition_node(int pos){

//   std::string start_label(dbg.node_label(pos));
//
//    if(start_label.compare(0,restriction_enzyme.length(),restriction_enzyme)==0){
//            return 1;
//    }
//    else
//        return 0;

    if(rest_nodes.find(pos)==rest_nodes.end())
        return 0;
    else
        return 1;


}

int found_res_node=0;

int search_index=0;
int origin;

struct triplet{
    int source;
    int target;
    int distance;

};

class trip_sp{
    public:
    int u,v,d;

    trip_sp(int a, int b, int c){
        u=a;
        v=b;
        d=c;
    }
};

//trip_sp asas(1,2,3);

vector< triplet > Medges;

vector<triplet> all_edges;

int factoring_loops;

int nodes_visited;

void DFS(debruijn_graph<> dbg, int start_node)
{
//    ofstream ofile;
//    ofile.open("traversal_result");
//    int cnt=0;
    if(visited[start_node]==1)
        return;
    else
        visited[start_node]++;

    search_index++;
    string res="AACTGAT";
    cout<<start_node<<"     "<<dbg.node_label(start_node)<<"     ";
    cout<<"outdegree = "<<dbg.outdegree(start_node)<<" "<<distance_from_start[start_node]<<"  visited"<<visited[start_node]<<endl;
//    cout<<"next node:"<<endl;

    if(check_restrcition_node(start_node)){
        cout<< " Found res node!! "<<endl;
        triplet t;
        t.target=start_node;
        t.distance=distance_from_start[start_node];
        t.source= origin;
        Medges.push_back(t);
        return;
    }
    else if(dbg.outdegree(start_node)==0){
        cout<<origin<< " nothing beyond!! "<< nodes_visited<<" ";
        return;
    }

    for (unsigned long x = 1; x < dbg.sigma + 1; x++) { // iterate through the alphabet of outgoing edges from node i
                // follow each strand or supernode

                ssize_t edge = dbg.outgoing_edge(start_node, x);
                if (edge == -1)
                    continue;
                else{
                    ssize_t pos = dbg._edge_to_node(edge);
                    distance_from_start[pos]=distance_from_start[start_node]+1;
                    DFS(dbg,pos);
                }
    }

//    cout<<endl;
}

void print_path(int node){
    cout<< origin << "--->" << node << " length : " << distance_from_start[node]<< endl;
    cout<< node << " ";
    while(before[node]!=0){
        cout<< before[node] <<" ";
        node=before[node];
    }
    cout<< endl;

}

void DFS_iterative(debruijn_graph<> dbg, int start_node)
{
//    ofstream ofile;
//    ofile.open("traversal_result");
//    int cnt=0;
    dfs_stack.push(start_node);

    while(!dfs_stack.empty()){
        int next_node=dfs_stack.top();
        dfs_stack.pop();

        if(check_restrcition_node(next_node)){
//        cout<< " Found res node!! ";
        print_path(next_node);
        triplet t;
        t.target=next_node;
        t.distance=distance_from_start[next_node];
        t.source= origin;
        Medges.push_back(t);
        continue;
        }
        else if(dbg.outdegree(next_node)==0){
//        cout<<origin<< " nothing beyond!! "<< nodes_visited<<" ";
        continue;
        }

        if(visited[next_node]==0){
            nodes_visited++;
            visited[next_node]++;

//        cout<<next_node<<"     "<<dbg.node_label(next_node)<<"     ";
//        cout<<"outdegree = "<<dbg.outdegree(next_node)<<" "<<distance_from_start[next_node]<<"  visited"<<visited[next_node]<<endl;

        for (unsigned long x = 1; x < dbg.sigma + 1; x++) { // iterate through the alphabet of outgoing edges from node i
                // follow each strand or supernode

                ssize_t edge = dbg.outgoing_edge(next_node, x);
                if (edge == -1)
                    continue;
                else{
                    ssize_t pos = dbg._edge_to_node(edge);
                    before[pos]=next_node;
                    distance_from_start[pos]=distance_from_start[next_node]+1;
                    dfs_stack.push(pos);
                }
        }
        }
    }


//    cout<<"next node:"<<endl;

//    cout<<endl;
}

typedef pair< int, int > iPair;

void find_min_dis(debruijn_graph<> dbg, int start_node)
{
    int max_lim=numeric_limits<int>::max();
    vector< int > dist;
    dist.resize(dbg.num_nodes(),max_lim);

    origin=start_node;

    vector < int > before;
    before.resize(dbg.num_nodes(),-1);
    before[origin]=0;

    priority_queue< iPair, vector <iPair> , greater<iPair> > pq;

    dist[start_node]=0;
    pq.push(make_pair(0,start_node));

    while(!pq.empty()){

        int top_node=pq.top().second;
        pq.pop();
        nodes_visited++;

        if(check_restrcition_node(top_node) && top_node!=start_node){
//        cout<< " Found res node!! ";
//        print_path(top_node);
        triplet t;
        t.target=top_node;
        t.distance=dist[top_node];
        t.source= origin;
        Medges.push_back(t);
        continue;
        }

        else if(dbg.outdegree(top_node)==0){
//        cout<< " nothing beyond!! ";
        continue;
        }

        for (unsigned long x = 1; x < dbg.sigma + 1; x++) { // iterate through the alphabet of outgoing edges from node i
                // follow each strand or supernodfind_min_distances(dbg,*it);e

                ssize_t edge = dbg.outgoing_edge(top_node, x);
                if (edge == -1)
                    continue;
                else{
                    ssize_t pos = dbg._edge_to_node(edge);
                    if(dist[pos]>dist[top_node]+1){
                        dist[pos]=dist[top_node]+1;
                        before[pos]=top_node;
                        pq.push(make_pair(dist[pos],pos));
                    }
                    else if(dist[pos]<dist[top_node]+1){
                        int temp=top_node;
                        while(before[temp]!=0){
                            temp=before[temp];
                            if(temp==pos){
                                cout << "CYCLE INDEED!! ";
                                break;
                            }
                        }
                        int size_of_cycle=dist[top_node]+1-dist[pos];
                        cout<<"Possible cycle/bulge!! origin:"<<origin<< " size of cycle/bulge:"<<dist[top_node]+1-dist[pos]<<endl;
                        if(size_of_cycle>200)
                            return;
                    }
                }
        }

    }


//    cout<<endl;
}

vector < vector<int> > visited_dis;


pair <int,int> find_h(float x, float p){

    float u,b,y1,y2,l,h;

    x=x/1000;

    float x1,x2;

    if(x<2.4)
        {u=0.858181;b=0.180196;}
    else if(x>=2.4 && x<3.6)
        {u=0.980760;b=0.071176;}
    else if(x>=3.6 && x<4.8)
        {u=1.003354; b=0.052800;}
    else
        {u=1.00482;b=0.042428;}

    x1=u+b*log(2*b*p);
    x2=u-b*log(2*b*p);

//    l=x1*x;
//    h=x2*x;

    y1=u-b*sign(p-0.5)*log(1-2*abs(p-0.5));
    y2=u-b*sign(1-p-0.5)*log(1-2*abs(1-p-0.5));

    l=y1*x*1000;
    h=y2*x*1000;

    return(make_pair((int)l,(int)h));

}



void find_all_paths_dp(debruijn_graph_shifted<> dbg, int node, int d){

//    vector< vector < int > > dp_table;
//    dp_table.resize(dbg.num_nodes());
    cout<<" Fact loops:"<<factoring_loops<<endl;
    vector< vector < int > > nodes_at_dis;
    nodes_at_dis.resize(d+1);

    nodes_at_dis[0].push_back(node);
//    dp_table[node].push_back(0);

//    ofstream node_table;
//    node_table.open("node_table.txt");
//    ofstream dis_table;
//    dis_table.open("distance_table.txt");
//    ofstream edge_outfile;
//    edge_outfile.open("edges.txt");

//    dis_table<<"0"<<" : "<<node<<endl;

    int nodes_visit=0;

    int i;

    vector<int> latest_dis(dbg.num_nodes(),0);

    vector<int> reached_status;
    reached_status.resize(dbg.num_nodes(),0);

    vector <int> before(dbg.num_nodes(),0);
    visited_dis.resize(dbg.num_nodes());

    for(i=1;i<=d;i++){

        set <int> reached;

//        dis_table<<i<<" : ";
        vector<int> temp;

        for(int j=0;j<nodes_at_dis[i-1].size();j++){

//            if(temp.size()>1000){
//                break;
//            }
            nodes_visit++;

            int top_node=nodes_at_dis[i-1][j];

            for (unsigned long x = 1; x < dbg.sigma + 1; x++) { // iterate through the alphabet of outgoing edges from node i
                // follow each strand or supernodfind_min_distances(dbg,*it);e

                ssize_t edge = dbg.outgoing_edge(top_node, x);
                if (edge == -1)
                    continue;
                else{
                    ssize_t pos = dbg._edge_to_node(edge);
                    if(reached.find(pos)!=reached.end()){continue;}
                    if(pos==top_node || pos==node)
                        continue;

                    if(loopy_nodes.find(pos)!=loopy_nodes.end()){

                        if(reached_status[pos]==1){
//                            cout<<"LOOPY NODES RESOLVED"<<pos<<endl;
                            continue;
                        }

                        else
                            {
//                                cout<<"LOOPY NODES FOUND"<<pos<<endl;
                                reached_status[pos]=1;}


                    }

                    if(latest_dis[pos]>0 && i-latest_dis[pos]<factoring_loops){continue;}
                    reached.insert(pos);
                    latest_dis[pos]=i;
                    before[pos]=top_node;
                    visited_dis[pos].push_back(i);
                    if(rest_nodes.find(pos)==rest_nodes.end()){
                        temp.push_back(pos);

                    }
                    else{
                        int flag=0;
                        triplet t;
                        t.source=node;
                        t.target=pos;
                        t.distance=i;
                        all_edges.push_back(t);

//                        set <int> check;
//                        int this_rest_node=pos;
//                        while(before[pos]!=0){
//                            cout<<before[pos]<<"  ";
//                            if(check.find(before[pos])==check.end())
//                                check.insert(before[pos]);
//                            else{
//                              return;
//                                flag=1;
//                            }
//                            pos=before[pos];
//
//                        }
//                        if(flag){
//                            cout<<"TROUBLE!!  "<<this_rest_node<<" ";
//                            return;
//                        }

                    }


//                    dp_table[pos].push_back(d);
//                    dis_table<<pos<<" ";
                }
            }


        }

        if(temp.size()==0){
//            dis_table<<endl<<"No outgoing edge";
            break;
        }


        nodes_at_dis[i]=(temp);

//        dis_table<<endl;

//        if(i%100==0)
//            cout <<endl<< i<< " "<<temp.size()<<" "<<nodes_visit<<" "<<all_edges.size();


    }

    cout<< node<<" "<<i<<" "<<nodes_visit<<" "<<all_edges.size()<<endl;

//    for(int i=0;i<dp_table.size();i++){
//        node_table<<i<<" : ";
//        for(int j=0;j<dp_table[i].size();j++){
//            node_table<<dp_table[i][j]<<" ";
//        }
//        node_table<<endl;
//    }

}

map <int, int> rest_map;
vector < vector < pair <int, int> > > multi_graph;

void load_simple_paths(){

ifstream infile;
infile.open("all_edges50k64");

int u,v,d;

string str;
int lastmapped=0;
int edges=0;
while(getline(infile, str))
    {
        edges++;
        istringstream ss(str);

        ss >> u;
        ss >> v;
        ss >> d;
        if(u==v){
//            cout<<u<<" "<<v<<" "<<d<<endl;
        }
        if(rest_map.find(u)==rest_map.end()){
            rest_map[u]=lastmapped++;
//                cout<<"weird node : "<<u<<endl;
        }
        if(rest_map.find(v)==rest_map.end()){
            rest_map[v]=lastmapped++;
//                cout<<"weird node : "<<u<<endl;
        }
        multi_graph.resize(lastmapped);
        multi_graph[rest_map[u]].push_back(make_pair(rest_map[v],d));
    }

//cout<<"Simple paths read"<<endl;
//cout<< "Number of nodes : "<<lastmapped<<endl<<"Number of edges : "<<edges<<endl;

}

vector < vector <int> > paths;

vector <int> temp_kgram;
int count_seed;

vector < trip_sp > curr_path;

int found_seed;

float min_sz_err=1000;

int check_frag(int u, int k, int k_max, int relax){
//    if(temp_kgram[k]>50000)
//        return 1;

    pair< int, int > tem=find_h(temp_kgram[k],0.1);
    int l=tem.first,h=tem.second;
//    int l=temp_kgram[k]-relax,h=temp_kgram[k]+relax;
    if(found_seed==1)
                return 1;
    if(k>k_max)
        return 0;
    for(int i=0;i<multi_graph[u].size();i++){
        if(multi_graph[u][i].second>=l && multi_graph[u][i].second<=h){

            if(k_max>k){
                trip_sp var(u, multi_graph[u][i].first, multi_graph[u][i].second);
                curr_path[k]=var;
                check_frag(multi_graph[u][i].first, k+1, k_max,relax);
            }
            else if(k_max==k){
                int tot_siz_err=0;
//                found_seed=1;
                trip_sp var(u, multi_graph[u][i].first, multi_graph[u][i].second);
                curr_path[k]=var;
//                for(int cp=0;cp<curr_path.size();cp++){
//                    tot_siz_err+=abs(curr_path[cp].d-temp_kgram[cp]);
//                    cout<<"("<<curr_path[cp].u<<"--"<<curr_path[cp].d<<"-->"<<curr_path[cp].v<<")";
//
//                }
                if((float)tot_siz_err/curr_path.size()<min_sz_err)
                    min_sz_err=(float)tot_siz_err/curr_path.size();
//                cout<< "Average sizing error : "<<(float)tot_siz_err/curr_path.size();
//                cout<<endl;
                count_seed++;
//                cout << k <<" "<<temp_kgram[k]<< " "<<multi_graph[u][i].second <<" "<<count_seed<<endl;
                return 1;
            }

        }
    }


}

int check_kgram(vector <int> kgram, int relax){
    int flag=0;
    pair< int, int > tem=find_h(kgram[0],0.1);
    int l=tem.first,h=tem.second;
//    int l=kgram[0]-relax,h=kgram[0]+relax;

    temp_kgram=kgram;
    count_seed=0;
    found_seed=0;
    curr_path.resize(kgram.size(),trip_sp(0,0,0));
//    for(int i=0;i<kgram.size();i++){
//        cout << kgram[i]<<" ";
//    }
//    cout<<endl;
    for(int i=0;i<multi_graph.size();i++){
        for(int j=0;j<multi_graph[i].size();j++){
            if(found_seed==1)
                return 1;
            if(multi_graph[i][j].second>=l && multi_graph[i][j].second<=h){
                trip_sp var(i, multi_graph[i][j].first, multi_graph[i][j].second);
                curr_path[0]=var;
                check_frag(multi_graph[i][j].first,1,kgram.size()-1,relax);

            }
        }

    }
    flag=count_seed;
//    cout<<flag<< endl;

    return flag;

//    for(int i=0;i<graph_kmers.size();i++){
//        if(abs(graph_kmers[i][0]-a)<relax && abs(graph_kmers[i][1]-b)<relax && abs(graph_kmers[i][2]-c)<relax){flag=1;
//        int total_dev= abs(graph_kmers[i][0]-a)+abs(graph_kmers[i][1]-b)+abs(graph_kmers[i][2]-c);
//        cout<<total_dev<<" || "<<graph_kmers[i][0]<<" "<<graph_kmers[i][1]<<" "<<graph_kmers[i][2]<< " || ";
//        break;}
//    }
//multi_graph[i][j].second
//    return flag;
}

void fill_gap(int u, int v, int d1, int d2, int hops){
    vector <vector < trip_sp > > active_nodes;
    active_nodes.resize(hops+1);
    active_nodes[0].push_back(trip_sp(0,u,0));


    for(int i=1;i<hops+1;i++){
        vector< trip_sp > temp_nodes;
        for(int j=0;j<active_nodes[i-1].size();j++){
            int th_node=active_nodes[i-1][j].v;
            int th_dis=active_nodes[i-1][j].d;
            for(int k=0;k<multi_graph[th_node].size();k++){
                    if(multi_graph[th_node][k].second+th_dis>d2)
                        continue;
                    else if(multi_graph[th_node][k].first==v && multi_graph[th_node][k].second+th_dis>d1){
//                        Track back path till u
                        cout << "("<<v<<") <-- "<<multi_graph[th_node][k].second+th_dis<<" --";
                        vector<int> temp;
                        temp.push_back(multi_graph[th_node][k].second+th_dis);
                        int before=th_node;
                        int index=j;
                        int t;
                        for(t=i-1;t>0;t--){
                            cout <<"("<< active_nodes[t][index].v<<") <-- "<<active_nodes[t][index].d<<" --";
                            temp.push_back(active_nodes[t][index].d);
                            index=active_nodes[t][index].u;

                        }
                        cout <<"("<< active_nodes[t][index].v<<")";
                        cout<<endl;
                        vector <int> flip;
                        int last=0;
                        for(int s=temp.size()-1;s>=0;s--){
                            cout<<temp[s]-last<<" ";
                            flip.push_back(temp[s]-last);
                            last=temp[s];
                        }
                        cout<<endl;
                        continue;
                    }
                    temp_nodes.push_back(trip_sp(j,multi_graph[th_node][k].first,multi_graph[th_node][k].second+th_dis));
            }
        }

        active_nodes[i]=(temp_nodes);

//        cout<<endl<< "Hop : "<<i<<" "<<temp_nodes.size()<<endl;
//        for(int ii=0; ii<active_nodes[i].size();ii++){
//            cout<<active_nodes[i][ii].u<<" "<<active_nodes[i][ii].v<<" "<<active_nodes[i][ii].d<<endl;
//        }

    }



}

float find_opt_alignment(vector <int> tar, vector <int> que){

    for(int i=0;i<que.size();i++){
        cout << "  "<<que[i];
    }

    for(int i=0;i<tar.size();i++){
        cout << endl<<tar[i];
    }

    cout<<endl<<endl;

    int p_max=2;
    int q_max=2;
    float penalty_for_indel=1;
    vector< vector <int> > dp_table(tar.size()+1);
    vector< vector <pair<int,int> > > back_track(tar.size());
    for(int i=0;i<tar.size()+1;i++){
        dp_table[i].resize(que.size()+1);
        back_track[i].resize(que.size());
    }
    for(int i=0;i<dp_table.size();i++){
        dp_table[i][0]=100000;
    }
    for(int i=0;i<dp_table[0].size();i++){
        dp_table[0][i]=100000;
    }

    dp_table[0][0]=0;

    for(int i=1;i<dp_table.size();i++){
        for(int j=1;j<dp_table[0].size();j++){
            pair <int,int> tem_pair(0,0);
            float min_sc=10000;
            for(int p=0;p<=p_max;p++){
                for(int q=0;q<=q_max;q++){
                    if(i-p-1<0 || j-q-1<0)
                        continue;
                    float sum_tar=0,sum_que=0;
                    for(int it_tar=i-p-1;it_tar<=i-1;it_tar++)
                        sum_tar+=tar[it_tar];
                    for(int it_que=j-q-1;it_que<=j-1;it_que++)
                        sum_que+=que[it_que];
                    float score=dp_table[i-p-1][j-q-1]+abs(p+q)*penalty_for_indel+abs(sum_tar-sum_que);
                    if(score<min_sc){
                        min_sc=score;
                        tem_pair=make_pair(i-p-1,j-q-1);
                    }

                }
            }
            back_track[i-1][j-1]=tem_pair;
            dp_table[i][j]= min_sc;
        }
    }

for(int i=0; i<tar.size()+1; i++){
    for(int j=0;j<que.size()+1;j++){
        cout<< dp_table[i][j]<<"  ";
    }
    cout<<endl;
}

for(int i=1; i<tar.size()+1; i++){
    for(int j=1;j<que.size()+1;j++){
        cout<<"("<< back_track[i-1][j-1].first<<","<<back_track[i-1][j-1].second<<") ";
    }
    cout<<endl;
}

cout<<endl;

int a=tar.size()-1,b=que.size()-1;
vector <pair<int,int>> alignments;
alignments.push_back(make_pair(tar.size(),que.size()));
while(back_track[a][b].first!=0 && back_track[a][b].second!=0){
    alignments.push_back(make_pair(back_track[a][b].first,back_track[a][b].second));
    cout<<"("<<back_track[a][b].first<<","<<back_track[a][b].second<<") ";
    int ta=back_track[a][b].first-1;
    int tb=back_track[a][b].second-1;
    a=ta;b=tb;
}

a=1,b=1;

cout<<endl;
for(int i=alignments.size()-1;i>=0;i--){
    for(int j1=a;j1<=alignments[i].first;j1++){
        cout<<tar[j1-1]<<" ";
    }
    cout<<" ---> ";

    for(int j2=b;j2<=alignments[i].second;j2++){
        cout<<que[j2-1]<<" ";
    }
    cout<<endl;
    a=alignments[i].first+1;
    b=alignments[i].second+1;

}
cout<<endl;

return(dp_table[tar.size()][que.size()]);


}


void print_multi_graph_stats(){

    ofstream ofile;
    ofile.open("multi_graph_stats");
    for(int i=0;i<multi_graph.size();i++){
//        for(int j=0;j<multi_graph[i].size();j++){
            ofile<<multi_graph[i].size()<<" ";
//        }
        ofile<<endl;
    }

}

int main(int argc, char* argv[]) {
//  parameters_t p;
//  parse_arguments(argc, argv, p);

  cout<<"Usage : big_d || factoring_loops || packed_file_path || restriction_nodes_file_path"<<endl;

  time_t ostart = clock();

  string nu2=argv[1];
  float num2=atof(argv[1]);

  string nu3=argv[2];
  factoring_loops=atoi(argv[2]);

  string packed_file=argv[3];
  string rest_node_path=argv[4];


  int big_d=num2;
    cout<< " D ="<<big_d<< " Factoring loops = "<<factoring_loops<<endl;

  string input_filename=packed_file;

  cerr << "loading dbg" << std::endl;
  debruijn_graph_shifted<> dbg;
  load_from_file(dbg, input_filename);


//  return 1;

//  string restriction_enzyme = p.restriction_seq;

  cout<<dbg.num_nodes()<<endl;
  cout<<dbg.num_edges()<<endl;

  ifstream infile;
  string tem = "restriction_nodes"+nu3;
  infile.open(rest_node_path.c_str());
  int new_node;

  string str;

  vector <int> rest_nodes_vec;

  while(getline(infile, str))
    {
        istringstream ss(str);

        ss >> new_node;
        rest_nodes.insert(new_node);
        rest_nodes_vec.push_back(new_node);

    }

//    cout<<"node under review:"<<num-1<<" "<<rest_nodes_vec[num-1]<<endl;


//  cout<< "num rest nodes:"<<rest_nodes.size()<<endl;
    infile.close();

    int run=0;

//    for (it=rest_nodes.begin(); it!=rest_nodes.end(); ++it){
//        cout<<"-"<<++run<<"-"<<endl;importdata('C:\Users\King\Desktop\Rmaps\rmaps2\rmap_relations_stats.txt');
//        find_all_paths_dp(dbg, 17320  , 25000);
        for(int i=0;i<rest_nodes_vec.size();i++){
            find_all_paths_dp(dbg, rest_nodes_vec[i] , big_d);
        }
//        find_all_paths_dp(dbg, rest_nodes_vec[num-1] , big_d);

   }
//    int sd=0;
//    int cnt=0;
//    int maxb=0;
//    ofstream explore;
//    explore.open("explore");
//    for(int i=0;i<visited_dis.size();i++){
//        if(visited_dis[i].size()==0)
//            continue;
//        if(visited_dis[i].back()-visited_dis[i][0]>maxb)
//            maxb=visited_dis[i].back()-visited_dis[i][0];
//        sd+=visited_dis[i].back()-visited_dis[i][0];
//        cnt++;
////        explore<<i<<"\t";
////        for(int j=0;j<visited_dis[i].size();j++)
////            explore<<visited_dis[i][j]<<" ";
//
////        explore<<endl;
//    }

//    cout<<"Average bulge size: " <<sd/cnt<<endl;
//    cout<<"Max bulge size: " <<maxb<<endl;

    ofstream ofile;
    string f_name="edges";
    ofile.open(f_name.c_str());

    for(int i=0;i<all_edges.size();i++){
            ofile<<all_edges[i].source<<" "<<all_edges[i].target<<" "<<all_edges[i].distance<<endl;

    }

    time_t oend = clock();

    ofstream time_file;
    string t_file="time";
    time_file.open(t_file.c_str());
    time_file<<double(oend - ostart) / CLOCKS_PER_SEC<<endl;
    time_file.close();

    return 1;

    ofile.close();
    ofile.open("mgraph.txt");

    ofstream ofile2;
    ofile2.open("tracker.txt");

    run=0;

    for (it=rest_nodes.begin(); it!=rest_nodes.end(); ++it){
            run++;
//            if (run<9)
//                continue;
            if (run>10000)
                break;
          nodes_visited=0;
          ofile2<<"origin:"<<*it<<" ";
//          find_min_dis(dbg,*it);
          ofile2<<nodes_visited<<" "<<Medges.size()<<endl;


    }

     for(int i=0;i<Medges.size();i++){
            ofile<<Medges[i].source<<" "<<Medges[i].target<<" "<<Medges[i].distance<<endl<<endl;

    }

//  find_restriction_kmers(dbg, restriction_enzyme);



    ofile.close();

}
