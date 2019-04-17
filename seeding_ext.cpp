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
#include "debruijn_graph.hpp"
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

bool debug=false;

vector<string> rmap_name;

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
        rmap_name.push_back(map_name);

//        cout<<tem<<endl;
    }

//    cout<< "total number of rmaps = "<<rmap.size()<<endl;
    file.close();

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

int nodes_visited;

typedef pair< int, int > iPair;

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



void find_all_paths_dp(debruijn_graph<> dbg, int node, int d){

//    vector< vector < int > > dp_table;
//    dp_table.resize(dbg.num_nodes());

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
                    if(pos==top_node)
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

                    if(latest_dis[pos]>0 && i-latest_dis[pos]<100){continue;}
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

vector< int > reverse_map_nodes;

void load_simple_paths(char* filename){


ifstream infile(filename, ifstream::in);

int u,v,d;


string str;
int lastmapped=0;
int edges=0;
while(getline(infile, str))
    {

        istringstream ss(str);

        ss >> u;
        ss >> v;
        ss >> d;
        if(d>100000)
            continue;

        edges++;
        if(u==v){
//            cout<<u<<" "<<v<<" "<<d<<endl;
        }
        if(rest_map.find(u)==rest_map.end()){
            rest_map[u]=lastmapped++;
            reverse_map_nodes.push_back(u);
//                cout<<"weird node : "<<u<<endl;
        }
        if(rest_map.find(v)==rest_map.end()){
            rest_map[v]=lastmapped++;
            reverse_map_nodes.push_back(v);
//                cout<<"weird node : "<<u<<endl;
        }
        multi_graph.resize(lastmapped);
        multi_graph[rest_map[u]].push_back(make_pair(rest_map[v],d));
    }

cout<<"Simple paths read"<<endl;
cout<< "Number of nodes : "<<lastmapped<<endl<<"Number of edges : "<<edges<<endl;



}

vector < vector <int> > paths;

vector <int> temp_kgram;
int count_seed;

vector < trip_sp > curr_path;

vector <vector < trip_sp > >seed;

struct all_seeds{
    vector< int > num_seeds;
    vector <vector < trip_sp > >rmap_seed;
    vector < pair <int, int> > seed_start_end;
    vector < set < int > > begins;
    vector < set < int > > ends;
    vector < pair<int,int> > scores;
    set <int> final_seeds;

};

vector < all_seeds > seed_set;

vector < vector < trip_sp > > best_seeds;

int found_seed;

float min_sz_err=1000;
int max_kk;

int found;

int kgram_size;

int check_frag(int u, int k, int k_max, int relax){
//    if(temp_kgram[k]>50000)
//        re    turn 1;
//    if(found>0)
//        return 1;
    int nxt_true=0;
    if (k==temp_kgram.size() && k>kgram_size-1){

        found++;
        if(max_kk<k){
            max_kk=k;
            seed.clear();
            vector<trip_sp> temp_sp=curr_path;
            temp_sp.resize(max_kk, trip_sp(0,0,0));
            seed.push_back(temp_sp);
//            seed.push_back(curr_path);
//            seed.back().resize(max_kk);
            return k;
        }
        else if(max_kk==k){
            vector<trip_sp> temp_sp=curr_path;
            temp_sp.resize(max_kk, trip_sp(0,0,0));
            seed.push_back(temp_sp);
//            seed.push_back(curr_path);
//            seed.back().resize(max_kk);
            return k;

        }
        else
            return k;
    }
    else if (k==temp_kgram.size() && k<=3)
        return k;


    pair< int, int > tem=find_h(temp_kgram[k],0.1);
//    int l=tem.first,h=tem.second;
    int l=temp_kgram[k]-relax,h=temp_kgram[k]+relax;

    for(int i=0;i<multi_graph[u].size();i++){
        if(multi_graph[u][i].second>=l && multi_graph[u][i].second<=h){

                nxt_true++;
                trip_sp var(u, multi_graph[u][i].first, multi_graph[u][i].second);
                curr_path[k]=var;
                check_frag(multi_graph[u][i].first, k+1, k_max,relax);

//                    else if(k_max==k){
//                        int tot_siz_err=0;
//        //                found_seed=1;
//                        trip_sp var(u, multi_graph[u][i].first, multi_graph[u][i].second);
//                        curr_path[k]=var;
//        //                for(int cp=0;cp<curr_path.size();cp++){
//        //                    tot_siz_err+=abs(curr_path[cp].d-temp_kgram[cp]);
//        //                    cout<<"("<<curr_path[cp].u<<"--"<<curr_path[cp].d<<"-->"<<curr_path[cp].v<<")";
//        //
//        //                }
//                        if((float)tot_siz_err/curr_path.size()<min_sz_err)
//                            min_sz_err=(float)tot_siz_err/curr_path.size();
//        //                cout<< "Average sizing error : "<<(float)tot_siz_err/curr_path.size();
//        //                cout<<endl;
//                        count_seed++;
//        //                cout << k <<" "<<temp_kgram[k]<< " "<<multi_graph[u][i].second <<" "<<count_seed<<endl;
//                        return 1;
//                    }

        }
    }

//    if(nxt_true==0 && k>3){
//        found++;
//        if(max_kk<k){
//            max_kk=k;
//            seed.clear();
//            vector<trip_sp> temp_sp=currcheck_kgram(kgram,300);_path;
//            temp_sp.resize(max_kk, trip_sp(0,0,0));
//            seed.push_back(temp_sp);
////            seed.push_back(curr_path);
////            seed.back().resize(max_kk);
//            return k;
//        }
//        else if(max_kk==k){
//            vector<trip_sp> temp_sp=curr_path;
//            temp_sp.resize(max_kk, trip_sp(0,0,0));
//            seed.push_back(temp_sp);
////            seed.push_back(curr_path);
////            seed.back().resize(max_kk);
//            return k;
//
//        }
//        else
//            return k;
//    }

}

int check_kgram(vector <int> kgram, int relax){
    int flag=0;
    found=0;
    pair< int, int > tem=find_h(kgram[0],0.1);
//    int l=tem.first,h=tem.second;
    int l=kgram[0]-relax,h=kgram[0]+relax;
    max_kk=1;
    seed.clear();
    temp_kgram=kgram;
    count_seed=0;
    found_seed=0;
    curr_path.clear();
    curr_path.resize(kgram.size(),trip_sp(0,0,0));
    for(int i=0;i<kgram.size();i++){
//        cout << kgram[i]<<" ";
    }
//    cout<<endl;
    for(int i=0;i<multi_graph.size();i++){
        for(int j=0;j<multi_graph[i].size();j++){
//            if(found_seed==1)
//                return 1;
            if(multi_graph[i][j].second>=l && multi_graph[i][j].second<=h){
                trip_sp var(i, multi_graph[i][j].first, multi_graph[i][j].second);
                curr_path[0]=var;
                check_frag(multi_graph[i][j].first,1,kgram.size()-1,relax);

            }
        }

    }
    flag=count_seed;
//    cout<<flag<< endl;



//    for(int i=0;i<max_kk;i++){
//        cout<<kgram[i]<<" ";
//    }
//    cout<<endl<<endl<<endl;
    if(debug==true){
//    for(int i=0;i<seed.size();i++){
//        for(int j=0;j<max_kk;j++){
//            cout<<"("<<seed[i][j].u<<"--"<<seed[i][j].d<<"-->"<<seed[i][j].v<<")";
//        }
//
//        cout<<endl;
//    }
    }
//    for(int i=0;i<graph_kmers.size();i++){
//        if(abs(graph_kmers[i][0]-a)<relax && abs(graph_kmers[i][1]-b)<relax && abs(graph_kmers[i][2]-c)<relax){flag=1;
//        int total_dev= abs(graph_kmers[i][0]-a)+abs(graph_kmers[i][1]-b)+abs(graph_kmers[i][2]-c);
//        cout<<total_dev<<" || "<<graph_kmers[i][0]<<" "<<graph_kmers[i][1]<<" "<<graph_kmers[i][2]<< " || ";
//        break;}
//    }
//multi_graph[i][j].second
//    return flag;

    if(seed.size()>1)
        seed.resize(1);

    return seed.size();
}

void trim_seeds(){
    cout<<endl<<"now trimming seeds"<<endl;

            ofstream ofile;
    ofile.open("gap_start_end_before_trim");

    for(int i=0;i<seed_set.size();i++){
        ofile<<"Rmap: "<<i<<" ";
        for(int j=0;j<seed_set[i].seed_start_end.size();j++){
                ofile<<"("<<seed_set[i].seed_start_end[j].first<<","<<seed_set[i].seed_start_end[j].second<<")->"<<seed_set[i].num_seeds[j]<<" ";

        }
        ofile<<endl;
    }
    ofile.close();

    for(int i=0;i<seed_set.size();i++){
        if(seed_set[i].seed_start_end.size()<2){
            cout<< i<<" ";
            continue;
        }
        for(int j=0;j<seed_set[i].seed_start_end.size()-1;j++){
            if(seed_set[i].seed_start_end[j+1].first-seed_set[i].seed_start_end[j].second<2){
                seed_set[i].seed_start_end[j].second--;seed_set[i].seed_start_end[j+1].first++;
                int bef=0;
                for(int k=0;k<j;k++)
                    bef+=seed_set[i].num_seeds[k];

                for(int k=bef;k<bef+seed_set[i].num_seeds[j];k++){
                    seed_set[i].rmap_seed[k].pop_back();
                }
                for(int k=bef+seed_set[i].num_seeds[j];k<bef+seed_set[i].num_seeds[j]+seed_set[i].num_seeds[j+1];k++){
                    seed_set[i].rmap_seed[k].erase(seed_set[i].rmap_seed[k].begin());
                }
            }
        }
    }

    ofile.open("gap_start_end_after_trim");

    for(int i=0;i<seed_set.size();i++){
        ofile<<"Rmap: "<<i<<" ";
        for(int j=0;j<seed_set[i].seed_start_end.size();j++){
                ofile<<"("<<seed_set[i].seed_start_end[j].first<<","<<seed_set[i].seed_start_end[j].second<<")->"<<seed_set[i].num_seeds[j]<<" ";

        }
        ofile<<endl;
    }

    for(int i=0;i<seed_set.size();i++){
        if(seed_set[i].seed_start_end.size()<2){
            continue;
        }
        seed_set[i].ends.resize(seed_set[i].seed_start_end.size());
        seed_set[i].begins.resize(seed_set[i].seed_start_end.size());
        for(int j=0;j<seed_set[i].seed_start_end.size()-1;j++){
                int bef=0;
                for(int k=0;k<j;k++)
                    bef+=seed_set[i].num_seeds[k];

                for(int k=bef;k<bef+seed_set[i].num_seeds[j];k++){
                    seed_set[i].ends[j].insert(seed_set[i].rmap_seed[k].back().v);
                }
                for(int k=bef+seed_set[i].num_seeds[j];k<bef+seed_set[i].num_seeds[j]+seed_set[i].num_seeds[j+1];k++){
                    seed_set[i].begins[j+1].insert(seed_set[i].rmap_seed[k].front().u);
//                    cout << " first node: "<<i<<" "<<j<<" "<< seed_set[i].rmap_seed[k].front().u;
                }

        }
    }

    cout<<endl<<"trimming seeds finished"<<endl;

}



vector < vector < trip_sp > > fill_gap(int u, int v, int d1, int d2, int hops){
//    cout<<endl<<" u:"<<u<<" v:"<<v<<" d1:"<<d1<<" d2:"<<d2<<" hops:"<<hops<<endl;
    vector <vector < trip_sp > > active_nodes;
    active_nodes.resize(hops+1);
    active_nodes[0].push_back(trip_sp(0,u,0));
    vector < vector < trip_sp > > all_paths_dis;

    for(int i=1;i<hops+1;i++){
//        cout<<"hops="<<i<<endl;
        vector < int > node_found;
        node_found.resize(multi_graph.size(),0);
        set <int> this_hop_nodes;
        vector< trip_sp > temp_nodes;
        if(active_nodes[i-1].size()==0)
            break;
        for(int j=0;j<active_nodes[i-1].size();j++){
            if(active_nodes[i-1][j].v==active_nodes[i-1][j].u)
                continue;
            int th_node=active_nodes[i-1][j].v;
            int th_dis=active_nodes[i-1][j].d;

            if(this_hop_nodes.find(th_node)!=this_hop_nodes.end()){
                continue;
            }
            else
                this_hop_nodes.insert(th_node);
//            cout<<"this_node:"<<th_node<<"this_dis"<<th_dis<<endl;
            if(th_dis>0 && node_found[th_node]>0 && abs(th_dis-node_found[th_node])<2000){
                    node_found[th_node]=th_dis;
                    cout<<"check!!!!!"<<th_dis<<" "<<node_found[th_node]<<" ";
                    continue;
            }
            node_found[th_node]=th_dis;
            for(int k=0;k<multi_graph[th_node].size();k++){
//                    cout<<th_node<<","<<multi_graph[th_node][k].first<<endl;
                    if(multi_graph[th_node][k].second+th_dis>d2)
                        continue;
                    else if(multi_graph[th_node][k].first==v && multi_graph[th_node][k].second+th_dis>d1){
//                        Track back path till u
//                        cout << "("<<v<<") <-- "<<multi_graph[th_node][k].second+th_dis<<" --";
                        vector<int> temp;
                        vector< int > path_nodes;
                        path_nodes.push_back(v);
                        temp.push_back(multi_graph[th_node][k].second+th_dis);
                        int before=th_node;
                        int index=j;
                        int t;
                        for(t=i-1;t>0;t--){
//                            cout <<"("<< active_nodes[t][index].v<<") <-- "<<active_nodes[t][index].d<<" --";
                            path_nodes.push_back(active_nodes[t][index].v);
                            temp.push_back(active_nodes[t][index].d);
                            index=active_nodes[t][index].u;

                        }
//                        cout <<"("<< active_nodes[t][index].v<<")";
                        path_nodes.push_back(active_nodes[t][index].v);
//                        cout<<endl;
                        vector <int> flip;
                        vector < trip_sp > paths_dis;
                        int last=0;
                        for(int s=temp.size()-1;s>=0;s--){
//                            cout<<temp[s]-last<<" ";
                            flip.push_back(temp[s]-last);
                            paths_dis.push_back(trip_sp(path_nodes[s+1],path_nodes[s],temp[s]-last));
                            last=temp[s];
                        }
//                        cout<<endl;
                        all_paths_dis.push_back(paths_dis);
                        continue;
                    }
                    temp_nodes.push_back(trip_sp(j,multi_graph[th_node][k].first,multi_graph[th_node][k].second+th_dis));
//                    cout<<j<<" "<<multi_graph[th_node][k].first<<" "<<multi_graph[th_node][k].second+th_dis<<endl;
            }
        }
//        cout<<"active node size:"<<temp_nodes.size()<<endl;
        active_nodes[i]=(temp_nodes);

//        if(i==1){
//            ofstream ofile;
//            ofile.open("see_nodes");
//            for(int m=0;m<temp_nodes.size();m++){
//                ofile<<temp_nodes[m].u<<","<<temp_nodes[m].v<<"->"<<temp_nodes[m].d<<endl;
//            }
//        }



//        cout<<endl<< "Hop : "<<i<<" "<<temp_nodes.size()<<endl;
//        for(int ii=0; ii<active_nodes[i].size();ii++){
//            cout<<active_nodes[i][ii].u<<" "<<active_nodes[i][ii].v<<" "<<active_nodes[i][ii].d<<endl;
//        }

//    cout<<"Active nodes size: "<<active_nodes[i].size()<<" for i="<<i<<endl;

    }

//    cout <<endl;
//for(int i=0;i<all_paths_dis.size();i++){
//    for(int j=0;j<all_paths_dis[i].size();j++){
//        cout<<"("<<all_paths_dis[i][j].u<<","<<all_paths_dis[i][j].d<<","<<all_paths_dis[i][j].v<<")";
//    }
//    cout<<endl;
//}


//if(all_paths_dis.size()==0)
//    cout<<" No paths found that satisfies condition ";
//else
//    cout<<"Number of paths found:"<<all_paths_dis.size();

return(all_paths_dis);

}

//vector < pair<float , float> > optimized_overlap_alignment(vector< float >& rmap1, vector< float >& rmap2, int p1, int p2);
//
//
//float find_opt_alignment(vector <int> tar, vector <int> que){
//    vector < float > rmap1,rmap2;
//    for(int i=0;i<que.size();i++){
//        rmap1.push_back((float)que[i]/1000);
////        cout<<rmap1[i]<<" ";
//    }
//    for(int i=0;i<tar.size();i++){
//        rmap2.push_back((float)tar[i]/1000);
////        cout<<rmap2[i]<<" ";
//    }
//
//    vector < pair<float ,float> > alignment_sc=optimized_overlap_alignment(rmap1,rmap2,0,0);
//
//    return(alignment_sc[0].first);
//
//    cout<<endl<<endl;
//
//    int p_max=3;
//    int q_max=3;
//    float penalty_for_indel=1;
//    vector< vector <int> > dp_table(tar.size()+1);
//    vector< vector <pair<int,int> > > back_track(tar.size());
//    for(int i=0;i<tar.size()+1;i++){
//        dp_table[i].resize(que.size()+1);
//        back_track[i].resize(que.size());
//    }
//    for(int i=0;i<dp_table.size();i++){
//        dp_table[i][0]=100000;
//    }
//    for(int i=0;i<dp_table[0].size();i++){
//        dp_table[0][i]=100000;
//    }
//
//    dp_table[0][0]=0;
//
//
//    cout<<"table initialized";
//
//    for(int i=1;i<dp_table.size();i++){
//        for(int j=1;j<dp_table[0].size();j++){
//            pair <int,int> tem_pair(0,0);
//            float min_sc=10000;
//            for(int p=0;p<=p_max;p++){
//                for(int q=0;q<=q_max;q++){
//                    if(i-p-1<0 || j-q-1<0)
//                        continue;
//                    float sum_tar=0,sum_que=0;
//                    for(int it_tar=i-p-1;it_tar<=i-1;it_tar++)
//                        sum_tar+=tar[it_tar];
//                    for(int it_que=j-q-1;it_que<=j-1;it_que++)
//                        sum_que+=que[it_que];
//                    float score=dp_table[i-p-1][j-q-1]+abs(p+q)*penalty_for_indel+abs(sum_tar-sum_que);
//                    if(score<min_sc){
//                        min_sc=score;
//                        tem_pair=make_pair(i-p-1,j-q-1);
//                    }
//
//                }
//            }
//            back_track[i-1][j-1]=tem_pair;
//            dp_table[i][j]= min_sc;
//        }
//    }
//
////for(int i=0; i<tar.size()+1; i++){
////    for(int j=0;j<que.size()+1;j++){
////        cout<< dp_table[i][j]<<"  ";
////    }
////    cout<<endl;
////}
////
////for(int i=1; i<tar.size()+1; i++){
////    for(int j=1;j<que.size()+1;j++){
////        cout<<"("<< back_track[i-1][j-1].first<<","<<back_track[i-1][j-1].second<<") ";
////    }
////    cout<<endl;
////}
////
////cout<<endl;
//
//int a=tar.size()-1,b=que.size()-1;
//vector <pair<int,int>> alignments;
//alignments.push_back(make_pair(tar.size(),que.size()));
//while(back_track[a][b].first!=0 && back_track[a][b].second!=0){
//    alignments.push_back(make_pair(back_track[a][b].first,back_track[a][b].second));
////    cout<<"("<<back_track[a][b].first<<","<<back_track[a][b].second<<") ";
//    int ta=back_track[a][b].first-1;
//    int tb=back_track[a][b].second-1;
//    a=ta;b=tb;
//}
//
//a=1,b=1;
//
//cout<<endl;
//for(int i=alignments.size()-1;i>=0;i--){
//    for(int j1=a;j1<=alignments[i].first;j1++){
//        cout<<tar[j1-1]<<" ";
//    }
//    cout<<" ---> ";
//
//    for(int j2=b;j2<=alignments[i].second;j2++){
//        cout<<que[j2-1]<<" ";
//    }
//    cout<<endl;
//    a=alignments[i].first+1;
//    b=alignments[i].second+1;
//
//}
//cout<<endl;
//
//return(dp_table[tar.size()][que.size()]);
//
//
//}


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

  time_t ostart = clock();

  cout<<"Usage : paths_file || om-gram_size ||  map_file"<<endl;

  string nu=argv[1];
  float num=atof(argv[1]);
  int top_seeds=atoi(argv[2]);

  string nu2=argv[2];
  float num2=atof(argv[2]);
  int play=atoi(argv[2]);
  play=5;

  string simple_paths_path=argv[1];

  ofstream output("rmap_alignements.txt");

  ofstream output_lengths("rmap_alignements_lengths.txt");


//  cout<< num<<" "<<num2<<endl;
//
//  pair<float,float> range;
//  range=find_h(num,num2);
//
//  cout<<range.first<<" "<<range.second;
//
//  return 1;

//  vector < int > tar,que;
//
//  tar.push_back(12);
//  tar.push_back(5);
//  tar.push_back(23);
//  tar.push_back(8);
//  tar.push_back(10);
//
//  que.push_back(11);
//  que.push_back(7);
//  que.push_back(7);
//  que.push_back(6);
//  que.push_back(13);
//  que.push_back(10);
//  que.push_back(9);
//
//  find_opt_alignment(tar,que);
//
//  return 1;
//
//
  debug=true;
  read_rmaps(argv[3]);
  load_simple_paths(argv[1]);

//  fill_gap(rest_map.find(19556996 )->second,rest_map.find(15569344)->second,11665,17498,4);
//
//
//
//
//  return 1;
  seed_set.resize(rmap.size());

//  print_multi_graph_stats();

//  return 1;


//  fill_gap(186,179,13000,16000,3);
////
//  return 1;

  kgram_size=top_seeds;

  best_seeds.resize(rmap.size());

  ofstream num_seeds;
  num_seeds.open("num_seeds");

//  ofstream seed_rep;
//  seed_rep.open("seed_report");

  vector< float > all_scores;
  vector< float > all_scores_scaled;

  vector<int> kgram;
  int cnt_1=0;
  int cnt_0=0;
  int rmaps_left=0;
  int seeds=0,no_seeds=0;
  int oor=0;
  time_t start_seed=clock();

  int tot_kgrams=0;

  int tot_num_seeds=0;

  vector < float > seed_coverage;

  int total_coverage=0;

//  ofstream seed_cov;
//  seed_cov.open("seed_coverage");

  ofstream seed_det;
  seed_det.open("seed_details");

//  ofstream per_alig;
//  per_alig.open("percent_aligned");
//
//  ofstream seed_stats;
//  seed_stats.open("seed_stats");

//  ofstream best_paths_view;
//  best_paths_view.open("best_paths_view");

//  vector< vector < trip_sp > > pts = fill_gap(563,257,0,15000,5);
//
//  cout<<"paths="<<pts.size();
//  return 1;


  for(int ol=0;ol<rmap.size();ol++){

    int flag=0;
    num=rmap[ol].size();
    cout<<endl<<"Rmap : "<<ol<<" Size:"<<num<<endl;
    int start=0;

//    kgram=rmap[ol];
//
//    check_kgram(kgram,300);

//    for(int j=0;j<rmap[ol].size()-kgram_size+1;j++){
//        vector < int > temp;
//        tot_kgrams++;
//        for(int jj=0;jj<kgram_size;jj++)
//            temp.push_back(rmap[ol][j+jj]);
//
//        kgram=temp;
//
//        check_kgram(kgram,30);
//    }

    int first_seed_start=0;
    int last_seed_finish=0;

    int seeds_for_this=0;

    vector < vector < pair <int,int> > > seeds_of_rmap;

    set < pair < int, int > > om_grams;

    set < pair < int, int > > ::iterator it;

    int lastmapped=0;

    for( int j=0; j<rmap[ol].size()-kgram_size+1; j++){

            int sum=0;

            for (int l=0; l<kgram_size; l++){
                sum+= rmap[ol][j+l];
            }

            om_grams.insert(make_pair(sum,lastmapped++));

    }

    int found=0;
    it=om_grams.end();

    vector < pair <int,int> > seed_locs;

    vector <vector < trip_sp > > rmap_seeds;


    while(found<4){

        int index=(--it)->second;

        kgram.clear();
        int totkgram=0;
        for(int i=index;i<index+kgram_size;i++){
            kgram.push_back(rmap[ol][i]);
            totkgram+=kgram.back();
        }
        tot_kgrams++;
        int numofseeds = check_kgram(kgram,500);

        if(numofseeds==0)
            continue;
        else
            found++;


        seed_locs.push_back(make_pair(index,rmap_seeds.size()));
        rmap_seeds.push_back(seed[0]);

        for(int ii=0;ii<seed.size();ii++){
            for(int jj=0;jj<kgram_size;jj++){
                seed_det<<"("<<seed[ii][jj].u<<"--"<<seed[ii][jj].d<<"-->"<<seed[ii][jj].v<<")";
            }

            seed_det<<endl;
        }

//        cout<<ol<< " " << index << " " <<numofseeds<<endl;

        seed_det<<endl;

//        if(numofseeds>0){
//            seeds_for_this++;
//            tot_num_seeds+=numofseeds;
//            if(first_seed_start=0)
//                first_seed_start=start;
//            last_seed_finish=start+kgram_size-1;
//        }
//
//
////        check_kgram(kgram,50);
//        seed_rep<<"("<<start<<","<<numofseeds<<")";
////        ","<<in_out.size()<<","<<totkgram/kgram_size<<") ";
//        seed_stats<<numofseeds<<" "<<totkgram/kgram_size<<endl;
//        start=start+1;
//        if(seed.size()<1)
//            no_seeds++;

    }


    vector < vector <int> > final_align, final_align_lengths;
    vector <int> rmap_aligned,rmap_path_lengths;
    rmap_aligned.resize(rmap[ol].size()+1,0);
    rmap_path_lengths.resize(rmap[ol].size()+1,0);

    for(int ii=0;ii<seed_locs.size();ii++){
        cout<<seed_locs[ii].first<<" and "<<seed_locs[ii].second<<endl;
        int pos=0;
        int jj=0;
        for(;jj<rmap_seeds[seed_locs[ii].second].size();jj++){
            rmap_aligned[seed_locs[ii].first+pos]=rmap_seeds[seed_locs[ii].second][jj].u;

            rmap_path_lengths[seed_locs[ii].first+pos]=rmap_seeds[seed_locs[ii].second][jj].d;
            pos++;
        }
        rmap_aligned[seed_locs[ii].first+pos]=rmap_seeds[seed_locs[ii].second][jj-1].v;

    }

    cout<<"Path in graph:"<<endl;
    for(int ii=0;ii<rmap_aligned.size();ii++){
        if(rmap_aligned[ii]>0)
            cout<<reverse_map_nodes[rmap_aligned[ii]]<<" ";
        else
            cout<<"0 ";
    }
    cout<<endl;


    cout<<"Path lengths in graph:"<<endl;
    for(int ii=0;ii<rmap_path_lengths.size();ii++){
        cout<<rmap_path_lengths[ii]<<" ";
    }
    cout<<endl;

    vector < pair<int,int> > start_end;

    int pos1=0,pos2=0;

    while(pos1<rmap_aligned.size() && pos2<rmap_aligned.size()){

        while(pos1<rmap_aligned.size() && rmap_aligned[pos1]==0){
            pos1++;
        }


        if(pos1>=rmap_aligned.size()-2)
            break;

        pos2=pos1;
        while(pos2<rmap_aligned.size() && rmap_aligned[pos2]!=0){
            pos2++;
        }

        start_end.push_back(make_pair(pos1,pos2));

        pos1=pos2;

    }

    cout<<endl<<endl;

    vector <vector < trip_sp > > extension_paths;

    for(int ii=1;ii<start_end.size();ii++){

//        cout<<start_end[ii-1].second<<" "<<start_end[ii].first<<endl;

        int sum=0;

        for(int jj=start_end[ii-1].second-1;jj<=start_end[ii].first-1;jj++){
//            cout<<rmap[ol][jj]<<endl;
            sum+=rmap[ol][jj];
        }

        cout<<endl<<sum<<" "<<reverse_map_nodes[rmap_aligned[start_end[ii-1].second-1]]<<" "<<reverse_map_nodes[rmap_aligned[start_end[ii].first]]<<endl;

        extension_paths = fill_gap(rmap_aligned[start_end[ii-1].second-1],rmap_aligned[start_end[ii].first],sum/2,sum*2,10);

        if(extension_paths.size()>0){

//                cout<<"EXTENSION SUCCESSFUL!"<<endl;

                vector<int> temp;
                vector<int> tem_lengths;

                for(int jj=start_end[ii-1].first;jj<start_end[ii-1].second;jj++){
                    temp.push_back(rmap_aligned[jj]);
                    if(jj<rmap_path_lengths.size())
                        tem_lengths.push_back(rmap_path_lengths[jj]);
                }

                if(extension_paths[0].size()>0){
                    tem_lengths.push_back(extension_paths[0][0].d);
                }

                for(int jj=1;jj<extension_paths[0].size();jj++){
                    temp.push_back(extension_paths[0][jj].u);
                    tem_lengths.push_back(extension_paths[0][jj].d);
                }

//                temp.push_back(extension_paths[0].back().v);

                for(int jj=start_end[ii].first;jj<start_end[ii].second;jj++){
                    temp.push_back(rmap_aligned[jj]);
                    if(jj<rmap_path_lengths.size())
                        tem_lengths.push_back(rmap_path_lengths[jj]);
                }

                final_align.push_back(temp);
                final_align_lengths.push_back(tem_lengths);

        }
        else{
            if(start_end[ii-1].second-start_end[ii-1].first > start_end[ii].second-start_end[ii].first){
                vector<int> temp;
                vector<int> tem_lengths;
                for(int jj=start_end[ii-1].first;jj<start_end[ii-1].second;jj++){
                    temp.push_back(rmap_aligned[jj]);
                    if(jj<rmap_path_lengths.size())
                        tem_lengths.push_back(rmap_path_lengths[jj]);
                }

                final_align.push_back(temp);
                final_align_lengths.push_back(tem_lengths);

            }
            else{
                vector<int> temp;
                vector<int> tem_lengths;
                for(int jj=start_end[ii].first;jj<start_end[ii].second;jj++){
                    temp.push_back(rmap_aligned[jj]);
                    if(jj<rmap_path_lengths.size())
                        tem_lengths.push_back(rmap_path_lengths[jj]);
                }

                final_align.push_back(temp);
                final_align_lengths.push_back(tem_lengths);
            }
        }

    }

    int max_size=0;
    int max_loc=0;

    for(int ii=0;ii<final_align.size();ii++){
        for(int jj=0;jj<final_align[ii].size();jj++){
            cout<<reverse_map_nodes[final_align[ii][jj]]<<" ";
        }
        cout<<endl;

        if(final_align.size()>max_size){
            max_size=final_align.size();
            max_loc=ii;
        }
    }


    if(max_size>0){
        output<<rmap_name[ol]<<"\t";
        output_lengths<<rmap_name[ol]<<"\t";
        for(int jj=0;jj<final_align[max_loc].size();jj++){
                output<<reverse_map_nodes[final_align[max_loc][jj]]<<" ";
//                if(jj<final_align_lengths[max_loc].size())
//                    output_lengths<<final_align_lengths[max_loc][jj]<<" ";
        }

        for(int jj=0;jj<final_align_lengths[max_loc].size()-1;jj++){
//                output<<reverse_map_nodes[final_align[max_loc][jj]]<<" ";
//                if(jj<final_align_lengths[max_loc].size())
                    output_lengths<<final_align_lengths[max_loc][jj]<<" ";
        }

        output<<endl;
        output_lengths<<endl;
    }




//    while(start < rmap[ol].size()-kgram_size+1){
////        cout<<"st"<<start<<" ";
//        kgram.clear();
//        int totkgram=0;
//        for(int i=start;i<start+kgram_size;i++){
//            kgram.push_back(rmap[ol][i]);
//            totkgram+=kgram.back();
//        }
//        tot_kgrams++;
//        int numofseeds = check_kgram(kgram,500);
//
//        for(int ii=0;ii<seed.size();ii++){
//            for(int jj=0;jj<kgram_size;jj++){
//                seed_det<<"("<<seed[ii][jj].u<<"--"<<seed[ii][jj].d<<"-->"<<seed[ii][jj].v<<")";
//            }
//
//            seed_det<<endl;
//        }
//
//        set <int> in_out;
//
//        vector < pair <int,int> > temp;
//
//        for(int ii=0;ii<seed.size();ii++){
//
//            int sum=seed[ii][0].u + seed[ii].back().v;
//            if(in_out.find(sum)==in_out.end()){
//                in_out.insert(sum);
//                seed_det<<seed[ii][0].u<<","<<seed[ii].back().v<<" || ";
//                temp.push_back(make_pair(seed[ii][0].u,seed[ii].back().v));
//            }
//
//        }
//
//        seeds_of_rmap.push_back(temp);
//
////        cout<<"seeds_of_rmap"<<seeds_of_rmap.size()<<endl;
//
//        seed_det<<endl;
//
//        if(numofseeds>0){
//            seeds_for_this++;
//            tot_num_seeds+=numofseeds;
//            if(first_seed_start=0)
//                first_seed_start=start;
//            last_seed_finish=start+kgram_size-1;
//        }
//
//
////        check_kgram(kgram,50);
//        seed_rep<<"("<<start<<","<<numofseeds<<")";
////        ","<<in_out.size()<<","<<totkgram/kgram_size<<") ";
//        seed_stats<<numofseeds<<" "<<totkgram/kgram_size<<endl;
//        start=start+1;
//        if(seed.size()<1)
//            no_seeds++;
////        for(int j=0;j<seed.size();j++){
////            seed_set[ol].rmap_seed.push_back(seed[j]);
////        }
//
////        vector < pair <float,float> > temp;
////        int cnt_good_seeds=0;
////        for(int m=0;m<seed.size();m++){
////            float this_score=0;
////            int tot_frag_length=0;
////            for(int n=0;n<seed[m].size();n++){
////                tot_frag_length+=rmap[ol][old_st+n];
////                this_score+=abs(seed[m][n].d-rmap[ol][old_st+n]);
////            }
////            this_score=(float)this_score/seed[m].size();
////            if(this_score<=10)
////                cnt_good_seeds++;
////            seed_set[ol].scores.push_back(make_pair(this_score,old_st));
////            temp.push_back(make_pair(this_score,(float)this_score/tot_frag_length));
////
//        }

//        int beg=0;
//        vector <int> rmap_aligned;
//        rmap_aligned.resize(rmap[ol].size(),0);
//
//        while(beg<seeds_of_rmap.size()){
//
//            if(seeds_of_rmap[beg].size()==0){
//                beg++;
//                continue;
//            }
//
//            if(beg+kgram_size+4>=seeds_of_rmap.size())
//                break;
//
//            seed_rep<<"seed under consideration "<<beg << " ";
//                        cout<<"seed under consideration "<<beg << " ";
//
//            int a;
//
//            for(a=beg+kgram_size+2;a<seeds_of_rmap.size();a++)
//                if(seeds_of_rmap[a].size()>0)
//                    break;
//
//            if(a==seeds_of_rmap.size()){
//                beg=a;
//                continue;
//            }
//
//
//            seed_rep<<"next seed"<< a<<endl;
//            cout<<"next seed"<< a<<endl;
//
//
//            int d2=0;
//            int d1=0;
//            int d;
//            int gap_len=0;
//
//            vector < int > gap;
//
//            for(int i=beg+kgram_size;i<a;i++){
//                cout << rmap[ol][i]<< " ";
//
//                d2+=rmap[ol][i];
//                gap.push_back(rmap[ol][i]);
//                gap_len++;
//            }
//
//            if(gap_len>6){beg++;continue;}
//            d1=d2*0.8;
//            d=
////            cout<<"frags in between:"<<gap_len;
//            d2=d2*1.2;
//
//
//
//            set <int> tem_set1;
//            set <int> tem_set2;
//            for(int i=0;i<seeds_of_rmap[beg].size();i++){
//                if(tem_set1.find(seeds_of_rmap[beg][i].second)==tem_set1.end()){
//                    tem_set1.insert(seeds_of_rmap[beg][i].second);
//                }
//            }
//            for(int ni=0;ni<seeds_of_rmap[a].size();ni++){
//                if(tem_set2.find(seeds_of_rmap[a][ni].first)==tem_set2.end()){
//                    tem_set2.insert(seeds_of_rmap[a][ni].first);
//                }
//            }
//
//            set < int >::iterator it1;
//            set < int >::iterator it2;
//            int num_paths_found=0;
//            vector<int> bst_que;
//            int max_sc=0;
//
//            int flag=0;
//
//            int fills_attempt=0;
//            vector < int >  best_path;
//            float best_score=0;
//
//            for(it1=tem_set1.begin();it1!=tem_set1.end();++it1){
//                for(it2=tem_set2.begin();it2!=tem_set2.end();++it2){
//                vector< vector < trip_sp > > paths_found;
////                seed_rep<<"start node: "<<*it1<<" end node: "<<*it2<<endl;
//
//                int relax=gap_len*500;
//                paths_found = fill_gap(*it1,*it2,d-relax,d+relax,gap_len+3);
//                if(paths_found.size()>0)
//                    seed_rep<<endl;
//                for(int i=0;i<paths_found.size();i++){
//
//                    vector<int> path;
//                    for(int j=0;j<paths_found[i].size();j++){
//                        path.push_back(paths_found[i][j].d);
//                        seed_rep<< paths_found[i][j].u<<" --> "<<paths_found[i][j].v<<" ("<<paths_found[i][j].d<<") ";
//                    }
//
//                    seed_rep<<endl;
//
//                    float oal = find_opt_alignment(gap,path);
//                    if(oal>best_score){
//                        best_score=oal;
//                        best_path=path;
//                    }
//                }
//                fills_attempt++;
//                if(paths_found.size()>0 && best_score>0)
//                    flag++;
////                cout<<" fill_gap("<<*it1<<","<<*it2<<","<<d2<<","<<gap_len+3<<") ="<<paths_found.size() ;
//                num_paths_found=num_paths_found+paths_found.size();
//                }
//            }
//
//
//
//            cout<<"fills attempted:"<<fills_attempt<<endl;
//
//            if(flag>0){
//                score_file<<(float)best_score/gap_len<<endl;
//                best_paths_view<<"Rmap:"<<ol<<" Seeds:"<<beg<<","<<a<< " Best score:"<<best_score<<endl;
//                for(int i=0;i<gap.size();i++){
//                    best_paths_view<<gap[i]<<" ";
//                }
//                best_paths_view<<endl;
//
//                for(int i=0;i<best_path.size();i++){
//                    best_paths_view<<best_path[i]<<" ";
//                }
//                best_paths_view<<endl;
//
//
//                for(int x=beg;x<a+kgram_size;x++){
//                    rmap_aligned[x]++;
//                }
//                int x=beg;
//                while(x>=0 && seeds_of_rmap[x].size()>0){
//                    rmap_aligned[x--]++;
//                }
//                x=a;
//                while(x<seeds_of_rmap.size() && seeds_of_rmap[x].size()>0){
//                    rmap_aligned[kgram_size+x-1]++;
//                    x++;
//                }
//
//                cout<<"gap filled. alignment length:"<<gap_len+2*kgram_size<<endl;
//
//                beg=a;
//
//
//
//            }
//            else
//                beg++;
//
//
//
//
//        }
//
//
//
//
//
//        cout<<"rmap_aligned.size():"<<rmap_aligned.size()<<" seeds_of_rmap.size():"<<seeds_of_rmap.size()<<endl;
//
//        int max_ali=0;
//        cout<<endl;
//        for(int x=0;x<rmap_aligned.size();x++){
//            cout<<rmap_aligned[x]<<" ";
//
//            if(rmap_aligned[x]==0)
//                continue;
//            int y;
//            for(y=x;y<rmap_aligned.size();y++){
//                if(rmap_aligned[y]==0)
//                    break;
//            }
//
//            if(max_ali<(y-x))
//                max_ali=y-x;
//
//        }
//        per_alig<<max_ali*100/rmap[ol].size()<<endl;
//
//        if(max_ali*100/rmap[ol].size()>=60)
//            percent_aligned++;
//        cout<<"percent of rmap aligned:"<<max_ali*100/rmap[ol].size();
//
//        seed_coverage.push_back((last_seed_finish-first_seed_start)*100/rmap[ol].size());
//
//        total_coverage+=seed_coverage.back();
//
//        seed_cov<<seeds_for_this<<" "<<((last_seed_finish-first_seed_start)*100/rmap[ol].size())<<endl;
//        seed_cov<<"Total kgrams:"<<tot_kgrams<< " No seeds: "<<no_seeds<<" Average num of seeds:"<<tot_num_seeds/tot_kgrams;
//




//        sort(temp.begin(),temp.end());
//        int max_score=0;
//        for(int te=0;te<10;te++){
//            if(te>=temp.size())
//                break;
//            all_scores.push_back(temp[te].first);
//            all_scores_scaled.push_back(temp[te].second);
//        }
//        max_score=all_scores.back();
//
//        seed_set[ol].num_seeds.push_back(seed.size());
//        seed_set[ol].seed_start_end.push_back(make_pair(old_st,start-1));
//
//        cout<<"("<<old_st<<","<<old_st+kgram_size-1<<" ->"<<max_score<< ")";
//        if(best_seeds[ol].size()==0 || old_st-best_seeds[ol].back().v>2){
//            best_seeds[ol].push_back(trip_sp(old_st,old_st+kgram_size-1,max_score));
//        }
//
//        if(cnt_good_seeds>0)
//            seed_rep<<"("<<old_st<<","<<old_st+kgram_size-1<<" ->"<<cnt_good_seeds<< ")";
//    }
//
//    sort(seed_set[ol].scores.begin(), seed_set[ol].scores.end());
//    set<int> temp;
//    for(int m=0;m<top_seeds;m++){
//        if(m==seed_set[ol].scores.size())
//            break;
//        seed_set[ol].final_seeds.insert(seed_set[ol].scores[m].second);
//    }
//    cout<<endl;
//    set<int>::iterator it;
//    vector< pair<int,int> > tem;
//    for(it=seed_set[ol].final_seeds.begin();it!=seed_set[ol].final_seeds.end();++it){
//        int st=seed_set[ol].seed_start_end[*it].first;
//        int en=seed_set[ol].seed_start_end[*it].second;
//        cout<<"("<<seed_set[ol].seed_start_end[*it].first<<","<<seed_set[ol].seed_start_end[*it].second<< ")";
//        tem.push_back(make_pair(st,en));
//    }
//    seed_file<<tem.size()<<endl;
//    if(tem.size()>1){
//
//            for(int m=0;m<tem.size()-1;m++){
//            int gap_length=0;
//            int st=tem[m].second+1;
//            while(st<tem[m+1].first){
//                gap_length+=rmap[ol][st++];
//            }
//            hist_file<<ol<<" "<<gap_length<<endl;
//        }
//    }
//
//
//    for(int m=0;m<seed_set[ol].scores.size();m++){
//
//
//        if(seed_set[ol].scores[m].first<=10){
//            cout<<seed_set[ol].scores[m].first<<" "<<seed_set[ol].scores[m].second<<endl;
//            seed_set[ol].final_seeds.insert(seed_set[ol].scores[m].second);
//        }
//
//    }


    cout<<endl;
}

//cout<<"Total kgrams:"<<tot_kgrams<< " No seeds: "<<no_seeds<<" Average num of seeds:"<<tot_num_seeds/tot_kgrams;
//
//cout<< endl<< "Average coverage"<<total_coverage/rmap.size();
//
//cout<<"Percent of Rmaps aligned : "<<percent_aligned*100/rmap.size();

return 1;

ofstream allsc;
allsc.open("all_scores");
ofstream allscscaled;
allscscaled.open("scores_scaled");

cout<<"Total kgrams"<<tot_kgrams;

for(int i=0;i<all_scores.size();i++){
        allsc<<all_scores[i]<<endl;
        allscscaled<<all_scores_scaled[i]<<endl;
}

//for(int m=0;m<seed_set[0].scores.size();m++){
//    cout<<seed_set[0].scores[m].first<<" ("<<seed_set[0].seed_start_end[seed_set[0].scores[m].second].first<<","<<seed_set[0].seed_start_end[seed_set[0].scores[m].second].second<<")"<<endl;
//}

ofstream bstseed;
bstseed.open("best_seeds");

ofstream gap_file;
gap_file.open("gap_lengths");

for(int i=0;i<rmap.size();i++){
    bstseed<<"Rmap "<<i<<endl;
    int gap_len=0;
    for(int j=0;j<best_seeds[i].size();j++){
        if(j>0){
            for(int jj=best_seeds[i][j-1].v;jj<best_seeds[i][j].u;jj++)
                gap_len+=rmap[i][jj];
        gap_file<<gap_len<<endl;
        }

        bstseed<<"("<<best_seeds[i][j].u<<","<<best_seeds[i][j].v<<")";
    }


    bstseed<<endl;


}

time_t end_seed=clock();

cout<<"seeding time : "<<double( end_seed-start_seed ) / CLOCKS_PER_SEC<<endl;


//return 1;

trim_seeds();

time_t end_trim=clock();

cout<<"trimming time : "<<double( end_trim-end_seed) / CLOCKS_PER_SEC<<endl;

//for(int i=0;i<seed_set[0].seed_start_end.size();i++){
//    cout<<seed_set[0].seed_start_end[i].first<<" "<<seed_set[0].seed_start_end[i].second<<endl;
//}
//
//for(int i=0;i<seed_set[0].begins.size();i++){
//
//    cout<<seed_set[0].begins[i].size()<<" "<<seed_set[0].ends[i].size()<<endl;
//}

//Extension now

int gaps=0, no_paths=0;
vector<int> path_best_scores;

int total_paths_found=0;
ofstream scores;
scores.open("scores");
ofstream offile;
offile.open("extension_stats");

int alignment_called=0;
int total_paths_filled=0;
for(int i=0;i<rmap.size();i++){
    offile<<"Rmap: "<<i;
    if(seed_set[i].seed_start_end.size()<2){
            offile<<" Not enough seeds!"<<endl;
            continue;
    }

    int no_fill_gap=0;
    int paths_found_this_rmap=0;
    for(int j=0;j<seed_set[i].seed_start_end.size()-1;j++){
        gaps++;
        int gap_start_in=seed_set[i].seed_start_end[j].second;
        int gap_end_in=seed_set[i].seed_start_end[j+1].first;
        if(gap_end_in-gap_start_in>8){
            no_paths++;
            continue;

        }

        int d=0;
        vector <int> tar;
        for(int k=gap_start_in;k<=gap_end_in;k++){
            d+=rmap[i][k];
            tar.push_back(rmap[i][k]);
        }

        pair< int, int > tem=find_h(d,0.1);
        int d1=d*0.1,d2=d*3;

        set < int >::iterator it1;
        set < int >::iterator it2;
        int num_paths_found=0;
        vector<int> bst_que;
        int max_sc=0;

        for(it1=seed_set[i].ends[j].begin();it1!=seed_set[i].ends[j].end();++it1){
            for(it2=seed_set[i].begins[j+1].begin();it2!=seed_set[i].begins[j+1].end();++it2){
                vector< vector < trip_sp > > paths_found;
//                cout<<"start node: "<<*it1<<" end node: "<<*it2<<endl;
                paths_found = fill_gap(*it1,*it2,0,d2,gap_end_in-gap_start_in+3);
                num_paths_found=num_paths_found+paths_found.size();
                total_paths_found=total_paths_found+paths_found.size();
//                cout<< total_paths_found <<" ";


                if(paths_found.size()>0){
                        paths_found_this_rmap++;
                        total_paths_filled++;
//                    for(int k=0;k<paths_found.size();k++){
//                        vector < int > que;
//                        for(int l=0;l<paths_found[k].size();l++){
//                            que.push_back(paths_found[k][l].d);
//                        }
//                        float al_sc=find_opt_alignment(tar,que);
//                        alignment_called++;
//                        if(al_sc>max_sc){
//                            bst_que=que;
//                            max_sc=al_sc;
//                        }
//
//
//                    }
                }



            }
        }
        if(num_paths_found>0){
//            cout << endl<<"paths found : "<<num_paths_found<<" dis:"<<d<<" d1: "<<d1<<" d2:"<<d2<<" hops: "<<gap_end_in-gap_start_in<<endl;
//            cout<<"highest scoring path :"<<max_sc<<endl;
            path_best_scores.push_back(max_sc);
            scores<<max_sc<<endl;
            no_fill_gap++;
//            for(int i=0;i<bst_que.size();i++){
//                cout<<bst_que[i]<<" ";
//            }
//            cout<<endl;
//
//            for(int i=0;i<tar.size();i++){
//                cout<<tar[i]<<" ";
//            }


        }
        else
            no_paths++;
    }
    offile<<" Gaps: "<<seed_set[i].seed_start_end.size()-1<<" Path found: "<<paths_found_this_rmap<<endl;
//    cout<<"***********Rmap:"<<i<<endl;
//    <<" Total gaps: "<<gaps<<" Gaps for which no paths found: "<<no_paths<<endl;
}

time_t end_extension=clock();

cout<<"extension time : "<<double( end_extension-end_trim) / CLOCKS_PER_SEC<<endl;


cout<<" Total gaps: "<<gaps<<endl<<" Gaps for which no paths found: "<<total_paths_filled<<endl<<" Total number of paths evaluated:"<<alignment_called<<endl;
//trim_seeds();

  return 1;

//  ifstream loop_file;
//  string line;
//  loop_file.open("loop_output/all_simple_loops_500");
//  int ln;
//  while(getline(loop_file,line)){
//    istringstream ss(line);
//    ss>>ln;
//    loopy_nodes.insert(ln);
////    cout<<ln<<endl;
//  }
//  string input_filename="/ufrc/boucher/kingdgp/cosmo/cosmo/experiments/32.list.packed";
//
//  ifstream input(input_filename, ios::in|ios::binary|ios::ate);
//  // Can add this to save a couple seconds off traversal - not really worth it.
//  //vector<size_t> minus_positions;
//  debruijn_graph<> dbg = debruijn_graph<>::load_from_packed_edges(input, "$ACGT"/*, &minus_positions*/);
//  input.close();
//
//
////  return 1;
//
////  string restriction_enzyme = p.restriction_seq;
//
//  cout<<dbg.num_nodes()<<endl;
//  cout<<dbg.num_edges()<<endl;
//
//
//  ifstream infile;
//  infile.open("restriction_nodes_32");
//  int new_node;
//
//  vector <int> rest_nodes_vec;
//
//  string str;
//
//  while(getline(infile, str))
//    {
//        istringstream ss(str);
//
//        ss >> new_node;
//        rest_nodes.insert(new_node);
//        rest_nodes_vec.push_back(new_node);
//
//    }
////    ofstream offile;
////    offile.open("DBG_dist_32");
////    for(int i=0;i<rest_nodes_vec.size();i++){
////        offile<<dbg.outdegree(rest_nodes_vec[i])<<endl;
////    }
////    return 1;
//    cout<<"node under review:"<<num-1<<" "<<rest_nodes_vec[num-1]<<endl;
//
//
////  cout<< "num rest nodes:"<<rest_nodes.size()<<endl;
//    infile.close();
//
//    int run=0;
//
////    for (it=rest_nodes.begin(); it!=rest_nodes.end(); ++it){
////        cout<<"-"<<++run<<"-"<<endl;importdata('C:\Users\King\Desktop\Rmaps\rmaps2\rmap_relations_stats.txt');
//////        find_all_paths_dp(dbg, 17320  , 25000);
////        for(int i=0;i<rest_nodes_vec.size();i++){
////            find_all_paths_dp(dbg, rest_nodes_vec[i] , 50000);
////        }
//        find_all_paths_dp(dbg, rest_nodes_vec[num-1] , 50000);
//
////   }
////    int sd=0;
////    int cnt=0;
////    int maxb=0;
////    ofstream explore;
////    explore.open("explore");
////    for(int i=0;i<visited_dis.size();i++){
////        if(visited_dis[i].size()==0)
////            continue;
////        if(visited_dis[i].back()-visited_dis[i][0]>maxb)
////            maxb=visited_dis[i].back()-visited_dis[i][0];
////        sd+=visited_dis[i].back()-visited_dis[i][0];
////        cnt++;
//////        explore<<i<<"\t";
//////        for(int j=0;j<visited_dis[i].size();j++)
//////            explore<<visited_dis[i][j]<<" ";
////
//////        explore<<endl;
////    }
//
////    cout<<"Average bulge size: " <<sd/cnt<<endl;
////    cout<<"Max bulge size: " <<maxb<<endl;
//
//    ofstream ofile;
//    string f_name="edges/"+nu+".edges";
//    ofile.open(f_name.c_str());
//
//    for(int i=0;i<all_edges.size();i++){
//            ofile<<all_edges[i].source<<" "<<all_edges[i].target<<" "<<all_edges[i].distance<<endl;
//
//    }
//
//    time_t oend = clock();
//
//    ofstream time_file;
//    string t_file="edges/"+nu+".time";
//    time_file.open(t_file.c_str());
//    time_file<<double(oend - ostart) / CLOCKS_PER_SEC<<endl;
//    time_file.close();
//
//    return 1;
//
//    ofile.close();
//    ofile.open("mgraph.txt");
//
//    ofstream ofile2;
//    ofile2.open("tracker.txt");
//
//    run=0;
//
//    for (it=rest_nodes.begin(); it!=rest_nodes.end(); ++it){
//            run++;
////            if (run<9)
////                continue;
//            if (run>10000)
//                break;
//          nodes_visited=0;
//          ofile2<<"origin:"<<*it<<" ";
//
//
//    }
//
//     for(int i=0;i<Medges.size();i++){
//            ofile<<Medges[i].source<<" "<<Medges[i].target<<" "<<Medges[i].distance<<endl<<endl;
//
//    }
//
////  find_restriction_kmers(dbg, restriction_enzyme);
//
//
//
//    ofile.close();

}
