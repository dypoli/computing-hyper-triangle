#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <set>
#include <map>
#include <ctime>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <cassert>
#include <bitset>
#include <list>
#include <cstdint>
//#include "..\console_utf8.h"
#include "read_data.cpp"
#include "HelperFunc.h"

using namespace std;

vector<vector<int>> ForTest;


int main(int argc, char *argv[])
{
    //init_console_utf8();
	clock_t start;
	clock_t run_start;
	int progress;

    
    string dataName = "unique-email-Enron.txt";

    string graphFile =dataName;
    // ==============================


	// Read data
	start = clock();
	vector< vector<int> > node2hyperedge;
	vector< vector<int> > hyperedge2node;
	vector< unordered_set<int> > hyperedge2node_set;
    vector< unordered_set<int> > node2hyperedge_set;
	read_data(graphFile, node2hyperedge, hyperedge2node, hyperedge2node_set);

	


	int V = node2hyperedge.size(), E = hyperedge2node.size();
	cout << "# of nodes: " << V << '\n';
	cout << "# of hyperedges: " << E << '\n';
	cout << "Reading original data done: "
		<< (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "------------------------------------------" << endl << endl;



	// Make adjacency list
	start = clock(); run_start = clock();
	hyperedge2node.resize(E); hyperedge2node_set.resize(E);
	vector<vector<int>> hyperedge_adj;

	hyperedge_adj.resize(E);
	vector<long long> upd_time(E, -1LL);


    // new structure, only for this version
    
    vector<int> Index_signal(E,-1);
    vector<vector<vector<int>>> Edge_node_Edge;
    Edge_node_Edge.resize(E);


    vector<vector<vector<int>>> Edge_node_Edge_count;
    Edge_node_Edge_count.resize(E);



    for(int i=0; i<E; i++){
        Edge_node_Edge[i].resize(hyperedge2node[i].size());
        Edge_node_Edge_count[i].resize(hyperedge2node[i].size());

        for(int j=0; j<hyperedge2node[i].size(); j++){
            Edge_node_Edge[i][j].resize(node2hyperedge[hyperedge2node[i][j]].size(),0);
            Edge_node_Edge_count[i][j].resize(node2hyperedge[hyperedge2node[i][j]].size(),0);

        }
    } 

	

    vector<unordered_map<int, uint32_t>> hyperedge_bitmap(E);


    


	for (int hyperedge_a = 0; hyperedge_a < E; hyperedge_a++){
		long long l_hyperedge_a = (long long)hyperedge_a;
		for (int i=0; i<hyperedge2node[hyperedge_a].size(); i++){
            int node=hyperedge2node[hyperedge_a][i];
			for (int j=0; j<node2hyperedge[node].size(); j++){
                int hyperedge_b=node2hyperedge[node][j];
				if (hyperedge_b == hyperedge_a) continue;
				if ((upd_time[hyperedge_b] >> 31) ^ hyperedge_a){
					hyperedge_adj[hyperedge_b].push_back(hyperedge_a);
					Edge_node_Edge[hyperedge_a][i][j]=hyperedge_b;
                }
			}
		}
	}

   
    

	// and convert the adjacent list to set for faster query
	removeDuplicatesFromHyperedgeAdj(hyperedge_adj);
    convertToVectorWithIndex(hyperedge_adj,Index_signal);
    vector<vector<int>> hyperedge_TriangleNum;
    hyperedge_TriangleNum.resize(E);
    for(int i=0; i<E; i++) hyperedge_TriangleNum[i].resize(hyperedge_adj[i].size(),0);

    vector<vector<int>> hyperedge_TriangleNum2;
	hyperedge_TriangleNum2.resize(E);
    for(int i=0; i<E; i++) hyperedge_TriangleNum2[i].resize(hyperedge_adj[i].size(),0);



    cout << "The raw data is read to obtain the original adjacency list, and this part is considered input in the algorithm.\n"
    << (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
    cout << "------------------------------------------" << endl << endl;



    string TrianglNum= "TriangleNumber.txt";
	  vector<int> h_triangle=ReadTriangleFile(TrianglNum2); 
    // =====================================================
    // Load deleted hyperedge IDs
    // =====================================================

    string graphFile3 ="DeleteHyperedge.txt";
  

    vector<int> deleted_edges;

    ifstream fin_del(graphFile3);

    int del_id;

    while (fin_del >> del_id) {
        deleted_edges.push_back(del_id);
    }

    fin_del.close();
    sort(deleted_edges.begin(), deleted_edges.end());


    // =====================================================
    // Process deletions one by one
    // =====================================================

    for (int delete_e : deleted_edges) {

        // invalid
        if (delete_e < 0 || delete_e >= hyperedge2node.size())
            continue;

        // already deleted
        if (hyperedge2node[delete_e].size() == 0)
            continue;

        // =================================================
        // 1. Find triangles BEFORE deletion
        // =================================================

        for (int e1 : hyperedge_adj[delete_e]) {

            // skip deleted neighbor
            if (hyperedge2node[e1].size() == 0)
                continue;

            for (int e2 : hyperedge_adj[e1]) {

                // enforce ordering
                if (e2 <= e1)
                    continue;

                // skip deleted neighbor
                if (hyperedge2node[e2].size() == 0)
                    continue;

                // check adjacency
                if (binary_search(
                        hyperedge_adj[delete_e].begin(),
                        hyperedge_adj[delete_e].end(),
                        e2
                    )) {

                    vector<int> patternNum =
                        get_specificPattern(
                            hyperedge2node[delete_e],
                            hyperedge2node[e1],
                            hyperedge2node[e2]
                        );

                    if (patternNum[0] > 0) {

                        h_triangle[patternNum[0] - 1] -= 1;
                    }
                }
            }
        }

        // =================================================
        // 2. Remove delete_e from neighbors' adjacency
        // =================================================

        for (int nb : hyperedge_adj[delete_e]) {

            auto &adj = hyperedge_adj[nb];

            adj.erase(
                remove(adj.begin(), adj.end(), delete_e),
                adj.end()
            );
        }

        // =================================================
        // 3. Remove delete_e from node2hyperedge
        // =================================================

        for (int v : hyperedge2node[delete_e]) {

            auto &incident = node2hyperedge[v];

            incident.erase(
                remove(incident.begin(), incident.end(), delete_e),
                incident.end()
            );
        }

        // =================================================
        // 4. Clear deleted hyperedge
        // =================================================

        hyperedge_adj[delete_e]={};

        hyperedge2node[delete_e]={};
    }

    

    cout << "update time.\n"
    << (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
    cout << "------------------------------------------" << endl << endl;
	return 0;
}
