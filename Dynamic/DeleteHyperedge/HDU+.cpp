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
    // Build deleted mark
    // =====================================================

    vector<bool> is_deleted(E, false);

    for (int e : deleted_edges) {

        if (e >= 0 && e < E)
            is_deleted[e] = true;
    }

    // =====================================================
    // Build:
    // 1. hyperedge_adj_delete
    // 2. hyperedge_adj_adjust
    // 3. updated hyperedge_adj
    // =====================================================

    vector<vector<int>> hyperedge_adj_delete(E);

    vector<vector<int>> hyperedge_adj_adjust(E);

    vector<vector<int>> hyperedge_adj_new(E);


    vector<vector<int>> hyperedge_adj_old=hyperedge_adj;

    start = clock();

    // =====================================================
    // Split adjacency
    // =====================================================

    for (int e = 0; e < E; e++) {

        for (int nb : hyperedge_adj[e]) {

            // avoid duplicate processing
            if (nb < e)
                continue;

            // ============================================
            // Case 1:
            // delete - delete
            // ============================================

            if (is_deleted[e] && is_deleted[nb]) {

                hyperedge_adj_delete[e].push_back(nb);
                hyperedge_adj_delete[nb].push_back(e);
            }

            // ============================================
            // Case 2:
            // delete - remain
            // ============================================

            else if (is_deleted[e] ^ is_deleted[nb]) {

                hyperedge_adj_adjust[e].push_back(nb);
                hyperedge_adj_adjust[nb].push_back(e);
            }

            // ============================================
            // Case 3:
            // remain - remain
            // → new hyperedge_adj
            // ============================================

            else {

                hyperedge_adj_new[e].push_back(nb);
                hyperedge_adj_new[nb].push_back(e);
            }
        }
    }

    // =====================================================
    // Replace old adjacency
    // =====================================================

    hyperedge_adj = hyperedge_adj_new;

    // =====================================================
    // Sort
    // =====================================================
    /*
    for (int i = 0; i < E; i++) {

        sort(
            hyperedge_adj_delete[i].begin(),
            hyperedge_adj_delete[i].end()
        );

        sort(
            hyperedge_adj_adjust[i].begin(),
            hyperedge_adj_adjust[i].end()
        );

        sort(
            hyperedge_adj[i].begin(),
            hyperedge_adj[i].end()
        );
    }*/


    ///////////////////////////////// type-1/ ////////////////////////////////////////////////
    for(int i=0; i<hyperedge_adj_old.size(); i++){
        int hyperedge_1=i;
        if (!is_deleted[hyperedge_1])  continue;


		for (int x = 0; x < (int)hyperedge_adj_adjust[i].size(); x++) {
            int hyperedge_2 = hyperedge_adj_adjust[i][x];




            int p = 0; 
            int q = 0;      
            
            while (p < hyperedge_adj_adjust[hyperedge_1].size() && q < hyperedge_adj_old[hyperedge_2].size()) {
                int a = hyperedge_adj_adjust[hyperedge_1][p];
                int b = hyperedge_adj_old[hyperedge_2][q];

                if (a == b ) {
                    int hyperedge_3 = a;
                    if(hyperedge_3>hyperedge_2){
                        vector<int> patternNum=get_specificPattern(hyperedge2node[hyperedge_1],hyperedge2node[hyperedge_2],hyperedge2node[hyperedge_3]);
                        if(patternNum[0]>0){
                            h_triangle[patternNum[0]-1]-=1;
                        }
                    }
                    p++; q++;
                }
                else if (a < b) {
                    p++;
                } else {
                    q++;
                }
            }
        }

    }



    ///////////////////////////////// type-2/ ////////////////////////////////////////////////
    for(int i=0; i<hyperedge_adj_old.size(); i++){
        int hyperedge_1=i;
        if (is_deleted[hyperedge_1])  continue;


		for (int x = 0; x < (int)hyperedge_adj_adjust[i].size(); x++) {
            int hyperedge_2 = hyperedge_adj_adjust[i][x];




            int p = 0;  
            int q = 0;     
            
            while (p < hyperedge_adj_adjust[hyperedge_1].size() && q < hyperedge_adj_delete[hyperedge_2].size()) {
                int a = hyperedge_adj_adjust[hyperedge_1][p];
                int b = hyperedge_adj_delete[hyperedge_2][q];

                if (a == b ) {
                    int hyperedge_3 = a;
                    if(hyperedge_3>hyperedge_2){
                        vector<int> patternNum=get_specificPattern(hyperedge2node[hyperedge_1],hyperedge2node[hyperedge_2],hyperedge2node[hyperedge_3]);
                        if(patternNum[0]>0){
                            h_triangle[patternNum[0]-1]-=1;
                        }
                    }
                    p++; q++;
                }
                else if (a < b) {
                    p++;
                } else {
                    q++;
                }
            }
        }

    }


    ///////////////////////////////// type-3/ ////////////////////////////////////////////////
    for(int i=0; i<hyperedge_adj_old.size(); i++){
        int hyperedge_1=i;
        if (!is_deleted[hyperedge_1])  continue;


		for (int x = 0; x < (int)hyperedge_adj_delete[i].size(); x++) {
            int hyperedge_2 = hyperedge_adj_delete[i][x];

            if(hyperedge_1>=hyperedge_2) continue;


            int p = 0;  
            int q = 0;      
            
            while (p < hyperedge_adj_delete[hyperedge_1].size() && q < hyperedge_adj_delete[hyperedge_2].size()) {
                int a = hyperedge_adj_delete[hyperedge_1][p];
                int b = hyperedge_adj_delete[hyperedge_2][q];

                if (a == b ) {
                    int hyperedge_3 = a;
                    if(hyperedge_3>hyperedge_2){
                        vector<int> patternNum=get_specificPattern(hyperedge2node[hyperedge_1],hyperedge2node[hyperedge_2],hyperedge2node[hyperedge_3]);
                        if(patternNum[0]>0){
                            h_triangle[patternNum[0]-1]-=1;
                        }
                    }
                    p++; q++;
                }
                else if (a < b) {
                    p++;
                } else {
                    q++;
                }
            }
        }

    }

  
   



    cout << "update time.\n"
    << (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
    cout << "------------------------------------------" << endl << endl;
	return 0;


}
