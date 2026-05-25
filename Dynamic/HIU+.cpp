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

    
    string dataName = "unique-high-primary-school.txt";

    string graphFile =
        "C:\\AllCode\\Triangle Counting\\HyperedgeAdd\\Updated_Dataset\\"
        + dataName;
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

	vector< vector<int> > hyperedge_new(E);
	string updateVertex = "C:\\AllCode\\Triangle Counting\\HyperedgeAdd\\Input_Data\\AffectedHE_"+dataName;
	Load_New_Hyperedges(updateVertex, hyperedge_new);




	string TrianglNum2= "C:\\AllCode\\Triangle Counting\\HyperedgeAdd\\Input_Data\\TriangleNum2_"+dataName;
	vector<int> h_triangle=ReadTriangleFile(TrianglNum2); 




    int total_hyperedges = hyperedge_new.size();

    hyperedge2node.resize(total_hyperedges);

    hyperedge_adj.resize(total_hyperedges);

    // =====================================================
    // New adjacency structures
    // =====================================================

    // Connections among new hyperedges only
    vector<vector<int>> hyperedge_adj_new(total_hyperedges);

    // Connections involving at least one new hyperedge
    vector<vector<int>> hyperedge_adj_update(total_hyperedges);

    vector<vector<int>> hyperedge_adj_old=hyperedge_adj;

    for (int new_e = 0; new_e < hyperedge_new.size(); new_e++) {

        // Skip empty hyperedges
        if (hyperedge_new[new_e].empty())
            continue;

        // =====================================================
        // 1. Insert into hyperedge2node
        // =====================================================

        hyperedge2node[new_e] = hyperedge_new[new_e];

        // =====================================================
        // 2. Find neighboring hyperedges
        // =====================================================

        unordered_set<int> neighbor_set;

        for (int v : hyperedge_new[new_e]) {

            for (int old_e : node2hyperedge[v]) {

                if (old_e == new_e)
                    continue;

                neighbor_set.insert(old_e);
            }
        }

        // =====================================================
        // 3. Build adjacency list
        // =====================================================

        vector<int> neighbors(
            neighbor_set.begin(),
            neighbor_set.end()
        );

        sort(
            neighbors.begin(),
            neighbors.end()
        );

        hyperedge_adj[new_e] = neighbors;

        // =====================================================
        // 4. Update extra adjacency structures
        // =====================================================

        for (int nb : neighbors) {

            // =============================================
            // A. New-New adjacency
            // Both hyperedges are new
            // =============================================

            if (!hyperedge_new[nb].empty()) {

                hyperedge_adj_new[new_e].push_back(nb);

                hyperedge_adj_new[nb].push_back(new_e);
            }

            // =============================================
            // B. Old-New adjacency
            // Exactly one hyperedge is new
            // =============================================

            if (hyperedge_new[nb].empty()) {

                hyperedge_adj_update[new_e].push_back(nb);

                hyperedge_adj_update[nb].push_back(new_e);
            }
            // =============================================
            // Original adjacency
            // =============================================

            hyperedge_adj[nb].push_back(new_e);
        }

        // =====================================================
        // 5. Update node2hyperedge
        // =====================================================

        for (int v : hyperedge_new[new_e]) {

            node2hyperedge[v].push_back(new_e);
        }
    }

    // =====================================================
    // 6. Sort all new adjacency structures
    // =====================================================

    for (int e = 0; e < total_hyperedges; e++) {

        sort(
            hyperedge_adj_new[e].begin(),
            hyperedge_adj_new[e].end()
        );

        sort(
            hyperedge_adj_update[e].begin(),
            hyperedge_adj_update[e].end()
        );
    }


 	///////////////////////// find new hyper-triangle ///////////////////////////////////



    ///////////////////////////////// type-1/ ////////////////////////////////////////////////
    for(int i=0; i<hyperedge_new.size(); i++){
        int hyperedge_1=i;
        if(hyperedge_new[i].size()==0) continue;


		for (int x = 0; x < (int)hyperedge_adj_update[i].size(); x++) {
            int hyperedge_2 = hyperedge_adj_update[i][x];




            int p = 0; 
            int q = 0;   
            
            while (p < hyperedge_adj_update[hyperedge_1].size() && q < hyperedge_adj_old[hyperedge_2].size()) {
                int a = hyperedge_adj_update[hyperedge_1][p];
                int b = hyperedge_adj_old[hyperedge_2][q];

                if (a == b ) {
                    int hyperedge_3 = a;
                    if(hyperedge_3>hyperedge_2){
                        vector<int> patternNum=get_specificPattern(hyperedge2node[hyperedge_1],hyperedge2node[hyperedge_2],hyperedge2node[hyperedge_3]);
                        if(patternNum[0]>0){
                            h_triangle[patternNum[0]-1]+=1;
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
    for(int i=0; i<hyperedge_new.size(); i++){
        int hyperedge_1=i;
        if(hyperedge_new[i].size()!=0) continue;


		for (int x = 0; x < (int)hyperedge_adj_update[i].size(); x++) {
            int hyperedge_2 = hyperedge_adj_update[i][x];




            int p = 0;  
            int q = 0;      
            
            while (p < hyperedge_adj_update[hyperedge_1].size() && q < hyperedge_adj_new[hyperedge_2].size()) {
                int a = hyperedge_adj_update[hyperedge_1][p];
                int b = hyperedge_adj_new[hyperedge_2][q];

                if (a == b ) {
                    int hyperedge_3 = a;
                    if(hyperedge_3>hyperedge_2){
                        vector<int> patternNum=get_specificPattern(hyperedge2node[hyperedge_1],hyperedge2node[hyperedge_2],hyperedge2node[hyperedge_3]);
                        if(patternNum[0]>0){
                            h_triangle[patternNum[0]-1]+=1;
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
    for(int i=0; i<hyperedge_new.size(); i++){
        int hyperedge_1=i;
        if(hyperedge_new[i].size()==0) continue;


		for (int x = 0; x < (int)hyperedge_adj_new[i].size(); x++) {
            int hyperedge_2 = hyperedge_adj_new[i][x];

            if(hyperedge_1>=hyperedge_2) continue;


            int p = 0;  
            int q = 0;      
            
            while (p < hyperedge_adj_new[hyperedge_1].size() && q < hyperedge_adj_new[hyperedge_2].size()) {
                int a = hyperedge_adj_new[hyperedge_1][p];
                int b = hyperedge_adj_new[hyperedge_2][q];

                if (a == b ) {
                    int hyperedge_3 = a;
                    if(hyperedge_3>hyperedge_2){
                        vector<int> patternNum=get_specificPattern(hyperedge2node[hyperedge_1],hyperedge2node[hyperedge_2],hyperedge2node[hyperedge_3]);
                        if(patternNum[0]>0){
                            h_triangle[patternNum[0]-1]+=1;
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

  
   

	cout<<"h_triangle:"<<endl;
	for(int i=0; i<h_triangle.size(); i++){
		cout<<h_triangle[i]<<" ";
	}
	cout<<endl;


    cout << "update time.\n"
    << (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
    cout << "------------------------------------------" << endl << endl;
	return 0;
}








