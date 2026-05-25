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

    
    string dataName="unique-email-Enron.txt";
	string graphFile = dataName;


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

	vector< vector<int> > hyperedge2node_new(E);
	string updateVertex = "NewVertices.txt";
	loadAffectedData(updateVertex, hyperedge2node_new,node2hyperedge);

	//string TrianglNum= "C:\\AllCode\\Triangle Counting\\DynamicCode\\ResultData\\TriangleNum_"+dataName;
	//std::vector<std::vector<int>> TriangleNum_record = ReadFileToVector(TrianglNum);


	string TrianglNum2= "TriangleNumber.txt";
	vector<int> h_triangle=ReadTriangleFile(TrianglNum2); 



	set<int> Updated_hyperedge_set;
	for(int i=0; i<E; i++){
		if(hyperedge2node_new[i].size()!=0) Updated_hyperedge_set.insert(i);
	}

	// An adjacency list for hyperedge i in which all hyperedges in hyperedge[i] have neither added nor removed any vertices.
	vector< vector<int> > hyperedge_adj_old(E); 
	
	// An adjacency list for hyperedge i in which all hyperedges in hyperedge[i] have either added nor removed any vertices.
	vector< vector<int> > hyperedge_adj_adjust(E);
	
	// new or removed connections after adding or deleting vertices
	vector< vector<int> > hyperedge_adj_new(E); 
	
	for(int i=0; i<E; i++){
		for(int j=0; j<hyperedge_adj[i].size(); j++){
			if(Updated_hyperedge_set.find(hyperedge_adj[i][j])!=Updated_hyperedge_set.end()){
				hyperedge_adj_adjust[i].push_back(hyperedge_adj[i][j]);
			}else{
				hyperedge_adj_old[i].push_back(hyperedge_adj[i][j]);
			}
		}
	}
	
	hyperedge_adj_new= build_hyperedge_adjacency(hyperedge2node_new,node2hyperedge);
    removeDuplicatesFromHyperedgeAdj(hyperedge_adj_new);
    for(int i = 0; i < E; i++){
        sort(hyperedge_adj_new[i].begin(), hyperedge_adj_new[i].end());
    }

    

	remove_existing_adjacency(hyperedge_adj_new,hyperedge_adj);
	vector<vector<int>> hyperedge2node_update(E);

	for(int i=0; i<hyperedge2node_new.size(); i++){
		if(hyperedge2node_new[i].size()==0){
			hyperedge2node_update[i]=hyperedge2node[i];
		}else{
			hyperedge2node_update[i]=concat_vectors(hyperedge2node[i],hyperedge2node_new[i]); 
		}
	}


start = clock();

///////////////////////////////////////  update hyper-triangle ///////////////////////////////////////




	for(int i=0; i<hyperedge2node_new.size(); i++){
		int hyperedge_1=i;
		if(hyperedge2node_new[i].size()==0) continue;
		//std::cout << "i:" << i << std::endl;



        for (int x = 0; x < (int)hyperedge_adj_old[i].size(); x++) {
            int hyperedge_2 = hyperedge_adj_old[i][x];
			


            // ---------- 双指针求交集 ----------
            int p = x + 1;  // i 的邻居中，j 之后的部分
            int q = 0;      // j 的邻居开头

            while (p < hyperedge_adj_old[i].size() && q < hyperedge_adj_old[hyperedge_2].size()) {
                int a = hyperedge_adj_old[i][p];
                int b = hyperedge_adj_old[hyperedge_2][q];

                if (a == b) {
                    int hyperedge_3 = a;
                    if (hyperedge_3 > hyperedge_2) {   
						vector<int> Previous_patternNum=get_specificPattern(hyperedge2node[hyperedge_1],hyperedge2node[hyperedge_2],hyperedge2node[hyperedge_3]);

						if (Previous_patternNum[0]>0){
							h_triangle[Previous_patternNum[0]-1]-=1;
						}

                        

                        

						vector<int> patternNum=get_specificPattern(hyperedge2node_update[hyperedge_1],hyperedge2node_update[hyperedge_2],hyperedge2node_update[hyperedge_3]);
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


    

    for(int i=0; i<hyperedge2node_new.size(); i++){
            int hyperedge_1=i;
            if(hyperedge2node_new[i].size()==0) continue;
            //std::cout << "i:" << i << std::endl;
		
        for (int x = 0; x < (int)hyperedge_adj_adjust[i].size(); x++) {
            int hyperedge_2 = hyperedge_adj_adjust[i][x];
			if(hyperedge_1>=hyperedge_2) continue;


            // ---------- 双指针求交集 ----------
            int p = 0;  // i 的邻居中，j 之后的部分
            int q = 0;      // j 的邻居开头

		

            while (p < hyperedge_adj_old[hyperedge_1].size() && q < hyperedge_adj_old[hyperedge_2].size()) {
                int a = hyperedge_adj_old[hyperedge_1][p];
                int b = hyperedge_adj_old[hyperedge_2][q];

                if (a == b) {
                    int hyperedge_3 = a;
                    
					vector<int> Previous_patternNum=get_specificPattern(hyperedge2node[hyperedge_1],hyperedge2node[hyperedge_2],hyperedge2node[hyperedge_3]);

					if(Previous_patternNum[0]>0){
						h_triangle[Previous_patternNum[0]-1]-=1;
					}
					


					//vector<int> hyperedge_1_node=concat_vectors(hyperedge2node[hyperedge_1],hyperedge2node_new[hyperedge_1]); // 保证 i < j < k，避免重复
					vector<int> patternNum=get_specificPattern(hyperedge2node_update[hyperedge_1],hyperedge2node_update[hyperedge_2],hyperedge2node_update[hyperedge_3]);
					if(patternNum[0]>0){
						h_triangle[patternNum[0]-1]+=1;
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

    
    
    


    for(int i=0; i<hyperedge2node_new.size(); i++){
            int hyperedge_1=i;
            if(hyperedge2node_new[i].size()==0) continue;

		for (int x = 0; x < (int)hyperedge_adj_adjust[i].size(); x++) {
            int hyperedge_2 = hyperedge_adj_adjust[i][x];
			if(hyperedge_1>=hyperedge_2) continue;


            // ---------- 双指针求交集 ----------
            int p = x+1;  // i 的邻居中，j 之后的部分
            int q = 0;      // j 的邻居开头
		

            while (p < hyperedge_adj_adjust[hyperedge_1].size() && q < hyperedge_adj_adjust[hyperedge_2].size()) {
                int a = hyperedge_adj_adjust[hyperedge_1][p];
                int b = hyperedge_adj_adjust[hyperedge_2][q];

                if (a == b) {
                    int hyperedge_3 = a;
                    if(hyperedge_3>hyperedge_2){

						vector<int> Previous_patternNum=get_specificPattern(hyperedge2node[hyperedge_1],hyperedge2node[hyperedge_2],hyperedge2node[hyperedge_3]);



						if(Previous_patternNum[0]>0){
							h_triangle[Previous_patternNum[0]-1]-=1;
						}



						//vector<int> hyperedge_1_node=concat_vectors(hyperedge2node[hyperedge_1],hyperedge2node_new[hyperedge_1]); // 保证 i < j < k，避免重复
						vector<int> patternNum=get_specificPattern(hyperedge2node_update[hyperedge_1],hyperedge2node_update[hyperedge_2],hyperedge2node_update[hyperedge_3]);
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



	///////////////////////// new hyper-triangle ///////////////////////////////////




    for(int i=0; i<hyperedge2node_new.size(); i++){
        int hyperedge_1=i;
        if(hyperedge2node_new[i].size()==0) continue;


		for (int x = 0; x < (int)hyperedge_adj_new[i].size(); x++) {
            int hyperedge_2 = hyperedge_adj_new[i][x];
			if(hyperedge_1>=hyperedge_2) continue;



            // ---------- 双指针求交集 ----------
            int p = 0;  // i 的邻居中，j 之后的部分
            int q = 0;      // j 的邻居开头

            while (p < hyperedge_adj[hyperedge_1].size() && q < hyperedge_adj[hyperedge_2].size()) {
                int a = hyperedge_adj[hyperedge_1][p];
                int b = hyperedge_adj[hyperedge_2][q];

                if (a == b) {
                    int hyperedge_3 = a;
                    
					//vector<int> hyperedge_1_node=concat_vectors(hyperedge2node[hyperedge_1],hyperedge2node_new[hyperedge_1]); // 保证 i < j < k，避免重复
					vector<int> patternNum=get_specificPattern(hyperedge2node_update[hyperedge_1],hyperedge2node_update[hyperedge_2],hyperedge2node_update[hyperedge_3]);
					if(patternNum[0]>0){
						h_triangle[patternNum[0]-1]+=1;
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



    for(int i=0; i<E; i++){
        int hyperedge_1=i;
        if(hyperedge2node_new[i].size()==0) continue;

		for (int x = 0; x < (int)hyperedge_adj_new[i].size(); x++) {
            int hyperedge_2 = hyperedge_adj_new[i][x];
			if(hyperedge_1>=hyperedge_2) continue;

            
            // ---------- 双指针求交集 ----------
            int p = 0;  // i 的邻居中，j 之后的部分
            int q = 0;      // j 的邻居开头

            while (p < hyperedge_adj_new[hyperedge_1].size() && q < hyperedge_adj_new[hyperedge_2].size()) {
                int a = hyperedge_adj_new[hyperedge_1][p];
                int b = hyperedge_adj_new[hyperedge_2][q];



                if (a == b) {
                    int hyperedge_3 = a;
                    if(hyperedge_3>hyperedge_2){
                        
						vector<int> patternNum=get_specificPattern(hyperedge2node_update[hyperedge_1],hyperedge2node_update[hyperedge_2],hyperedge2node_update[hyperedge_3]);
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
    
    

    for(int i=0; i<E; i++){
        int hyperedge_1=i;
        if(hyperedge2node_new[i].size()==0) continue; 

		for (int x = 0; x < hyperedge_adj[i].size(); x++) {
            int hyperedge_2 = hyperedge_adj[i][x];
            if(hyperedge_1>=hyperedge_2) continue;


            // ---------- 双指针求交集 ----------
            int p = 0;  // i 的邻居中，j 之后的部分
            int q = 0;      // j 的邻居开头

            while (p < hyperedge_adj_new[hyperedge_1].size() && q < hyperedge_adj_new[hyperedge_2].size()) {
                int a = hyperedge_adj_new[hyperedge_1][p];
                int b = hyperedge_adj_new[hyperedge_2][q];

                if (a == b) {
                    int hyperedge_3 = a;
                    
                    vector<int> patternNum=get_specificPattern(hyperedge2node_update[hyperedge_1],hyperedge2node_update[hyperedge_2],hyperedge2node_update[hyperedge_3]);
                    if(patternNum[0]>0){
                        h_triangle[patternNum[0]-1]+=1;
                        
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







