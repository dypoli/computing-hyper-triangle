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
#include <utility>
#include <bitset>
#include <chrono>

#include "read_data.cpp"
#include "motif_id.cpp"
#include <random>

#include <omp.h>
#include <atomic>
#include <thread>

using namespace std;




////////////////////////////////////////////////////      helper function   from line     //////////////////////////////////////////////////////////////

struct Compare {
    const std::vector<std::vector<int>>& hyperedge2node;

    Compare(const std::vector<std::vector<int>>& hyperedge2node)
        : hyperedge2node(hyperedge2node) {}

    bool operator()(int lhs, int rhs) const {
        if (hyperedge2node[lhs].size() == hyperedge2node[rhs].size()) {
            return lhs < rhs;
        }
        return hyperedge2node[lhs].size() < hyperedge2node[rhs].size();
    }
};




map<int,vector<vector<int>>> result;
vector<long long> h_triangle;
vector<vector<long long>> par_h_triangle;
vector< vector<int> > node2hyperedge;
vector< vector<int> > hyperedge2node;
vector< unordered_set<int> > hyperedge2node_set;
vector< unordered_map<int,int>> PairTypeInc; //<<id1, id2>,intersection>, <<id1, id2>, inclusion>.... note that we assume the pair.first < pair.second
vector< unordered_map<int,vector<int>> > PairTypeInt; //<<id1, id2>,intersection>, <<id1, id2>, inclusion>.... note that we assume the pair.first < pair.second
vector<vector<int>> HyperEdgeConnectionInc; //<<id2>.....> HyperEdge id1 and id2 are connect through inclusion type
vector<vector<int>> HyperEdgeConnectionInt; //<<id2>.....> HyperEdge id1 and id2 are connect through intersection type
vector< unordered_set<int> >  TotalSame;
vector< vector<int> > HyperEdgeConnectionRefinedInc;
vector< vector<int>> IntersectionList;
//vector< unordered_set<int>> IntersectionSet;
vector<int> sortedHyperEdge;
//vector< unordered_map<int,int>> PairCheck;
vector<unordered_map<int,vector<bool>>> PairCount;
vector<unordered_map<int,int>> PairId;
int TotalPair=0;
vector<bool> checkBefore; 

vector< unordered_map<int,vector<int>> > IntersecPair; //<<id1, id2>,intersection>, <<id1, id2>, inclusion>.... note that we assume the pair.first < pair.second



int binarySearch(const std::vector<int>& list, int id) {
    int left = 0;
    int right = list.size() - 1;

    while (left <= right) {
        int mid = left + (right - left) / 2;

        if (list[mid] == id) {
            return mid; // 找到元素，返回其索引
        } else if (list[mid] < id) {
            left = mid + 1; // 搜索右半部分
        } else {
            right = mid - 1; // 搜索左半部分
        }
    }

    return -1; // 未找到元素，返回-1
}


bool compareHyperEdge(const std::pair<int, std::vector<int>>& a, const std::pair<int, std::vector<int>>& b) {
    if (a.second.size() != b.second.size()) {
        return a.second.size() < b.second.size();
    }
    return a.first < b.first;
}

bool AltCompareHyperEdge(int a, int b) {
    if (hyperedge2node[a].size() == hyperedge2node[b].size()) {
        return a < b;
    }
    return hyperedge2node[a].size() < hyperedge2node[b].size();
}

void mergeVectors(const std::vector<std::vector<std::vector<int>>>& par_IntersectionList) {
    
    if (!par_IntersectionList.empty()) {
        size_t numInnerVectors = par_IntersectionList[0].size();
        
        // 遍历par_IntersectionList中的每个元素
        
            // 合并相同位置的vector<int>
		#pragma omp parallel for
		for (size_t i = 0; i < numInnerVectors; ++i) {
			for (const std::vector<std::vector<int>>& subList : par_IntersectionList) {
				IntersectionList[i].insert(IntersectionList[i].end(), subList[i].begin(), subList[i].end());
			}
        }
    }
    return;
}

void mergeMaps(const std::vector<std::vector<std::unordered_map<int, std::vector<int>>>>& par_IntersecPair) {


	if (!par_IntersecPair.empty()) {
        size_t numInnerMaps = par_IntersecPair[0].size();


        // 遍历par_IntersecPair中的每个元素
        

            // 合并相同位置的unordered_map<int, std::vector<int>>
		#pragma omp parallel for
		for (size_t i = 0; i < numInnerMaps; ++i) {
			for (const std::vector<std::unordered_map<int, std::vector<int>>>& subList : par_IntersecPair) {
				for (const auto& kvp : subList[i]) {
					int key = kvp.first;
					const std::vector<int>& values = kvp.second;
					IntersecPair[i][key].insert(IntersecPair[i][key].end(), values.begin(), values.end());
				}
			}
        }
    }

    return;
}

void mergeHyperEdgeConnections(const std::vector<std::vector<std::vector<int>>>& par_HyperEdgeConnectionInt) {


    if (!par_HyperEdgeConnectionInt.empty()) {
        size_t numInnerVectors = par_HyperEdgeConnectionInt[0].size();

        // 遍历 par_HyperEdgeConnectionInt 中的每个元素

            // 合并相同位置的 vector<int>
		#pragma omp parallel for
		for (size_t i = 0; i < numInnerVectors; ++i) {
			for (const std::vector<std::vector<int>>& subList : par_HyperEdgeConnectionInt) {
				HyperEdgeConnectionInt[i].insert(HyperEdgeConnectionInt[i].end(), subList[i].begin(), subList[i].end());
			}
        }
    }

    return;
}

void mergePairTypeInt(const std::vector<std::vector<std::unordered_map<int, std::vector<int>>>>& par_PairTypeInt) {


    if (!par_PairTypeInt.empty()) {
        size_t numInnerMaps = par_PairTypeInt[0].size();

        // 遍历 par_PairTypeInt 中的每个元素
        

            // 合并相同位置的 unordered_map<int, std::vector<int>>
		#pragma omp parallel for
		for (size_t i = 0; i < numInnerMaps; ++i) {
			for (const std::vector<std::unordered_map<int, std::vector<int>>>& subList : par_PairTypeInt) {
                for (const auto& kvp : subList[i]) {
                    int key = kvp.first;
                    const std::vector<int>& values = kvp.second;
                    PairTypeInt[i][key].insert(PairTypeInt[i][key].end(), values.begin(), values.end());
                }
            }
        }
    }

    return;
}

void mergePairTypeInc(const std::vector<std::vector<std::unordered_map<int, int>>>& par_PairTypeInc) {


    if (!par_PairTypeInc.empty()) {
        size_t numInnerMaps = par_PairTypeInc[0].size();
        

        // 遍历 par_PairTypeInc 中的每个元素
        

            // 合并相同位置的 unordered_map<int, int>
		#pragma omp parallel for
		for (size_t i = 0; i < numInnerMaps; ++i) {
			for (const std::vector<std::unordered_map<int, int>>& subList : par_PairTypeInc) {
                for (const auto& kvp : subList[i]) {
                    int key = kvp.first;
                    int value = kvp.second;
                    PairTypeInc[i][key] += value;  // 合并方式可以根据需求修改
                }
            }
        }
    }

    return;
}

void mergeHyperEdgeConnectionRefinedInc(const std::vector<std::vector<std::vector<int>>>& par_HyperEdgeConnectionRefinedInc) {

    if (!par_HyperEdgeConnectionRefinedInc.empty()) {
        size_t numInnerVectors = par_HyperEdgeConnectionRefinedInc[0].size();


        // 遍历 par_HyperEdgeConnectionRefinedInc 中的每个元素
        
            // 合并相同位置的 vector<int>
		#pragma omp parallel for
		for (size_t i = 0; i < numInnerVectors; ++i) {
			for (const std::vector<std::vector<int>>& subList : par_HyperEdgeConnectionRefinedInc) {
                HyperEdgeConnectionRefinedInc[i].insert(HyperEdgeConnectionRefinedInc[i].end(), subList[i].begin(), subList[i].end());
            }
        }
    }

    return;
}

void mergeHyperEdgeConnectionInc(const std::vector<std::vector<std::vector<int>>>& par_HyperEdgeConnectionInc) {


    if (!par_HyperEdgeConnectionInc.empty()) {
        size_t numInnerVectors = par_HyperEdgeConnectionInc[0].size();
        


        // 遍历 par_HyperEdgeConnectionInc 中的每个元素
        

            // 合并相同位置的 vector<int>
		#pragma omp parallel for
		for (size_t i = 0; i < numInnerVectors; ++i) {
			for (const std::vector<std::vector<int>>& subList : par_HyperEdgeConnectionInc) {
                HyperEdgeConnectionInc[i].insert(HyperEdgeConnectionInc[i].end(), subList[i].begin(), subList[i].end());
            }
        }
    }

    return;
}

std::vector<int> sortHyperedge2Node(const std::vector<std::vector<int>>& hyperedge2node) {
    std::vector<std::pair<int, std::vector<int>>> indexedHyperedge2Node;
    
    // 将每个vector<int>与其index存储到一个pair中
    for (int i = 0; i < hyperedge2node.size(); ++i) {
        indexedHyperedge2Node.emplace_back(i, hyperedge2node[i]);
    }
    
    // 使用compare函数进行排序
    std::sort(indexedHyperedge2Node.begin(), indexedHyperedge2Node.end(), compareHyperEdge);
    
    // 构建排序后的结果
    std::vector<int> result;
    for (const auto& pair : indexedHyperedge2Node) {
        result.push_back(pair.first);
    }
    
    return result;
}


vector<int>* findPair(vector< unordered_map<int,vector<int>> >& PairTypeInt, int id1, int id2){
	if(id1<id2){
		auto it = PairTypeInt[id1].find(id2);
		if(it!=PairTypeInt[id1].end()){
			return &it->second;
		}
		
	}else{
		auto it = PairTypeInt[id2].find(id1);
		if(it!=PairTypeInt[id2].end()){
			return &it->second;
		}
		
	}
	return NULL;
}

int findPairSize(vector< unordered_map<int,vector<int>> >& PairTypeInt, int id1, int id2){
	if(id1<id2){
		auto it = PairTypeInt[id1].find(id2);
		if(it!=PairTypeInt[id1].end()){
			return it->second.size();
		}
		
	}else{
		auto it = PairTypeInt[id2].find(id1);
		if(it!=PairTypeInt[id2].end()){
			return it->second.size();
		}
		
	}
	return 0;
}

vector<int>* findPair2(vector< unordered_map<int,int >>& PairTypeInc, vector< vector<int> >& hyperedge2node, int id1, int id2){
	if(id1<id2){
		auto it = PairTypeInc[id1].find(id2);
		if(it==PairTypeInc[id1].end()){
			return NULL;
		}
		int includId=it->second;
		return &hyperedge2node[includId];
	}else{
		auto it = PairTypeInc[id2].find(id1);
		if(it==PairTypeInc[id2].end()){
			return NULL;
		}
		int includId=it->second;
		return &hyperedge2node[includId];
	}
}

int findPairSize2(vector< unordered_map<int,int >>& PairTypeInc, vector< vector<int> >& hyperedge2node, int id1, int id2){
	if(id1<id2){
		auto it = PairTypeInc[id1].find(id2);
		if(it!=PairTypeInc[id1].end()){
			int includId=it->second;
			return hyperedge2node[includId].size();
		}
		
	}else{
		auto it = PairTypeInc[id2].find(id1);
		if(it!=PairTypeInc[id2].end()){
			int includId=it->second;
			return hyperedge2node[includId].size();
		}
		
	}
		return 0;
	
}

int DistinguishAAA(int part1, int part2, int part3, int part12, int part13, int part23, int part123){
	int count=0;
	if(part12!=0){
		count+=1;
	}
	if(part23!=0){
		count+=1;
	}
	if(part13!=0){
		count+=1;
	}

	if(part1!=0 && part2!=0 && part3!=0){//motif 2 6 12 16 26
		if(part123==0){
			return 26;
		}
		if(count==0){
			return 2;
		}else if(count==1){
			return 6;
		}else if(count==2){
			return 12;
		}else if(count==3){
			return 16;
		}
	}else if((part1==0 && part2!=0 && part3!=0) || (part1!=0 && part2==0 && part3!=0) ||(part1!=0 && part2!=0 && part3==0) ){//11 15 25
		if(part123==0){
			return 25;
		}
		if(count==2){
			return 11;
		}else if(count==3){
			return 15;
		}
	}else if(part1==0 && part2==0 && part3==0){//13 23
		if(part123==0){
			return 23;
		}else{
			return 13;
		}
	}else{// 14 24
		if(part123==0){
			return 24;
		}else{
			return 14;
		}
	}
	return -1;
}

bool compareByNum(const int& a, const int& b, vector< vector<int> >& hyperedge2node) {
    return hyperedge2node[a].size() < hyperedge2node[b].size();
}

bool CheckSame(vector< unordered_set<int> >& TotalSame, int id1, int id2){
	if(id1<id2){
		auto it = TotalSame[id1].find(id2);
		if(it==TotalSame[id1].end()){
			return true;
		}
		return false;
	}else{
		auto it = TotalSame[id2].find(id1);
		if(it==TotalSame[id2].end()){
			return true;
		}
		return false;
	}

}

vector<long long> TransferPattern(vector<long long> h_triangle_final){
	h_triangle_final[0]=h_triangle[2];
	h_triangle_final[1]=h_triangle[0];
	h_triangle_final[2]=h_triangle[3];
	h_triangle_final[3]=h_triangle[6];
	h_triangle_final[4]=h_triangle[7];
	h_triangle_final[5]=h_triangle[4];
	h_triangle_final[6]=h_triangle[8];
	h_triangle_final[7]=h_triangle[9];
	h_triangle_final[8]=h_triangle[1];
	h_triangle_final[9]=h_triangle[5];
	h_triangle_final[10]=h_triangle[10];
	h_triangle_final[11]=h_triangle[11];
	h_triangle_final[12]=h_triangle[12];

	h_triangle_final[13]=h_triangle[13];
	h_triangle_final[14]=h_triangle[14];
	h_triangle_final[15]=h_triangle[15];
	h_triangle_final[16]=h_triangle[22];
	h_triangle_final[17]=h_triangle[23];
	h_triangle_final[18]=h_triangle[24];
	h_triangle_final[19]=h_triangle[25];
	return h_triangle_final;
}

////////////////////////////////////////////////////   end of helper function     //////////////////////////////////////////////////////////////



void preprocess(int currentId, vector< unordered_set<int> >& TotalSame, vector< vector<int>>& IntersectionList, vector< unordered_map<int,vector<int>> >& IntersecPair,  vector< vector<int>>& HyperEdgeConnectionRefinedInc, vector< vector<int> >& HyperEdgeConnectionInc, vector< vector<int> >& HyperEdgeConnectionInt, vector< unordered_map<int,int> >& PairTypeInc, vector< unordered_map<int,vector<int>> >& PairTypeInt, vector<vector<int>>& Hyperedge2Node2Hyperedge, vector<int>& NodeIdList){
	while(true){
		int id=-1;
		bool check=false;
		for(size_t i=0; i<Hyperedge2Node2Hyperedge.size();i++){
			if(Hyperedge2Node2Hyperedge[i].size()>0){
				check=true;
				if(id==-1){
					id=Hyperedge2Node2Hyperedge[i][0];
				}
				if(id>Hyperedge2Node2Hyperedge[i][0]){
					id=Hyperedge2Node2Hyperedge[i][0];
				}
			}
		}
		if(!check){
			break;
		}
		int TtoalNumber=0;
		for(size_t i=0; i<Hyperedge2Node2Hyperedge.size();i++){
			if(Hyperedge2Node2Hyperedge[i].size()>0){
				if(id==Hyperedge2Node2Hyperedge[i][0]){
					//Hyperedge2Node2Hyperedge[i].erase(Hyperedge2Node2Hyperedge[i].begin());
					TtoalNumber+=1;
				}
			}
		}

		int MinimumId=id;
		int EdgeNodeNum1=hyperedge2node[currentId].size();
		int EdgeNodeNum2=hyperedge2node[MinimumId].size();
		int InclusionNum=TtoalNumber;
		
		if((EdgeNodeNum1>EdgeNodeNum2 && EdgeNodeNum2==InclusionNum) || (EdgeNodeNum1<EdgeNodeNum2 && EdgeNodeNum1==InclusionNum) ){//hyperedge2 is included in hyperedge1 or hyperedge1 is included in hyperedge2
			pair<int, int> pair2 = make_pair(currentId, MinimumId);
			//PairTypeInc.insert(make_pair(pair2,pair1.second));
			//StoreEdgeType(MinimumId, currentId, HyperEdgeType, "inclusion");
			for(size_t i=0; i<Hyperedge2Node2Hyperedge.size();i++){
				if(Hyperedge2Node2Hyperedge[i].size()>0){
					if(id==Hyperedge2Node2Hyperedge[i][0]){
						Hyperedge2Node2Hyperedge[i].erase(Hyperedge2Node2Hyperedge[i].begin());
					}
				}
			}	
			if(EdgeNodeNum1>EdgeNodeNum2){
				PairTypeInc[currentId].insert(make_pair(MinimumId,MinimumId));
			}else{
				PairTypeInc[currentId].insert(make_pair(MinimumId,currentId));
			}
			HyperEdgeConnectionInc[currentId].push_back(MinimumId);
			HyperEdgeConnectionInc[MinimumId].push_back(currentId);
			HyperEdgeConnectionRefinedInc[currentId].push_back(MinimumId);
		}else if(InclusionNum!=EdgeNodeNum2 && InclusionNum!=EdgeNodeNum1){//intersection
			PairTypeInt[currentId].insert(make_pair(MinimumId,vector<int>()));
			TotalPair+=1;

			if(AltCompareHyperEdge(currentId, MinimumId)){
				IntersectionList[currentId].push_back(MinimumId);
				IntersecPair[currentId].insert(make_pair(MinimumId,vector<int>()));

				IntersecPair[currentId][MinimumId].push_back(TotalPair-1);
			}else{
				IntersectionList[MinimumId].push_back(currentId);
				IntersecPair[MinimumId].insert(make_pair(currentId,vector<int>()));

				IntersecPair[MinimumId][currentId].push_back(TotalPair-1);
			}
			
			for(size_t i=0; i<Hyperedge2Node2Hyperedge.size();i++){
				if(Hyperedge2Node2Hyperedge[i].size()>0){
					if(id==Hyperedge2Node2Hyperedge[i][0]){
						Hyperedge2Node2Hyperedge[i].erase(Hyperedge2Node2Hyperedge[i].begin());
						PairTypeInt[currentId][MinimumId].push_back(NodeIdList[i]);
						if(AltCompareHyperEdge(currentId, MinimumId)){
							IntersecPair[currentId][MinimumId].push_back(NodeIdList[i]);
						}else{
							IntersecPair[MinimumId][currentId].push_back(NodeIdList[i]);
						}
					}
				}
			}


			
			//StoreEdgeType(MinimumId, currentId, HyperEdgeType, "intersection");
			HyperEdgeConnectionInt[currentId].push_back(MinimumId);
			HyperEdgeConnectionInt[MinimumId].push_back(currentId);
			
		}else if(InclusionNum==EdgeNodeNum2 && InclusionNum==EdgeNodeNum1){

			for(size_t i=0; i<Hyperedge2Node2Hyperedge.size();i++){
				if(Hyperedge2Node2Hyperedge[i].size()>0){
					if(id==Hyperedge2Node2Hyperedge[i][0]){
						Hyperedge2Node2Hyperedge[i].erase(Hyperedge2Node2Hyperedge[i].begin());
					}
				}
			}	
			TotalSame[currentId].insert(MinimumId);
		}else{
			for(size_t i=0; i<Hyperedge2Node2Hyperedge.size();i++){
				if(Hyperedge2Node2Hyperedge[i].size()>0){
					if(id==Hyperedge2Node2Hyperedge[i][0]){
						Hyperedge2Node2Hyperedge[i].erase(Hyperedge2Node2Hyperedge[i].begin());
					}
				}
			}	
		}
	}

	//store the type of hyperEdge currentId "A"/"B"/"AB"
}

int par_CountTTT(){

	#pragma omp parallel for
	for(int i=0; i<sortedHyperEdge.size();i++){
		int tid = omp_get_thread_num();
		int currentEdge=sortedHyperEdge[i];
		
		for(int j=0;j<IntersectionList[currentEdge].size();j++){
			int secondEdge=IntersectionList[currentEdge][j];
			vector<int>* CommonVer=&IntersecPair[currentEdge][secondEdge];
			int FirstId=(*CommonVer)[0];


			for(int k=0;k<IntersectionList[secondEdge].size();k++){
				int thirdEdge=IntersectionList[secondEdge][k];
				if(AltCompareHyperEdge(IntersectionList[currentEdge][IntersectionList[currentEdge].size()-1], thirdEdge)){
					break;
				}
				

				auto it2 = IntersecPair[currentEdge].find(thirdEdge);
				if(it2== IntersecPair[currentEdge].end() ){
					continue;
				}
			
				vector<int>* firstCommon=&it2->second;

				auto it3 = IntersecPair[secondEdge].find(thirdEdge);
				vector<int>* secondCommon=&it3->second;
				
				int TripleCommon=0;
				int ii = 1, jj = 1;
				int n1 = firstCommon->size(), n2 = secondCommon->size();
				while (ii < n1 && jj < n2 ) {
					if ((*firstCommon)[ii] == (*secondCommon)[jj] ) {
						TripleCommon+=1;
						ii++;
						jj++;
						
					}
					else if ((*firstCommon)[ii] < (*secondCommon)[jj]) {
						ii++;
					}
					else {
						jj++;
					}
					
				}

				
				int part123=TripleCommon;
				int part12= CommonVer->size()-1-part123;
				int part13=firstCommon->size()-1-part123;
				int part23=secondCommon->size()-1-part123;
				int part1=hyperedge2node[currentEdge].size()-part12-part13-part123;
				int part2=hyperedge2node[secondEdge].size()-part12-part23-part123;
				int part3=hyperedge2node[thirdEdge].size()-part23-part13-part123;
				int motifNum=DistinguishAAA(part1,part2, part3, part12, part13, part23, part123);
				par_h_triangle[tid][motifNum-1]+=1;

			}
		}
	}
	return 1;
}

int par_CountTTC(){

	#pragma omp parallel for
	for (int ii=0; ii<PairTypeInc.size();ii++) {
		int tid = omp_get_thread_num();
		for(const auto& pair : PairTypeInc[ii]){
			
			int firstId=ii;
			int secondId=pair.first;
			int includeId=pair.second;
			vector<int> CommonVer = hyperedge2node[includeId];
			if(HyperEdgeConnectionInt[firstId].size()==0 && HyperEdgeConnectionInt[secondId].size()==0 ){//this two hyperedges only contian in inclusion type
				continue;
			}

			int i = 0;  
			int j = 0;  
			
			while (i < HyperEdgeConnectionInt[firstId].size() || j < HyperEdgeConnectionInt[secondId].size()) {
                if(i == HyperEdgeConnectionInt[firstId].size()){
                    int thirdId=HyperEdgeConnectionInt[secondId][j];
                    ++j;
					
                    if(!CheckSame(TotalSame, firstId, thirdId)){
						continue;
					}
					if(thirdId==firstId){
						continue;
					}
                    int Common2=findPairSize2(PairTypeInc, hyperedge2node, thirdId, firstId);
                    if(Common2!=0){
                        continue;
                    }


					int firstCommonSize=findPairSize2(PairTypeInc,hyperedge2node, firstId, secondId);
					int secondCommonSize=findPairSize(PairTypeInt, thirdId, secondId);
					if(firstCommonSize+secondCommonSize==hyperedge2node[secondId].size()){
						//result[19-1].push_back(motifEdges);
						par_h_triangle[tid][19-1]+=1;
					}else{
						//result[20-1].push_back(motifEdges);
						par_h_triangle[tid][20-1]+=1;
					}

                }else if(j == HyperEdgeConnectionInt[secondId].size()){
                    int thirdId=HyperEdgeConnectionInt[firstId][i];
                    ++i;
                    if(thirdId==secondId){
                        continue;
                    }
					
                    if(!CheckSame(TotalSame, secondId, thirdId)){
                        continue;
                    }

                    int Common2=findPairSize2(PairTypeInc, hyperedge2node, thirdId, secondId);
                    if(Common2!=0){
                        continue;
                    }
   
                    int firstCommonSize=findPairSize2(PairTypeInc,hyperedge2node, firstId, secondId);
                    int secondCommonSize=findPairSize(PairTypeInt, thirdId, firstId);
                    if(firstCommonSize+secondCommonSize==hyperedge2node[firstId].size()){
                        //result[19-1].push_back(motifEdges);
                        par_h_triangle[tid][19-1]+=1;
                    }else{
                        //result[20-1].push_back(motifEdges);
                        par_h_triangle[tid][20-1]+=1;
                    }

                }else{
                    if (HyperEdgeConnectionInt[firstId][i] < HyperEdgeConnectionInt[secondId][j]) {
                        int thirdId=HyperEdgeConnectionInt[firstId][i];
                        ++i;
                        if(thirdId==secondId){
							continue;
						}
						
						if(!CheckSame(TotalSame, secondId, thirdId)){
							continue;
						}
                        int Common2=findPairSize2(PairTypeInc, hyperedge2node, thirdId, secondId);
                        if(Common2!=0){
                            continue;
                        }
   
                        int firstCommonSize=findPairSize2(PairTypeInc,hyperedge2node, firstId, secondId);
                        int secondCommonSize=findPairSize(PairTypeInt, thirdId, firstId);
                        if(firstCommonSize+secondCommonSize==hyperedge2node[firstId].size()){
                            //result[19-1].push_back(motifEdges);
                            par_h_triangle[tid][19-1]+=1;
                        }else{
                            //result[20-1].push_back(motifEdges);
                            par_h_triangle[tid][20-1]+=1;
                        }

                    } else if (HyperEdgeConnectionInt[firstId][i] > HyperEdgeConnectionInt[secondId][j]) {
                        int thirdId=HyperEdgeConnectionInt[secondId][j];
                        ++j;
						
                        if(!CheckSame(TotalSame, firstId, thirdId)){
                            continue;
                        }
                        if(thirdId==firstId){
                            continue;
                        }

                        int Common2=findPairSize2(PairTypeInc, hyperedge2node, thirdId, firstId);
                        if(Common2!=0){
                            continue;
                        }

                        int firstCommonSize=findPairSize2(PairTypeInc,hyperedge2node, firstId, secondId);
                        int secondCommonSize=findPairSize(PairTypeInt, thirdId, secondId);
                        if(firstCommonSize+secondCommonSize==hyperedge2node[secondId].size()){
                            //result[19-1].push_back(motifEdges);
                            par_h_triangle[tid][19-1]+=1;
                        }else{
                            //result[20-1].push_back(motifEdges);
                            par_h_triangle[tid][20-1]+=1;
                        }

                    } else {
                        int commonHyperedgeId=HyperEdgeConnectionInt[firstId][i];
                        vector<int>* firstCommon=findPair(PairTypeInt, commonHyperedgeId, firstId);
                        vector<int>* secondCommon=findPair(PairTypeInt, commonHyperedgeId, secondId);
                        
                        if(firstCommon->size()!=secondCommon->size()){//belong to motif type 9, 10
                            
                            int ABC;
                            int BC;
                            int tempId;
                            if(firstCommon->size()<secondCommon->size()){
                                ABC=firstCommon->size();
                                tempId=secondId;
                                BC= secondCommon->size();
                            }else{
                                ABC=secondCommon->size();
                                tempId=firstId;
                                BC=firstCommon->size();
                            }
                            if(BC+CommonVer.size()-ABC==hyperedge2node[tempId].size()){//belong to motif type 9
                                //result[8].push_back(motifEdges);
                                par_h_triangle[tid][8]+=1;
                            }else{//belong to motif type 10
                                //result[9].push_back(motifEdges);
                                par_h_triangle[tid][9]+=1;
                            }
                            
                            
                        }else {
                            //result[4].push_back(motifEdges);
                            par_h_triangle[tid][4]+=1;
                        }
                        ++i;
                        ++j;
                    }
                }



				
			}
		}
    }
	return 0;
}

int par_CountCCC(){
	#pragma omp parallel for
	for (int ii=0; ii<PairTypeInc.size();ii++) {
		int tid = omp_get_thread_num();
		for(auto& pair : PairTypeInc[ii]){
			int firstId=ii;
			int secondId=pair.first;
			for(int i=0; i<HyperEdgeConnectionRefinedInc[secondId].size(); i++){
				int commonHyperedgeId=HyperEdgeConnectionRefinedInc[secondId][i];
				if(commonHyperedgeId>HyperEdgeConnectionRefinedInc[firstId][HyperEdgeConnectionRefinedInc[firstId].size()-1]){
					break;
				}
				auto it = PairTypeInc[firstId].find(commonHyperedgeId);
				
				if(it!=PairTypeInc[firstId].end()){//firstId are connect to commonHyperedgeId
					
					//vector<int> motifEdges;
					//motifEdges.push_back(firstId);
					//motifEdges.push_back(secondId);
					//motifEdges.push_back(commonHyperedgeId);
					//result[motifNum-1].push_back(motifEdges);
					par_h_triangle[tid][2]+=1;
				}
			}
		}
    }
	return 0;
}

int par_CountTCC(){
	#pragma omp parallel for
	for (int ii=0; ii<PairTypeInt.size();ii++) {
		int tid = omp_get_thread_num();
		for(const auto& pair : PairTypeInt[ii]){
			int firstId=ii;
			int secondId=pair.first;
			vector<int> CommonVer = pair.second;
			if(HyperEdgeConnectionInc[firstId].size()==0 || HyperEdgeConnectionInc[secondId].size()==0 ){//this two hyperedges only contian in inclusion type
				continue;
			}

			if(HyperEdgeConnectionInc[firstId][HyperEdgeConnectionInc[firstId].size()-1]<HyperEdgeConnectionInc[secondId][0] || HyperEdgeConnectionInc[firstId][0]>HyperEdgeConnectionInc[secondId][HyperEdgeConnectionInc[secondId].size()-1]){
				continue;
			}
			
			int i = 0;  
			int j = 0;  

			while (i < HyperEdgeConnectionInc[firstId].size() && j < HyperEdgeConnectionInc[secondId].size()) {
				if (HyperEdgeConnectionInc[firstId][i] < HyperEdgeConnectionInc[secondId][j]) {
					++i;
				} else if (HyperEdgeConnectionInc[firstId][i] > HyperEdgeConnectionInc[secondId][j]) {
					++j;
				} else {
					int commonHyperedgeId=HyperEdgeConnectionInc[firstId][i];
					vector<int>* firstCommon=findPair2(PairTypeInc, hyperedge2node, commonHyperedgeId, firstId);
					vector<int>* secondCommon=findPair2(PairTypeInc, hyperedge2node,commonHyperedgeId, secondId);
					
					if(hyperedge2node[firstId].size()>hyperedge2node[commonHyperedgeId].size() && hyperedge2node[secondId].size()>hyperedge2node[commonHyperedgeId].size()){//h-motif 1 or 4
						if(CommonVer.size()==firstCommon->size()){
							//result[0].push_back(motifEdges);
							par_h_triangle[tid][0]+=1;
						}else{
							//result[3].push_back(motifEdges);
							par_h_triangle[tid][3]+=1;
						}
					}else{//h-motif 7 or 8
						int countNum=firstCommon->size()+secondCommon->size()-CommonVer.size();
						if(countNum==hyperedge2node[commonHyperedgeId].size()){
							//result[6].push_back(motifEdges);
							par_h_triangle[tid][6]+=1;
						}else{
							//result[7].push_back(motifEdges);
							par_h_triangle[tid][7]+=1;
						}
					}
					
					++i;
					++j;
				}
			}
		}
    }
	return 0;
}

int par_ExactOptimal(){
	
    par_CountTTT();
	
    par_CountTTC();

    par_CountTCC();

    par_CountCCC();
    return 1;
}



inline long long convert_id(int hyperedge_a, int hyperedge_b){
	return hyperedge_a * (1LL << 31) + hyperedge_b;
}

int main(int argc, char *argv[])
{
	chrono::system_clock::time_point start;
	chrono::system_clock::time_point run_start;
	chrono::duration<double> dur;
	int progress;

	//for test
    //int num_threads = stoi(argv[1]);
	int num_threads =  stoi(argv[1]);
	omp_set_num_threads(num_threads);

	for(int i=0; i<num_threads; i++){
		par_h_triangle.push_back(vector<long long>());
		par_h_triangle[i].resize(26,0);
	}


	for(int i=0; i<26; i++){
		vector<vector<int>> temp;
		result.insert(make_pair(i,temp));
		h_triangle.push_back(0);
	}

    

	string dataName="dblp_graph.txt";
	string graphFile = "Dataset/"+dataName;
	// Read data
	start = chrono::system_clock::now();
	
	read_data(graphFile, node2hyperedge, hyperedge2node, hyperedge2node_set);
	int V = node2hyperedge.size(), E = hyperedge2node.size();
	cout << "# of nodes: " << V << '\n';
	cout << "# of hyperedges: " << E << '\n';
    dur = std::chrono::system_clock::now() - start;
	cout << "Reading data done: "
		<< dur.count() << " sec" << endl;
	cout << "------------------------------------------" << endl << endl;

	//exact count

	

	
	
	run_start = chrono::system_clock::now();

	vector<vector<int>> copyVec=node2hyperedge;
	
	for(int i=0; i<hyperedge2node.size();i++){
		PairTypeInc.push_back(unordered_map<int,int>());
		PairTypeInt.push_back(unordered_map<int,vector<int>> ());
		HyperEdgeConnectionInc.push_back(vector<int>());
		HyperEdgeConnectionInt.push_back(vector<int>());
		TotalSame.push_back(unordered_set<int> ());
		HyperEdgeConnectionRefinedInc.push_back(vector<int>());
        IntersectionList.push_back(vector<int>());
		PairCount.push_back(unordered_map<int,vector<bool>>());
		IntersecPair.push_back(unordered_map<int,vector<int>> ());
		PairId.push_back(unordered_map<int,int>());
	}
	
	//////////////////////////////////////////////     only for parallel algorithm only //////////////////////////////////////////////////////
	vector<vector< unordered_map<int,int>>> par_PairTypeInc=vector<vector< unordered_map<int,int>>>(num_threads); 
	vector<vector< unordered_map<int,vector<int>>>> par_PairTypeInt=vector<vector< unordered_map<int,vector<int>>>>(num_threads); 
	vector<vector<vector<int>>> par_HyperEdgeConnectionInc=vector<vector<vector<int>>>(num_threads); 
	vector<vector<vector<int>>> par_HyperEdgeConnectionInt=vector<vector<vector<int>>>(num_threads); 
	vector<vector< vector<int>>> par_HyperEdgeConnectionRefinedInc=vector<vector< vector<int>>>(num_threads);
	vector<vector< vector<int>>> par_IntersectionList=vector<vector< vector<int>>>(num_threads);
	vector<vector< unordered_map<int,vector<int>> >> par_IntersecPair=vector<vector< unordered_map<int,vector<int>> >>(num_threads);
	
	for(int i=0; i<num_threads; i++){
		for(int j=0; j<hyperedge2node.size();j++){
			par_PairTypeInc[i].push_back(unordered_map<int,int>());
			par_PairTypeInt[i].push_back(unordered_map<int,vector<int>> ());
			par_HyperEdgeConnectionInc[i].push_back(vector<int>());
			par_HyperEdgeConnectionInt[i].push_back(vector<int>());
			par_HyperEdgeConnectionRefinedInc[i].push_back(vector<int>());
			par_IntersectionList[i].push_back(vector<int>());
			par_IntersecPair[i].push_back(unordered_map<int,vector<int>> ());

		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    start = chrono::system_clock::now();
	

	
    #pragma omp parallel for
	for(int id=0; id<hyperedge2node.size();id++){
		int tid = omp_get_thread_num();
		vector<int> connectNode=hyperedge2node[id];
		vector<vector<int>> Hyperedge2Node2Hyperedge;
		vector<int> NodeIdList;
		//first delete hyperedge in connectNode
		for(size_t i=0; i<connectNode.size();i++){
			int NodeId=connectNode[i];
			int NodeIndex=binarySearch(copyVec[NodeId], id);
			if(copyVec[NodeId].size()-1<=NodeIndex){//this node should be removed			
				//copyVec[NodeId][0]=-1;
			}else{
				//copyVec[NodeId].erase(copyVec[NodeId].begin());
				Hyperedge2Node2Hyperedge.push_back(vector<int>());
				int curSize=Hyperedge2Node2Hyperedge.size();
				Hyperedge2Node2Hyperedge[curSize-1].insert(Hyperedge2Node2Hyperedge[curSize-1].end(), copyVec[NodeId].begin()+NodeIndex+1,  copyVec[NodeId].end());	
				NodeIdList.push_back(NodeId);		
			}			
		}
	
		preprocess(id,TotalSame, par_IntersectionList[tid], par_IntersecPair[tid], par_HyperEdgeConnectionRefinedInc[tid], par_HyperEdgeConnectionInc[tid], par_HyperEdgeConnectionInt[tid], par_PairTypeInc[tid], par_PairTypeInt[tid], Hyperedge2Node2Hyperedge, NodeIdList);
	}



	//////////////////////////////////////////         combine vector         ///////////////////////////////////////////////////////////

	mergeVectors(par_IntersectionList);
	mergeMaps(par_IntersecPair);
	mergeHyperEdgeConnections(par_HyperEdgeConnectionInt);
	mergePairTypeInc(par_PairTypeInc);
	mergePairTypeInt(par_PairTypeInt);
	mergeHyperEdgeConnectionRefinedInc(par_HyperEdgeConnectionRefinedInc);
	mergeHyperEdgeConnectionInc(par_HyperEdgeConnectionInc);
	
 	#pragma omp parallel for
	for (int i=0; i<HyperEdgeConnectionInc.size(); i++) {
        sort(HyperEdgeConnectionInc[i].begin(), HyperEdgeConnectionInc[i].end());
    }
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	
	dur = std::chrono::system_clock::now() - start;
	double time1=dur.count();
	cout << "Adjacency list construction done: "
		<< (double)(time1)<< " sec" << endl;
	cout << "------------------------------------------" << endl << endl;

	chrono::system_clock::time_point start2 = std::chrono::system_clock::now();

	////////////////////////////////////     build     hyperEdge    network      /////////////////////////
	sortedHyperEdge=sortHyperedge2Node(hyperedge2node);

	Compare comp(hyperedge2node);
    for (auto& vec : IntersectionList) {
        std::sort(vec.begin(), vec.end(), comp);
    }

	
	

	dur = std::chrono::system_clock::now() - start2;
	double time2=dur.count();
	
    cout << "build tree done: "
		<< (double)(time2) << " sec" << endl;
	cout<<'\n';

	///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// /////////////////////////
	
	start = std::chrono::system_clock::now();

 	par_ExactOptimal();
	////////////////////////////////////////////////////  add all of them into one vector<long long>   /////////////////////////////////////////////////////
	for(int i=0; i<num_threads; i++){
		for(int j=0; j<26; j++){
			h_triangle[j]+=par_h_triangle[i][j];
		}
	}
	//////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// /////////////////////////

    dur = std::chrono::system_clock::now() - run_start;
	double runtime = dur.count();

    dur = std::chrono::system_clock::now() - start;
	double time3=dur.count();
	vector<long long> h_triangle_final(20, 0);
	h_triangle_final=TransferPattern(h_triangle_final);
	long long triangle=0;
	int index = 0;
	for (int i = 1; i <= 20; i++){
		cout << fixed << "h-triangle " << ++index << ": " << fixed << h_triangle_final[i-1] << endl;
		
		triangle+=h_triangle_final[i-1];
		
	}
	cout << "counting motif data done: "
		<< (double)(time3) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "Total runtime: " << runtime << endl;
	cout<<"triangle num:  "<<triangle<<"\n";

	return 0;
}
