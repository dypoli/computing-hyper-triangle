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
#include "read_data.cpp"
#include "motif_id.cpp"
#include <random>

using namespace std;

vector< unordered_map<int,vector<int>> > IntersecPair; //<<id1, id2>,intersection>, <<id1, id2>, inclusion>.... note that we assume the pair.first < pair.second
std::vector<std::vector<int>> OneLevelInclusionList;
std::vector<std::vector<int>> OneLevelIntersectionList;
int TotalPair=0;
vector< vector<int> > hyperedge2node;
vector< vector<int>> IntersectionList;
vector<vector<int>> Vertex2Pair;
map<int,vector<vector<int>>> result;
vector<long long> h_triangle;
vector<int> sortedHyperEdge;
vector<pair<int,int>> PairId2Edge;
int maxNode=0;
//vector< unordered_set<int>> IntersectionSet;
vector<int> signalList;
vector<int> VertexSignal;
vector<int> deleteList;
int IdSize=0;

/////////////////////////////           use only one time        /////////////////////////
vector< vector<int>> EnclosureList;
vector< vector<int> > node2hyperedge;
vector< unordered_set<int> > hyperedge2node_set;
vector< vector<int>> InclusionList;

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



void AddMotif(int currentNode, int secondNode, int thirdNode, int type){
	vector<int> temp;
	temp.push_back(currentNode);
	temp.push_back(secondNode);
	temp.push_back(thirdNode);
	sort(temp.begin(),temp.end());
	vector<int> motifEdges;
	motifEdges.push_back(temp[0]);
	motifEdges.push_back(temp[1]);
	motifEdges.push_back(temp[2]);
	result[type-1].push_back(motifEdges);
		
			
}

bool compareHyperEdge(const std::pair<int, std::vector<int>>& a, const std::pair<int, std::vector<int>>& b) {
    if (a.second.size() != b.second.size()) {
        return a.second.size() < b.second.size();
    }
    return a.first < b.first;
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

vector<int>* createIntVector(int num_elements) {

    size_t element_size = sizeof(int);
    size_t total_bytes = num_elements * element_size;

    int* intArray = static_cast<int*>(calloc(num_elements, element_size));

    if (intArray != nullptr) {
        // 将分配的内存用于创建vector<int>对象
        std::vector<int>* intVector = new std::vector<int>(intArray, intArray + num_elements);

        // 返回指向vector<int>对象的指针
        return intVector;
    } else {
        // 处理内存分配失败的情况
        std::cerr << "Memory allocation failed!" << std::endl;
        return nullptr;
    }
}

bool AltCompareHyperEdge(int a, int b) {
    if (hyperedge2node[a].size() == hyperedge2node[b].size()) {
        return a < b;
    }
    return hyperedge2node[a].size() < hyperedge2node[b].size();
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


void preprocess(int currentId, vector< unordered_set<int> >& TotalSame, vector<vector<int>>& Hyperedge2Node2Hyperedge, vector<int>& NodeIdList){
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
			for(size_t i=0; i<Hyperedge2Node2Hyperedge.size();i++){
				if(Hyperedge2Node2Hyperedge[i].size()>0){
					if(id==Hyperedge2Node2Hyperedge[i][0]){
						Hyperedge2Node2Hyperedge[i].erase(Hyperedge2Node2Hyperedge[i].begin());
					}
				}
			}	
			
			if(EdgeNodeNum1>EdgeNodeNum2){
				InclusionList[MinimumId].push_back(currentId);
                EnclosureList[currentId].push_back(MinimumId);
		
			}else{
				InclusionList[currentId].push_back(MinimumId);
                EnclosureList[MinimumId].push_back(currentId);

			}
		
		}else if(InclusionNum!=EdgeNodeNum2 && InclusionNum!=EdgeNodeNum1){//intersection
			TotalPair+=1;
			
			if(AltCompareHyperEdge(currentId, MinimumId)){
				IntersectionList[currentId].push_back(MinimumId);
				IntersecPair[currentId].insert(make_pair(MinimumId,vector<int>()));
			}else{
				IntersectionList[MinimumId].push_back(currentId);
				IntersecPair[MinimumId].insert(make_pair(currentId,vector<int>()));
			}
			
			for(size_t i=0; i<Hyperedge2Node2Hyperedge.size();i++){
				if(Hyperedge2Node2Hyperedge[i].size()>0){
					if(id==Hyperedge2Node2Hyperedge[i][0]){
						Hyperedge2Node2Hyperedge[i].erase(Hyperedge2Node2Hyperedge[i].begin());
						
						if(AltCompareHyperEdge(currentId, MinimumId)){
							IntersecPair[currentId][MinimumId].push_back(NodeIdList[i]);
						}else{
							IntersecPair[MinimumId][currentId].push_back(NodeIdList[i]);
						}
					}
				}
			}
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

}

void CountDenseTTT(){
    vector<int> pointer=vector<int>(maxNode);
    pair<int, int> nextPair;
	int lastPair=-1;
	vector<int> IdList=vector<int>(IdSize,0);
    for(int i=0; i<sortedHyperEdge.size();i++){
		int currentEdge=sortedHyperEdge[i];
		for(int j=0;j<IntersectionList[currentEdge].size();j++){
			int simi=0;
			//vector<int> IdList=vector<int>();
			int secondEdge=IntersectionList[currentEdge][j];
			//IdList.push_back(secondEdge);
			vector<int>* CommonVer=&IntersecPair[currentEdge][secondEdge];
			
			if(signalList[secondEdge]==currentEdge){
				for(int k=1; k<(*CommonVer).size();k++){
					int curVer=(*CommonVer)[k];
					VertexSignal[curVer]+=1;
				}
				continue;
			}
			vector<int>* secondCommon;
			int FirstId=(*CommonVer)[0];
            int round=0;
            while(true){
                nextPair.first=-1;
                for(int k=1; k<(*CommonVer).size();k++){
					
                    int curIndex=k-1;
                    int curVer=(*CommonVer)[k];
			
                    if(round==0 && nextPair.first==-1){
                        //Vertex2Pair[curVer].erase(Vertex2Pair[curVer].begin());
	
						VertexSignal[curVer]+=1;
                        int Tindex=VertexSignal[curVer];
                        while(Tindex<Vertex2Pair[curVer].size()){
                            int pairId=Vertex2Pair[curVer][Tindex];
							int thirdEdge=PairId2Edge[pairId].second;
                            if(PairId2Edge[pairId].first!=currentEdge){
								pointer[curIndex]=Vertex2Pair[curVer].size();
                                break;
                            }
                            pointer[curIndex]=Tindex;
                            nextPair.first=thirdEdge;
                            nextPair.second=1;
		
                            break;
                        }

                    }else if(round==0 ){
                        //Vertex2Pair[curVer].erase(Vertex2Pair[curVer].begin());
	
                        VertexSignal[curVer]+=1;
                        int Tindex=VertexSignal[curVer];
						
                        while(Tindex<Vertex2Pair[curVer].size()){
                            int pairId=Vertex2Pair[curVer][Tindex];
							int thirdEdge=PairId2Edge[pairId].second;
                            if(PairId2Edge[pairId].first!=currentEdge){
								pointer[curIndex]=Vertex2Pair[curVer].size();
                                break;
                            }
							if(thirdEdge==nextPair.first){
								pointer[curIndex]=Tindex;
								nextPair.second+=1;
								break;
							}else if(!AltCompareHyperEdge(thirdEdge, nextPair.first)){
								pointer[curIndex]=Tindex;
								break;
							}else if(AltCompareHyperEdge(thirdEdge, nextPair.first)){
								pointer[curIndex]=Tindex;
								nextPair.first=thirdEdge;
                            	nextPair.second=1;
							}
                            break;
                        }
                    }else if(nextPair.first==-1){
						
						int Tindex=pointer[curIndex];
                        while(Tindex<Vertex2Pair[curVer].size()){
                            int pairId=Vertex2Pair[curVer][Tindex];
							int thirdEdge=PairId2Edge[pairId].second;
                            if(PairId2Edge[pairId].first!=currentEdge){
								pointer[curIndex]=Vertex2Pair[curVer].size();
                                break;
                            }
							if(thirdEdge==lastPair){
								Tindex+=1;
								if(Tindex==Vertex2Pair[curVer].size()){
									pointer[curIndex]=Vertex2Pair[curVer].size();
								}
                                continue;
							}
                            pointer[curIndex]=Tindex;
                            nextPair.first=thirdEdge;
                            nextPair.second=1;
                            break;
                        }
			
					}else{
						
                        int Tindex=pointer[curIndex];
                        while(Tindex<Vertex2Pair[curVer].size()){
                            int pairId=Vertex2Pair[curVer][Tindex];
							int thirdEdge=PairId2Edge[pairId].second;
                            if(PairId2Edge[pairId].first!=currentEdge){
								pointer[curIndex]=Vertex2Pair[curVer].size();
                                break;
                            }
							if(thirdEdge==lastPair){
								Tindex+=1;
								if(Tindex==Vertex2Pair[curVer].size()){
									pointer[curIndex]=Vertex2Pair[curVer].size();
								}
                                continue;
							}
							if(thirdEdge==nextPair.first){
								pointer[curIndex]=Tindex;
								nextPair.second+=1;
								break;
							}else if(!AltCompareHyperEdge(thirdEdge, nextPair.first)){
								pointer[curIndex]=Tindex;
								break;
							}
	
							if(AltCompareHyperEdge(thirdEdge, nextPair.first)){
								pointer[curIndex]=Tindex;
								nextPair.first=thirdEdge;
                            	nextPair.second=1;
							}
                            break;
                        }
                    }
                }
                round+=1;
                if(nextPair.first==-1){
                    break;
                }
				
				lastPair=nextPair.first;
				int thirdEdge=nextPair.first;
				auto it3 = IntersecPair[currentEdge].find(thirdEdge);
				vector<int>* firstCommon=&it3->second;



				auto it2 = IntersecPair[secondEdge].find(thirdEdge);
				
				if(it2!= IntersecPair[secondEdge].end() ){
					secondCommon=&it2->second;
			
					auto it3 = IntersecPair[currentEdge].find(thirdEdge);
					vector<int>* firstCommon=&it3->second;

				
					int TripleCommon=nextPair.second;
					
					int part123=TripleCommon;
					int part12= CommonVer->size()-1-part123;
					int part13=firstCommon->size()-1-part123;
					int part23=secondCommon->size()-1-part123;		
					int part1=hyperedge2node[currentEdge].size()-part12-part13-part123;
					int part2=hyperedge2node[secondEdge].size()-part12-part23-part123;
					int part3=hyperedge2node[thirdEdge].size()-part23-part13-part123;
					int motifNum=DistinguishAAA(part1,part2, part3, part12, part13, part23, part123);
					h_triangle[motifNum-1]+=1;

				}
				
				for(int p=0; p<simi; p++){
					
					int secondEdge2=IdList[p];

					auto it2 = IntersecPair[secondEdge2].find(thirdEdge);
				
					if(it2== IntersecPair[secondEdge2].end() ){
						continue;
					}
					vector<int>* secondCommon2=&it2->second;

					int TripleCommon=nextPair.second;
					
					int part123=TripleCommon;
					int part12= CommonVer->size()-1-part123;
					int part13=firstCommon->size()-1-part123;
					int part23=secondCommon2->size()-1-part123;		
					int part1=hyperedge2node[currentEdge].size()-part12-part13-part123;
					int part2=hyperedge2node[secondEdge2].size()-part12-part23-part123;
					int part3=hyperedge2node[thirdEdge].size()-part23-part13-part123;
					int motifNum=DistinguishAAA(part1,part2, part3, part12, part13, part23, part123);
					h_triangle[motifNum-1]+=1;
			

				}
				
				if(firstCommon->size()==CommonVer->size() && nextPair.second==CommonVer->size()-1){
					signalList[thirdEdge]=currentEdge;
					IdList[simi]=thirdEdge;
					simi+=1;

				}
				
            }
            

		}
	}

}



inline long long convert_id(int hyperedge_a, int hyperedge_b){
	return hyperedge_a * (1LL << 31) + hyperedge_b;
}

int main(int argc, char *argv[])
{
	clock_t start;
	clock_t run_start;
	int progress;
	

	for(int i=0; i<26; i++){
		vector<vector<int>> temp;
		result.insert(make_pair(i,temp));
		h_triangle.push_back(0);
	}
	string dataName="unique-high-primary-school.txt";
	string graphFile = "Dataset/"+dataName;
	// Read data
	start = clock();
	
	read_data(graphFile, node2hyperedge, hyperedge2node, hyperedge2node_set);

    
	int V = node2hyperedge.size(), E = hyperedge2node.size();
	cout << "# of nodes: " << V << '\n';
	cout << "# of hyperedges: " << E << '\n';
	cout << "Reading data done: "
		<< (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "------------------------------------------" << endl << endl;

	//exact count
	run_start = clock();

	
	vector< unordered_set<int> >  TotalSame;
	vector<vector<int>> copyVec=node2hyperedge;
	vector<int> NodeIdList;

	for(int i=0; i<hyperedge2node.size();i++){
		TotalSame.push_back(unordered_set<int> ());

		InclusionList.push_back(vector<int>());
        EnclosureList.push_back(vector<int>());
        IntersectionList.push_back(vector<int>());
		IntersecPair.push_back(unordered_map<int,vector<int>> ());

	}
	start = clock();

	for(int id=0; id<hyperedge2node.size();id++){
		vector<int> connectNode=hyperedge2node[id];
		vector<vector<int>> Hyperedge2Node2Hyperedge;

		//first delete hyperedge in connectNode
		for(size_t i=0; i<connectNode.size();i++){
			int NodeId=connectNode[i];
			if(copyVec[NodeId].size()==1){//this node should be removed			
				copyVec[NodeId][0]=-1;
			}else{
				copyVec[NodeId].erase(copyVec[NodeId].begin());
				Hyperedge2Node2Hyperedge.push_back(copyVec[NodeId]);	
				NodeIdList.push_back(NodeId);		
			}			
		}
	
		preprocess(id,TotalSame, Hyperedge2Node2Hyperedge, NodeIdList);
		Hyperedge2Node2Hyperedge.clear();
		NodeIdList.clear();
	}
    clock_t time1=clock() - start;
	cout << "Adjacency list construction done: "
		<< (double)(time1) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "------------------------------------------" << endl << endl;
	clock_t start2 = clock();

	////////////////////////////////////     build     hyperEdge    network      /////////////////////////
    for(int i=0; i<hyperedge2node.size();i++){
        if(maxNode<hyperedge2node[i].size()){
            maxNode=hyperedge2node[i].size();
        }
    }

	sortedHyperEdge=sortHyperedge2Node(hyperedge2node);
    Compare comp(hyperedge2node);
    for (auto& vec : IntersectionList) {
        std::sort(vec.begin(), vec.end(), comp);
    }
	for (auto& vec : InclusionList) {
        std::sort(vec.begin(), vec.end(), comp);
    }
    Vertex2Pair.resize(V);
    PairId2Edge.resize(TotalPair);
    signalList=vector<int>(E,-1);
	VertexSignal=vector<int>(V,0);
    int id=0;
    for(int i=0; i<sortedHyperEdge.size(); i++){
		int currentEdge=sortedHyperEdge[i];
        for(int j=0; j<IntersectionList[currentEdge].size();j++){
            int secondEdge=IntersectionList[currentEdge][j];
            vector<int>* CommonVer=&IntersecPair[currentEdge][secondEdge];
            for(int k: (*CommonVer)){
                Vertex2Pair[k].push_back(id);
            }
            (*CommonVer).insert((*CommonVer).begin(), id);
            PairId2Edge[id]=make_pair(currentEdge, secondEdge);
            id+=1;
        }
    }
	for(int i=0; i<node2hyperedge.size();i++){
		if(IdSize<node2hyperedge[i].size()){
			IdSize=node2hyperedge[i].size();
		}
	}
    clock_t time2=clock() - start2;
    cout << "build tree done: "
		<< (double)(time2) / CLOCKS_PER_SEC << " sec" << endl;
	cout<<'\n';

	///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// /////////////////////////


	start = clock();

	CountDenseTTT();

	cout << "counting normal AAA done: "
		<< (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
	cout<<'\n';

	double runtime = (double)(clock() - run_start) / CLOCKS_PER_SEC;
	clock_t time3=clock() - start;
	vector<long long> h_triangle_final(20, 0);
	h_triangle_final=TransferPattern(h_triangle_final);
	int index = 0;
	for (int i = 1; i <= 20; i++){
		cout << fixed << "h-triangle " << ++index << ": " << fixed << h_triangle_final[i-1] << endl;
		
	}
	cout << "counting motif data done: "
		<< (double)(time3) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "Total runtime: " << runtime << endl;

	
	return 0;
}
