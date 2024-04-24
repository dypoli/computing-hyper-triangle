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
#include <cstdlib>
#include <chrono>
#include "read_data.cpp"
#include <random>

using namespace std;




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
vector< vector<int> > node2hyperedge;
vector< vector<int> > hyperedge2node;
vector< unordered_set<int> > hyperedge2node_set;
vector< unordered_map<int,int>> PairTypeInc; //<<id1, id2>,intersection>, <<id1, id2>, inclusion>.... note that we assume the pair.first < pair.second
vector< unordered_map<int,vector<int>> > PairTypeInt; //<<id1, id2>,intersection>, <<id1, id2>, inclusion>.... note that we assume the pair.first < pair.second
vector<vector<int>> AdjacentList;
vector< vector<int>> InclusionList;
vector< unordered_set<int> >  TotalSame;
vector< vector<int>> OrderInclusionList;
vector< vector<int>> IntersectionList;
vector< vector<int>> OrderIntersectionList;
//vector< unordered_set<int>> IntersectionSet;
vector<int> sortedHyperEdge;
//vector< unordered_map<int,int>> PairCheck;
vector<unordered_map<int,vector<bool>>> PairCount;
vector<unordered_map<int,int>> PairId;
int TotalPair=0;
int TotalInclude=0;
int TotalIntersec=0;
vector<pair<int,int>> PairId2Edge;
vector<pair<int,int>> Include2Edge;
vector<pair<int,int>> Intersec2Edge;

long long partialAAB=0;
vector< unordered_set<int>> IncludePair; //<<id1, id2>,intersection>, <<id1, id2>, inclusion>.... note that we assume the pair.first < pair.second
vector< unordered_map<int,vector<int>> > IntersecPair; //<<id1, id2>,intersection>, <<id1, id2>, inclusion>.... note that we assume the pair.first < pair.second

void generateRandomNumbers(vector<int> &result, int n, int m) {
 
    // 设置随机种子
    std::srand(std::time(0));

    // 创建一个包含1到n的所有数字的数组
    std::vector<int> allNumbers;
    for (int i = 0; i <= n; ++i) {
        allNumbers.push_back(i);
    }

    // 从数组中随机选择m个数字
    for (int i = 0; i < m; ++i) {
        // 随机生成一个索引
        int randomIndex = std::rand() % allNumbers.size();

        // 将选定的数字添加到结果中
        result.push_back(allNumbers[randomIndex]);

        // 从数组中移除已选定的数字，以确保不重复选择
        allNumbers.erase(allNumbers.begin() + randomIndex);
    }

}



// 随机洗牌id数组
std::vector<int> ShuffleIds( long long TotalInclude) {
	std::vector<int> ids(TotalInclude);
    for (int i = 0; i < TotalInclude; ++i) {
        ids[i] = i;
    }
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(ids.begin(), ids.end(), gen);
	return ids;
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

vector<int>* GetIntersectPair(int firstID, int secondID){
    if(AltCompareHyperEdge(firstID, secondID)){
        auto it3 = IntersecPair[firstID].find(secondID);
        if(it3==IntersecPair[firstID].end()){
            return NULL;
        }
		return &it3->second;
    }else{
        auto it3 = IntersecPair[secondID].find(firstID);
        if(it3==IntersecPair[secondID].end()){
            return NULL;
        }
		return &it3->second;
    }
}

vector<int>* GetIncludePair(int firstID, int secondID){
    if(hyperedge2node[firstID].size()<hyperedge2node[secondID].size()){
        auto it3 = IncludePair[firstID].find(secondID);
        if(it3!=IncludePair[firstID].end()){
            return &hyperedge2node[firstID];
        }

    }else if(hyperedge2node[firstID].size()>hyperedge2node[secondID].size()){
        auto it3 = IncludePair[secondID].find(firstID);
        if(it3!=IncludePair[secondID].end()){
            return &hyperedge2node[secondID];
        }
	
    }
    return NULL;
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

void TimeSet(int sampling_Time, double & A_sampling_Time_1, double & A_sampling_Time_2, double & B_sampling_Time_1, double & B_sampling_Time_2){

	double B_sampling_Time=(double)sampling_Time*(double)TotalInclude/(double)TotalPair;
	B_sampling_Time_2=(double)B_sampling_Time*(double)TotalInclude/(double)TotalPair;
	B_sampling_Time_2=B_sampling_Time_2>1.0?B_sampling_Time_2:1.0;


	B_sampling_Time_1=(double)B_sampling_Time-B_sampling_Time_2;
    B_sampling_Time_1=B_sampling_Time_1>1.0?B_sampling_Time_1:1.0;

	double A_sampling_Time=(double)sampling_Time-(double)B_sampling_Time;

	if(A_sampling_Time<1.0){
		A_sampling_Time_1=1.0;
		A_sampling_Time_2=1.0;
	}else{
		A_sampling_Time_2=(double)A_sampling_Time*(double)TotalInclude/(double)TotalPair;
		A_sampling_Time_2=A_sampling_Time_2>1.0?A_sampling_Time_2:1.0;
		A_sampling_Time_1=(double)A_sampling_Time-A_sampling_Time_2;
		A_sampling_Time=A_sampling_Time>1.0?A_sampling_Time:1.0;
	}

}


void preprocess(int currentId,vector< unordered_set<int> >& TotalSame, vector< vector<int>>& InclusionList,  vector< vector<int> >& hyperedge2node, vector< unordered_map<int,int> >& PairTypeInc, vector< unordered_map<int,vector<int>> >& PairTypeInt, vector<vector<int>>& Hyperedge2Node2Hyperedge, vector<int>& NodeIdList){
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
			TotalInclude+=1;
            TotalPair+=1;
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
				IncludePair[MinimumId].insert(currentId);
                //Include2Edge[TotalInclude-1]=make_pair(MinimumId, currentId);
                OrderInclusionList[MinimumId].push_back(currentId);
			}else{
				IncludePair[currentId].insert(MinimumId);
                //Include2Edge[TotalInclude-1]=make_pair(currentId, MinimumId);
                OrderInclusionList[currentId].push_back(MinimumId);
			}
            AdjacentList[MinimumId].push_back(currentId);
            AdjacentList[currentId].push_back(MinimumId);
            InclusionList[MinimumId].push_back(currentId);
            InclusionList[currentId].push_back(MinimumId);
		}else if(InclusionNum!=EdgeNodeNum2 && InclusionNum!=EdgeNodeNum1){//intersection
			TotalIntersec+=1;
            TotalPair+=1;
			if(AltCompareHyperEdge(currentId, MinimumId)){
				IntersecPair[currentId].insert(make_pair(MinimumId,vector<int>()));
                //Intersec2Edge[TotalIntersec-1]=make_pair(currentId, MinimumId);
                OrderIntersectionList[currentId].push_back(MinimumId);
			}else{
				IntersecPair[MinimumId].insert(make_pair(currentId,vector<int>()));
                //Intersec2Edge[TotalIntersec-1]=make_pair(MinimumId, currentId);
                OrderIntersectionList[MinimumId].push_back(currentId);
				
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
            AdjacentList[MinimumId].push_back(currentId);
            AdjacentList[currentId].push_back(MinimumId);
			IntersectionList[currentId].push_back(MinimumId);
            IntersectionList[MinimumId].push_back(currentId);
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


void ApproOptimal(int sampling_Time){
   

    double A_sampling_Time_1, A_sampling_Time_2, B_sampling_Time_1, B_sampling_Time_2;
   	TimeSet(sampling_Time, A_sampling_Time_1, A_sampling_Time_2, B_sampling_Time_1, B_sampling_Time_2);
	/////////////////////   approxiamtely count AAA //////////////////////////////

	
    long long A_sampling_time=0;
    long long B_sampling_time=0;
    vector<long long> temp_motif;
    temp_motif.resize(26, 0);
	vector<int> randomList=ShuffleIds(TotalIntersec);
	int index=0;
    auto startTime = std::chrono::system_clock::now();
    while(true){
        
        std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - startTime;

        if (elapsed_seconds.count() >= A_sampling_Time_1) {
            break;
        }
		if(index>=TotalIntersec){
			break;
		}

        A_sampling_time+=1;
        int PairId = randomList[index];
		index+=1;
        int currentEdge=Intersec2Edge[PairId].first;;
        int secondEdge=Intersec2Edge[PairId].second;
		if(OrderIntersectionList[currentEdge].size()==0){
			continue;
		}
        
        vector<int>* CommonVer=&IntersecPair[currentEdge][secondEdge];

        for(int j=0;j<OrderIntersectionList[secondEdge].size(); j++){
   
            
            int thirdEdge=OrderIntersectionList[secondEdge][j];

            if(OrderIntersectionList[currentEdge][OrderIntersectionList[currentEdge].size()-1]< thirdEdge){
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
            int ii = 0, jj = 0;
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
            int part12= CommonVer->size()-part123;
            int part13=firstCommon->size()-part123;
            int part23=secondCommon->size()-part123;
            int part1=hyperedge2node[currentEdge].size()-part12-part13-part123;
            int part2=hyperedge2node[secondEdge].size()-part12-part23-part123;
            int part3=hyperedge2node[thirdEdge].size()-part23-part13-part123;
            int motifNum=DistinguishAAA(part1,part2, part3, part12, part13, part23, part123);
            temp_motif[motifNum-1]+=1;
            
        }
            
    }
    

    h_triangle[1]=temp_motif[1]*TotalIntersec/(A_sampling_time);
    h_triangle[5]=temp_motif[5]*TotalIntersec/(A_sampling_time);
    for(int i=10; i<16; i++){
        h_triangle[i]=temp_motif[i]*TotalIntersec/(A_sampling_time);
    }
    for(int i=22; i<26; i++){
        h_triangle[i]=temp_motif[i]*TotalIntersec/(A_sampling_time);
    }


	/////////////////////   approxiamtely count ABB //////////////////////////////
    
	A_sampling_time=0;
	
	randomList=ShuffleIds(TotalIntersec);
	index=0;
    startTime = std::chrono::system_clock::now();
	while(true){
        
        std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - startTime;

        if (elapsed_seconds.count() >= A_sampling_Time_2) {
            break;
        }
		if(index>=TotalIntersec){
			break;
		}

        A_sampling_time+=1;
        int PairId = randomList[index];
		index+=1;
        int currentEdge=Intersec2Edge[PairId].first;;
        int secondEdge=Intersec2Edge[PairId].second;
        if(InclusionList[currentEdge].size()==0){
			continue;
		}
        vector<int>* CommonVer=&IntersecPair[currentEdge][secondEdge];

        for(int j=0;j<InclusionList[secondEdge].size(); j++){
   
            
            int thirdEdge=InclusionList[secondEdge][j];
			
			if(InclusionList[currentEdge][InclusionList[currentEdge].size()-1]<thirdEdge){
				break;
			}
			
            vector<int>* firstCommon=GetIncludePair(currentEdge, thirdEdge);
			if(firstCommon==NULL ){
                continue;
            }
			
			vector<int>* secondCommon=GetIncludePair(secondEdge, thirdEdge);

           	int TripleCommon=secondCommon->size()<firstCommon->size()?secondCommon->size():firstCommon->size();
			TripleCommon=TripleCommon<CommonVer->size()?TripleCommon:CommonVer->size();
			int part123=TripleCommon;
			int part12= CommonVer->size()-part123;
			int part13=firstCommon->size()-part123;
			int part23=secondCommon->size()-part123;
			int part1=hyperedge2node[currentEdge].size()-part12-part13-part123;
			int part2=hyperedge2node[secondEdge].size()-part12-part23-part123;
			int part3=hyperedge2node[thirdEdge].size()-part23-part13-part123;
			
			if(part12==0 && part23==0 && part13==0){
				temp_motif[0]+=1;
			}else if(part12*part23 + part13*part23 + part12*part13==0 ){
				temp_motif[3]+=1;
			}else if(part1+part2+part3==0 ){
				temp_motif[6]+=1;
			}else{
				temp_motif[7]+=1;
			}

		}
    }

    h_triangle[0]=temp_motif[0]*TotalIntersec/(A_sampling_time);
    h_triangle[3]=temp_motif[3]*TotalIntersec/(A_sampling_time);
    h_triangle[6]=temp_motif[6]*TotalIntersec/(A_sampling_time);
    h_triangle[7]=temp_motif[7]*TotalIntersec/(A_sampling_time);
    
	


	/////////////////////   approxiamtely count AAB //////////////////////////////

	randomList=ShuffleIds(TotalInclude);
	index=0;
	B_sampling_time=0;
	startTime = std::chrono::system_clock::now();
	while(true){
        
        std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - startTime;

        if (elapsed_seconds.count() >= B_sampling_Time_1) {
            break;
        }
		if(index>=TotalInclude){
			break;
		}

        B_sampling_time+=1;
        int PairId = randomList[index];
		index+=1;
        int currentEdge=Include2Edge[PairId].first;
        int secondEdge=Include2Edge[PairId].second;
        if(IntersectionList[currentEdge].size()==0){
			continue;
		}
        vector<int>* CommonVer=GetIncludePair(currentEdge, secondEdge);

        for(int j=0;j<IntersectionList[secondEdge].size(); j++){
   
            
            int thirdEdge=IntersectionList[secondEdge][j];
			
			if(IntersectionList[currentEdge][IntersectionList[currentEdge].size()-1]<thirdEdge){
				break;
			}
			vector<int>* firstCommon=GetIntersectPair(currentEdge, thirdEdge);
            if(firstCommon==NULL ){
                continue;
            }

			vector<int>* secondCommon=GetIntersectPair(secondEdge, thirdEdge);


           	int TripleCommon=CommonVer->size()<firstCommon->size()?CommonVer->size():firstCommon->size();
			TripleCommon=TripleCommon<secondCommon->size()?TripleCommon:secondCommon->size();
			int part123=TripleCommon;
			int part12= CommonVer->size()-part123;
			int part13=firstCommon->size()-part123;
			int part23=secondCommon->size()-part123;
			int part1=hyperedge2node[currentEdge].size()-part12-part13-part123;
			int part2=hyperedge2node[secondEdge].size()-part12-part23-part123;
			int part3=hyperedge2node[thirdEdge].size()-part23-part13-part123;
			
			if(part12*part23 + part13*part23 + part12*part13==0 ){
				temp_motif[4]+=1;
			}else if(part1*part2 + part1*part3 + part2*part3==0 ){
				temp_motif[8]+=1;
			}else{
				temp_motif[9]+=1;
			}

		}
    }
	h_triangle[4]=temp_motif[4]*TotalInclude/(1*B_sampling_time);
    h_triangle[8]=temp_motif[8]*TotalInclude/(1*B_sampling_time);
    h_triangle[9]=temp_motif[9]*TotalInclude/(1*B_sampling_time);

	/////////////////////   approxiamtely count BBB //////////////////////////////
	B_sampling_time=0;
	
	randomList=ShuffleIds(TotalInclude);
	index=0;
    startTime = std::chrono::system_clock::now();
	while(true){
        
        std::chrono::duration<double> elapsed_seconds = std::chrono::system_clock::now() - startTime;

        if (elapsed_seconds.count() >= B_sampling_Time_2) {
            break;
        }
		if(index>=TotalInclude){
			break;
		}

        B_sampling_time+=1;
        int PairId = randomList[index];
		index+=1;
        int currentEdge=Include2Edge[PairId].first;
        int secondEdge=Include2Edge[PairId].second;
        if(OrderInclusionList[currentEdge].size()==0){
			continue;
		}
        vector<int>* CommonVer=GetIncludePair(currentEdge, secondEdge);

        for(int j=0;j<OrderInclusionList[secondEdge].size(); j++){
   
            
            int thirdEdge=OrderInclusionList[secondEdge][j];
			if(thirdEdge==currentEdge){
				continue;
			}
			
			if(OrderInclusionList[currentEdge][OrderInclusionList[currentEdge].size()-1]<thirdEdge){
				break;
			}
			vector<int>* firstCommon=GetIncludePair(currentEdge, thirdEdge);
            if(firstCommon==NULL ){
                continue;
            }
			temp_motif[2]+=1;
			

		}
    }
	h_triangle[2]=temp_motif[2]*TotalInclude/(1*B_sampling_time);
}



inline long long convert_id(int hyperedge_a, int hyperedge_b){
	return hyperedge_a * (1LL << 31) + hyperedge_b;
}

int main(int argc, char *argv[])
{
	clock_t start;
	clock_t run_start;
	int progress;
	int sampling_time = stoi(argv[1]);
	//for test
	map<int,vector<vector<int>>> result2;
	vector<long long> h_triangle2;

	for(int i=0; i<26; i++){
		vector<vector<int>> temp;
		result.insert(make_pair(i,temp));
		h_triangle.push_back(0);

		vector<vector<int>> temp2;
		result2.insert(make_pair(i,temp2));
		h_triangle2.push_back(0);
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

	
	run_start = clock();
	vector<vector<int>> copyVec=node2hyperedge;
	vector<int> NodeIdList;

	
	for(int i=0; i<hyperedge2node.size();i++){
		PairTypeInc.push_back(unordered_map<int,int>());
		PairTypeInt.push_back(unordered_map<int,vector<int>> ());

		TotalSame.push_back(unordered_set<int> ());
		InclusionList.push_back(vector<int>());
        IntersectionList.push_back(vector<int>());
        OrderInclusionList.push_back(vector<int>());
        OrderIntersectionList.push_back(vector<int>());
        AdjacentList.push_back(vector<int>());
        IntersectionList.push_back(vector<int>());

		PairCount.push_back(unordered_map<int,vector<bool>>());

		IncludePair.push_back(unordered_set<int>());
		IntersecPair.push_back(unordered_map<int,vector<int>> ());

		PairId.push_back(unordered_map<int,int>());
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
	
		preprocess(id,TotalSame,InclusionList,  hyperedge2node, PairTypeInc, PairTypeInt, Hyperedge2Node2Hyperedge, NodeIdList);
		Hyperedge2Node2Hyperedge.clear();
		NodeIdList.clear();
	}
    Intersec2Edge.resize(TotalIntersec);
    int nn=0;
    for(int i=0;  i<IntersectionList.size(); i++){
        for(int j=0; j<IntersectionList[i].size();j++){
            int a= i;
            int b=IntersectionList[i][j];
            if(AltCompareHyperEdge(a,b)){
                Intersec2Edge[nn]=make_pair(a, b);
                nn+=1;
            }
        }
    }

    Include2Edge.resize(TotalInclude);
    nn=0;
    for(int i=0;  i<InclusionList.size(); i++){
        for(int j=0; j<InclusionList[i].size();j++){
            int a= i;
            int b=InclusionList[i][j];
            if(AltCompareHyperEdge(a,b)){
                Include2Edge[nn]=make_pair(a, b);
                nn+=1;
            }
        }
    }


	clock_t time1=clock() - start;
	cout << "Adjacency list construction done: "
		<< (double)(time1) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "------------------------------------------" << endl << endl;
	clock_t start2 = clock();

	////////////////////////////////////     build     hyperEdge    network      /////////////////////////


	clock_t time2=clock() - start2;

	///////////////////////// ///////////////////////// ///////////////////////// ///////////////////////// /////////////////////////

	start = clock();

    ApproOptimal(sampling_time);
	

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
