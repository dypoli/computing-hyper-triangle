
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
#include <list>
#include <climits>
#include <sstream>
#include <string>
using namespace std;



// struct used to store the hyperedge

struct Hyperedge { 
    int ID;
    Hyperedge* prev;
    Hyperedge* next;
    Hyperedge(int id) : ID(id), prev(nullptr), next(nullptr) {}
};


struct Triangle { 
    int sigal;
    int node_1;
    int node_2;
    int node_3;
    Triangle(int a, int b, int c) : sigal(1), node_1(a), node_2(b), node_3(c){}
};

struct BloomNode; 
struct EdgeInfo;
struct ChildNode;

struct ChildNode { 
    int Id;
    vector<BloomNode*> Bloom_List;
    ChildNode(int id) : Id(id), Bloom_List(vector<BloomNode*>()) {}
};

struct BloomNode { 
    int sigal;
    int node_1;
    int node_2;
    vector<ChildNode*> Child_List;
    BloomNode(int a, int b) : sigal(1), node_1(a), node_2(b), Child_List(vector<ChildNode*>()) {}
};

struct EdgeInfo {
    vector<BloomNode*> Bloom_List;
    ChildNode* Child;
    EdgeInfo(): Bloom_List(vector<BloomNode*>()), Child(nullptr) {}
};

struct SupportLayer {
    int support;
    unordered_set<int> hyperedges;
    SupportLayer(int supportnum): support(supportnum) {}
};

struct SupportData {
    vector<SupportLayer*> Support_List;
    vector<SupportLayer*> HyperedgeToSupport_List;
};



inline long long convert_id(int hyperedge_a, int hyperedge_b){
	return hyperedge_a * (1LL << 31) + hyperedge_b;
}

/*       help function, used to check the list       */


vector<int> get_sorted_hyperedges_by_degree(const vector<vector<int>>& hyperedge_adj) {
    int n = hyperedge_adj.size();
    vector<pair<int, int>> degrees;  // pair<degree, id>

    for (int i = 0; i < n; ++i) {
        degrees.emplace_back(hyperedge_adj[i].size(), i);
    }


    sort(degrees.begin(), degrees.end(), [](const pair<int, int>& a, const pair<int, int>& b) {
        if (a.first != b.first) return a.first > b.first;
        return a.second > b.second;
    });

    vector<int> sorted_ids;
    for (const auto& p : degrees) {
        sorted_ids.push_back(p.second);
    }

    return sorted_ids;
}

void sort_adj_lists_by_degree_desc(vector<vector<int>>& hyperedge_adj) {
    int n = hyperedge_adj.size();
    vector<int> degree(n);
    for (int i = 0; i < n; ++i) {
        degree[i] = hyperedge_adj[i].size();
    }

    for (int i = 0; i < n; ++i) {
        sort(hyperedge_adj[i].begin(), hyperedge_adj[i].end(), [&](int a, int b) {
            if (degree[a] != degree[b]) return degree[a] > degree[b];
            return a > b;
        });
    }
}

void filter_and_sort_adj_lists_by_degree(vector<vector<int>>& hyperedge_adj) {
    int n = hyperedge_adj.size();
    vector<int> degree(n);

    for (int i = 0; i < n; ++i) {
        degree[i] = hyperedge_adj[i].size();
    }

    for (int i = 0; i < n; ++i) {
        vector<int> filtered;

        for (int j : hyperedge_adj[i]) {
            if (degree[j] < degree[i] || (degree[j] == degree[i] && j < i)) {
                filtered.push_back(j);
            }
        }

        sort(filtered.begin(), filtered.end(), [&](int a, int b) {
            if (degree[a] != degree[b]) return degree[a] > degree[b];
            return a > b;
        });

        hyperedge_adj[i] = std::move(filtered);
    }
}


std::vector<Hyperedge*> buildHyperedgePtrMap(Hyperedge* head, int total_hyperedges) {
    std::vector<Hyperedge*> hyperedge_ptrs(total_hyperedges, nullptr);
    for (Hyperedge* curr = head; curr != nullptr; curr = curr->next) {
        hyperedge_ptrs[curr->ID] = curr;
    }
    return hyperedge_ptrs;
}

Hyperedge* buildHyperedgeList(const std::vector<int>& hyperedge_sup) {
    std::vector<Hyperedge*> nodes;
    for (int i = 0; i < hyperedge_sup.size(); ++i) {
        nodes.push_back(new Hyperedge(i));
    }

   
    std::sort(nodes.begin(), nodes.end(), [&](Hyperedge* a, Hyperedge* b) {
        if (hyperedge_sup[a->ID] != hyperedge_sup[b->ID])
            return hyperedge_sup[a->ID] < hyperedge_sup[b->ID];
        return a->ID < b->ID;
    });

    
    for (int i = 0; i < nodes.size(); ++i) {
        if (i > 0) nodes[i]->prev = nodes[i - 1];
        if (i < nodes.size() - 1) nodes[i]->next = nodes[i + 1];
    }

    return nodes.empty() ? nullptr : nodes[0]; 
}

void deleteHyperedgeList(Hyperedge* head) {
    while (head) {
        Hyperedge* temp = head;
        head = head->next;
        delete temp;
    }
}

vector<int> buildSortedHyperedgesBySup(const vector<int>& hyperedge_sup) {
    int n = hyperedge_sup.size();
    vector<int> SortSup(n);

    
    for (int i = 0; i < n; ++i) {
        SortSup[i] = i;
    }

    
    sort(SortSup.begin(), SortSup.end(), [&](int a, int b) {
        if (hyperedge_sup[a] != hyperedge_sup[b])
            return hyperedge_sup[a] < hyperedge_sup[b];
        return a < b;
    });

    return SortSup;
}


vector<vector<int>> buildHyperedgeAdjSup(const vector<vector<int>>& hyperedge_adj, const vector<int>& hyperedge_sup) {
    vector<vector<int>> hyperedge_adjSup(hyperedge_adj.size());

    for (size_t i = 0; i < hyperedge_adj.size(); ++i) {
        for (int j : hyperedge_adj[i]) {
            if (
                hyperedge_sup[j] > hyperedge_sup[i] ||
                (hyperedge_sup[j] == hyperedge_sup[i] && j > static_cast<int>(i))
            ) {
                hyperedge_adjSup[i].push_back(j);
            }
        }

        
        sort(hyperedge_adjSup[i].begin(), hyperedge_adjSup[i].end(),
            [&](int a, int b) {
                if (hyperedge_sup[a] != hyperedge_sup[b])
                    return hyperedge_sup[a] < hyperedge_sup[b];
                return a < b;
            }
        );
    }

    return hyperedge_adjSup;
}

map<int, vector<int>> groupIndicesBySupport(const vector<int>& hyperedge_sup) {
    map<int, vector<int>> hyperedge_supLayer;

    for (size_t i = 0; i < hyperedge_sup.size(); ++i) {
        int sup = hyperedge_sup[i];
        hyperedge_supLayer[sup].push_back(i);
    }

    return hyperedge_supLayer;
}

void writeHyperedgeAdjToFile(const vector<vector<int>>& hyperedge_adj, const string& filename) {
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        return;
    }

    for (const auto& edge : hyperedge_adj) {
        for (size_t i = 0; i < edge.size(); ++i) {
            outfile << edge[i];
            if (i + 1 < edge.size()) outfile << " ";
        }
        outfile << "\n";
    }

    outfile.close();
    cout << "Data written to " << filename << endl;
}

void removeDuplicatesFromHyperedgeAdj(vector<vector<int>>& hyperedge_adj) {
    for (auto& edge : hyperedge_adj) {
        edge.erase(unique(edge.begin(), edge.end()), edge.end());
    }
}


void completeSymmetricValues(
    const std::vector<std::vector<int>>& hyperedge_adj,
    std::vector<std::vector<int>>& hyperedge_TriangleNum
) {
    int n = hyperedge_adj.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < hyperedge_adj[i].size(); ++j) {
            int nei = hyperedge_adj[i][j];
            if (i < nei) {
                
                auto& nei_list = hyperedge_adj[nei];
                auto it = std::find(nei_list.begin(), nei_list.end(), i);
                if (it != nei_list.end()) {
                    int pos = it - nei_list.begin();
                   
                    hyperedge_TriangleNum[nei][pos] = hyperedge_TriangleNum[i][j];
                }
            }
        }
    }
}



vector<unordered_set<int>> convertToSet(const vector<vector<int>>& hyperedge_adj) {
    vector<unordered_set<int>> hyperedge_adj_set;
    hyperedge_adj_set.reserve(hyperedge_adj.size());

    for (const auto& vec : hyperedge_adj) {
        unordered_set<int> s(vec.begin(), vec.end());
        hyperedge_adj_set.push_back(std::move(s));
    }

    return hyperedge_adj_set;
}

vector<unordered_set<int>> convertToSetFiltered(
    const vector<vector<int>>& hyperedge_adj,
    const vector<int>& hyperedge2order
) {
    vector<unordered_set<int>> hyperedge_adj_set;
    hyperedge_adj_set.reserve(hyperedge_adj.size());

    for (int i = 0; i < hyperedge_adj.size(); ++i) {
        unordered_set<int> s;
        int current_order = hyperedge2order[i];

        for (int neighbor : hyperedge_adj[i]) {
            if (hyperedge2order[neighbor] > current_order) {
                s.insert(neighbor);
            }
        }

        hyperedge_adj_set.push_back(std::move(s));
    }

    return hyperedge_adj_set;
}





vector<unordered_set<int>> convertToSetFilteredByIndex(const vector<vector<int>>& hyperedge_adj) {
    vector<unordered_set<int>> hyperedge_adj_set;
    hyperedge_adj_set.reserve(hyperedge_adj.size());

    for (int i = 0; i < hyperedge_adj.size(); ++i) {
        unordered_set<int> s;
        for (int neighbor : hyperedge_adj[i]) {
            if (neighbor > i) {
                s.insert(neighbor);
            }
        }
        hyperedge_adj_set.push_back(std::move(s));
    }

    return hyperedge_adj_set;
}


void WriteZeroPatternToTxt(const string &filename,
                           const vector<vector<int>> &ZeroPatternList) {
    ofstream fout(filename);
    if (!fout.is_open()) {
        cerr << "Unable to open file: " << filename << endl;
        return;
    }

    for (const auto &triple : ZeroPatternList) {
        fout << triple[0] << "," << triple[1] << "," << triple[2] << "\n";
    }

    fout.close();
}



vector<vector<int>> convertToVectorFiltered(
    const vector<vector<int>>& hyperedge_adj,
    const vector<int>& hyperedge2order
) {
    vector<vector<int>> hyperedge_adj_set;
    hyperedge_adj_set.resize(hyperedge_adj.size());

    for (int i = 0; i < hyperedge_adj.size(); ++i) {
        vector<int> s;
        int current_order = hyperedge2order[i];

        for (int neighbor : hyperedge_adj[i]) {
            if (hyperedge2order[neighbor] > current_order) {
                s.push_back(neighbor);
            }
        }

        hyperedge_adj_set[i]=s;
    }

    return hyperedge_adj_set;
}

vector<vector<int>> convertToVectorWithIndex(
    const vector<vector<int>>& hyperedge_adj,
    vector<int>& Index_signal
) {
    vector<vector<int>> hyperedge_adj_set;
    hyperedge_adj_set.resize(hyperedge_adj.size());

    for (int i = 0; i < hyperedge_adj.size(); ++i) {
        vector<int> s;

        for (int j=0;j<hyperedge_adj[i].size();++j) {
            int neighbor=hyperedge_adj[i][j];
            if (neighbor > i) {
                if(Index_signal[i]==-1) Index_signal[i]=j;
            }
        }

    }

    return hyperedge_adj_set;
}

vector<vector<int>> convertToVectorFilteredByIndex(
    const vector<vector<int>>& hyperedge_adj
) {
    vector<vector<int>> hyperedge_adj_set;
    hyperedge_adj_set.resize(hyperedge_adj.size());

    for (int i = 0; i < hyperedge_adj.size(); ++i) {
        vector<int> s;

        for (int neighbor : hyperedge_adj[i]) {
            if (neighbor > i) {
                s.push_back(neighbor);
            }
        }

        hyperedge_adj_set[i]=s;
    }

    return hyperedge_adj_set;
}

void removeNode(Hyperedge*& head, Hyperedge* node) {
    if (node->prev) node->prev->next = node->next;
    else head = node->next; 

    if (node->next) node->next->prev = node->prev;

    node->prev = node->next = nullptr;
}

void insertNodeSorted(Hyperedge*& head, Hyperedge* node, const std::vector<int>& hyperedge_sup) {
    if (!head) {
        head = node;
        return;
    }

    Hyperedge* curr = head;
    while (curr) {
        int support_curr = hyperedge_sup[curr->ID];
        int support_node = hyperedge_sup[node->ID];

        if (support_node < support_curr || 
            (support_node == support_curr && node->ID < curr->ID)) {
           
            node->next = curr;
            node->prev = curr->prev;
            if (curr->prev) curr->prev->next = node;
            else head = node;
            curr->prev = node;
            return;
        }

        if (!curr->next) break;
        curr = curr->next;
    }


    curr->next = node;
    node->prev = curr;
    node->next = nullptr;
}

void RecurseFunc(int depth, vector<int>& ValidHyperedge, vector<int>& hyperedge_sup, vector<int>& hyperedge_adj_i, 
    vector<unordered_set<int>>& hyperedge_adjset, unordered_map<int, int>& hyperwedge_adj_i,
    vector<Hyperedge*>& hyperedge_ptrs, Hyperedge* iter) {

    int hyperedge_2=hyperedge_adj_i[depth];
    if (depth == hyperedge_adj_i.size()-1) {
        if(hyperedge_sup[hyperedge_2]>0) {
            ValidHyperedge.insert(ValidHyperedge.begin(), hyperedge_2);
        }
        return;
    }else{
        RecurseFunc(depth+1, ValidHyperedge, hyperedge_sup, hyperedge_adj_i, hyperedge_adjset, hyperwedge_adj_i, hyperedge_ptrs, iter);  
        if(hyperedge_sup[hyperedge_2]<=0) return;
        ValidHyperedge.insert(ValidHyperedge.begin(), hyperedge_2);
		if(hyperwedge_adj_i[hyperedge_2]==0) return;

        for(int j=1; j<ValidHyperedge.size(); j++){
            int hyperedge_3=ValidHyperedge[j];
            if(!hyperedge_adjset[hyperedge_2].count(hyperedge_3)) continue;

            hyperedge_sup[hyperedge_2]-=1;
            hyperedge_sup[hyperedge_3]-=1;
            removeNode(iter, hyperedge_ptrs[hyperedge_2]);
            removeNode(iter, hyperedge_ptrs[hyperedge_3]);

            
            insertNodeSorted(iter, hyperedge_ptrs[hyperedge_2], hyperedge_sup);
            insertNodeSorted(iter, hyperedge_ptrs[hyperedge_3], hyperedge_sup);

        }
    }

}


SupportData BuildSupportData(const vector<int>& hyperedge_sup) {
  
    int max_support = 0;
    for (int s : hyperedge_sup) {
        if (s > max_support) max_support = s;
    }


    vector<SupportLayer*> Support_List;
    for (int s = 0; s <= max_support; ++s) {
        Support_List.push_back(new SupportLayer(s));
    }


    vector<SupportLayer*> HyperedgeToSupport_List(hyperedge_sup.size());

    for (int hyperedge_id = 0; hyperedge_id < hyperedge_sup.size(); ++hyperedge_id) {
        int support = hyperedge_sup[hyperedge_id];
        SupportLayer* layer = Support_List[support];
        layer->hyperedges.insert(hyperedge_id);
        HyperedgeToSupport_List[hyperedge_id] = layer;
    }

    return {Support_List, HyperedgeToSupport_List};
}

queue<int> SetToQueue(const unordered_set<int>& input_set) {
    queue<int> q;
    for (int val : input_set) {
        q.push(val);
    }
    return q;
}

queue<int> VectorToQueue(const vector<int>& input_set) {
    queue<int> q;
    for (int i=0; i<input_set.size(); i++) {
        q.push(input_set[i]);
    }
    return q;
}




std::vector<std::unordered_set<int>> BuildSupportList(const std::vector<int>& hyperedge_sup) {
 
    int max_support = 0;
    for (int sup : hyperedge_sup) {
        max_support = std::max(max_support, sup);
    }


    std::vector<std::unordered_set<int>> Support_List(max_support + 1);


    for (int hyperedge = 0; hyperedge < hyperedge_sup.size(); ++hyperedge) {
        int support = hyperedge_sup[hyperedge];
        Support_List[support].insert(hyperedge);
    }

    return Support_List;
}

void sortHyperedges(
    const vector<int>& hyperedge_NodeDegree,
    vector<vector<int>>& hyperedge_adj,
    vector<int>& hyperedge_order
) {
    int E = hyperedge_NodeDegree.size();

  
    for (int i = 0; i < E; i++) {
        sort(hyperedge_adj[i].begin(), hyperedge_adj[i].end(), [&](int a, int b) {
            if (hyperedge_NodeDegree[a] != hyperedge_NodeDegree[b])
                return hyperedge_NodeDegree[a] < hyperedge_NodeDegree[b];
            return a < b; // tie-breaker: smaller id first
        });
    }


    hyperedge_order.resize(E);
    for (int i = 0; i < E; i++) hyperedge_order[i] = i;

    sort(hyperedge_order.begin(), hyperedge_order.end(), [&](int a, int b) {
        if (hyperedge_NodeDegree[a] != hyperedge_NodeDegree[b])
            return hyperedge_NodeDegree[a] < hyperedge_NodeDegree[b];
        return a < b; // tie-breaker
    });
}

bool OrderCheck(const vector<int>& hyperedge_NodeDegree, int a, int b){
    if (hyperedge_NodeDegree[a] != hyperedge_NodeDegree[b])
            return hyperedge_NodeDegree[a] < hyperedge_NodeDegree[b];
        return a < b; // tie-breaker: smaller id first
}


vector<unordered_set<int>> generateHyperedgeNeighborNodeSet(
    const vector<unordered_set<int>>& hyperedge_adjset,
    const vector<vector<int>>& hyperedge2node
) {
    int E = hyperedge_adjset.size();
    vector<unordered_set<int>> hyperedge_NeighborNodeSet(E);

    for (int i = 0; i < E; ++i) {
        for (int neighbor : hyperedge_adjset[i]) {
            for (int node : hyperedge2node[neighbor]) {
                hyperedge_NeighborNodeSet[i].insert(node);
            }
        }
    }

    return hyperedge_NeighborNodeSet;
}


std::vector<std::vector<int>> readBitmapIndices(const std::string& filename) {
    std::vector<std::vector<int>> result;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "OPEN FAIL: " << filename << std::endl;
        return result;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::vector<int> indices;
        std::istringstream iss(line);
        int num;
        while (iss >> num) {
            indices.push_back(num);
        }
        result.push_back(indices);
    }

    file.close();
    return result;
}


std::vector<std::list<int>> convertToListVector(const std::vector<std::vector<int>>& vecVec) {
    std::vector<std::list<int>> result;
    result.reserve(vecVec.size()); // 

    for (const auto& innerVec : vecVec) {
        result.emplace_back(innerVec.begin(), innerVec.end());
    }

    return result;
}
/*       end of the help function         */


void WriteDegreeToFile(const vector<int>& Degree_hyperedges, const string& filename) {
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return;
    }

    for (int degree : Degree_hyperedges) {
        outfile << degree << '\n';
    }

    outfile.close();
}


void WriteAdjacencyToFile(const vector<vector<int>>& hyperedge_adj, const string& filename) {
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return;
    }

    for (const auto& row : hyperedge_adj) {
        for (size_t i = 0; i < row.size(); ++i) {
            outfile << row[i];
            if (i != row.size() - 1) outfile << ' ';  
        }
        outfile << '\n';  
    }

    outfile.close();
}

void loadAffectedData(
        const std::string &filename,
        std::vector<std::vector<int>> &hyperedge_adj_new,
        std::vector<std::vector<int>> &node2hyperedge_new) 
{
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        std::stringstream ss(line);
        int hyperedge_id;
        char colon;

        ss >> hyperedge_id >> colon;  // e.g. read "7498:"

        if (!ss || colon != ':') continue; // Malformed line, skip

        int node;
        while (ss >> node) {
            // 1) Insert into hyperedge_adj_new
            if (hyperedge_id >= hyperedge_adj_new.size())
                hyperedge_adj_new.resize(hyperedge_id + 1);

            hyperedge_adj_new[hyperedge_id].push_back(node);

            // 2) Insert into node2hyperedge_new
            if (node >= node2hyperedge_new.size())
                node2hyperedge_new.resize(node + 1);

            node2hyperedge_new[node].push_back(hyperedge_id);
        }
    }
}


set<int> loadDeleteVertex(const string &filePath) {
    set<int> deletedVertices;
    ifstream fin(filePath);
    if (!fin.is_open()) {
        cerr << "Error: Cannot open file: " << filePath << endl;
        return deletedVertices;
    }

    int x;
    while (fin >> x) {      // Read integers one by one
        deletedVertices.insert(x);
    }

    fin.close();
    return deletedVertices;
}


void removeDeletedVertices(
        const vector<vector<int>> &hyperedge2node,
        const set<int> &deletedSet,
        vector<vector<int>> &hyperedge2node_new,
        vector<vector<int>> &hyperedge2node_delete,
        set<int> &DeletedHyperedge)
{
    int E = hyperedge2node.size();
    hyperedge2node_new.clear();
    hyperedge2node_new.resize(E);
    DeletedHyperedge.clear();

    for (int e = 0; e < E; ++e) {

        for (int v : hyperedge2node[e]) {
            if (deletedSet.count(v) == 0) {
                hyperedge2node_new[e].push_back(v);
            }else{
                hyperedge2node_delete[e].push_back(v);
            }
        }

        // If the hyperedge becomes empty after deletion, record it
        if (hyperedge2node_new[e].empty()) {
            DeletedHyperedge.insert(e);
        }
    }
}

void removeDeletedHyperedges(
        const vector<vector<int>> &node2hyperedge,
        const set<int> &deletedSet,
        vector<vector<int>> &node2hyperedge_new,
        set<int> &DeletedHyperedge)
{
    for(int i=0; i<node2hyperedge.size(); i++){
        if(deletedSet.count(i)) continue;
        for(int j=0;j<node2hyperedge[i].size();j++){
            if(DeletedHyperedge.count(node2hyperedge[i][j])) continue;
            node2hyperedge_new[i].push_back(node2hyperedge[i][j]);
        }
    }

}


void classifyHyperedgeAdj(  
        const vector<vector<int>> &hyperedge2node_delete,  // Deleted vertices for each hyperedge
        const vector<vector<int>> &hyperedge_adj_update,   // Adjacency before deletion
        vector<vector<int>> &hyperedge_adj_old,            // Unaffected neighbors
        vector<vector<int>> &hyperedge_adj_adjust)         // Neighbors affected by deletion
{
    int E = hyperedge_adj_update.size();

    // Initialize outputs
    hyperedge_adj_old.assign(E, {});
    hyperedge_adj_adjust.assign(E, {});

    for (int e = 0; e < E; ++e) {

        // Iterate its original neighbors
        for (int nb : hyperedge_adj_update[e]) {

            // Neighbor hyperedge was affected by vertex deletions
            if (hyperedge2node_delete[nb].size()!=0) {
                hyperedge_adj_adjust[e].push_back(nb);
            }

            // Neighbor hyperedge was not affected
            else {
                hyperedge_adj_old[e].push_back(nb);
            }
        }
    }
}


std::vector<std::vector<int>> build_hyperedge_adjacency(
        const std::vector<std::vector<int>> &hyperedge2node_new,
        const std::vector<std::vector<int>> &node2hyperedge)
{
    int E = hyperedge2node_new.size();
    std::vector<std::vector<int>> hyperedge_adj_new(E);

    // Visited flags for fast de-duplication
    std::vector<char> visited(E, 0);

    for (int h = 0; h < E; h++) {
        const auto &nodes = hyperedge2node_new[h];

        for (int v : nodes) {
            if (v >= node2hyperedge.size()) continue;

            // Visit all hyperedges containing node v
            for (int h2 : node2hyperedge[v]) {
                if (h2 == h) continue;        // Skip itself
                if (!visited[h2]) {          // De-dup
                    visited[h2] = 1;
                    hyperedge_adj_new[h].push_back(h2);
                }
            }
        }

        // Clear visited flags
        for (int h2 : hyperedge_adj_new[h]) {
            visited[h2] = 0;
        }
    }

    return hyperedge_adj_new;
}



void buildnewHyperedgeAdj(
        const vector<vector<int>> &hyperedge2node_new,
        const vector<vector<int>> &node2hyperedge,
        vector<vector<int>> &hyperedge_adj_update)
{
    int E = hyperedge2node_new.size();
    hyperedge_adj_update.assign(E, {});

    for (int e = 0; e < E; ++e) {

        vector<int> &adj = hyperedge_adj_update[e];

        // For each node v in hyperedge e
        for (int v : hyperedge2node_new[e]) {

            // Hyperedges incident to v (sorted)
            const vector<int> &list = node2hyperedge[v];

            // Add all candidate neighboring hyperedges
            adj.insert(adj.end(), list.begin(), list.end());
        }

        // De-dup + sort
        sort(adj.begin(), adj.end());
        adj.erase(unique(adj.begin(), adj.end()), adj.end());

        // Remove itself (e)
        adj.erase(remove(adj.begin(), adj.end(), e), adj.end());
    }
}

void getDeletedAdj(
        const vector<vector<int>> &hyperedge_adj,       // Old adjacency (sorted)
        const vector<vector<int>> &hyperedge_adj_new,   // New adjacency (sorted)
        vector<vector<int>> &hyperedge_adj_delete)      // Output: removed neighbors
{
    int E = hyperedge_adj.size();
    hyperedge_adj_delete.assign(E, {});

    for (int e = 0; e < E; ++e) {

        const vector<int> &oldList = hyperedge_adj[e];
        const vector<int> &newList = hyperedge_adj_new[e];

        vector<int> &delList = hyperedge_adj_delete[e];

        int i = 0, j = 0;

        while (i < (int)oldList.size() && j < (int)newList.size()) {

            if (oldList[i] == newList[j]) {
                // Present in new -> keep
                i++;
                j++;
            }
            else if (oldList[i] < newList[j]) {
                // In old but not in new -> removed
                delList.push_back(oldList[i]);
                i++;
            }
            else { // oldList[i] > newList[j]
                j++;  // Extra entries in newList; advance
            }
        }

        // Remaining items in oldList are removed
        while (i < (int)oldList.size()) {
            delList.push_back(oldList[i]);
            i++;
        }
    }
}

void remove_existing_adjacency(
    std::vector<std::vector<int>> &hyperedge_adj_new,
    const std::vector<std::vector<int>> &hyperedge_adj)
{
    int E = hyperedge_adj_new.size();

    for (int h = 0; h < E; ++h) {

        const std::vector<int> &old_list = hyperedge_adj[h];
        const std::vector<int> &new_list = hyperedge_adj_new[h];

        std::vector<int> filtered;
        filtered.reserve(new_list.size());

        int i = 0, j = 0;
        while (i < (int)new_list.size() && j < (int)old_list.size()) {
            if (new_list[i] < old_list[j]) {
                filtered.push_back(new_list[i]);
                ++i;
            } else if (new_list[i] > old_list[j]) {
                ++j;
            } else {
                // Equal -> already exists -> skip new_list[i]
                ++i;
                ++j;
            }
        }

        // Copy remaining new_list entries (not in old_list)
        while (i < (int)new_list.size()) {
            filtered.push_back(new_list[i]);
            ++i;
        }

        hyperedge_adj_new[h] = std::move(filtered);
    }
}


vector<int> concat_vectors(const vector<int>& a, const vector<int>& b)
{
    vector<int> result;
    result.reserve(a.size() + b.size()); // Reserve capacity for efficiency

    result.insert(result.end(), a.begin(), a.end());
    result.insert(result.end(), b.begin(), b.end());

    return result;
}




std::vector<std::vector<int>> ReadFileToVector(const std::string &filename)
{
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return {};
    }

    std::vector<std::vector<int>> result;
    std::string line;

    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        std::vector<int> nums;
        std::stringstream ss(line);
        std::string tmp;

        while (std::getline(ss, tmp, ',')) {
            if (!tmp.empty()) {
                nums.push_back(std::stoi(tmp));
            }
        }

        result.push_back(nums);
    }

    fin.close();
    return result;
}


std::vector<int> ReadTriangleFile(const std::string &filename)
{
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Unable to open file: " << filename << std::endl;
        return {};
    }

    std::vector<int> h_triangle;
    std::string line;

    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        // Find colon position
        size_t pos = line.find(':');
        if (pos == std::string::npos) continue;

        // Number after colon, e.g. " 2456"
        std::string num_str = line.substr(pos + 1);

        // Trim whitespace
        num_str.erase(0, num_str.find_first_not_of(" \t"));
        num_str.erase(num_str.find_last_not_of(" \t") + 1);

        // Convert to int
        int value = std::stoi(num_str);

        h_triangle.push_back(value);
    }

    fin.close();
    return h_triangle;
}

struct Regions {
    int x_only = 0;
    int y_only = 0;
    int z_only = 0;
    int xy_only = 0;
    int xz_only = 0;
    int yz_only = 0;
    int xyz = 0;
};



int DistinguishAAA(Regions r){

    int part1=r.x_only;
    int part2=r.y_only;
    int part3=r.z_only;
    int part12=r.xy_only;
    int part13=r.xz_only;
    int part23=r.yz_only;
    int part123=r.xyz;

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

	if(part1!=0 && part2!=0 && part3!=0){
		if(part123==0){
			return 20;
		}
		if(count==0){
			return 9;
		}else if(count==1){
			return 10;
		}else if(count==2){
			return 12;
		}else if(count==3){
			return 16;
		}
	}else if((part1==0 && part2!=0 && part3!=0) || (part1!=0 && part2==0 && part3!=0) ||(part1!=0 && part2!=0 && part3==0) ){//11 15 25
		if(part123==0){
			return 19;
		}
		if(count==2){
			return 11;
		}else if(count==3){
			return 15;
		}
	}else if(part1==0 && part2==0 && part3==0){
		if(part123==0){
			return 17;
		}else{
			return 13;
		}
	}else{
		if(part123==0){
			return 18;
		}else{
			return 14;
		}
	}
	return -1;
}


// ====================================================
// compute_regions: derive 7 regions from sorted x,y,z
// ====================================================
Regions compute_regions(const vector<int>& x,
                        const vector<int>& y,
                        const vector<int>& z)
{
    size_t i = 0, j = 0, k = 0;
    Regions r;

    while (i < x.size() || j < y.size() || k < z.size())
    {
        int vx = (i < x.size() ? x[i] : INT_MAX);
        int vy = (j < y.size() ? y[j] : INT_MAX);
        int vz = (k < z.size() ? z[k] : INT_MAX);

        int mn = min(vx, min(vy, vz));

        bool inX = (i < x.size() && x[i] == mn);
        bool inY = (j < y.size() && y[j] == mn);
        bool inZ = (k < z.size() && z[k] == mn);

        if (inX && inY && inZ) r.xyz++;
        else if (inX && inY)   r.xy_only++;
        else if (inX && inZ)   r.xz_only++;
        else if (inY && inZ)   r.yz_only++;
        else if (inX)          r.x_only++;
        else if (inY)          r.y_only++;
        else if (inZ)          r.z_only++;

        if (inX) i++;
        if (inY) j++;
        if (inZ) k++;
    }

    //std::cout << "x_only: " << r.x_only << " y_only: " << r.y_only << " z_only: " << r.z_only << " xy_only: " << r.xy_only << " xz_only: " << r.xz_only << " yz_only: " << r.yz_only << " xyz: " << r.xyz << std::endl;

    return r;
}

// ====================================================
// Classification result structure
// ====================================================
struct XYZResult {
    int a;               // Class 1~4
    vector<int> x, y, z; // x,y,z ordered per rule
    Regions reg;         // Final 7 regions (computed once)
    int check;
    int check2;
    int type;
    vector<int> empty_sig={0,0,0,0};
};

// ====================================================
// Determine containment / intersection using the 7 regions
// ====================================================
bool contain_xy(const Regions& r){ return (r.y_only==0 && r.yz_only==0); }
bool contain_yx(const Regions& r){ return (r.x_only==0 && r.xz_only==0); }

bool contain_xz(const Regions& r){ return (r.z_only==0 && r.yz_only==0 ); }
bool contain_zx(const Regions& r){ return (r.x_only==0 && r.xy_only==0 ); }

bool contain_yz(const Regions& r){ return (r.z_only==0 && r.xz_only==0 ); }
bool contain_zy(const Regions& r){ return (r.y_only==0 && r.xy_only==0 ); }


// ====================================================
// classify_xyz_regions: derive a,x,y,z + reg
// ====================================================
XYZResult classify_xyz_regions(const vector<int>& A,
                               const vector<int>& B,
                               const vector<int>& C)
{
    XYZResult res;

    vector<vector<int>> E = {A,B,C};
    vector<int> idx = {0,1,2};

    
    const vector<int>& x = E[idx[0]];
    const vector<int>& y = E[idx[1]];
    const vector<int>& z = E[idx[2]];

    Regions r = compute_regions(x,y,z);

    bool c_xy = contain_xy(r);
    bool c_yx = contain_yx(r);
    bool c_xz = contain_xz(r);
    bool c_zx = contain_zx(r);
    bool c_yz = contain_yz(r);
    bool c_zy = contain_zy(r);

    res.empty_sig={r.x_only,r.y_only,r.z_only};

    if((c_xy && c_yx) || (c_xz && c_zx) || (c_yz && c_zy)){
        res.a = -1;
        return res;
    }

    

    int cc = (c_xy||c_yx) + (c_xz||c_zx) + (c_yz||c_zy);

    // Case 1: three containments
    if (cc == 3) {
        res.a = 1;
        res.x = x; res.y = y; res.z = z;
        res.reg = r;
        //res.empty_sig={r.x_only,r.y_only,r.z_only};
        //return res;
    }

    // Case 2: two containments + one intersection (z is contained)
    if (cc == 2 ) {
        if(c_xy && c_xz){
            res.a = 3;
            res.x = y; res.y = z; res.z = x;
            res.reg = r;
            res.check=r.x_only;
            //res.empty_sig={r.y_only,r.z_only,r.x_only};
            //return res;
        }else if(c_yx && c_yz){
            res.a = 3;
            res.x = z; res.y = x; res.z = y;
            res.reg = r;
            res.check=r.y_only;
            //res.empty_sig={r.z_only,r.x_only,r.y_only};
            //return res;
        }else if(c_zx && c_zy){
            res.a = 3;
            res.x = x; res.y = y; res.z = z;
            res.reg = r; 
            res.check=r.z_only;
            //res.empty_sig={r.x_only,r.y_only,r.z_only};
            //return res; 
        }else if(c_xz && c_yz){
            res.a = 2;
            res.x = x; res.y = y; res.z = z;
            res.reg = r;
            res.check=r.xy_only;
            //res.empty_sig={r.x_only,r.y_only,r.z_only};
            //return res;
        }else if(c_zx && c_yx){
            res.a = 2;
            res.x = z; res.y = y; res.z = x;
            res.reg = r;
            res.check=r.yz_only;
            //res.empty_sig={r.z_only,r.y_only,r.x_only};
            //return res;
        }else if(c_xy && c_zy){
            res.a = 2;
            res.x = x; res.y = z; res.z = y;
            res.reg = r;
            res.check=r.xz_only;
            //res.empty_sig={r.x_only,r.z_only,r.y_only};
            //return res;
        }
    }

    // Case 4: one containment + two intersections
    if (cc == 1 ) {
        if (c_xz) {     // y ⊇ z
            res.a = 4;
            res.x = y; res.y = x; res.z = z;
            res.reg = r;
            res.check=r.xy_only;
            res.check2=r.x_only;
            //res.empty_sig={r.y_only,r.x_only,r.z_only}; 
            //return res;
        }else if(c_yz){
            res.a = 4;
            res.x = x; res.y = y; res.z = z;
            res.reg = r;
            res.check=r.xy_only;
            res.check2=r.y_only;
            //res.empty_sig={r.x_only,r.y_only,r.z_only};
            //return res;
        }else if (c_xy){
            res.a = 4;
            res.x = z; res.y = x; res.z = y;
            res.reg = r;
            res.check=r.xz_only;
            res.check2=r.x_only;
            //res.empty_sig={r.z_only,r.x_only,r.y_only};
            //return res;
        }else if (c_zy){
            res.a = 4;
            res.x = x; res.y = z; res.z = y;
            res.reg = r;
            res.check=r.xz_only;
            res.check2=r.z_only;
            //res.empty_sig={r.x_only,r.z_only,r.y_only};
            //return res;
        }else if (c_zx){
            res.a = 4;
            res.x = y; res.y = z; res.z = x;
            res.reg = r;
            res.check=r.yz_only;
            res.check2=r.z_only;
            //res.empty_sig={r.y_only,r.z_only,r.x_only};
            //return res;
        }else if (c_yx){
            res.a = 4;
            res.x = z; res.y = y; res.z = x;
            res.reg = r;
            res.check=r.yz_only;
            res.check2=r.y_only;
            //res.empty_sig={r.z_only,r.y_only,r.x_only};
            //return res;
        }
    }


    if(cc==0){
        res.a = DistinguishAAA(r);
    }


    return res;
}

// ====================================================
// Final pattern (1–8) identification function
// ====================================================
vector<int> get_pattern(const vector<int>& A,
                const vector<int>& B,
                const vector<int>& C)
{
    XYZResult r = classify_xyz_regions(A,B,C);
    //if (r.a == -1) return -1;

    const Regions& reg = r.reg;




    // -------- Pattern 1 --------
    if (r.a == 1) return concat_vectors({1},r.empty_sig);

    // -------- Pattern 2 / 3 --------
    if (r.a == 2) {
        if (r.check == 0) return concat_vectors({2},r.empty_sig);
        else return concat_vectors({3},r.empty_sig);
    }

    // -------- Pattern 4 / 5 --------
    if (r.a == 3) {
        if (r.check == 0) return concat_vectors({4},r.empty_sig);
        else return concat_vectors({5},r.empty_sig);
    }

    // -------- Pattern 6 / 7 / 8 --------
    if (r.a == 4) {
        if (r.check == 0) return concat_vectors({6},r.empty_sig);
        else {
            if (r.check2 == 0) return concat_vectors({7},r.empty_sig);
            else return concat_vectors({8},r.empty_sig);
        }
    }
    //std::cout << "/* message */" << std::endl;

    return concat_vectors({r.a},r.empty_sig);
}


vector<int> DistinguishAAA2(Regions r){

    int part1=r.x_only;
    int part2=r.y_only;
    int part3=r.z_only;
    int part12=r.xy_only;
    int part13=r.xz_only;
    int part23=r.yz_only;
    int part123=r.xyz;

	int count=0;
    int pat10=-1;
    int pat12=-1;
	if(part12!=0){
		count+=1;
        pat10=3;
	}else{
        pat12=3;
    }
	if(part23!=0){
		count+=1;
        pat10=1;
	}else{
        pat12=1;
    }
	if(part13!=0){
		count+=1;
        pat10=2;
	}else{
        pat12=2;
    }

	if(part1!=0 && part2!=0 && part3!=0){
		if(part123==0){
			return {20,0};
		}
		if(count==0){
			return {9,0};
		}else if(count==1){
			return {10,pat10};
		}else if(count==2){
			return {12,pat12};
		}else if(count==3){
			return {16,0};
		}
	}else if((part1==0 && part2!=0 && part3!=0) || (part1!=0 && part2==0 && part3!=0) ||(part1!=0 && part2!=0 && part3==0) ){//11 15 25
		if(part123==0){
			return {19,0};
		}
		if(count==2){
			return {11,0};
		}else if(count==3){
			return {15,0};
		}
	}else if(part1==0 && part2==0 && part3==0){
		if(part123==0){
			return {17,0};
		}else{
			return {13,0};
		}
	}else{
		if(part123==0){
			return {18,0};
		}else{
			return {14,0};
		}
	}
	return {-1,-1};
}

XYZResult classify_xyz_regions2(const vector<int>& A,
                               const vector<int>& B,
                               const vector<int>& C)
{
    XYZResult res;

    vector<vector<int>> E = {A,B,C};
    vector<int> idx = {0,1,2};

    
    const vector<int>& x = E[idx[0]];
    const vector<int>& y = E[idx[1]];
    const vector<int>& z = E[idx[2]];

    Regions r = compute_regions(x,y,z);

    bool c_xy = contain_xy(r);
    bool c_yx = contain_yx(r);
    bool c_xz = contain_xz(r);
    bool c_zx = contain_zx(r);
    bool c_yz = contain_yz(r);
    bool c_zy = contain_zy(r);

    res.empty_sig={r.x_only,r.y_only,r.z_only,-1};

    if((c_xy && c_yx) || (c_xz && c_zx) || (c_yz && c_zy)){
        res.a = -1;
        return res;
    }

    

    int cc = (c_xy||c_yx) + (c_xz||c_zx) + (c_yz||c_zy);

    // Case 1: three containments
    if (cc == 3) {
        res.a = 1;
        res.x = x; res.y = y; res.z = z;
        res.reg = r;

        if(res.empty_sig[0]>0){
            if(B.size() < C.size() ){
                res.empty_sig[1]=-1;

            }else{
                res.empty_sig[2]=-1;

            }
        }else if(res.empty_sig[1]>0){
            if(A.size() < C.size() ){
                res.empty_sig[0]=-1;

            }else{
                res.empty_sig[2]=-1;

            }
        }else{
            if(A.size() < B.size() ){
                res.empty_sig[0]=-1;

            }else{
                res.empty_sig[1]=-1;
            }
        }
        //res.empty_sig={r.x_only,r.y_only,r.z_only};
        //return res;
    }

    // Case 2: two containments + one intersection (z is contained)
    if (cc == 2 ) {
        if(c_xy && c_xz){
            res.a = 3;
            res.x = y; res.y = z; res.z = x;
            res.reg = r;
            res.check=r.x_only;
            if(r.x_only==0){
                res.empty_sig[0]=-1;
            }
            //res.empty_sig={r.y_only,r.z_only,r.x_only};
            //return res;
        }else if(c_yx && c_yz){
            res.a = 3;
            res.x = z; res.y = x; res.z = y;
            res.reg = r;
            res.check=r.y_only;
            if(r.y_only==0){
                res.empty_sig[1]=-1;
            }
            //res.empty_sig={r.z_only,r.x_only,r.y_only};
            //return res;
        }else if(c_zx && c_zy){
            res.a = 3;
            res.x = x; res.y = y; res.z = z;
            res.reg = r; 
            res.check=r.z_only;
            if(r.z_only==0){
                res.empty_sig[2]=-1;
            }
            //res.empty_sig={r.x_only,r.y_only,r.z_only};
            //return res; 
        }else if(c_xz && c_yz){
            res.a = 2;
            res.x = x; res.y = y; res.z = z;
            res.reg = r;
            res.check=r.xy_only;
            //res.empty_sig={r.x_only,r.y_only,r.z_only};
            //return res;
        }else if(c_zx && c_yx){
            res.a = 2;
            res.x = z; res.y = y; res.z = x;
            res.reg = r;
            res.check=r.yz_only;
            //res.empty_sig={r.z_only,r.y_only,r.x_only};
            //return res;
        }else if(c_xy && c_zy){
            res.a = 2;
            res.x = x; res.y = z; res.z = y;
            res.reg = r;
            res.check=r.xz_only;
            //res.empty_sig={r.x_only,r.z_only,r.y_only};
            //return res;
        }
    }

    // Case 4: one containment + two intersections
    if (cc == 1 ) {
        if (c_xz) {     // y ⊇ z
            res.a = 4;
            res.x = y; res.y = x; res.z = z;
            res.reg = r;
            res.check=r.xy_only;
            res.check2=r.x_only;
            if(r.x_only==0){
                res.empty_sig[0]=-1;
            }
            if(r.xy_only==0){
                res.empty_sig[3]=2;
            }else{
                res.empty_sig[3]=1;
            }
            //res.empty_sig={r.y_only,r.x_only,r.z_only}; 
            //return res;
        }else if(c_yz){
            res.a = 4;
            res.x = x; res.y = y; res.z = z;
            res.reg = r;
            res.check=r.xy_only;
            res.check2=r.y_only;
            if(r.y_only==0){
                res.empty_sig[1]=-1;
            }
            if(r.xy_only==0){
                res.empty_sig[3]=1;
            }else{
                res.empty_sig[3]=2;
            }
            //res.empty_sig={r.x_only,r.y_only,r.z_only};
            //return res;
        }else if (c_xy){
            res.a = 4;
            res.x = z; res.y = x; res.z = y;
            res.reg = r;
            res.check=r.xz_only;
            res.check2=r.x_only;
            if(r.x_only==0){
                res.empty_sig[0]=-1;
            }
            if(r.xz_only==0){
                res.empty_sig[3]=3;
            }else{
                res.empty_sig[3]=1;
            }
            //res.empty_sig={r.z_only,r.x_only,r.y_only};
            //return res;
        }else if (c_zy){
            res.a = 4;
            res.x = x; res.y = z; res.z = y;
            res.reg = r;
            res.check=r.xz_only;
            res.check2=r.z_only;
            if(r.z_only==0){
                res.empty_sig[2]=-1;
            }
            if(r.xz_only==0){
                res.empty_sig[3]=1;
            }else{
                res.empty_sig[3]=3;
            }
            //res.empty_sig={r.x_only,r.z_only,r.y_only};
            //return res;
        }else if (c_zx){
            res.a = 4;
            res.x = y; res.y = z; res.z = x;
            res.reg = r;
            res.check=r.yz_only;
            res.check2=r.z_only;
            if(r.z_only==0){
                res.empty_sig[2]=-1;
            }
            if(r.yz_only==0){
                res.empty_sig[3]=2;
            }else{
                res.empty_sig[3]=3;
            }
            //res.empty_sig={r.y_only,r.z_only,r.x_only};
            //return res;
        }else if (c_yx){
            res.a = 4;
            res.x = z; res.y = y; res.z = x;
            res.reg = r;
            res.check=r.yz_only;
            res.check2=r.y_only;
            if(r.y_only==0){
                res.empty_sig[1]=-1;
            }
            if(r.yz_only==0){
                res.empty_sig[3]=3;
            }else{
                res.empty_sig[3]=2;
            }
            //res.empty_sig={r.z_only,r.y_only,r.x_only};
            //return res;
        }
    }


    if(cc==0){
        vector<int> result=DistinguishAAA2(r);
        res.a = result[0];
        if(res.a==10 || res.a==12){
            res.empty_sig[3]=result[1];
        }
    }


    return res;
}

vector<int> get_specificPattern(const vector<int>& A,
                const vector<int>& B,
                const vector<int>& C)
{
    XYZResult r = classify_xyz_regions2(A,B,C);
    //if (r.a == -1) return -1;

    const Regions& reg = r.reg;




    // -------- Pattern 1 --------
    if (r.a == 1){
        return concat_vectors({1},r.empty_sig);
    } 

    // -------- Pattern 2 / 3 --------
    if (r.a == 2) {
        if (r.check == 0) return concat_vectors({2},r.empty_sig);
        else return concat_vectors({3},r.empty_sig);
    }

    // -------- Pattern 4 / 5 --------
    if (r.a == 3) {
        if (r.check == 0){
            return concat_vectors({4},r.empty_sig);
        } else{
            return concat_vectors({5},r.empty_sig);
        } 
    }

    // -------- Pattern 6 / 7 / 8 --------
    if (r.a == 4) {
        if (r.check == 0) return concat_vectors({6},r.empty_sig);
        else {
            if (r.check2 == 0) return concat_vectors({7},r.empty_sig);
            else return concat_vectors({8},r.empty_sig);
        }
    }
    //std::cout << "/* message */" << std::endl;

    return concat_vectors({r.a},r.empty_sig);
}


vector<vector<int>> build_same_node_adj_direct(
    const vector<vector<int>>& hyperedge2node,
    const set<int>& Updated_hyperedge_set)   // <-- add this parameter
{
    int E = hyperedge2node.size();
    vector<vector<int>> hyperedge_adj_same(E);

    for (int i = 0; i < E; i++) {
        for (int j = i + 1; j < E; j++) {

            if (hyperedge2node[i].size() != hyperedge2node[j].size())
                continue;

            bool same = true;
            const vector<int>& A = hyperedge2node[i];
            const vector<int>& B = hyperedge2node[j];

            for (int t = 0; t < A.size(); t++) {
                if (A[t] != B[t]) {
                    same = false;
                    break;
                }
            }

            if (same) {

                // Filter: skip neighbors that are in Updated_hyperedge_set
                if (Updated_hyperedge_set.count(j) == 0)
                    hyperedge_adj_same[i].push_back(j);

                if (Updated_hyperedge_set.count(i) == 0)
                    hyperedge_adj_same[j].push_back(i);
            }
        }
    }

    return hyperedge_adj_same;
}



void update_pattern(int hyperedge_1, int hyperedge_2, int hyperedge_3, std::vector<std::vector<int>> &Empty_TriangleNum_record,
    vector<int> &Previous_patternNum, bool is_add){
    int sign=1;
    if(!is_add) sign=-1;
    if(Previous_patternNum[0]==1){
        if(Previous_patternNum[1]==-1){
            Empty_TriangleNum_record[hyperedge_1][20]+=sign;
        }
        if(Previous_patternNum[2]==-1){
            Empty_TriangleNum_record[hyperedge_2][20]+=sign;
        }
        if(Previous_patternNum[3]==-1){
            Empty_TriangleNum_record[hyperedge_3][20]+=sign;
        }
    }else if(Previous_patternNum[0]==4){
        if(Previous_patternNum[1]==-1){
            Empty_TriangleNum_record[hyperedge_1][21]+=sign;
        }
        if(Previous_patternNum[2]==-1){
            Empty_TriangleNum_record[hyperedge_2][21]+=sign;
        }
        if(Previous_patternNum[3]==-1){
            Empty_TriangleNum_record[hyperedge_3][21]+=sign;
        }
    }else if(Previous_patternNum[0]==7){
        if(Previous_patternNum[1]==-1){
            Empty_TriangleNum_record[hyperedge_1][22]+=sign;
        }
        if(Previous_patternNum[2]==-1){
            Empty_TriangleNum_record[hyperedge_2][22]+=sign;
        }
        if(Previous_patternNum[3]==-1){
            Empty_TriangleNum_record[hyperedge_3][22]+=sign;
        }
    }

    if(Previous_patternNum[1]<=0){
        Empty_TriangleNum_record[hyperedge_1][Previous_patternNum[0]-1]+=sign;
    }
    if(Previous_patternNum[2]<=0){
        Empty_TriangleNum_record[hyperedge_2][Previous_patternNum[0]-1]+=sign;
    }
    if(Previous_patternNum[3]<=0){
        Empty_TriangleNum_record[hyperedge_3][Previous_patternNum[0]-1]+=sign;
    }
}


void update_twice_pattern(int hyperedge_1, int hyperedge_2, int hyperedge_3, std::vector<std::vector<int>> &Empty_TriangleNum_record, std::vector<std::vector<int>> &Empty_TriangleNum_record_copy,
    vector<int> &Previous_patternNum, bool is_add){
    int sign=1;
    if(!is_add) sign=-1;
    if(Previous_patternNum[0]==1){
        if(Previous_patternNum[1]==-1){
            Empty_TriangleNum_record[hyperedge_1][20]+=sign;
            Empty_TriangleNum_record_copy[hyperedge_1][20]+=sign;
        }
        if(Previous_patternNum[2]==-1){
            Empty_TriangleNum_record[hyperedge_2][20]+=sign;
            Empty_TriangleNum_record_copy[hyperedge_2][20]+=sign;
        }
        if(Previous_patternNum[3]==-1){
            Empty_TriangleNum_record[hyperedge_3][20]+=sign;
            Empty_TriangleNum_record_copy[hyperedge_3][20]+=sign;
        }
    }else if(Previous_patternNum[0]==4){
        if(Previous_patternNum[1]==-1){
            Empty_TriangleNum_record[hyperedge_1][21]+=sign;
            Empty_TriangleNum_record_copy[hyperedge_1][21]+=sign;
        }
        if(Previous_patternNum[2]==-1){
            Empty_TriangleNum_record[hyperedge_2][21]+=sign;
            Empty_TriangleNum_record_copy[hyperedge_2][21]+=sign;
        }
        if(Previous_patternNum[3]==-1){
            Empty_TriangleNum_record[hyperedge_3][21]+=sign;
            Empty_TriangleNum_record_copy[hyperedge_3][21]+=sign;
        }
    }else if(Previous_patternNum[0]==7){
        if(Previous_patternNum[1]==-1){
            Empty_TriangleNum_record[hyperedge_1][22]+=sign;
            Empty_TriangleNum_record_copy[hyperedge_1][22]+=sign;
        }
        if(Previous_patternNum[2]==-1){
            Empty_TriangleNum_record[hyperedge_2][22]+=sign;
            Empty_TriangleNum_record_copy[hyperedge_2][22]+=sign;
        }
        if(Previous_patternNum[3]==-1){
            Empty_TriangleNum_record[hyperedge_3][22]+=sign;
            Empty_TriangleNum_record_copy[hyperedge_3][22]+=sign;
        }
    }

    if(Previous_patternNum[1]<=0){
        Empty_TriangleNum_record[hyperedge_1][Previous_patternNum[0]-1]+=sign;
        Empty_TriangleNum_record_copy[hyperedge_1][Previous_patternNum[0]-1]+=sign;
    }
    if(Previous_patternNum[2]<=0){
        Empty_TriangleNum_record[hyperedge_2][Previous_patternNum[0]-1]+=sign;
        Empty_TriangleNum_record_copy[hyperedge_2][Previous_patternNum[0]-1]+=sign;
    }
    if(Previous_patternNum[3]<=0){
        Empty_TriangleNum_record[hyperedge_3][Previous_patternNum[0]-1]+=sign;
        Empty_TriangleNum_record_copy[hyperedge_3][Previous_patternNum[0]-1]+=sign;
    }
}



void TransformAndWriteSpecificPattern(
    const std::vector<std::vector<int>> &TriangleNum_record,
    const std::string &filename)
{
    // First try opening the file
    std::ofstream fout(filename);

    // If opening fails, try creating an empty file
    if (!fout.is_open()) {
        std::ofstream create_file(filename);
        if (create_file.is_open()) {
            create_file.close();
            std::cout << "File did not exist. Created an empty file: " << filename << std::endl;
        } else {
            std::cerr << "Failed to create file: " << filename << std::endl;
            return;
        }

        // Try opening again
        fout.open(filename);
        if (!fout.is_open()) {
            std::cerr << "Failed to open file after creation: " << filename << std::endl;
            return;
        }
    }

    int R = TriangleNum_record.size();

    for (int i = 0; i < R; ++i) {
        const auto &h_triangle_record = TriangleNum_record[i];

        if (h_triangle_record.size() < 23) {
            std::cerr << "Row " << i << " has less than 23 elements, skipped." << std::endl;
            continue;
        }

        std::vector<int> out(23);

        for (int x = 0; x < 23; x++) {
            out[x] = h_triangle_record[x];
        }

        for (int j = 0; j < 23; ++j) {
            if (j > 0) fout << ",";
            fout << out[j];
        }
        fout << "\n";
    }

    fout.close();
    std::cout << "Transformation completed and written to file: " << filename << std::endl;
}

void write_set_to_file(const vector<int> &adjust_hyperedge, const string &path) {
    ofstream fout(path);
    if (!fout) {
        cout << "Unable to open file: " << path << endl;
        return;
    }

    for (int x=0; x<adjust_hyperedge.size(); x++) {
        if(adjust_hyperedge[x]==1) fout << x << "\n";   // One number per line
    }

    fout.close();
}


void load_adjust_hyperedge(vector<int> &adjust_hyperedge, const string &path) {
    ifstream fin(path);
    if (!fin.is_open()) {
        // If the file does not exist, create an empty file
        ofstream fout(path);
        if (fout.is_open()) {
            fout.close();
        }
        return;   // Do not modify adjust_hyperedge
    }

    int x;
    while (fin >> x) {
        if (x >= 0 && x < (int)adjust_hyperedge.size()) {
            adjust_hyperedge[x] = 1;   // Mark positions that appear
        }
    }

    fin.close();
}




void write_hyperedge2node_to_file(const vector<vector<int>> &hyperedge2node,
                                  const string &graphFile) {
    ofstream fout(graphFile);
    if (!fout.is_open()) {
        cout << "Unable to open file: " << graphFile << endl;
        return;
    }

    for (const auto &row : hyperedge2node) {
        for (int i = 0; i < (int)row.size(); i++) {
            fout << row[i];
            if (i + 1 < (int)row.size()) fout << ",";  // Add comma
        }
        fout << "\n";
    }

    fout.close();
}


void write_h_triangle_to_file(const string &filename, const vector<long long> &h_triangle_final) {
    ofstream fout(filename);
    if (!fout.is_open()) {
        cerr << "Unable to open file!" << endl;
        return;
    }

    long long triangle = 0;
    int index = 0;

    for (int i = 1; i <= (int)h_triangle_final.size(); i++) {

        fout << "h-triangle " << ++index << ": "
             << h_triangle_final[i - 1] << endl;

        triangle += h_triangle_final[i - 1];
    }

    fout.close();
}


void write_h_triangle_to_file_v2(const string &filename, const vector<int> &h_triangle_final) {
    ofstream fout(filename);
    if (!fout.is_open()) {
        cerr << "Unable to open file!" << endl;
        return;
    }

    long long triangle = 0;
    int index = 0;

    for (int i = 1; i <= (int)h_triangle_final.size(); i++) {

        fout << "h-triangle " << ++index << ": "
             << h_triangle_final[i - 1] << endl;

        triangle += h_triangle_final[i - 1];
    }

    fout.close();
}

void copy_file_content(const string &source_path, const string &dest_path) {
    ifstream fin(source_path);
    if (!fin.is_open()) {
        cerr << "Unable to open source file: " << source_path << endl;
        return;
    }

    ofstream fout(dest_path);
    if (!fout.is_open()) {
        cerr << "Unable to open destination file: " << dest_path << endl;
        fin.close();
        return;
    }

    string line;
    while (getline(fin, line)) {
        fout << line << "\n";
    }

    fin.close();
    fout.close();
}

void write_h_triangle(const vector<int>& h_triangle, const string& filepath) {
    ofstream out(filepath);
    
    if (!out.is_open()) {
        cerr << "Error: cannot open file " << filepath << endl;
        return;
    }

    for (size_t i = 0; i < h_triangle.size(); ++i) {
        out << "h-triangle " << (i + 1) << ": " << h_triangle[i] << '\n';
    }

    out.close();
}

void write_file_with_content(const string &dest_path, const string &source_path) {
    ifstream fin(source_path);
    if (!fin.is_open()) {
        cerr << "Error: Cannot open source file: " << source_path << endl;
        return;
    }

    ofstream fout(dest_path);
    if (!fout.is_open()) {
        cerr << "Error: Cannot open destination file: " << dest_path << endl;
        fin.close();
        return;
    }

    string line;
    while (getline(fin, line)) {
        fout << line << "\n";
    }

    fin.close();
    fout.close();
}


void UpdateHTriangle(
    std::vector<int> &h_triangle,
    const std::vector<std::vector<int>> &Empty_TriangleNum_record,
    int i)
{
    // Basic safety check (prevent out-of-range access)
    if (i >= Empty_TriangleNum_record.size()) {
        std::cerr << "Index i out of range in Empty_TriangleNum_record\n";
        return;
    }

    const auto &rec = Empty_TriangleNum_record[i];

    if (rec.size() < 23 || h_triangle.size() < 20) {
        std::cerr << "Input size is insufficient.\n";
        return;
    }

    // === Original logic ===
    h_triangle[2] += rec[0] - rec[20];
    h_triangle[5] += rec[20];

    h_triangle[4] += rec[21];
    h_triangle[6] += rec[3] - rec[21];

    h_triangle[7] += rec[22];
    h_triangle[10] += rec[6] - rec[22];

    h_triangle[8] += rec[1];
    h_triangle[9] += rec[2];

    h_triangle[7] += rec[4];
    h_triangle[9] += rec[5];

    h_triangle[11] += rec[7];
    h_triangle[11] += rec[10];

    h_triangle[13] += rec[12];
    h_triangle[14] += rec[13];
    h_triangle[15] += rec[14];

    h_triangle[17] += rec[16];
    h_triangle[18] += rec[17];
    h_triangle[19] += rec[18];
}






void Load_New_Hyperedges(
    const string &updateVertex,
    vector<vector<int>> &hyperedge2node_new
) {

    ifstream fin(updateVertex);

    if (!fin.is_open()) {

        cerr << "Unable to open file: "
             << updateVertex << endl;

        return;
    }

    string line;

    while (getline(fin, line)) {

        if (line.empty())
            continue;

        stringstream ss(line);

        int hid;

        char colon;

        // Read:
        // 1562:
        ss >> hid >> colon;

        vector<int> vertices;

        int v;

        while (ss >> v) {

            vertices.push_back(v);
        }

        // =====================================
        // Sort ascending
        // =====================================
        sort(
            vertices.begin(),
            vertices.end()
        );

        // Ensure the vector is large enough
        if ((int)hyperedge2node_new.size() <= hid) {

            hyperedge2node_new.resize(hid + 1);
        }

        hyperedge2node_new[hid] = vertices;
    }

    fin.close();
}
