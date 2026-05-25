


#include <iostream>
#include <fstream>
#include <unordered_set>
#include <algorithm>
using namespace std;

void read_data(string path, vector< vector<int> >& node2hyperedge,
               vector< vector<int> >& hyperedge2node,
               vector< unordered_set<int> >& hyperedge2node_set)
{
    ifstream graphFile(path.c_str());
    string line;
    int num_hyperedges = 0;

    while (getline(graphFile, line))
    {
        unordered_set<int> tokens_set;
        vector<int> tokens;

        // --------------------- 新增：空行处理 ---------------------
        if (line.empty()) {
            // 空 hyperedge：vector 和 set 都是空的
            hyperedge2node.push_back(tokens);          // 空 vector
            hyperedge2node_set.push_back(tokens_set);  // 空 set
            num_hyperedges++;
            continue;   // 不进行后续解析
        }
        // --------------------------------------------------------

        bool EOL = false;
        int idx;

        while (!EOL)
        {
            int pos = line.find(",");
            if (pos == string::npos) {
                pos = line.size();
                EOL = true;
                idx = stoi(line);
            } else {
                idx = stoi(line.substr(0, pos));
                line.erase(0, pos + 1);
            }

            while (idx >= (int)node2hyperedge.size()) {
                node2hyperedge.push_back(vector<int>());
            }

            if (node2hyperedge[idx].empty() ||
                node2hyperedge[idx].back() != num_hyperedges)
            {
                node2hyperedge[idx].push_back(num_hyperedges);
                tokens.push_back(idx);
                tokens_set.insert(idx);
            }
        }

        sort(tokens.begin(), tokens.end());
        hyperedge2node.push_back(tokens);
        hyperedge2node_set.push_back(tokens_set);
        num_hyperedges++;
    }

    graphFile.close();
}
