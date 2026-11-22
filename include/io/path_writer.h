#ifndef path_writer_h
#define path_writer_h

#include <iostream>
#include <fstream>
#include <string>
#include <sdbg/sdbg.h>
#include <vector>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <chrono>
#include <iomanip>
using namespace std;

class PathWriter{
    #ifdef DEBUG
    private:
        SDBG& sdbg;
        string genome_id;
    public:
        PathWriter(string mode, SDBG& sdbg, const vector<uint64_t>& path, string genome_id, string type);
        string FetchNodeLabel(size_t node);
        void CollectPathsIntoStringStreams(const vector<uint64_t>& path);
    #endif
};

#endif
