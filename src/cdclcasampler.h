#pragma once

#include "../minisat_ext/BlackBoxSolver.h"
#include "../minisat_ext/Ext.h"
#include <vector>
#include <queue>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <unistd.h>
#include <sstream>
#include <cassert>
#include <cstring>
#include <chrono>
#include <filesystem>
#include <set>

using namespace std;

struct hash_pair {
    size_t operator()(const pair<int, int>& p) const {
        return hash<int>()(p.first) ^ hash<int>()(p.second);
    }
};

class DSU {

private:
    unordered_map<int, int> parent;
public:
    int find(int x) {
        if (parent.find(x) == parent.end()) {
            parent[x] = x; 
        }
        if (parent[x] != x) {
            parent[x] = find(parent[x]); 
        }
        return parent[x];
    }

    void unionSets(int x, int y) {
        int rootX = find(x);
        int rootY = find(y);
        if (rootX != rootY) {
            parent[rootX] = rootY;
        }
    }
};

class SimpleBitSet {
public:
    SimpleBitSet(): siz(0), cap(0), cnt(0ll){
        arr = nullptr;
    }
    SimpleBitSet(size_t _siz): siz(_siz), cnt(0ll){
        cap = (_siz + 31) >> 5;
        if (cap == 0) cap = 1;
        arr = new uint32_t[cap]{0};
    }
    SimpleBitSet(const SimpleBitSet& bs): siz(bs.siz), cap(bs.cap), cnt(bs.cnt){
        arr = new uint32_t[bs.cap]{0};
        memcpy(arr, bs.arr, sizeof(uint32_t) * cap);
    }
    SimpleBitSet(SimpleBitSet&& bs): siz(bs.siz), cap(bs.cap), cnt(bs.cnt){
        arr = bs.arr;
        bs.arr = nullptr;
        bs.siz = 0, bs.cap = 0, bs.cnt = 0;
    }
    ~SimpleBitSet(){ delete [] arr; }
    SimpleBitSet& operator=(const SimpleBitSet& bs){
        if (this == &bs) return *this;
        delete [] arr;
        siz = bs.siz;
        cap = bs.cap;
        cnt = bs.cnt;
        arr = new uint32_t[cap]{0};
        memcpy(arr, bs.arr, sizeof(uint32_t) * cap);
        return *this;
    }
    SimpleBitSet& operator=(SimpleBitSet&& bs){
        if (this == &bs) return *this;
        delete [] arr;
        siz = bs.siz;
        cap = bs.cap;
        cnt = bs.cnt;
        arr = bs.arr;
        bs.arr = nullptr;
        bs.siz = 0, bs.cap = 0, bs.cnt = 0;
        return *this;
    }
    bool set(size_t idx){
        uint32_t& t1 = arr[idx >> 5];
        uint32_t cg = (1u << (idx & 31));
        if ((t1 & cg) == 0u){
            ++cnt;
            t1 |= cg;
            return true;
        }
        return false;
    }
    bool unset(size_t idx){
        uint32_t& t1 = arr[idx >> 5];
        uint32_t cg = (1u << (idx & 31));
        if ((t1 & cg) != 0u){
            --cnt;
            t1 ^= cg;
            return true;
        }
        return false;
    }
    bool get(size_t idx) const {
        return (arr[idx >> 5] & (1 << (idx & 31))) != 0u;
    }
    long long count() const {
        return cnt;
    }
    void difference_(const SimpleBitSet& bs){
        for (size_t i = 0; i < cap; ++i){
            uint32_t tmp = arr[i] & bs.arr[i];
            if (tmp != 0u){
                arr[i] ^= tmp;
                cnt -= __builtin_popcount(tmp);
            }
        }
    }
    bool init_iter(){
        iter_curpos = iter_curbit = 0;
        while (iter_curpos < cap && arr[iter_curpos] == 0u){
            ++iter_curpos;
        }
        if (iter_curpos < cap){
            iter_curbit = __builtin_ctz(arr[iter_curpos]);
            iter_curpos_msb = 31 - __builtin_clz(arr[iter_curpos]);
            return true;
        }
        return false;
    }
    size_t iter_get() const {
        return (iter_curpos << 5) | iter_curbit;
    }
    bool iter_next(){
        if (iter_curpos == cap) return false;

        if (iter_curbit < iter_curpos_msb){
            uint32_t tmp = arr[iter_curpos];
            do { ++iter_curbit; } while (((tmp >> iter_curbit) & 1u) == 0u);
            return true;
        }
        ++iter_curpos, iter_curbit = 0u;
        while (iter_curpos < cap && arr[iter_curpos] == 0u){
            ++iter_curpos;
        }
        if (iter_curpos == cap){
            return false;
        }
        iter_curbit = __builtin_ctz(arr[iter_curpos]);
        iter_curpos_msb = 31 - __builtin_clz(arr[iter_curpos]);
        return true;
    }
    void bitwise_or(SimpleBitSet& rhs) {
        if (rhs.cap != cap) return;     
        for (size_t w = 0; w < cap; ++w) {
            uint32_t old = arr[w];
            arr[w] |= rhs.arr[w];
            cnt += __builtin_popcount(arr[w] ^ old);
        }
    }
    void bitwise_xor(const SimpleBitSet& rhs) {
        if (rhs.cap != cap) return; 
        for (size_t w = 0; w < cap; ++w) {
                uint32_t old = arr[w];   
                arr[w] ^= rhs.arr[w];    
                cnt += __builtin_popcount(arr[w]) - __builtin_popcount(old);
        }
    }
    SimpleBitSet bitwise_xor_inplace(const SimpleBitSet& rhs) const {
        if (rhs.cap != cap) {
            return SimpleBitSet();
        }
        SimpleBitSet ret(*this);    
        ret.bitwise_xor(rhs);         
        return ret;                  
    }
    void clear() {
        fill(arr, arr + cap, 0u);
        cnt = 0;
    }
private:
    uint32_t* arr;
    size_t siz;
    size_t cap;
    long long cnt;
    
    size_t iter_curpos;
    uint32_t iter_curbit;
    int iter_curpos_msb;
};

class CDCLCASampler
{
public:
    CDCLCASampler(string cnf_file_path, int seed);
    ~CDCLCASampler();

    void SetTransferSize(int transfer_set_size);
    void SetCandidateSetSize(int candidate_set_size); 
    inline void SetTestcaseSetSavePath(string testcase_set_path) { testcase_set_save_path_ = testcase_set_path; } 
    void SetDefaultPara();
    void GenerateInitTestcaseCDCL();
    void GenerateTestcase();
    void GenerateCandidateTestcaseSet();
    int SelectTestcaseFromCandidateSetByTupleNum();
    long long Get2TupleMapIndex(long i, long v_i, long j, long v_j);
    void Init2TupleInfo();
    void Update2TupleInfo();
    void InitSampleWeightByAppearance();
    void UpdateSampleWeightByAppearance();
    void ReduceCNF();
    void Init();
    void GenerateCoveringArray();
    void SaveTestcaseSet(vector<vector<int>>,string);
    void remove_temp_files();
    void find_uncovered_tuples(bool simplify);
    void get_cadical_solution(vector<int>& tc);
    void read_cnf_header(ifstream& ifs, int& nvar, int& nclauses);
    bool read_cnf();
    bool check_no_clauses();
    void clear_pq();
    void clear_final();
    void new_choose_final_plain();
    void Update2TupleInfo1();
    void Reduce_redundancy(bool last);
    void processbackbones(vector<vector<int>>& clauses, const vector<int>& backbones);
    bool extract_backbones(const string& cnf_path, vector<int>& nums) ;
    vector<pair<int, int> > get_prob_in(bool init);
    SimpleBitSet get_gain(const vector<int>& asgn);

private:
    DSU dsu;
    int original_num_var_;
    int transfer_count;
    int seed_;
    int t_wise_;
    int candidate_set_size_;
    int num_var_ = 0;
    int num_clauses_ = 0;
    int reduce_size = 150 ;
    int num_generated_testcase_;
    int selected_candidate_index_;

    double cpu_time_;
    long long num_tuple_;
    long long num_combination_all_possible_;
    long long num_tuple_all_possible_;
    long long num_tuple_all_exact_;
        
    string tmp_cnf_file;
    string cnf_instance_name_;
    string cnf_file_path_;
    string reduced_cnf_file_path_;
    string testcase_set_save_path_;
    string get_random_prefix();  
    mt19937_64 rnd_file_id_gen_;

    SimpleBitSet selected_candidate_bitset_;
    SimpleBitSet already_t_wise;

    vector<int> pq_idx;
    vector<vector<int>> uncovered_tuples;
    vector<vector<int>> clauses;
    vector<vector<int>> testcase_set_;
    vector<vector<int>> candidate_testcase_set_;
    vector<vector<int>> candidate_sample_init_solution_set_;
    vector<int> var_positive_appearance_count_;
    vector<double> var_positive_sample_weight_;
    vector<pair<vector<int>, SimpleBitSet> > pq;
    vector<int> backbones;
    vector<pair<float, float>> sample_proportion;
    vector<pair<float, float>> new_get_prob(bool init);
    vector<SimpleBitSet> testcase_SimpleBitset;
    vector<vector<int>> tri_clauses;

    unordered_set<pair<int, int>,hash_pair> double_clause_Set;
    unordered_set<pair<int, int>,hash_pair> equivalence_set;
    unordered_set<int> equivalence_number_set;
    unordered_map<int, int> varSignMap;
    unordered_map<int, int> To_original_varSignMap;
    unordered_map<int, vector<int>> equivalence_map;
    unordered_map<int, int> difference_(const unordered_map<int, int>& map1, const unordered_map<int, int>& map2); 
    
    CDCLSolver::Solver * cdcl_solver;
    ExtMinisat::SamplingSolver *cdcl_sampler;
};