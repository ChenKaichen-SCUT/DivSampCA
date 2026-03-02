#pragma once
#include "../minisat_ext/BlackBoxSolver.h"
#include "../minisat_ext/Ext.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <unordered_set>
#include <stack>
#include <vector>
#include <cstdint>
#include <cstring>
#include <immintrin.h>
#include <algorithm>
using namespace std;

class SignedDSU {
private:
    unordered_map<int,int> parent;
    unordered_map<int,int> parity; 
public:
    int find(int x) {
        if (!parent.count(x)) {
            parent[x] = x;
            parity[x] = 1;
            return x;
        }
        if (parent[x] == x) return x;

        int p = parent[x];
        int root = find(p);
        parity[x] *= parity[p];
        parent[x] = root;
        return root;
    }
    void unite(int x, int y, bool same) {
        int rx = find(x);
        int ry = find(y);
        int px = parity[x];
        int py = parity[y];

        if (rx == ry) return;
        parent[rx] = ry;
        parity[rx] = px * (same ? 1 : -1) * py;
    }

    int sign(int x) {
        find(x);
        return parity[x];
    }
};

class SimpleBitSet {
public:
    SimpleBitSet(): siz(0), cap(0), cnt(0ll) {
        arr = nullptr;
    }

    SimpleBitSet(size_t _siz): siz(_siz), cnt(0ll) {
        cap = (_siz + 63) >> 6;
        if (cap == 0) cap = 1;
        arr = new uint64_t[cap]{0};
    }

    SimpleBitSet(const SimpleBitSet& bs): siz(bs.siz), cap(bs.cap), cnt(bs.cnt) {
        arr = new uint64_t[cap]{0};
        memcpy(arr, bs.arr, sizeof(uint64_t) * cap);
    }

    SimpleBitSet(SimpleBitSet&& bs): siz(bs.siz), cap(bs.cap), cnt(bs.cnt) {
        arr = bs.arr;
        bs.arr = nullptr;
        bs.siz = 0; bs.cap = 0; bs.cnt = 0;
    }

    ~SimpleBitSet() { delete [] arr; }

    SimpleBitSet& operator=(const SimpleBitSet& bs) {
        if (this == &bs) return *this;
        delete [] arr;
        siz = bs.siz;
        cap = bs.cap;
        cnt = bs.cnt;
        arr = new uint64_t[cap]{0};
        memcpy(arr, bs.arr, sizeof(uint64_t) * cap);
        return *this;
    }

    SimpleBitSet& operator=(SimpleBitSet&& bs) {
        if (this == &bs) return *this;
        delete [] arr;
        siz = bs.siz;
        cap = bs.cap;
        cnt = bs.cnt;
        arr = bs.arr;
        bs.arr = nullptr;
        bs.siz = 0; bs.cap = 0; bs.cnt = 0;
        return *this;
    }

    bool set(size_t idx) {
        uint64_t& t1 = arr[idx >> 6];
        uint64_t mask = 1ULL << (idx & 63);
        if ((t1 & mask) == 0) {
            ++cnt;
            t1 |= mask;
            return true;
        }
        return false;
    }

    bool unset(size_t idx) {
        uint64_t& t1 = arr[idx >> 6];
        uint64_t mask = 1ULL << (idx & 63);
        if ((t1 & mask) != 0) {
            --cnt;
            t1 ^= mask;
            return true;
        }
        return false;
    }

    bool get(size_t idx) const {
        return (arr[idx >> 6] & (1ULL << (idx & 63))) != 0;
    }

    long long count() const {
        return cnt;
    }

    void difference_(const SimpleBitSet& bs) {
        for (size_t i = 0; i < cap; ++i) {
            uint64_t tmp = arr[i] & bs.arr[i];
            if (tmp != 0) {
                arr[i] ^= tmp;
                cnt -= __builtin_popcountll(tmp);
            }
        }
    }

    void bitwise_or(const SimpleBitSet& rhs) {
        if (rhs.cap != cap) return;
        for (size_t i = 0; i < cap; ++i) {
            uint64_t old = arr[i];
            arr[i] |= rhs.arr[i];
            cnt += __builtin_popcountll(arr[i] ^ old);
        }
    }

    void bitwise_xor(const SimpleBitSet& rhs) {
        if (rhs.cap != cap) return;
        for (size_t i = 0; i < cap; ++i) {
            uint64_t old = arr[i];
            arr[i] ^= rhs.arr[i];
            cnt += __builtin_popcountll(arr[i]) - __builtin_popcountll(old);
        }
    }

    SimpleBitSet bitwise_xor_inplace(const SimpleBitSet& rhs) const {
        if (rhs.cap != cap) return SimpleBitSet();
        SimpleBitSet ret(*this);
        ret.bitwise_xor(rhs);
        return ret;
    }

    void clear() {
        std::fill(arr, arr + cap, 0ULL);
        cnt = 0;
    }

    bool init_iter() {
        iter_curpos = iter_curbit = 0;
        while (iter_curpos < cap && arr[iter_curpos] == 0ULL) ++iter_curpos;
        if (iter_curpos < cap) {
            iter_curbit = __builtin_ctzll(arr[iter_curpos]);
            iter_curpos_msb = 63 - __builtin_clzll(arr[iter_curpos]);
            return true;
        }
        return false;
    }

    size_t iter_get() const {
        return (iter_curpos << 6) | iter_curbit;
    }

    bool iter_next() {
        if (iter_curpos == cap) return false;

        uint64_t tmp = arr[iter_curpos];
        if (iter_curbit < iter_curpos_msb) {
            do { ++iter_curbit; } while (((tmp >> iter_curbit) & 1ULL) == 0);
            return true;
        }

        ++iter_curpos;
        iter_curbit = 0;
        while (iter_curpos < cap && arr[iter_curpos] == 0ULL) ++iter_curpos;
        if (iter_curpos == cap) return false;

        tmp = arr[iter_curpos];
        iter_curbit = __builtin_ctzll(tmp);
        iter_curpos_msb = 63 - __builtin_clzll(tmp);
        return true;
    }

private:
    uint64_t* arr;
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

    void SetDefaultPara();
    void SetTransferSize(int transfer_size);
    void SetCandidateSetSize(int candidate_set_size); 
    inline void SetTestcaseSetSavePath(string testcase_set_path) { testcase_set_save_path_ = testcase_set_path; } 
    void GenerateCoveringArray();
    void Init();
    bool read_cnf();
    void ReduceCNF();
    void Init2TupleInfo();
    void GenerateInitTestcaseCDCL();
    void Update_already_t_wise();
    bool extract_backbones(const string& cnf_path, vector<int>& nums) ;
    void processbackbones(vector<vector<int>>& clauses, const vector<int>& backbones);
    void RecoverSolutions();
    void InitCombinationOffset();
    inline long long Get2TupleMapIndex(long i, long v_i, long j, long v_j) const noexcept {
        long long base = (v_i << 1 | v_j) * num_combination_all_possible_;
        long long pos = combination_offset[i] + j - i - 1;
        return base + pos;
    }
    void SaveTestcaseSet(vector<vector<int>>,string);
    void Reduce_redundancy();
    void remove_temp_files();
    void clear_final();
    void choose_final_plain();
    void find_uncovered_tuples();
    void Update2TupleInfo();
    void clear_pq();
    void GenerateTestcase();
    int SelectTestcaseFromCandidateSetByTupleNum();
    SimpleBitSet get_gain(const vector<int>& asgn);
    void GenerateCandidateTestcaseSet();
    vector<pair<float, float>> get_prob_in();

private:
    string cnf_file_path_;
    int seed_;
    string testcase_set_save_path_;
    int transfer_size_;
    int candidate_set_size_;
    vector<vector<int>> candidate_testcase_set_;
    string cnf_instance_name_;
    string reduced_cnf_file_path_;
    vector<vector<int>> clauses;
    void read_cnf_header(ifstream& ifs, int& nvar, int& nclauses);
    int num_var_;
    int num_clauses_;
    int original_num_var_;
    long long num_tuple_;
    long long num_combination_all_possible_;
    long long num_tuple_all_possible_;
    long long num_tuple_all_exact_;
    vector<pair<float, float>> sample_proportion;
    vector<vector<int>> testcase_set_;
    vector<vector<int>> final_testcase_set_;
    vector<int> backbones;
    SignedDSU dsu;
    vector<int> new2old;
    SimpleBitSet selected_candidate_bitset_;
    SimpleBitSet already_t_wise;
    int num_generated_testcase_;
    vector<long long> combination_offset;
    string tmp_cnf_file;
    vector<vector<int>> uncovered_tuples;
    double cpu_time_;
    vector<pair<vector<int>, SimpleBitSet> > pq;
    vector<int> pq_idx;
    int selected_candidate_index_;
    CDCLSolver::Solver * cdcl_solver;
    ExtMinisat::SamplingSolver *cdcl_sampler;
};